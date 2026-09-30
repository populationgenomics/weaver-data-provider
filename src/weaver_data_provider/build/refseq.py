"""Cut one RefSeq annotation release into gene bundles.

The inputs are the files NCBI publishes for an annotated assembly — the GFF3 annotation, the
transcript and protein sequence sets, and its alignments of every RefSeq transcript to the assembly —
plus HGNC's complete set, for gene identity, and the MANE summary, for ranking.

Bundled: every transcript the annotation names with a versioned accession, with each placement NCBI
publishes for it on any sequence of the assembly, the CDS the complete placements state, and its
protein. Left out, and counted by kind in the build's report: RNA features naming no transcript
accession. Refused, so that no bundle carries a placement or a CDS the builder made up: a coding
transcript with no alignment, an alignment that does not tile its record or lies outside the
annotation's span, a CDS bound on a genome-only base, placements whose CDS bounds disagree with the
record's protein confirming none of them, and an annotation whose features break the shape this
reader relies on. The reasons are in `docs/design/placements.md`.
"""

from __future__ import annotations

import collections
import csv
import dataclasses
import gzip
import pathlib
import re
import sys
import urllib.parse
from collections.abc import Iterable, Iterator

import pysam.libcalignedsegment
import pysam.libcalignmentfile

from weaver_data_provider import build, refget
from weaver_data_provider.v1 import bundle_pb2

_TAGS = {
    'MANE Select': bundle_pb2.TAG_MANE_SELECT,
    'MANE Plus Clinical': bundle_pb2.TAG_MANE_PLUS_CLINICAL,
    'RefSeq Select': bundle_pb2.TAG_REFSEQ_SELECT,
}
_ACCESSION = re.compile(r'^([A-Z]{2}_\d+)\.(\d+)$')
_CIGAR_OPS = {7: '=', 8: 'X', 1: 'I', 2: 'D', 4: 'I'}  # soft-clipped transcript bases have no genome counterpart
_BASES = 'TCAG'
_AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_CODONS = {
    a + b + c: _AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_HGNC_COLUMNS = frozenset({'hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id'})
_MANE_COLUMNS = frozenset({'RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'})


@dataclasses.dataclass(frozen=True)
class Release:
    """One RefSeq annotation release on one assembly, as local files."""

    assembly: bundle_pb2.Assembly
    release: str  # the publisher's identifier, "RS_2024_08"
    annotation: pathlib.Path  # GFF3, gzipped
    transcripts: pathlib.Path  # RNA FASTA, gzipped
    proteins: pathlib.Path  # protein FASTA, gzipped
    alignments: tuple[pathlib.Path, ...]  # known and model alignment BAMs, each with its .bai beside it
    hgnc: pathlib.Path  # HGNC complete set, TSV
    mane: pathlib.Path  # MANE summary, gzipped TSV

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [('annotation', self.annotation), ('transcripts', self.transcripts), ('proteins', self.proteins)]
        roles += [(f'alignments/{i}', path) for i, path in enumerate(self.alignments)]
        return [*roles, ('hgnc', self.hgnc), ('mane', self.mane)]


@dataclasses.dataclass
class _Gene:
    gene_id: str
    symbol: str
    dbxrefs: dict[str, str]
    synonyms: list[str]
    description: str

    @property
    def identity(self) -> str:
        """The GeneID, shared by every feature of one gene, including its copy on a patch or scaffold.

        Raises:
            build.BuildError: If the feature carries no GeneID.
        """
        gene_id = self.dbxrefs.get('GeneID')
        if gene_id is None:
            raise build.BuildError(f'{self.gene_id}: a gene feature with no GeneID in its Dbxref')
        return gene_id


@dataclasses.dataclass(frozen=True)
class _CdsSegment:
    """One CDS row: its 1-based closed span, and whether the annotation leaves either end open.

    An open end (`start_range=.,N` or `end_range=N,.`) is one the record's coding sequence runs past:
    the sequence lacks the bases.
    """

    start: int
    end: int
    start_open: bool
    end_open: bool


@dataclasses.dataclass
class _AnnotatedPlacement:
    """One placement of a transcript as the annotation states it: its span on one sequence, and its CDS there."""

    sequence: str
    span: tuple[int, int]  # 1-based closed
    cds: list[_CdsSegment] = dataclasses.field(default_factory=list)

    @property
    def determines_cds(self) -> bool:
        """Whether this placement states both bounds of the record's CDS: it has rows and its outer ends are closed."""
        if not self.cds:
            return False
        lowest = min(self.cds, key=lambda c: c.start)
        highest = max(self.cds, key=lambda c: c.end)
        return not (lowest.start_open or highest.end_open)


@dataclasses.dataclass
class _Transcript:
    accession: str
    version: int
    gene_id: str
    biotype: str
    tags: list[int]
    placements: dict[str, _AnnotatedPlacement]  # by sequence
    protein: tuple[str, int] | None = None

    @property
    def coding(self) -> bool:
        return any(p.cds for p in self.placements.values())

    def determining(self, aligned: Iterable[str]) -> list[str]:
        """The sequences among `aligned` whose placement states the whole CDS: the ones the CDS is projected through."""
        return [c for c in aligned if c in self.placements and self.placements[c].determines_cds]


@dataclasses.dataclass(frozen=True)
class _Placement:
    chromosome: str  # the assembly sequence: a chromosome, alternate locus or patch
    minus: bool
    exons: list[bundle_pb2.Exon]


def _sequence_kind(accession: str) -> str:
    """`chromosome`, `scaffold` (an alternate locus or an unplaced one) or `patch`, from the RefSeq prefix."""
    kind = build.SEQUENCE_KINDS.get(accession[:3])
    if kind is None:
        raise build.BuildError(f'{accession}: not a sequence of an NCBI assembly (NC_, NT_ or NW_)')
    return kind


def _placement_order(chromosome: str) -> tuple[int, str]:
    """The order a transcript's alignments are kept in, which the proto promises: by sequence kind, then accession."""
    return (list(build.SEQUENCE_KINDS.values()).index(_sequence_kind(chromosome)), chromosome)


def _split(accession: str) -> tuple[str, int]:
    match = _ACCESSION.match(accession)
    if match is None:
        raise build.BuildError(f'{accession!r} is not a versioned accession')
    return match.group(1), int(match.group(2))


def _attributes(field: str) -> dict[str, str]:
    return dict(kv.split('=', 1) for kv in field.split(';') if '=' in kv)


def _values(attribute: str) -> list[str]:
    """A GFF3 attribute's comma-separated values, each percent-decoded after the split."""
    return [urllib.parse.unquote(v) for v in attribute.split(',') if v]


def _dbxrefs(value: str) -> dict[str, str]:
    """`GeneID:672,HGNC:HGNC:1100,MIM:113705` as a map from database to id (`HGNC` -> `HGNC:1100`)."""
    out: dict[str, str] = {}
    for ref in _values(value):
        if ':' in ref:
            db, ident = ref.split(':', 1)
            out.setdefault(db, ident)
    return out


@dataclasses.dataclass
class _Annotation:
    """Genes and the transcripts placed on the assembly's sequences, as the GFF3's features are read in order."""

    genes: dict[str, _Gene] = dataclasses.field(default_factory=dict)
    transcripts: dict[str, _Transcript] = dataclasses.field(default_factory=dict)  # by versioned accession
    # RNA features naming no transcript accession, by feature kind
    without_accession: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)
    _by_rna_id: dict[str, tuple[_Transcript, _AnnotatedPlacement]] = dataclasses.field(default_factory=dict)

    def add_gene(self, attrs: dict[str, str]) -> None:
        """A gene feature, by its feature id; every placement of a gene is one."""
        self.genes[attrs['ID']] = _Gene(
            gene_id=attrs['ID'],
            symbol=attrs.get('Name', ''),
            dbxrefs=_dbxrefs(attrs.get('Dbxref', '')),
            synonyms=_values(attrs.get('gene_synonym', '')),
            description=urllib.parse.unquote(attrs.get('description', '')),
        )

    def add_rna(self, chromosome: str, kind: str, span: tuple[int, int], attrs: dict[str, str]) -> None:
        """A transcript's placement on one sequence; a further feature for a transcript already seen adds one.

        Raises:
            build.BuildError: If the transcript is placed twice on one sequence, which NCBI does not do.
        """
        if 'transcript_id' not in attrs:
            self.without_accession[kind] += 1
            return
        versioned = attrs['transcript_id']
        record = self.transcripts.get(versioned)
        if record is None:
            accession, version = _split(versioned)
            record = self.transcripts[versioned] = _Transcript(
                accession=accession,
                version=version,
                gene_id=attrs['Parent'],
                biotype=kind,
                tags=[_TAGS[t] for t in _values(attrs.get('tag', '')) if t in _TAGS],
                placements={},
            )
        if chromosome in record.placements:
            raise build.BuildError(f'{versioned}: placed twice on {chromosome}')
        placement = record.placements[chromosome] = _AnnotatedPlacement(chromosome, span)
        self._by_rna_id[attrs['ID']] = (record, placement)

    def add_cds(self, span: tuple[int, int], attrs: dict[str, str]) -> None:
        """A CDS segment of the placement its parent RNA feature is."""
        placed = self._by_rna_id.get(attrs.get('Parent', ''))
        if placed is None:
            return
        record, placement = placed
        placement.cds.append(
            _CdsSegment(span[0], span[1], start_open='start_range' in attrs, end_open='end_range' in attrs)
        )
        if record.protein is None and 'protein_id' in attrs:
            record.protein = _split(attrs['protein_id'])


def _read_annotation(path: pathlib.Path) -> _Annotation:
    """Genes and every transcript placement, with each placement's CDS spans and the protein id.

    A transcript placed more than once — on X and Y, or on a chromosome and a patch — is one record
    with a span per sequence; each further placement's features carry a suffixed id
    (`rna-NM_000451.4-2`), and its CDS features are the bounds on that sequence.
    """
    annotation = _Annotation()
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields == ['']:
                continue
            if len(fields) != 9:
                raise build.BuildError(f'{path}: a feature line with {len(fields)} fields, not 9: {line[:80]!r}')
            kind, span, attrs = fields[2], (int(fields[3]), int(fields[4])), _attributes(fields[8])
            if kind in ('gene', 'pseudogene'):
                annotation.add_gene(attrs)
            elif attrs.get('Parent', '').startswith('gene-') and kind != 'exon':
                annotation.add_rna(fields[0], kind, span, attrs)
            elif kind == 'CDS':
                annotation.add_cds(span, attrs)
    return annotation


def _read_fasta(path: pathlib.Path) -> dict[str, bytes]:
    out: dict[str, bytes] = {}
    name, chunks = None, []
    with gzip.open(path, 'rt', encoding='ascii') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(chunks).upper().encode('ascii')
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        out[name] = ''.join(chunks).upper().encode('ascii')
    return out


# An exon in genome order: genome start and inclusive end, query start and end, and its cigar ops.
_GenomeOrderExon = tuple[int, int, int, int, list[tuple[int, str]]]


def _genome_order_exons(read: pysam.libcalignedsegment.AlignedSegment) -> tuple[list[_GenomeOrderExon], int]:
    """A spliced alignment cut at its introns, left to right on the genome, with the query's length."""
    if read.cigartuples is None:
        raise build.BuildError(f'{read.query_name}: a primary alignment on {read.reference_name} with no cigar')
    exons: list[_GenomeOrderExon] = []
    ops: list[tuple[int, str]] = []
    g = read.reference_start
    q = 0
    exon_g0, exon_q0 = g, q
    for op, n in read.cigartuples:
        if op == 3:  # N: intron
            exons.append((exon_g0, g - 1, exon_q0, q, ops))
            g += n
            ops, exon_g0, exon_q0 = [], g, q
            continue
        if op not in _CIGAR_OPS:
            raise build.BuildError(f'{read.query_name}: cigar op {op} is not one an NCBI transcript alignment uses')
        ch = _CIGAR_OPS[op]
        if ch in '=XD':
            g += n
        if ch in '=XI':
            q += n
        if ops and ops[-1][1] == ch:
            ops[-1] = (ops[-1][0] + n, ch)
        else:
            ops.append((n, ch))
    exons.append((exon_g0, g - 1, exon_q0, q, ops))
    return exons, q


def _exons_of(read: pysam.libcalignedsegment.AlignedSegment) -> list[bundle_pb2.Exon]:
    """A spliced alignment as exons with cigars in transcript orientation.

    The BAM cigar runs left to right on the genome; for a minus-strand transcript the exon order and
    each exon's ops are reversed so transcript coordinates ascend 5' to 3'.
    """
    exons, tx_len = _genome_order_exons(read)
    out = []
    for g0, g1, q0, q1, genome_order in exons:
        if read.is_reverse:
            t0, t1, exon_ops = tx_len - q1, tx_len - q0, list(reversed(genome_order))
        else:
            t0, t1, exon_ops = q0, q1, genome_order
        out.append(
            bundle_pb2.Exon(
                transcript_start=t0,
                transcript_end=t1,
                genome_start=g0,
                genome_end_inclusive=g1,
                cigar=''.join(f'{n}{ch}' for n, ch in exon_ops),
            )
        )
    out.sort(key=lambda e: e.transcript_start)
    return out


def _read_alignments(
    paths: Iterable[pathlib.Path], transcripts: dict[str, _Transcript]
) -> dict[tuple[str, str], _Placement]:
    """Each transcript's alignment on each sequence the annotation places it on, first seen wins."""
    placements: dict[tuple[str, str], _Placement] = {}
    for path in paths:
        with pysam.libcalignmentfile.AlignmentFile(str(path)) as bam:
            for read in bam:
                name = read.query_name
                if name is None or read.is_secondary or read.is_supplementary:
                    continue
                record = transcripts.get(name)
                chromosome = read.reference_name
                if record is None or chromosome not in record.placements or (name, chromosome) in placements:
                    continue
                placements[name, chromosome] = _Placement(chromosome, read.is_reverse, _exons_of(read))
    return placements


def _transcript_index(placement: _Placement, genome_position: int, *, what: str) -> int:
    """The 0-based transcript index a 0-based genome position projects to, through the exon cigars.

    Raises:
        build.BuildError: If the position is in no exon, or is a genome-only base (a `D` run), which no
            transcript index names; NCBI never puts a CDS bound on one.
    """
    for exon in placement.exons:
        if not exon.genome_start <= genome_position <= exon.genome_end_inclusive:
            continue
        offset = exon.genome_end_inclusive - genome_position if placement.minus else genome_position - exon.genome_start
        t = g = 0
        for count, op in re.findall(r'(\d+)([=XID])', exon.cigar):
            n = int(count)
            if op in '=X':
                if g + n > offset:
                    return exon.transcript_start + t + (offset - g)
                g += n
                t += n
            elif op == 'D':
                if g + n > offset:
                    raise build.BuildError(
                        f'{what}: genome position {genome_position} is a genome-only base in the alignment'
                    )
                g += n
            else:
                t += n
    raise build.BuildError(f'{what}: genome position {genome_position} lies in no exon of the alignment')


def _columns(reader: csv.DictReader[str], needed: frozenset[str], *, named: str) -> None:
    if missing := needed - set(reader.fieldnames or ()):
        raise build.BuildError(f'{named}: columns {sorted(missing)} are not in the table; its schema has changed')


def _read_hgnc(path: pathlib.Path) -> dict[str, dict[str, str]]:
    """HGNC's complete set keyed by NCBI GeneID."""
    with path.open(encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        _columns(rows, _HGNC_COLUMNS, named='hgnc')
        return {row['entrez_id']: row for row in rows if row['entrez_id']}


def _read_mane(path: pathlib.Path) -> dict[str, tuple[str, str]]:
    """RefSeq transcript accession.version -> (Ensembl partner accession.version, MANE status)."""
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        _columns(rows, _MANE_COLUMNS, named='mane')
        return {row['RefSeq_nuc']: (row['Ensembl_nuc'], row['MANE_status']) for row in rows}


def _preferred(hgnc: dict[str, str] | None, column: str, annotation: str) -> str:
    """HGNC's value where it has the gene and states one; the annotation's otherwise."""
    if hgnc is not None and hgnc[column]:
        return hgnc[column]
    return annotation


def _gene_message(gene: _Gene, hgnc: dict[str, str] | None) -> bundle_pb2.Gene:
    message = bundle_pb2.Gene(
        symbol=_preferred(hgnc, 'symbol', gene.symbol),
        hgnc_id=_preferred(hgnc, 'hgnc_id', gene.dbxrefs.get('HGNC', '')),
        ncbi_gene_id=gene.dbxrefs.get('GeneID', ''),
        ensembl_gene_id=_preferred(hgnc, 'ensembl_gene_id', gene.dbxrefs.get('Ensembl', '')),
        name=_preferred(hgnc, 'name', gene.description),
    )
    if hgnc is not None:
        message.previous_symbols.extend(s for s in hgnc['prev_symbol'].split('|') if s)
        message.alias_symbols.extend(s for s in hgnc['alias_symbol'].split('|') if s)
    for synonym in gene.synonyms:
        if synonym not in message.alias_symbols and synonym != message.symbol:
            message.alias_symbols.append(synonym)
    if not message.symbol:
        raise build.BuildError(f'{gene.gene_id}: no symbol in the annotation or HGNC')
    return message


def _sequence(digests: dict[str, bundle_pb2.Sequence], residues: bytes, alphabet: bundle_pb2.Alphabet) -> str:
    """The digest of `residues`, adding the sequence to the bundle's set on first sight."""
    digest = refget.digest(residues)
    if digest not in digests:
        digests[digest] = bundle_pb2.Sequence(digest=digest, alphabet=alphabet, residues=residues)
    return digest


def _check_placement(versioned: str, record: _Transcript, placement: _Placement, length: int) -> None:
    """An alignment has to tile the record's sequence and lie inside the annotation's span for it."""
    exons = placement.exons
    if exons[0].transcript_start != 0 or exons[-1].transcript_end != length:
        raise build.BuildError(
            f'{versioned}: the alignment on {placement.chromosome} covers transcript bases '
            f'{exons[0].transcript_start}..{exons[-1].transcript_end} of a {length}-base record'
        )
    low, high = record.placements[placement.chromosome].span
    if min(e.genome_start for e in exons) + 1 < low or max(e.genome_end_inclusive for e in exons) + 1 > high:
        raise build.BuildError(
            f'{versioned}: the alignment on {placement.chromosome} lies outside the annotation span {low}..{high}'
        )


def _alignment(assembly: bundle_pb2.Assembly, placement: _Placement) -> bundle_pb2.Alignment:
    return bundle_pb2.Alignment(
        assembly=assembly,
        chromosome=placement.chromosome,
        strand=bundle_pb2.STRAND_MINUS if placement.minus else bundle_pb2.STRAND_PLUS,
        source=bundle_pb2.ALIGNMENT_SOURCE_NCBI_BAM,
        exons=placement.exons,
    )


def _encodes(cds: bytes, protein: bytes) -> bool:
    """Whether `cds` is whole codons reading, in the standard genetic code, as `protein` and then a stop.

    Strict on purpose: a selenocysteine, a non-AUG initiator or an ambiguous base reads differently,
    and a record with one of those is arbitrated by nothing here.
    """
    text = cds.decode('ascii')
    translated = ''.join(_CODONS.get(text[i : i + 3], 'X') for i in range(0, len(text), 3))
    return translated == protein.decode('ascii') + '*'  # a trailing partial codon reads X, never the stop


def _cds(versioned: str, annotated: _AnnotatedPlacement, placement: _Placement) -> bundle_pb2.Cds:
    """The CDS as transcript indices, projected through this placement from the bounds stated on its sequence."""
    low = min(c.start for c in annotated.cds) - 1
    high = max(c.end for c in annotated.cds) - 1
    first, last = (high, low) if placement.minus else (low, high)
    return bundle_pb2.Cds(
        start_index=_transcript_index(placement, first, what=f'{versioned} CDS start on {placement.chromosome}'),
        end_index_inclusive=_transcript_index(placement, last, what=f'{versioned} CDS end on {placement.chromosome}'),
    )


@dataclasses.dataclass(frozen=True)
class _Loaded:
    """A release's inputs, read."""

    assembly: bundle_pb2.Assembly
    annotation: _Annotation
    sequences: dict[str, bytes]  # transcript records by versioned accession
    proteins: dict[str, bytes]
    placements: dict[tuple[str, str], _Placement]  # by (versioned accession, sequence)
    hgnc: dict[str, dict[str, str]]  # by NCBI GeneID
    mane: dict[str, tuple[str, str]]
    # placements whose CDS projected differently from the one the record's protein confirms, as bundling finds them
    disagreeing: list[str] = dataclasses.field(default_factory=list)


def _counted(counts: collections.Counter[str]) -> str:
    """`12 (tRNA 8, rRNA 4)`, most common first; `0` when empty."""
    if not counts:
        return '0'
    return f'{counts.total()} (' + ', '.join(f'{kind} {n}' for kind, n in counts.most_common()) + ')'


def _report(release: Release, loaded: _Loaded) -> None:
    annotation = loaded.annotation
    aligned: dict[str, list[str]] = collections.defaultdict(list)
    for versioned, chromosome in loaded.placements:
        aligned[versioned].append(chromosome)
    unplaced = collections.Counter(r.biotype for t, r in annotation.transcripts.items() if t not in aligned)
    by_kind = collections.Counter(_sequence_kind(chromosome) for _, chromosome in loaded.placements)
    undetermined = sum(r.coding and not r.determining(aligned[t]) for t, r in annotation.transcripts.items())
    print(
        f'{release.release}: {len(annotation.transcripts)} transcripts, '
        f'placements aligned: {_counted(by_kind)}; bundled without a placement: {_counted(unplaced)}; '
        f'coding transcripts with no placement stating their whole CDS: {undetermined}; '
        f'RNA features left out for naming no transcript accession: {_counted(annotation.without_accession)}',
        file=sys.stderr,
    )


def _load(release: Release) -> _Loaded:
    annotation = _read_annotation(release.annotation)
    loaded = _Loaded(
        assembly=release.assembly,
        annotation=annotation,
        sequences=_read_fasta(release.transcripts),
        proteins=_read_fasta(release.proteins),
        placements=_read_alignments(release.alignments, annotation.transcripts),
        hgnc=_read_hgnc(release.hgnc),
        mane=_read_mane(release.mane),
    )
    _report(release, loaded)
    return loaded


def _check_coverage(loaded: _Loaded) -> None:
    """Every coding transcript has an alignment, and every transcript a record."""
    transcripts = loaded.annotation.transcripts
    unaligned = sorted(
        t for t, r in transcripts.items() if r.coding and not any((t, c) in loaded.placements for c in r.placements)
    )
    if unaligned:
        raise build.BuildError(
            f'{len(unaligned)} coding transcripts the annotation places have no alignment, '
            f'e.g. {unaligned[:5]}; a bundle would have to invent one'
        )
    if missing := sorted(t for t in transcripts if t not in loaded.sequences):
        raise build.BuildError(f'{len(missing)} transcripts have no sequence in the transcript set, e.g. {missing[:5]}')


def _protein_residues(versioned: str, protein: tuple[str, int], loaded: _Loaded) -> bytes:
    """A transcript's protein from the protein set.

    Raises:
        build.BuildError: If the set lacks it.
    """
    residues = loaded.proteins.get(f'{protein[0]}.{protein[1]}')
    if residues is None:
        raise build.BuildError(f'{versioned}: its protein {protein[0]}.{protein[1]} is not in the protein set')
    return residues


def _add_protein(
    bundle: bundle_pb2.GeneBundle, digests: dict[str, bundle_pb2.Sequence], record: _Transcript, loaded: _Loaded
) -> None:
    if record.protein is None:
        return
    accession, version = record.protein
    residues = _protein_residues(f'{record.accession}.{record.version}', record.protein, loaded)
    bundle.proteins.add(
        accession=accession,
        version=version,
        sequence_digest=_sequence(digests, residues, bundle_pb2.ALPHABET_PROTEIN),
        transcript_accession=record.accession,
        transcript_version=record.version,
    )


def _add_transcript(
    bundle: bundle_pb2.GeneBundle, digests: dict[str, bundle_pb2.Sequence], record: _Transcript, loaded: _Loaded
) -> None:
    """One transcript into its gene's bundle: its record, MANE status, alignments, CDS and protein."""
    versioned = f'{record.accession}.{record.version}'
    residues = loaded.sequences[versioned]
    transcript = bundle.transcripts.add(
        accession=record.accession,
        version=record.version,
        publisher=bundle_pb2.PUBLISHER_REFSEQ,
        status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
        biotype=record.biotype,
        tags=record.tags,
        sequence_digest=_sequence(digests, residues, bundle_pb2.ALPHABET_NUCLEOTIDE),
    )
    partner = loaded.mane.get(versioned)
    if partner is not None:
        transcript.mane_partner = partner[0]
        tag = _TAGS.get(partner[1])
        if tag is not None and tag not in transcript.tags:
            transcript.tags.append(tag)
    projected: dict[tuple[int, int], list[str]] = {}  # each CDS the determining placements project to, by sequence
    for chromosome in sorted(record.placements, key=_placement_order):
        placement = loaded.placements.get((versioned, chromosome))
        if placement is None:
            continue
        _check_placement(versioned, record, placement, len(residues))
        transcript.alignments.append(_alignment(loaded.assembly, placement))
        annotated = record.placements[chromosome]
        if annotated.determines_cds:
            cds = _cds(versioned, annotated, placement)
            projected.setdefault((cds.start_index, cds.end_index_inclusive), []).append(chromosome)
    if projected:
        start, end = _agreed_cds(versioned, record, projected, residues, loaded)
        transcript.cds.start_index, transcript.cds.end_index_inclusive = start, end
    transcript.cds_undetermined = record.coding and not projected
    if record.protein is not None:
        transcript.protein_accession, transcript.protein_version = record.protein
    _add_protein(bundle, digests, record, loaded)


def _agreed_cds(
    versioned: str,
    record: _Transcript,
    projected: dict[tuple[int, int], list[str]],
    residues: bytes,
    loaded: _Loaded,
) -> tuple[int, int]:
    """The CDS the determining placements project to; where they differ, the one the record encodes its protein over.

    Raises:
        build.BuildError: If the placements disagree and the record names no protein, or encodes it over none of the
            projections or more than one.
    """
    if len(projected) == 1:
        return next(iter(projected))
    described = ', '.join(f'{a}..{b} on {", ".join(seqs)}' for (a, b), seqs in sorted(projected.items()))
    if record.protein is None:
        raise build.BuildError(
            f'{versioned}: the CDS projects to different transcript indices ({described}), and no protein is named '
            'to arbitrate'
        )
    protein = _protein_residues(versioned, record.protein, loaded)
    confirmed = [bounds for bounds in projected if _encodes(residues[bounds[0] : bounds[1] + 1], protein)]
    if len(confirmed) != 1:
        raise build.BuildError(
            f'{versioned}: the CDS projects to different transcript indices ({described}), and in the standard '
            f'genetic code the record encodes its protein over {"none" if not confirmed else "more than one"} of them'
        )
    loaded.disagreeing.extend(
        f'{versioned} on {seq}' for bounds, seqs in projected.items() if bounds != confirmed[0] for seq in seqs
    )
    return confirmed[0]


def _gene_bundle(gene: _Gene, members: Iterable[_Transcript], loaded: _Loaded) -> bundle_pb2.GeneBundle:
    bundle = bundle_pb2.GeneBundle(gene=_gene_message(gene, loaded.hgnc.get(gene.dbxrefs.get('GeneID', ''))))
    digests: dict[str, bundle_pb2.Sequence] = {}
    for record in sorted(members, key=lambda r: (r.accession, r.version)):
        _add_transcript(bundle, digests, record, loaded)
    bundle.sequences.extend(digests[d] for d in sorted(digests))
    return bundle


def bundles(release: Release) -> Iterator[bundle_pb2.GeneBundle]:
    """Every gene's bundle in the release, in a stable order: by symbol, then by GeneID.

    Raises:
        build.BuildError: On any of the refusals the module docstring lists.
    """
    loaded = _load(release)
    _check_coverage(loaded)
    genes: dict[str, _Gene] = {}  # by GeneID, the feature read first
    by_gene: dict[str, list[_Transcript]] = collections.defaultdict(list)
    for record in loaded.annotation.transcripts.values():
        gene = loaded.annotation.genes.get(record.gene_id)
        if gene is None:
            raise build.BuildError(f'{record.gene_id}: transcripts name a gene the annotation does not declare')
        genes.setdefault(gene.identity, gene)
        by_gene[gene.identity].append(record)
    out: list[tuple[str, str, bundle_pb2.GeneBundle]] = []
    for identity, members in by_gene.items():
        bundle = _gene_bundle(genes[identity], members, loaded)
        out.append((bundle.gene.symbol, identity, bundle))
    for _, _, bundle in sorted(out, key=lambda item: (item[0], item[1])):
        yield bundle
    if loaded.disagreeing:
        print(
            f"{release.release}: {len(loaded.disagreeing)} placements project a CDS the record's protein does not "
            f'confirm, e.g. {loaded.disagreeing[:5]}; their alignments are kept, their CDS is not',
            file=sys.stderr,
        )
