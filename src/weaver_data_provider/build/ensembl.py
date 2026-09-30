"""Cut one Ensembl release into gene bundles.

The inputs are the files Ensembl publishes for a release — the GFF3 annotation over every sequence of
the assembly, the GTF over the same, whose tags say whether a CDS is complete, the cDNA and ncRNA
sequence sets and the peptide set — plus the assembly's genome as `weaver-data-build genome` cut it,
NCBI's RefSeq alignment BAMs, HGNC's complete set and the MANE summary.

An Ensembl transcript is the genome spliced at its exons, so its placement is the annotation's own
coordinates, and the builder checks each record against the genome before it bundles it. A MANE
transcript Ensembl places on no chromosome — one annotated only on the fix patch that supplies its 5'
end — also gets its RefSeq partner's NCBI alignments to the chromosomes, after a check that the
sequence NCBI aligned is this record.

Bundled: every transcript with a published sequence, under its Ensembl gene, placed on the RefSeq
accession of the sequence the annotation names. Left out, and counted in the build's report:
transcripts Ensembl publishes no sequence for, by feature type, and transcripts on a sequence the
assembly no longer holds — a scaffold retired by a later patch release that Ensembl still annotates —
by sequence. Refused: a record that is not the genome spliced at its exons, a transcript on a sequence
the annotation declares no region for, a sequence with no RefSeq accession among its aliases, a
feature id or GTF row that repeats, a coding transcript with no sequence or one the GTF does not name,
a partner NCBI aligned with another sequence than this record, without a sequence, on a chromosome
the genome lacks or in a way that does not describe the record against the genome, and a protein
missing from the peptide set. A transcript whose CDS the GTF marks incomplete at either end is
bundled with `cds_undetermined` set. The reasons are in `docs/design/placements.md`.
"""

from __future__ import annotations

import collections
import dataclasses
import gzip
import pathlib
import re
import sys
import urllib.parse
from collections.abc import Iterable, Iterator

from weaver_data_provider import build
from weaver_data_provider import genome as genome_mod
from weaver_data_provider.build import common
from weaver_data_provider.v1 import bundle_pb2

_TAGS = {
    'MANE_Select': bundle_pb2.TAG_MANE_SELECT,
    'MANE_Plus_Clinical': bundle_pb2.TAG_MANE_PLUS_CLINICAL,
    'Ensembl_canonical': bundle_pb2.TAG_ENSEMBL_CANONICAL,
}
_REFSEQ_ALIAS = re.compile(r'^N[CTW]_\d+\.\d+$')
_REGION_KINDS = frozenset({'chromosome', 'scaffold', 'supercontig'})
_GENE_KINDS = frozenset({'gene', 'ncRNA_gene', 'pseudogene'})
_GTF_ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')


@dataclasses.dataclass(frozen=True)
class Release:
    """One Ensembl release on one assembly, as local files."""

    assembly: bundle_pb2.Assembly
    release: str  # Ensembl's release number, "116"
    annotation: pathlib.Path  # the GFF3 over chromosomes, patches, haplotypes and scaffolds, gzipped
    completeness: pathlib.Path  # the GTF over the same, gzipped; its cds_start_NF and cds_end_NF tags
    transcripts: tuple[pathlib.Path, ...]  # the cDNA and ncRNA FASTAs, gzipped
    proteins: pathlib.Path  # the peptide FASTA, gzipped
    genome: pathlib.Path  # the assembly as `weaver-data-build genome` cut it
    alignments: tuple[pathlib.Path, ...]  # NCBI's RefSeq alignment BAMs, each with its .bai beside it
    hgnc: pathlib.Path  # HGNC complete set, TSV
    mane: pathlib.Path  # MANE summary, gzipped TSV

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [('annotation', self.annotation), ('completeness', self.completeness)]
        roles += [(f'transcripts/{i}', path) for i, path in enumerate(self.transcripts)]
        roles += [('proteins', self.proteins), ('genome', self.genome / genome_mod.CATALOGUE)]
        roles += [(f'alignments/{i}', path) for i, path in enumerate(self.alignments)]
        return [*roles, ('hgnc', self.hgnc), ('mane', self.mane)]


@dataclasses.dataclass
class _Gene:
    gene_id: str  # "ENSG00000012048", the stable id
    symbol: str
    description: str


@dataclasses.dataclass
class _Transcript:
    accession: str  # "ENST00000357654", the stable id
    version: int
    gene_id: str
    kind: str  # the feature type: mRNA, lnc_RNA, ...
    biotype: str  # Ensembl's: protein_coding, nonsense_mediated_decay, ...
    tags: list[int]
    sequence: str  # the annotation's name for the sequence it is on
    minus: bool
    exons: list[tuple[int, int]] = dataclasses.field(default_factory=list)  # 1-based closed, as read
    cds: list[tuple[int, int]] = dataclasses.field(default_factory=list)  # 1-based closed
    protein: tuple[str, int] | None = None

    @property
    def versioned(self) -> str:
        return f'{self.accession}.{self.version}'


@dataclasses.dataclass
class _Annotation:
    """Sequences, genes and transcripts, as the GFF3's features are read in order."""

    # the annotation's own sequence name -> the RefSeq accession among its aliases
    aliases: dict[str, str] = dataclasses.field(default_factory=dict)
    genes: dict[str, _Gene] = dataclasses.field(default_factory=dict)  # by feature id, "gene:ENSG…"
    transcripts: dict[str, _Transcript] = dataclasses.field(default_factory=dict)  # by feature id, "transcript:ENST…"

    def add_region(self, name: str, attrs: dict[str, str]) -> None:
        """A sequence of the assembly: its RefSeq accession is among its aliases.

        Raises:
            build.BuildError: If no alias is a RefSeq accession.
        """
        refseq = [a for a in common.values(attrs.get('Alias', '')) if _REFSEQ_ALIAS.match(a)]
        if len(refseq) != 1:
            raise build.BuildError(
                f'{name}: {len(refseq)} RefSeq accessions among its aliases {attrs.get("Alias", "")!r}'
            )
        self.aliases[name] = refseq[0]

    def add_gene(self, attrs: dict[str, str]) -> None:
        """A gene feature, by its id.

        Raises:
            build.BuildError: If the id repeats.
        """
        if attrs['ID'] in self.genes:
            raise build.BuildError(f'{attrs["ID"]}: a gene feature id that repeats')
        self.genes[attrs['ID']] = _Gene(
            gene_id=attrs['gene_id'],
            symbol=attrs.get('Name', ''),
            description=urllib.parse.unquote(attrs.get('description', '')),
        )

    def add_transcript(self, sequence: str, kind: str, strand: str, attrs: dict[str, str]) -> None:
        """A transcript feature, by its id.

        Raises:
            build.BuildError: If the id repeats.
        """
        if attrs['ID'] in self.transcripts:
            raise build.BuildError(f'{attrs["ID"]}: a transcript feature id that repeats')
        self.transcripts[attrs['ID']] = _Transcript(
            accession=attrs['transcript_id'],
            version=int(attrs['version']),
            gene_id=attrs['Parent'],
            kind=kind,
            biotype=attrs.get('biotype', ''),
            tags=[_TAGS[t] for t in common.values(attrs.get('tag', '')) if t in _TAGS],
            sequence=sequence,
            minus=strand == '-',
        )

    def add_exon(self, span: tuple[int, int], attrs: dict[str, str]) -> None:
        self.transcripts[attrs['Parent']].exons.append(span)

    def add_cds(self, span: tuple[int, int], attrs: dict[str, str]) -> None:
        transcript = self.transcripts[attrs['Parent']]
        transcript.cds.append(span)
        if transcript.protein is None:
            transcript.protein = (attrs['protein_id'], int(attrs['version']))


def _read_annotation(path: pathlib.Path) -> _Annotation:
    """Every sequence, gene and transcript of the release, with each transcript's exons, CDS and protein.

    Raises:
        build.BuildError: If a feature line is malformed, or a sequence has no RefSeq accession.
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
            kind, span, attrs = fields[2], (int(fields[3]), int(fields[4])), common.attributes(fields[8])
            if kind in _REGION_KINDS:
                annotation.add_region(fields[0], attrs)
            elif kind in _GENE_KINDS:
                annotation.add_gene(attrs)
            elif attrs.get('ID', '').startswith('transcript:'):
                annotation.add_transcript(fields[0], kind, fields[6], attrs)
            elif kind == 'exon':
                annotation.add_exon(span, attrs)
            elif kind == 'CDS':
                annotation.add_cds(span, attrs)
    return annotation


def _read_completeness(path: pathlib.Path) -> dict[str, bool]:
    """Whether each transcript's CDS is whole, from the GTF's `cds_start_NF` and `cds_end_NF` tags.

    The GFF3 states a CDS's coordinates and phase but not whether the annotation reaches its real start
    or end; the GTF for the same release carries that as tags on the transcript row.

    Raises:
        build.BuildError: If a transcript has two rows.
    """
    whole: dict[str, bool] = {}
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] != 'transcript':
                continue
            attrs: dict[str, list[str]] = collections.defaultdict(list)
            for key, value in _GTF_ATTRIBUTE.findall(fields[8]):
                attrs[key].append(value)
            versioned = f'{attrs["transcript_id"][0]}.{attrs["transcript_version"][0]}'
            if versioned in whole:
                raise build.BuildError(f'{path}: {versioned} has two transcript rows')
            whole[versioned] = not ({'cds_start_NF', 'cds_end_NF'} & set(attrs['tag']))
    return whole


def _placement(transcript: _Transcript, sequence: str) -> common.Placement:
    """The transcript's exons as its placement: in transcript order, each a match of its length."""
    exons = sorted(transcript.exons, reverse=transcript.minus)
    out, t = [], 0
    for start, end in exons:
        n = end - start + 1
        out.append(
            bundle_pb2.Exon(
                transcript_start=t,
                transcript_end=t + n,
                genome_start=start - 1,
                genome_end_inclusive=end - 1,
                cigar=f'{n}=',
            )
        )
        t += n
    return common.Placement(sequence, transcript.minus, out)


def _spliced(placement: common.Placement, genome: genome_mod.Genome) -> bytes:
    """The genome at the placement's exons, in transcript orientation, with ambiguity codes read as N."""
    ranges = sorted((e.genome_start, e.genome_end_inclusive + 1) for e in placement.exons)
    joined = common.unambiguous(''.join(genome.fetch_many(placement.chromosome, ranges)).encode('ascii'))
    return common.reverse_complement(joined) if placement.minus else joined


def _fits(record: bytes, placement: common.Placement, genome: genome_mod.Genome) -> bool:
    """Whether an alignment describes this record against the genome: every `=` run is the record's bases there."""
    for exon in placement.exons:
        (genomic,) = genome.fetch_many(placement.chromosome, [(exon.genome_start, exon.genome_end_inclusive + 1)])
        bases = common.unambiguous(genomic.encode('ascii'))
        if placement.minus:
            bases = common.reverse_complement(bases)
        t = g = 0
        for count, op in re.findall(r'(\d+)([=XID])', exon.cigar):
            n = int(count)
            if op == '=' and record[exon.transcript_start + t : exon.transcript_start + t + n] != bases[g : g + n]:
                return False
            if op in '=XI':
                t += n
            if op in '=XD':
                g += n
        if exon.transcript_start + t != exon.transcript_end or g != len(bases):
            return False
    return True


def _cds_indices(transcript: _Transcript, placement: common.Placement) -> tuple[int, int]:
    """The CDS the annotation states, as transcript indices; Ensembl's CDS includes the stop codon."""
    low = min(c[0] for c in transcript.cds) - 1
    high = max(c[1] for c in transcript.cds) - 1
    start_g, end_g = (high, low) if transcript.minus else (low, high)
    return (
        placement.transcript_index(start_g, what=f'{transcript.versioned} CDS start'),
        placement.transcript_index(end_g, what=f'{transcript.versioned} CDS end'),
    )


@dataclasses.dataclass(frozen=True)
class _Loaded:
    """A release's inputs, read."""

    assembly: bundle_pb2.Assembly
    annotation: _Annotation
    sequences: dict[str, bytes]  # transcript records by versioned accession
    proteins: dict[str, bytes]
    genome: genome_mod.Genome
    whole_cds: dict[str, bool]  # by versioned accession, from the GTF
    hgnc: dict[str, dict[str, str]]  # by Ensembl gene id
    mane: dict[str, tuple[str, str]]  # Ensembl accession.version -> (RefSeq partner accession.version, MANE status)
    partner_placements: dict[str, list[common.Placement]]  # NCBI's chromosome alignments of each MANE partner
    without_sequence: collections.Counter[str]  # transcripts Ensembl publishes no sequence for, by feature type
    # transcripts on a sequence the genome does not hold, by that sequence's RefSeq accession
    off_assembly: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)


def _load(release: Release) -> _Loaded:
    genome = genome_mod.Genome(str(release.genome))
    if genome.assembly != release.assembly:
        name = bundle_pb2.Assembly.Name
        raise build.BuildError(
            f'{release.genome}: a genome on {name(genome.assembly)} for a {name(release.assembly)} release'
        )
    annotation = _read_annotation(release.annotation)
    sequences: dict[str, bytes] = {}
    for path in release.transcripts:
        sequences.update(common.read_fasta(path))
    mane = {row['Ensembl_nuc']: (row['RefSeq_nuc'], row['MANE_status']) for row in common.read_mane(release.mane)}
    partners = {mane[t.versioned][0] for t in annotation.transcripts.values() if t.versioned in mane}
    on_chromosomes = common.read_alignments(
        release.alignments,
        lambda name, chromosome: name in partners and common.sequence_kind(chromosome) == 'chromosome',
    )
    partner_placements: dict[str, list[common.Placement]] = collections.defaultdict(list)
    for (name, _), placement in on_chromosomes.items():
        partner_placements[name].append(placement)
    loaded = _Loaded(
        assembly=release.assembly,
        annotation=annotation,
        sequences=sequences,
        proteins=common.read_fasta(release.proteins),
        genome=genome,
        whole_cds=_read_completeness(release.completeness),
        hgnc=common.read_hgnc(release.hgnc, 'ensembl_gene_id'),
        mane=mane,
        partner_placements=dict(partner_placements),
        without_sequence=collections.Counter(
            t.kind for t in annotation.transcripts.values() if t.versioned not in sequences
        ),
    )
    _check_transcripts(loaded)
    return loaded


def _check_transcripts(loaded: _Loaded) -> None:
    """Every transcript is on a declared sequence; every coding one has a sequence and a completeness statement.

    Raises:
        build.BuildError: Where one has not.
    """
    for transcript in loaded.annotation.transcripts.values():
        accession = loaded.annotation.aliases.get(transcript.sequence)
        if accession is None:
            raise build.BuildError(
                f'{transcript.versioned}: on {transcript.sequence}, a sequence the annotation does not declare'
            )
        if accession not in loaded.genome:
            loaded.off_assembly[accession] += 1
        if transcript.cds and transcript.versioned not in loaded.sequences:
            raise build.BuildError(f'{transcript.versioned}: a coding transcript with no sequence in the sequence sets')
        if transcript.cds and transcript.versioned not in loaded.whole_cds:
            raise build.BuildError(f'{transcript.versioned}: not in the GTF, which says whether its CDS is complete')


def _report(release: Release, loaded: _Loaded, tally: _Tally) -> None:
    print(
        f'{release.release}: {len(loaded.annotation.transcripts)} transcripts, placements: '
        f'{common.counted(tally.placed)}; MANE partner alignments copied: {tally.copied}; MANE transcripts on no '
        f'chromosome whose partner NCBI aligns to none either: {tally.unpartnered}; coding transcripts whose CDS the '
        f'GTF marks incomplete: {tally.undetermined}; transcripts left out for having no published sequence: '
        f'{common.counted(loaded.without_sequence)}; left out for being on a sequence the assembly no longer holds: '
        f'{common.counted(loaded.off_assembly)}',
        file=sys.stderr,
    )


@dataclasses.dataclass
class _Tally:
    placed: collections.Counter[str] = dataclasses.field(default_factory=collections.Counter)
    copied: int = 0
    unpartnered: int = 0  # MANE transcripts on no chromosome whose partner NCBI aligns to no chromosome
    undetermined: int = 0


def _own_placement(transcript: _Transcript, residues: bytes, loaded: _Loaded) -> common.Placement:
    """The annotation's placement, checked: the record is the genome spliced at the exons.

    Raises:
        build.BuildError: If the record differs from the splice.
    """
    placement = _placement(transcript, loaded.annotation.aliases[transcript.sequence])
    spliced = _spliced(placement, loaded.genome)
    if spliced != residues:
        raise build.BuildError(
            f'{transcript.versioned}: its {len(residues)}-base record is not the genome spliced at its exons '
            f'({len(spliced)} bases); an edit the annotation does not state'
        )
    return placement


def _partner_alignments(
    transcript: _Transcript, residues: bytes, own: str, loaded: _Loaded, tally: _Tally
) -> list[bundle_pb2.Alignment]:
    """NCBI's chromosome alignments of the MANE partner, for a transcript Ensembl places on no chromosome.

    Raises:
        build.BuildError: If NCBI's alignment carries no sequence or another sequence than this record, does not
            describe it against the genome, or is on a chromosome not in the genome.
    """
    partner = loaded.mane.get(transcript.versioned)
    if partner is None or common.sequence_kind(own) == 'chromosome':
        return []
    placements = loaded.partner_placements.get(partner[0], [])
    if not placements:
        tally.unpartnered += 1
    out = []
    for placement in placements:
        if placement.query is None:
            raise build.BuildError(
                f"{transcript.versioned}: NCBI's alignment of its MANE partner {partner[0]} on {placement.chromosome} "
                'carries no sequence to check against the record'
            )
        if placement.query != residues:
            raise build.BuildError(
                f'{transcript.versioned}: the sequence NCBI aligned for its MANE partner {partner[0]} on '
                f'{placement.chromosome} is not this record; the pair is not identical'
            )
        if placement.chromosome not in loaded.genome:
            raise build.BuildError(
                f'{transcript.versioned}: its MANE partner is aligned on {placement.chromosome}, which the genome lacks'
            )
        placement.check_tiles(transcript.versioned, len(residues))
        if not _fits(residues, placement, loaded.genome):
            raise build.BuildError(
                f'{transcript.versioned}: its MANE partner {partner[0]} aligns on {placement.chromosome} in a way that '
                'does not describe this record against the genome'
            )
        out.append(placement.message(loaded.assembly, bundle_pb2.ALIGNMENT_SOURCE_MANE_PARTNER))
    return out


def _add_transcript(
    bundle: bundle_pb2.GeneBundle,
    digests: dict[str, bundle_pb2.Sequence],
    transcript: _Transcript,
    loaded: _Loaded,
    tally: _Tally,
) -> None:
    residues = loaded.sequences[transcript.versioned]
    own = _own_placement(transcript, residues, loaded)
    message = bundle.transcripts.add(
        accession=transcript.accession,
        version=transcript.version,
        publisher=bundle_pb2.PUBLISHER_ENSEMBL,
        status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
        biotype=transcript.biotype,
        tags=transcript.tags,
        sequence_digest=common.sequence(digests, residues, bundle_pb2.ALPHABET_NUCLEOTIDE),
    )
    partner = loaded.mane.get(transcript.versioned)
    if partner is not None:
        message.mane_partner = partner[0]
        tag = common.MANE_STATUS_TAGS.get(partner[1])
        if tag is not None and tag not in message.tags:
            message.tags.append(tag)
    alignments = [own.message(loaded.assembly, bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION)]
    alignments.extend(_partner_alignments(transcript, residues, own.chromosome, loaded, tally))
    tally.copied += len(alignments) - 1
    alignments.sort(key=lambda a: common.placement_order(a.chromosome))
    message.alignments.extend(alignments)
    tally.placed.update(common.sequence_kind(a.chromosome) for a in alignments)
    if transcript.cds:
        if loaded.whole_cds[transcript.versioned]:
            message.cds.start_index, message.cds.end_index_inclusive = _cds_indices(transcript, own)
        else:
            message.cds_undetermined = True
            tally.undetermined += 1
    if transcript.protein is not None:
        accession, version = transcript.protein
        protein = loaded.proteins.get(f'{accession}.{version}')
        if protein is None:
            raise build.BuildError(
                f'{transcript.versioned}: its protein {accession}.{version} is not in the peptide set'
            )
        message.protein_accession, message.protein_version = accession, version
        bundle.proteins.add(
            accession=accession,
            version=version,
            sequence_digest=common.sequence(digests, protein, bundle_pb2.ALPHABET_PROTEIN),
            transcript_accession=transcript.accession,
            transcript_version=transcript.version,
        )


def _gene_bundle(gene: _Gene, members: Iterable[_Transcript], loaded: _Loaded, tally: _Tally) -> bundle_pb2.GeneBundle:
    hgnc = loaded.hgnc.get(gene.gene_id)
    message = common.gene_message(
        feature_id=gene.gene_id,
        symbol=gene.symbol or gene.gene_id,
        hgnc_id='',
        ncbi_gene_id='',
        ensembl_gene_id=gene.gene_id,
        description=gene.description,
        synonyms=(),
        hgnc=hgnc,
    )
    bundle = bundle_pb2.GeneBundle(gene=message)
    digests: dict[str, bundle_pb2.Sequence] = {}
    for transcript in sorted(members, key=lambda t: (t.accession, t.version)):
        _add_transcript(bundle, digests, transcript, loaded, tally)
    bundle.sequences.extend(digests[d] for d in sorted(digests))
    return bundle


def bundles(release: Release) -> Iterator[bundle_pb2.GeneBundle]:
    """Every gene's bundle in the release, in a stable order: by symbol, then by Ensembl gene id.

    Raises:
        build.BuildError: On any of the refusals the module docstring lists.
    """
    loaded = _load(release)
    by_gene: dict[str, list[_Transcript]] = collections.defaultdict(list)
    for transcript in loaded.annotation.transcripts.values():
        if transcript.versioned not in loaded.sequences:
            continue
        if loaded.annotation.aliases[transcript.sequence] not in loaded.genome:
            continue
        if transcript.gene_id not in loaded.annotation.genes:
            raise build.BuildError(f'{transcript.gene_id}: transcripts name a gene the annotation does not declare')
        by_gene[transcript.gene_id].append(transcript)
    tally = _Tally()
    out: list[tuple[str, str, bundle_pb2.GeneBundle]] = []
    for gene_id, members in by_gene.items():
        gene = loaded.annotation.genes[gene_id]
        bundle = _gene_bundle(gene, members, loaded, tally)
        out.append((bundle.gene.symbol, gene.gene_id, bundle))
    for _, _, bundle in sorted(out, key=lambda item: (item[0], item[1])):
        yield bundle
    _report(release, loaded, tally)
