"""Cut one RefSeq annotation release, or NCBI's historical set of retired transcripts, into gene bundles.

The inputs are the files NCBI publishes for an annotated assembly — the GFF3 annotation, the
transcript and protein sequence sets, the transcripts' GenBank records, and its alignments of every
RefSeq transcript to the assembly — plus HGNC's complete set, for gene identity, and the MANE summary,
for ranking. NCBI's historical set — every `NM_`/`NR_` version it has ever aligned to GRCh38, about
half of them since replaced or suppressed, each with the alignment it last had — comes as the same
GFF3, BAM and GenBank files but no sequence sets; the records supply the sequences and the proteins,
and a table of each version's status in Entrez, which the records cannot state, says which are
retired. A retired version is bundled without tags or a MANE partner, since those name what the
publisher recommends today. A few of the set's oldest records state their CDS in a form current
records never do — a `complement` location, a closed start read from its third base — and are left
out and counted, where the same in a release fails the build: the set is an archive, and a record
the reader declines to interpret is better absent than guessed at.

Bundled: every transcript the annotation names with a versioned accession, with each placement NCBI
publishes for it on any sequence of the assembly, the CDS its GenBank record states, and its protein.
The annotation states the CDS again on each placement, in that sequence's coordinates; each placement
that states it whole is projected through its alignment as a cross-check, and one that lands elsewhere
is a `build.BuildWarning` naming both. An end of the CDS the record marks open — `<` or `>`, the CDS
running off the record — is bundled as open, so that no position is numbered from it. Left out, and
counted by kind in the build's report: RNA features naming no transcript accession. Refused, so that no
bundle carries a placement or a CDS the builder made up: a coding transcript with no alignment, one the
annotation gives a CDS and its record none, an alignment that does not tile its record or lies outside
the annotation's span, a record CDS outside its sequence, and an annotation or record file whose shape
breaks what this reader relies on. The reasons are in `docs/design/placements.md`.
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
import warnings
from collections.abc import Iterable, Iterator

from weaver_data_provider import build
from weaver_data_provider.build import common, flatfile
from weaver_data_provider.v1 import bundle_pb2

_TAGS = {**common.MANE_STATUS_TAGS, 'RefSeq Select': bundle_pb2.TAG_REFSEQ_SELECT}
# Entrez's status words for a RefSeq record (`build/entrez.py`); any other word fails the build.
_STATUSES = {
    'live': bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
    'replaced': bundle_pb2.TRANSCRIPT_STATUS_SUPERSEDED,
    'suppressed': bundle_pb2.TRANSCRIPT_STATUS_SUPPRESSED,
}
_STATUS_COLUMNS = frozenset({'accession_version', 'status'})
_ACCESSION = re.compile(r'^([A-Z]{2}_\d+)\.(\d+)$')


@dataclasses.dataclass(frozen=True)
class Release:
    """One RefSeq annotation release on one assembly, as local files."""

    assembly: bundle_pb2.Assembly
    release: str  # the publisher's identifier, "RS_2024_08"
    annotation: pathlib.Path  # GFF3, gzipped
    alignments: tuple[pathlib.Path, ...]  # known and model alignment BAMs, each with its .bai beside it
    hgnc: pathlib.Path  # HGNC complete set, TSV
    mane: pathlib.Path  # MANE summary, gzipped TSV
    records: pathlib.Path  # the transcripts' GenBank records, gzipped: each one's CDS in its own coordinates
    # An annotation release's sequence sets, RNA and protein FASTA, gzipped; the historical set has none,
    # and its sequences and proteins are read from the records.
    transcripts: pathlib.Path | None = None
    proteins: pathlib.Path | None = None
    # The historical set's table of each version's status in Entrez (`weaver-data-build status`);
    # an annotation release names only current versions and has none.
    status: pathlib.Path | None = None

    def __post_init__(self) -> None:
        """An annotation release has both sequence sets and no status table; the historical set the reverse.

        Raises:
            ValueError: If the inputs are any other combination.
        """
        sets = (self.transcripts is None, self.proteins is None)
        if sets == (False, False) and self.status is None:
            return
        if sets == (True, True) and self.status is not None:
            return
        raise ValueError(
            'a release has transcript and protein sequence sets and no status table; '
            'the historical set has a status table and no sequence sets'
        )

    @property
    def historical(self) -> bool:
        return self.status is not None

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [('annotation', self.annotation)]
        if self.transcripts is not None and self.proteins is not None:
            roles += [('transcripts', self.transcripts), ('proteins', self.proteins)]
        roles += [(f'alignments/{i}', path) for i, path in enumerate(self.alignments)]
        roles += [('hgnc', self.hgnc), ('mane', self.mane), ('records', self.records)]
        if self.status is not None:
            roles.append(('status', self.status))
        return roles


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
        """Whether this placement states both bounds of the CDS, to check the record's: its outer ends are closed."""
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


def _split(accession: str) -> tuple[str, int]:
    match = _ACCESSION.match(accession)
    if match is None:
        raise build.BuildError(f'{accession!r} is not a versioned accession')
    return match.group(1), int(match.group(2))


def _dbxrefs(value: str) -> dict[str, str]:
    """`GeneID:672,HGNC:HGNC:1100,MIM:113705` as a map from database to id (`HGNC` -> `HGNC:1100`)."""
    out: dict[str, str] = {}
    for ref in common.values(value):
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
            synonyms=common.values(attrs.get('gene_synonym', '')),
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
                tags=[_TAGS[t] for t in common.values(attrs.get('tag', '')) if t in _TAGS],
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
            kind, span, attrs = fields[2], (int(fields[3]), int(fields[4])), common.attributes(fields[8])
            if kind in ('gene', 'pseudogene'):
                annotation.add_gene(attrs)
            elif attrs.get('Parent', '').startswith('gene-') and kind != 'exon':
                annotation.add_rna(fields[0], kind, span, attrs)
            elif kind == 'CDS':
                annotation.add_cds(span, attrs)
    return annotation


def _read_alignments(
    paths: Iterable[pathlib.Path], transcripts: dict[str, _Transcript]
) -> dict[tuple[str, str], common.Placement]:
    """Each transcript's alignment on each sequence the annotation places it on, first seen wins."""
    return common.read_alignments(
        paths, lambda name, chromosome: name in transcripts and chromosome in transcripts[name].placements
    )


def _check_placement(versioned: str, record: _Transcript, placement: common.Placement, length: int) -> None:
    """An alignment has to tile the record's sequence and lie inside the annotation's span for it."""
    placement.check_tiles(versioned, length)
    exons = placement.exons
    low, high = record.placements[placement.chromosome].span
    if min(e.genome_start for e in exons) + 1 < low or max(e.genome_end_inclusive for e in exons) + 1 > high:
        raise build.BuildError(
            f'{versioned}: the alignment on {placement.chromosome} lies outside the annotation span {low}..{high}'
        )


def _cds(versioned: str, annotated: _AnnotatedPlacement, placement: common.Placement) -> bundle_pb2.Cds:
    """The CDS as transcript indices, projected through this placement from the bounds stated on its sequence."""
    low = min(c.start for c in annotated.cds) - 1
    high = max(c.end for c in annotated.cds) - 1
    first, last = (high, low) if placement.minus else (low, high)
    return bundle_pb2.Cds(
        start_index=placement.transcript_index(first, what=f'{versioned} CDS start on {placement.chromosome}'),
        end_index_inclusive=placement.transcript_index(last, what=f'{versioned} CDS end on {placement.chromosome}'),
    )


@dataclasses.dataclass(frozen=True)
class _Loaded:
    """A release's inputs, read."""

    assembly: bundle_pb2.Assembly
    annotation: _Annotation
    sequences: dict[str, bytes]  # transcript records by versioned accession
    proteins: dict[str, bytes]
    placements: dict[tuple[str, str], common.Placement]  # by (versioned accession, sequence)
    hgnc: dict[str, list[dict[str, str]]]  # by NCBI GeneID
    mane: dict[str, tuple[str, str]]
    records: dict[str, flatfile.Record]  # by versioned accession
    status: dict[str, bundle_pb2.TranscriptStatus] | None  # by versioned accession; None for an annotation release
    left_out: dict[str, flatfile.UninterpretableError]  # left out for a CDS statement the reader refuses
    other_protein: list[str]  # transcripts whose annotation names another protein than their record does

    def protein_of(self, record: _Transcript) -> tuple[str, int] | None:
        """The transcript's protein: the one its GenBank record names, else the one the annotation's CDS rows name.

        The record is the publisher's statement about the transcript; an annotation generated from alignments
        can name another version's protein (`docs/design/retired-versions.md`).
        """
        stated = self.records.get(f'{record.accession}.{record.version}')
        if stated is None or stated.protein_id is None:
            return record.protein
        return _split(stated.protein_id)

    def status_of(self, versioned: str) -> bundle_pb2.TranscriptStatus:
        """A transcript's status: current in an annotation release, as the table states in the historical set."""
        if self.status is None:
            return bundle_pb2.TRANSCRIPT_STATUS_CURRENT
        return self.status[versioned]


def _read_status(path: pathlib.Path) -> dict[str, bundle_pb2.TranscriptStatus]:
    """The status table: each version's status in Entrez, as `weaver-data-build status` writes it.

    Raises:
        build.BuildError: If the table lacks a column, names a version twice, or states a status word
            that is not one of Entrez's.
    """
    out: dict[str, bundle_pb2.TranscriptStatus] = {}
    with path.open(encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        if missing := _STATUS_COLUMNS - set(rows.fieldnames or ()):
            raise build.BuildError(f'{path}: columns {sorted(missing)} are not in the status table')
        for row in rows:
            versioned, word = row['accession_version'], row['status']
            if versioned in out:
                raise build.BuildError(f'{path}: {versioned} has two status rows')
            status = _STATUSES.get(word)
            if status is None:
                raise build.BuildError(f'{path}: {versioned} has status {word!r}, not one of {sorted(_STATUSES)}')
            out[versioned] = status
    return out


def _report(release: Release, loaded: _Loaded) -> None:
    annotation = loaded.annotation
    aligned: dict[str, list[str]] = collections.defaultdict(list)
    for versioned, chromosome in loaded.placements:
        aligned[versioned].append(chromosome)
    unplaced = collections.Counter(r.biotype for t, r in annotation.transcripts.items() if t not in aligned)
    by_kind = collections.Counter(common.sequence_kind(chromosome) for _, chromosome in loaded.placements)
    open_cds = collections.Counter(
        kind
        for t, record in loaded.records.items()
        if t in annotation.transcripts
        and record.cds is not None
        and (kind := common.open_ends(record.cds.start_open, record.cds.end_open)) is not None
    )
    left_out = ''
    if loaded.left_out:
        left_out = f'records left out for a CDS statement the reader refuses: {len(loaded.left_out)}; '
    other_protein = ''
    if loaded.other_protein:
        other_protein = (
            f"transcripts whose annotation names another protein than their record, the record's bundled: "
            f'{len(loaded.other_protein)} (e.g. {", ".join(loaded.other_protein[:3])}); '
        )
    by_status = ''
    if loaded.status is not None:
        statuses = collections.Counter(
            bundle_pb2.TranscriptStatus.Name(loaded.status[t]).removeprefix('TRANSCRIPT_STATUS_').lower()
            for t in annotation.transcripts
            if t in loaded.status
        )
        by_status = f'by status: {common.counted(statuses)}; '
    print(
        f'{release.release}: {len(annotation.transcripts)} transcripts, {by_status}{left_out}{other_protein}'
        f'placements aligned: {common.counted(by_kind)}; bundled without a placement: {common.counted(unplaced)}; '
        f'coding transcripts whose record leaves a CDS end open: {common.counted(open_cds)}; '
        f'RNA features left out for naming no transcript accession: {common.counted(annotation.without_accession)}',
        file=sys.stderr,
    )


def _sequence_sets(release: Release, records: dict[str, flatfile.Record]) -> tuple[dict[str, bytes], dict[str, bytes]]:
    """The transcript and protein sequences: the release's FASTA sets, or what the records carry.

    Raises:
        build.BuildError: If two records translate one protein version to different residues.
    """
    if release.transcripts is not None and release.proteins is not None:
        return common.read_fasta(release.transcripts), common.read_fasta(release.proteins)
    sequences = {v: r.sequence for v, r in records.items() if r.sequence is not None}
    proteins: dict[str, bytes] = {}
    for versioned, record in records.items():
        if record.protein is None:
            continue
        protein_id, residues = record.protein
        if proteins.setdefault(protein_id, residues) != residues:
            raise build.BuildError(f'{versioned}: translates {protein_id} to other residues than an earlier record')
    return sequences, proteins


def _other_protein(annotation: _Annotation, records: dict[str, flatfile.Record]) -> list[str]:
    """Transcripts whose annotation names another protein than their record does, for the report."""
    out = []
    for versioned, record in annotation.transcripts.items():
        stated = records.get(versioned)
        if record.protein is None or stated is None or stated.protein_id is None:
            continue
        annotated = f'{record.protein[0]}.{record.protein[1]}'
        if annotated != stated.protein_id:
            out.append(f'{versioned} ({annotated} for {stated.protein_id})')
    return out


def _load(release: Release) -> _Loaded:
    annotation = _read_annotation(release.annotation)
    left_out: dict[str, flatfile.UninterpretableError] = {}
    records = flatfile.read_records(
        release.records, sequences=release.transcripts is None, left_out=left_out if release.historical else None
    )
    for versioned in left_out:
        annotation.transcripts.pop(versioned, None)
    sequences, proteins = _sequence_sets(release, records)
    loaded = _Loaded(
        assembly=release.assembly,
        annotation=annotation,
        sequences=sequences,
        proteins=proteins,
        placements=_read_alignments(release.alignments, annotation.transcripts),
        hgnc=common.read_hgnc(release.hgnc, 'entrez_id'),
        mane={row['RefSeq_nuc']: (row['Ensembl_nuc'], row['MANE_status']) for row in common.read_mane(release.mane)},
        records=records,
        status=_read_status(release.status) if release.status is not None else None,
        left_out=left_out,
        other_protein=_other_protein(annotation, records),
    )
    _report(release, loaded)
    if left_out:
        shown = ', '.join(f'{v}: {reason}' for v, reason in list(left_out.items())[:5])
        print(f'{release.release}: records left out, e.g. {shown}', file=sys.stderr)
    return loaded


def _check_coverage(loaded: _Loaded) -> None:
    """Every coding transcript has an alignment and a GenBank CDS; every transcript a sequence and a status."""
    transcripts = loaded.annotation.transcripts
    stated_cds = {t for t, r in loaded.records.items() if r.cds is not None}
    coding = {t for t, r in transcripts.items() if r.coding or t in stated_cds}
    unaligned = sorted(t for t in coding if not any((t, c) in loaded.placements for c in transcripts[t].placements))
    if unaligned:
        raise build.BuildError(
            f'{len(unaligned)} coding transcripts the annotation places have no alignment, '
            f'e.g. {unaligned[:5]}; a bundle would have to invent one'
        )
    if unstated := sorted(t for t, r in transcripts.items() if r.coding and t not in stated_cds):
        raise build.BuildError(
            f'{len(unstated)} transcripts the annotation gives a CDS state none in their record, e.g. {unstated[:5]}'
        )
    if missing := sorted(t for t in transcripts if t not in loaded.sequences):
        raise build.BuildError(f'{len(missing)} transcripts have no sequence in the transcript set, e.g. {missing[:5]}')
    if loaded.status is not None and (unknown := sorted(t for t in transcripts if t not in loaded.status)):
        raise build.BuildError(f'{len(unknown)} transcripts have no row in the status table, e.g. {unknown[:5]}')


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
    bundle: bundle_pb2.GeneBundle,
    digests: dict[str, bundle_pb2.Sequence],
    record: _Transcript,
    protein: tuple[str, int] | None,
    loaded: _Loaded,
) -> None:
    if protein is None:
        return
    accession, version = protein
    residues = _protein_residues(f'{record.accession}.{record.version}', protein, loaded)
    bundle.proteins.add(
        accession=accession,
        version=version,
        sequence_digest=common.sequence(digests, residues, bundle_pb2.ALPHABET_PROTEIN),
        transcript_accession=record.accession,
        transcript_version=record.version,
    )


def _cross_check(
    versioned: str, record: _Transcript, aligned: list[common.Placement], stated: flatfile.RecordCds
) -> None:
    """Warn where a placement that states the whole CDS projects it elsewhere than the record states it.

    A placement whose annotation leaves an outer end of the CDS open — the start codon in bases that
    sequence lacks — states no bound to check there, and is passed over. A record stating a CDS the
    annotation gives none of is a warning too: the transcript is bundled with the record's CDS and the
    protein its record names.
    """
    if not record.coding:
        warnings.warn(
            f'{versioned}: CDS {stated.start}..{stated.end} taken from its record; the annotation gives it no CDS',
            build.BuildWarning,
            stacklevel=1,
        )
        return
    for placement in aligned:
        annotated = record.placements[placement.chromosome]
        if not annotated.determines_cds:
            continue
        try:
            projected = _cds(versioned, annotated, placement)
        except build.BuildError as error:
            disagreement = f'the annotation does not project: {error}'
        else:
            if (projected.start_index, projected.end_index_inclusive) == (stated.start, stated.end):
                continue
            disagreement = (
                f'the annotation projects {projected.start_index}..{projected.end_index_inclusive} '
                f'through the alignment on {placement.chromosome}'
            )
        warnings.warn(
            f'{versioned}: CDS {stated.start}..{stated.end} taken from its record; {disagreement}',
            build.BuildWarning,
            stacklevel=1,
        )


def _add_transcript(
    bundle: bundle_pb2.GeneBundle, digests: dict[str, bundle_pb2.Sequence], record: _Transcript, loaded: _Loaded
) -> None:
    """One transcript into its gene's bundle: its record, status, MANE status, alignments, CDS and protein.

    A retired version carries no tags and no MANE partner: those name what the publisher recommends
    today, and a tag the annotation gave the version while it was current no longer holds.
    """
    versioned = f'{record.accession}.{record.version}'
    residues = loaded.sequences[versioned]
    status = loaded.status_of(versioned)
    current = status == bundle_pb2.TRANSCRIPT_STATUS_CURRENT
    transcript = bundle.transcripts.add(
        accession=record.accession,
        version=record.version,
        publisher=bundle_pb2.PUBLISHER_REFSEQ,
        status=status,
        biotype=record.biotype,
        tags=record.tags if current else [],
        sequence_digest=common.sequence(digests, residues, bundle_pb2.ALPHABET_NUCLEOTIDE),
    )
    partner = loaded.mane.get(versioned) if current else None
    if partner is not None:
        transcript.mane_partner = partner[0]
        tag = _TAGS.get(partner[1])
        if tag is not None and tag not in transcript.tags:
            transcript.tags.append(tag)
    aligned: list[common.Placement] = []
    for chromosome in sorted(record.placements, key=common.placement_order):
        placement = loaded.placements.get((versioned, chromosome))
        if placement is None:
            continue
        _check_placement(versioned, record, placement, len(residues))
        transcript.alignments.append(placement.message(loaded.assembly, bundle_pb2.ALIGNMENT_SOURCE_NCBI_BAM))
        aligned.append(placement)
    stated = loaded.records[versioned].cds if versioned in loaded.records else None
    if stated is not None:
        transcript.cds.start_index, transcript.cds.end_index_inclusive = stated.start, stated.end
        transcript.cds.start_open, transcript.cds.end_open = stated.start_open, stated.end_open
        common.check_cds(versioned, transcript.cds, len(residues))
        _cross_check(versioned, record, aligned, stated)
    protein = loaded.protein_of(record)
    if protein is not None:
        transcript.protein_accession, transcript.protein_version = protein
    _add_protein(bundle, digests, record, protein, loaded)


def _gene_bundle(gene: _Gene, members: Iterable[_Transcript], loaded: _Loaded) -> bundle_pb2.GeneBundle:
    hgnc, others = common.hgnc_row(loaded.hgnc.get(gene.dbxrefs.get('GeneID', '')), gene.symbol)
    message = common.gene_message(
        feature_id=gene.gene_id,
        symbol=gene.symbol,
        hgnc_id=gene.dbxrefs.get('HGNC', ''),
        ncbi_gene_id=gene.dbxrefs.get('GeneID', ''),
        ensembl_gene_id=gene.dbxrefs.get('Ensembl', ''),
        description=gene.description,
        synonyms=[*gene.synonyms, *others],
        hgnc=hgnc,
    )
    bundle = bundle_pb2.GeneBundle(gene=message)
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
