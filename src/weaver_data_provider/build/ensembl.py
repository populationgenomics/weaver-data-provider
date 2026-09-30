"""Cut one Ensembl release into gene bundles.

The inputs are the files Ensembl publishes for a release — the primary-chromosome GFF3 annotation, the
cDNA and ncRNA sequence sets, and the peptide set — plus HGNC's complete set, for gene identity, the
MANE summary, for each MANE transcript's RefSeq partner, NCBI's assembly report, to name Ensembl's
chromosomes (`1`, `X`, `MT`) by the accessions the genome uses, and the genome itself, built for the
same assembly.

Ensembl publishes no alignments: a transcript is its annotated exons on the genome, and its record is
the genome spliced at them. The builder checks that record against the genome, exon by exon, and writes
each exon's cigar from the comparison — `=`, or `X` where the record's base is not the genome's; a
record whose length differs from its exons fails the build. The CDS is the annotation's, projected
through those exons.

What is and is not bundled. Every transcript the annotation places on a chromosome that has a record in
the sequence sets. A transcript with no record — Ensembl's TEC and artifact models — is left out, and a
coding one fails the build. Ensembl places some transcripts only on patches and alternate loci,
MANE-tagged ones among them; their records have no transcript in the primary-chromosome annotation and
are counted, with the MANE ones named, in the build's report. The pseudoautosomal regions are two
genes, one on X and one on Y, each with transcripts of its own. A gene the annotation names nothing is
named by its stable id.
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

from weaver_data_provider import build
from weaver_data_provider import genome as genome_mod
from weaver_data_provider.build import _common
from weaver_data_provider.v1 import bundle_pb2

_TAGS = {
    'MANE_Select': bundle_pb2.TAG_MANE_SELECT,
    'MANE_Plus_Clinical': bundle_pb2.TAG_MANE_PLUS_CLINICAL,
    'Ensembl_canonical': bundle_pb2.TAG_ENSEMBL_CANONICAL,
}
_SOURCE = re.compile(r'\s*\[Source:[^]]*\]$')
_COMPLEMENT = str.maketrans('ACGTRYKMBVDHN', 'TGCAYRMKVBHDN')
_REPORT_COLUMNS = frozenset({'Sequence-Name', 'RefSeq-Accn'})


@dataclasses.dataclass(frozen=True)
class Release:
    """One Ensembl release on one assembly, as local files."""

    assembly: bundle_pb2.Assembly
    release: str  # the publisher's release number, "116"
    annotation: pathlib.Path  # the primary-chromosome GFF3, gzipped
    transcripts: tuple[pathlib.Path, ...]  # cDNA and ncRNA FASTA, gzipped
    proteins: pathlib.Path  # peptide FASTA, gzipped
    hgnc: pathlib.Path  # HGNC complete set, TSV
    mane: pathlib.Path  # MANE summary, gzipped TSV
    assembly_report: pathlib.Path  # NCBI's assembly report, naming each sequence's RefSeq accession
    genome: pathlib.Path  # the genome built for this assembly, whose bases each record is checked against

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [('annotation', self.annotation)]
        roles += [(f'transcripts/{i}', path) for i, path in enumerate(self.transcripts)]
        roles += [('proteins', self.proteins), ('hgnc', self.hgnc), ('mane', self.mane)]
        # the genome by its catalogue, which records the blocks file's digest
        return [*roles, ('assembly_report', self.assembly_report), ('genome', self.genome / genome_mod.CATALOGUE)]


@dataclasses.dataclass
class _Gene:
    gene_id: str  # the stable id, "ENSG00000139618"
    symbol: str  # the annotation's Name; empty where it gives none
    description: str


@dataclasses.dataclass
class _Transcript:
    accession: str
    version: int
    gene_id: str
    biotype: str  # the annotation's feature type: mRNA, lnc_RNA, ...
    tags: list[int]
    seqid: str  # Ensembl's name for the chromosome, "13"
    minus: bool
    exons: list[tuple[int, int]] = dataclasses.field(default_factory=list)  # 1-based closed genome spans
    cds_spans: list[tuple[int, int]] = dataclasses.field(default_factory=list)
    protein: tuple[str, int] | None = None

    @property
    def versioned(self) -> str:
        return f'{self.accession}.{self.version}'


def _attributes(field: str) -> dict[str, str]:
    return {k: urllib.parse.unquote(v) for k, _, v in (kv.partition('=') for kv in field.split(';') if kv)}


def _add_feature(
    genes: dict[str, _Gene], transcripts: dict[str, _Transcript], fields: list[str], attrs: dict[str, str]
) -> None:
    ident = attrs.get('ID', '')
    if ident.startswith('gene:'):
        genes[ident] = _Gene(
            gene_id=attrs['gene_id'],
            symbol=attrs.get('Name', ''),
            description=_SOURCE.sub('', attrs.get('description', '')),
        )
    elif ident.startswith('transcript:'):
        transcripts[ident] = _Transcript(
            accession=attrs['transcript_id'],
            version=int(attrs['version']),
            gene_id=attrs['Parent'],
            biotype=fields[2],
            tags=[_TAGS[t] for t in attrs.get('tag', '').split(',') if t in _TAGS],
            seqid=fields[0],
            minus=fields[6] == '-',
        )
    elif fields[2] == 'exon':
        transcripts[attrs['Parent']].exons.append((int(fields[3]), int(fields[4])))
    elif fields[2] == 'CDS':
        record = transcripts[attrs['Parent']]
        record.cds_spans.append((int(fields[3]), int(fields[4])))
        record.protein = (attrs['protein_id'], int(attrs['version']))


def _read_annotation(path: pathlib.Path) -> tuple[dict[str, _Gene], dict[str, _Transcript]]:
    """Genes and their transcripts, with each transcript's exons, CDS spans and protein.

    Features follow the ones they belong to, as Ensembl writes them: a gene before its transcripts, a
    transcript before its exons and CDS.

    Raises:
        build.BuildError: If a feature line has other than nine fields, or names a transcript not yet
            declared.
    """
    genes: dict[str, _Gene] = {}
    transcripts: dict[str, _Transcript] = {}
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise build.BuildError(f'{path}: a feature line with {len(fields)} fields, not 9: {line[:80]!r}')
            try:
                _add_feature(genes, transcripts, fields, _attributes(fields[8]))
            except KeyError as error:
                raise build.BuildError(f'{path}: a {fields[2]} feature without {error}: {line[:120]!r}') from error
    return genes, transcripts


def _read_report(path: pathlib.Path) -> dict[str, str]:
    """NCBI's assembly report: each sequence's name, as Ensembl gives it, -> its RefSeq accession."""
    with path.open(encoding='utf-8') as fh:
        header: list[str] | None = None
        rows: list[str] = []
        for line in fh:
            if line.startswith('#'):
                if 'Sequence-Name' in line:
                    header = line.lstrip('#').strip().split('\t')
            else:
                rows.append(line)
    if header is None:
        raise build.BuildError(f'{path}: no column header; not an NCBI assembly report')
    table = csv.DictReader(rows, fieldnames=header, delimiter='\t')
    _common.columns(table, _REPORT_COLUMNS, named='assembly report')
    return {row['Sequence-Name']: row['RefSeq-Accn'] for row in table if row['RefSeq-Accn'] not in ('', 'na')}


def _cigar(record: str, genome: str) -> str:
    """`=` and `X` runs comparing two sequences of one length, base by base."""
    ops: list[list[int | str]] = []
    for a, b in zip(record, genome, strict=True):
        op = '=' if a == b else 'X'
        if ops and ops[-1][1] == op:
            ops[-1][0] = int(ops[-1][0]) + 1
        else:
            ops.append([1, op])
    return ''.join(f'{n}{op}' for n, op in ops)


def _placement(record: _Transcript, residues: bytes, chromosome: str, genome: genome_mod.Genome) -> _common.Placement:
    """The transcript's annotated exons on `chromosome`, each with the cigar its record and the genome give.

    Raises:
        build.BuildError: If the record's length is not its exons'.
    """
    spans = sorted(record.exons, reverse=record.minus)  # transcript order
    length = sum(e - s + 1 for s, e in spans)
    if length != len(residues):
        raise build.BuildError(
            f'{record.versioned}: a {len(residues)}-base record for {length} bases of annotated exons; '
            'with no published alignment, a record that is not its exons cannot be placed'
        )
    bases = genome.fetch_many(chromosome, [(s - 1, e) for s, e in spans])
    text = residues.decode('ascii')
    exons, t = [], 0
    for (s, e), genomic in zip(spans, bases, strict=True):
        oriented = genomic.translate(_COMPLEMENT)[::-1] if record.minus else genomic
        n = e - s + 1
        exons.append(
            bundle_pb2.Exon(
                transcript_start=t,
                transcript_end=t + n,
                genome_start=s - 1,
                genome_end_inclusive=e - 1,
                cigar=_cigar(text[t : t + n], oriented),
            )
        )
        t += n
    return _common.Placement(chromosome, record.minus, exons)


@dataclasses.dataclass(frozen=True)
class _Loaded:
    """A release's inputs, read."""

    assembly: bundle_pb2.Assembly
    genes: dict[str, _Gene]
    transcripts: dict[str, _Transcript]  # the ones with a record
    sequences: dict[str, bytes]  # records by versioned accession
    proteins: dict[str, bytes]
    hgnc: dict[str, dict[str, str]]  # by Ensembl gene id
    mane: dict[str, tuple[str, str]]  # by Ensembl accession.version
    chromosomes: dict[str, str]  # Ensembl's sequence names -> the genome's
    genome: genome_mod.Genome


def _chromosomes(
    transcripts: Iterable[_Transcript], report: dict[str, str], genome: genome_mod.Genome
) -> dict[str, str]:
    seqids = {t.seqid for t in transcripts}
    if unnamed := sorted(s for s in seqids if report.get(s, '') not in genome):
        raise build.BuildError(
            f'sequences {unnamed[:5]} have no accession in the assembly report that the genome holds'
        )
    return {s: report[s] for s in seqids}


def _load(release: Release) -> _Loaded:
    genes, annotated = _read_annotation(release.annotation)
    sequences: dict[str, bytes] = {}
    for path in release.transcripts:
        sequences |= _common.read_fasta(path)
    genome = genome_mod.Genome(str(release.genome))
    if genome.assembly != release.assembly:
        name = bundle_pb2.Assembly.Name
        raise build.BuildError(f'the genome is on {name(genome.assembly)}, the release on {name(release.assembly)}')
    without_record = [t for t in annotated.values() if t.versioned not in sequences]
    if coding := sorted(t.versioned for t in without_record if t.cds_spans):
        raise build.BuildError(
            f'{len(coding)} coding transcripts have no record in the sequence sets, e.g. {coding[:5]}'
        )
    transcripts = {t.versioned: t for t in annotated.values() if t.versioned in sequences}
    mane = _common.read_mane(release.mane, key='Ensembl_nuc', partner='RefSeq_nuc')
    elsewhere = sorted(v for v in sequences if v not in transcripts)
    elsewhere_mane = [v for v in elsewhere if v in mane]
    print(
        f'Ensembl {release.release}: {len(transcripts)} transcripts on chromosomes; left out for having no record: '
        f'{_common.counted(collections.Counter(t.biotype for t in without_record))}; records for transcripts '
        f'placed off these chromosomes, left out: {len(elsewhere)}, of them MANE {len(elsewhere_mane)}'
        + (f': {", ".join(elsewhere_mane)}' if elsewhere_mane else ''),
        file=sys.stderr,
    )
    return _Loaded(
        assembly=release.assembly,
        genes=genes,
        transcripts=transcripts,
        sequences=sequences,
        proteins=_common.read_fasta(release.proteins),
        hgnc=_common.read_hgnc(release.hgnc, key='ensembl_gene_id'),
        mane=mane,
        chromosomes=_chromosomes(transcripts.values(), _read_report(release.assembly_report), genome),
        genome=genome,
    )


def _gene_message(gene: _Gene, hgnc: dict[str, str] | None) -> bundle_pb2.Gene:
    return _common.gene_message(
        gene.gene_id,
        hgnc,
        symbol=_common.preferred(hgnc, 'symbol', gene.symbol or gene.gene_id),
        hgnc_id=_common.preferred(hgnc, 'hgnc_id', ''),
        ncbi_gene_id=_common.preferred(hgnc, 'entrez_id', ''),
        ensembl_gene_id=gene.gene_id,
        name=_common.preferred(hgnc, 'name', gene.description),
        synonyms=[],
    )


def _add_transcript(
    bundle: bundle_pb2.GeneBundle, digests: dict[str, bundle_pb2.Sequence], record: _Transcript, loaded: _Loaded
) -> None:
    """One transcript into its gene's bundle: its record, MANE partner, placement, CDS and protein."""
    residues = loaded.sequences[record.versioned]
    transcript = bundle.transcripts.add(
        accession=record.accession,
        version=record.version,
        publisher=bundle_pb2.PUBLISHER_ENSEMBL,
        status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
        biotype=record.biotype,
        tags=record.tags,
        sequence_digest=_common.sequence(digests, residues, bundle_pb2.ALPHABET_NUCLEOTIDE),
    )
    partner = loaded.mane.get(record.versioned)
    if partner is not None:
        transcript.mane_partner = partner[0]
        tag = _common.MANE_TAGS.get(partner[1])
        if tag is not None and tag not in transcript.tags:
            transcript.tags.append(tag)
    placement = _placement(record, residues, loaded.chromosomes[record.seqid], loaded.genome)
    transcript.alignments.append(_common.alignment(loaded.assembly, placement, bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION))
    if record.cds_spans:
        transcript.cds.CopyFrom(_common.project_cds(record.versioned, record.cds_spans, placement))
    if record.protein is not None:
        _common.add_protein(bundle, digests, transcript, record.protein, loaded.proteins)


def _gene_bundle(gene: _Gene, members: Iterable[_Transcript], loaded: _Loaded) -> bundle_pb2.GeneBundle:
    bundle = bundle_pb2.GeneBundle(gene=_gene_message(gene, loaded.hgnc.get(gene.gene_id)))
    digests: dict[str, bundle_pb2.Sequence] = {}
    for record in sorted(members, key=lambda r: (r.accession, r.version)):
        _add_transcript(bundle, digests, record, loaded)
    bundle.sequences.extend(digests[d] for d in sorted(digests))
    return bundle


def bundles(release: Release) -> Iterator[bundle_pb2.GeneBundle]:
    """Every gene's bundle in the release, in a stable order: by symbol, then by stable gene id.

    Raises:
        build.BuildError: On a coding transcript with no record, a record whose length is not its exons',
            a protein missing from the peptide set, a chromosome the assembly report or the genome does not
            name, or a genome on another assembly.
    """
    loaded = _load(release)
    by_gene: dict[str, list[_Transcript]] = collections.defaultdict(list)
    for record in loaded.transcripts.values():
        by_gene[record.gene_id].append(record)
    out: list[tuple[str, str, bundle_pb2.GeneBundle]] = []
    for gene_key, members in by_gene.items():
        gene = loaded.genes.get(gene_key)
        if gene is None:
            raise build.BuildError(f'{gene_key}: transcripts name a gene the annotation does not declare')
        bundle = _gene_bundle(gene, members, loaded)
        out.append((bundle.gene.symbol, gene.gene_id, bundle))
    for _, _, bundle in sorted(out, key=lambda item: (item[0], item[1])):
        yield bundle
