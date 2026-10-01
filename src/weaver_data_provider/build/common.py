"""What the RefSeq and Ensembl builders share: the publishers' file formats and the shape of a placement.

Each builder reads its publisher's annotation in its own module; the sequence sets, HGNC, the MANE
summary, NCBI's alignment BAMs and the per-exon cigar a placement is stored as are common.
"""

from __future__ import annotations

import collections
import csv
import dataclasses
import gzip
import pathlib
import re
import urllib.parse
from collections.abc import Callable, Iterable

import pysam.libcalignedsegment
import pysam.libcalignmentfile

from weaver_data_provider import build, refget
from weaver_data_provider.v1 import bundle_pb2

MANE_STATUS_TAGS = {'MANE Select': bundle_pb2.TAG_MANE_SELECT, 'MANE Plus Clinical': bundle_pb2.TAG_MANE_PLUS_CLINICAL}
HGNC_COLUMNS = frozenset({'hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id'})
MANE_COLUMNS = frozenset({'RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'})
_CIGAR_OPS = {7: '=', 8: 'X', 1: 'I', 2: 'D', 4: 'I'}  # soft-clipped transcript bases have no genome counterpart
_COMPLEMENT = bytes.maketrans(b'ACGTN', b'TGCAN')
_AMBIGUOUS = bytes.maketrans(b'RYKMSWBDHV', b'NNNNNNNNNN')


def reverse_complement(residues: bytes) -> bytes:
    return residues.translate(_COMPLEMENT)[::-1]


def unambiguous(residues: bytes) -> bytes:
    """The residues with every IUPAC ambiguity code read as N, which is how a published transcript writes one."""
    return residues.translate(_AMBIGUOUS)


def sequence_kind(accession: str) -> str:
    """`chromosome`, `scaffold` (an alternate locus or an unplaced one) or `patch`, from the RefSeq prefix.

    Raises:
        build.BuildError: If the accession is not one of an NCBI assembly's sequences.
    """
    kind = build.SEQUENCE_KINDS.get(accession[:3])
    if kind is None:
        raise build.BuildError(f'{accession}: not a sequence of an NCBI assembly (NC_, NT_ or NW_)')
    return kind


def placement_order(chromosome: str) -> tuple[int, str]:
    """The order a transcript's alignments are kept in, which the proto promises: by sequence kind, then accession."""
    return (list(build.SEQUENCE_KINDS.values()).index(sequence_kind(chromosome)), chromosome)


def attributes(field: str) -> dict[str, str]:
    """A GFF3 attributes column as a map; a repeated key keeps its last value."""
    return dict(kv.split('=', 1) for kv in field.split(';') if '=' in kv)


def values(attribute: str) -> list[str]:
    """A GFF3 attribute's comma-separated values, each percent-decoded after the split."""
    return [urllib.parse.unquote(v) for v in attribute.split(',') if v]


def read_fasta(path: pathlib.Path) -> dict[str, bytes]:
    """Every record of a gzipped FASTA, upper-cased, keyed by the first word of its header."""
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


def _columns(reader: csv.DictReader[str], needed: frozenset[str], *, named: str) -> None:
    if missing := needed - set(reader.fieldnames or ()):
        raise build.BuildError(f'{named}: columns {sorted(missing)} are not in the table; its schema has changed')


def read_hgnc(path: pathlib.Path, key: str) -> dict[str, list[dict[str, str]]]:
    """HGNC's complete set grouped by an id column, `entrez_id` or `ensembl_gene_id`; rows without one are dropped.

    A key can hold more than one row: HGNC maps two of its genes to one Ensembl gene in three cases.

    Raises:
        build.BuildError: If the table lacks a column the builders read.
    """
    out: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    with path.open(encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        _columns(rows, HGNC_COLUMNS, named='hgnc')
        for row in rows:
            if row[key]:
                out[row[key]].append(row)
    return dict(out)


def hgnc_row(rows: list[dict[str, str]] | None, symbol: str) -> tuple[dict[str, str] | None, list[str]]:
    """The HGNC row naming a gene, and every symbol of the other rows its key holds.

    The row is the one whose symbol is the annotation's, or else the lowest HGNC id. The other rows'
    symbols, previous and alias ones too, stay findable as the gene's aliases.
    """
    if not rows:
        return None, []
    named = [r for r in rows if r['symbol'] == symbol]
    row = named[0] if named else min(rows, key=lambda r: int(r['hgnc_id'].removeprefix('HGNC:')))
    others = [
        s
        for other in rows
        if other is not row
        for s in (other['symbol'], *other['prev_symbol'].split('|'), *other['alias_symbol'].split('|'))
        if s
    ]
    return row, others


def read_mane(path: pathlib.Path) -> list[dict[str, str]]:
    """The MANE summary's rows: each names a RefSeq and an Ensembl transcript, versioned, and the pair's status.

    Raises:
        build.BuildError: If the table lacks a column the builders read.
    """
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        _columns(rows, MANE_COLUMNS, named='mane')
        return list(rows)


def preferred(hgnc: dict[str, str] | None, column: str, annotation: str) -> str:
    """HGNC's value where it has the gene and states one; the annotation's otherwise."""
    if hgnc is not None and hgnc[column]:
        return hgnc[column]
    return annotation


def gene_message(
    *,
    feature_id: str,
    symbol: str,
    hgnc_id: str,
    ncbi_gene_id: str,
    ensembl_gene_id: str,
    description: str,
    synonyms: Iterable[str],
    hgnc: dict[str, str] | None,
) -> bundle_pb2.Gene:
    """A gene's identity across nomenclatures: HGNC's where it knows the gene, the annotation's otherwise.

    Raises:
        build.BuildError: If neither names a symbol.
    """
    message = bundle_pb2.Gene(
        symbol=preferred(hgnc, 'symbol', symbol),
        hgnc_id=preferred(hgnc, 'hgnc_id', hgnc_id),
        ncbi_gene_id=preferred(hgnc, 'entrez_id', ncbi_gene_id),
        ensembl_gene_id=preferred(hgnc, 'ensembl_gene_id', ensembl_gene_id),
        name=preferred(hgnc, 'name', description),
    )
    if hgnc is not None:
        message.previous_symbols.extend(s for s in hgnc['prev_symbol'].split('|') if s)
        message.alias_symbols.extend(s for s in hgnc['alias_symbol'].split('|') if s)
    for synonym in synonyms:
        if synonym not in (message.symbol, *message.previous_symbols, *message.alias_symbols):
            message.alias_symbols.append(synonym)
    if not message.symbol:
        raise build.BuildError(f'{feature_id}: no symbol in the annotation or HGNC')
    return message


def sequence(digests: dict[str, bundle_pb2.Sequence], residues: bytes, alphabet: bundle_pb2.Alphabet) -> str:
    """The digest of `residues`, adding the sequence to the bundle's set on first sight."""
    digest = refget.digest(residues)
    if digest not in digests:
        digests[digest] = bundle_pb2.Sequence(digest=digest, alphabet=alphabet, residues=residues)
    return digest


def counted(counts: collections.Counter[str]) -> str:
    """`12 (tRNA 8, rRNA 4)`, most common first; `0` when empty."""
    if not counts:
        return '0'
    return f'{counts.total()} (' + ', '.join(f'{kind} {n}' for kind, n in counts.most_common()) + ')'


# ---- alignments ------------------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Placement:
    """One transcript's alignment on one sequence, as the exons a bundle stores."""

    chromosome: str  # the assembly sequence: a chromosome, scaffold or patch
    minus: bool
    exons: list[bundle_pb2.Exon]
    query: bytes | None = None  # the aligned record as the BAM carries it, in transcript orientation

    def message(self, assembly: bundle_pb2.Assembly, source: bundle_pb2.AlignmentSource) -> bundle_pb2.Alignment:
        return bundle_pb2.Alignment(
            assembly=assembly,
            chromosome=self.chromosome,
            strand=bundle_pb2.STRAND_MINUS if self.minus else bundle_pb2.STRAND_PLUS,
            source=source,
            exons=self.exons,
        )

    def check_tiles(self, versioned: str, length: int) -> None:
        """The alignment has to cover the record's bases from its first to its last.

        Raises:
            build.BuildError: If it does not.
        """
        first, last = self.exons[0].transcript_start, self.exons[-1].transcript_end
        if first != 0 or last != length:
            raise build.BuildError(
                f'{versioned}: the alignment on {self.chromosome} covers transcript bases {first}..{last} of a '
                f'{length}-base record'
            )

    def transcript_index(self, genome_position: int, *, what: str) -> int:
        """The 0-based transcript index a 0-based genome position projects to, through the exon cigars.

        Raises:
            build.BuildError: If the position is in no exon, or is a genome-only base (a `D` run), which no
                transcript index names.
        """
        for exon in self.exons:
            if not exon.genome_start <= genome_position <= exon.genome_end_inclusive:
                continue
            offset = exon.genome_end_inclusive - genome_position if self.minus else genome_position - exon.genome_start
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


def read_alignments(
    paths: Iterable[pathlib.Path], wanted: Callable[[str, str], bool]
) -> dict[tuple[str, str], Placement]:
    """Each wanted transcript's primary alignment on each wanted sequence, first seen wins.

    Args:
        paths: NCBI's alignment BAMs.
        wanted: Whether an alignment of this transcript on this sequence is one to keep.
    """
    placements: dict[tuple[str, str], Placement] = {}
    for path in paths:
        with pysam.libcalignmentfile.AlignmentFile(str(path)) as bam:
            for read in bam:
                name = read.query_name
                if name is None or read.is_secondary or read.is_supplementary:
                    continue
                chromosome = read.reference_name
                if chromosome is None or (name, chromosome) in placements or not wanted(name, chromosome):
                    continue
                query = read.query_sequence.upper().encode('ascii') if read.query_sequence else None
                if query is not None and read.is_reverse:
                    query = reverse_complement(query)
                placements[name, chromosome] = Placement(chromosome, read.is_reverse, _exons_of(read), query)
    return placements
