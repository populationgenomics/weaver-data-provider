"""What the RefSeq and Ensembl builders share, whoever published the release.

Reading the tables beside a release — sequence sets, HGNC, MANE — and making the messages a bundle is
built from.
"""

from __future__ import annotations

import collections
import csv
import dataclasses
import gzip
import pathlib
import re

from weaver_data_provider import build, refget
from weaver_data_provider.v1 import bundle_pb2

HGNC_COLUMNS = frozenset({'hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id'})
MANE_COLUMNS = frozenset({'RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'})
MANE_TAGS = {'MANE Select': bundle_pb2.TAG_MANE_SELECT, 'MANE Plus Clinical': bundle_pb2.TAG_MANE_PLUS_CLINICAL}


@dataclasses.dataclass(frozen=True)
class Placement:
    """A transcript on one chromosome: its exons in transcript order, each cigar in transcript orientation."""

    chromosome: str
    minus: bool
    exons: list[bundle_pb2.Exon]


def read_fasta(path: pathlib.Path) -> dict[str, bytes]:
    """Each record's residues, upper-cased, keyed by the first word of its header."""
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


def columns(reader: csv.DictReader[str], needed: frozenset[str], *, named: str) -> None:
    """Refuse a table that lacks a column this builder reads, rather than read it wrong."""
    if missing := needed - set(reader.fieldnames or ()):
        raise build.BuildError(f'{named}: columns {sorted(missing)} are not in the table; its schema has changed')


def read_hgnc(path: pathlib.Path, *, key: str) -> dict[str, list[dict[str, str]]]:
    """HGNC's complete set grouped by `key`, `entrez_id` or `ensembl_gene_id`; rows without one are left out.

    A key can hold more than one row: three Ensembl genes each carry two HGNC genes.
    """
    out: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    with path.open(encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        columns(rows, HGNC_COLUMNS, named='hgnc')
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


def read_mane(path: pathlib.Path, *, key: str, partner: str) -> dict[str, tuple[str, str]]:
    """MANE's pairs: one publisher's accession.version (`key`) -> (the other's, `partner`, and the MANE status)."""
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        rows = csv.DictReader(fh, delimiter='\t')
        columns(rows, MANE_COLUMNS, named='mane')
        return {row[key]: (row[partner], row['MANE_status']) for row in rows}


def preferred(hgnc: dict[str, str] | None, column: str, annotation: str) -> str:
    """HGNC's value where it has the gene and states one; the annotation's otherwise."""
    if hgnc is not None and hgnc[column]:
        return hgnc[column]
    return annotation


def gene_message(
    gene_id: str,
    hgnc: dict[str, str] | None,
    *,
    symbol: str,
    hgnc_id: str,
    ncbi_gene_id: str,
    ensembl_gene_id: str,
    name: str,
    synonyms: list[str],
) -> bundle_pb2.Gene:
    """A gene with HGNC's previous and alias symbols, then the annotation's synonyms not already among them.

    Raises:
        build.BuildError: If the gene has no symbol.
    """
    message = bundle_pb2.Gene(
        symbol=symbol, hgnc_id=hgnc_id, ncbi_gene_id=ncbi_gene_id, ensembl_gene_id=ensembl_gene_id, name=name
    )
    if hgnc is not None:
        message.previous_symbols.extend(s for s in hgnc['prev_symbol'].split('|') if s)
        message.alias_symbols.extend(s for s in hgnc['alias_symbol'].split('|') if s)
    for synonym in synonyms:
        if synonym not in message.alias_symbols and synonym != message.symbol:
            message.alias_symbols.append(synonym)
    if not message.symbol:
        raise build.BuildError(f'{gene_id}: no symbol in the annotation or HGNC')
    return message


def sequence(digests: dict[str, bundle_pb2.Sequence], residues: bytes, alphabet: bundle_pb2.Alphabet) -> str:
    """The digest of `residues`, adding the sequence to the bundle's set on first sight."""
    digest = refget.digest(residues)
    if digest not in digests:
        digests[digest] = bundle_pb2.Sequence(digest=digest, alphabet=alphabet, residues=residues)
    return digest


def alignment(
    assembly: bundle_pb2.Assembly, placement: Placement, source: bundle_pb2.AlignmentSource
) -> bundle_pb2.Alignment:
    return bundle_pb2.Alignment(
        assembly=assembly,
        chromosome=placement.chromosome,
        strand=bundle_pb2.STRAND_MINUS if placement.minus else bundle_pb2.STRAND_PLUS,
        source=source,
        exons=placement.exons,
    )


def transcript_index(placement: Placement, genome_position: int, *, what: str) -> int:
    """The 0-based transcript index a 0-based genome position projects to, through the exon cigars.

    Raises:
        build.BuildError: If the position is in no exon, or is a genome-only base (a `D` run), which no
            transcript index names.
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


def project_cds(versioned: str, cds_spans: list[tuple[int, int]], placement: Placement) -> bundle_pb2.Cds:
    """The annotation's CDS, 1-based closed genome spans, projected onto the transcript through a placement."""
    low = min(s for s, _ in cds_spans) - 1
    high = max(e for _, e in cds_spans) - 1
    first, last = (high, low) if placement.minus else (low, high)
    return bundle_pb2.Cds(
        start_index=transcript_index(placement, first, what=f'{versioned} CDS start'),
        end_index_inclusive=transcript_index(placement, last, what=f'{versioned} CDS end'),
    )


def add_protein(
    bundle: bundle_pb2.GeneBundle,
    digests: dict[str, bundle_pb2.Sequence],
    transcript: bundle_pb2.Transcript,
    protein: tuple[str, int],
    proteins: dict[str, bytes],
) -> None:
    """A transcript's protein into its bundle, and onto the transcript.

    Raises:
        build.BuildError: If the protein is not in the protein set.
    """
    accession, version = protein
    transcript.protein_accession, transcript.protein_version = accession, version
    residues = proteins.get(f'{accession}.{version}')
    if residues is None:
        raise build.BuildError(
            f'{transcript.accession}.{transcript.version}: its protein {accession}.{version} is not in the protein set'
        )
    bundle.proteins.add(
        accession=accession,
        version=version,
        sequence_digest=sequence(digests, residues, bundle_pb2.ALPHABET_PROTEIN),
        transcript_accession=transcript.accession,
        transcript_version=transcript.version,
    )


def counted(counts: collections.Counter[str]) -> str:
    """`12 (tRNA 8, rRNA 4)`, most common first; `0` when empty."""
    if not counts:
        return '0'
    return f'{counts.total()} (' + ', '.join(f'{kind} {n}' for kind, n in counts.most_common()) + ')'
