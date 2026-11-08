"""Cut GenBank mRNA records into gene bundles: each record's sequence, CDS and protein, as it states them, unplaced.

The inputs are GenBank flat files — whichever records the caller fetched, from the nine reference
standards to every human mRNA — plus NCBI's `gene2accession`, which names each record's gene, and HGNC's
complete set, for the gene's identity. A record states no place on any genome, so none is bundled: the
provider serves the record on its own sequence and refuses what needs a genome
(`docs/design/genbank-and-uniprot.md`).

Bundled: every mRNA record with one CDS naming its protein, whose translation is the conceptual
translation of the CDS's bases, of a human gene `gene2accession` names. A record NCBI ties to several
genes goes into each one's bundle. Left out, and counted by reason in the build's report, rather than
bundled with a CDS or protein the record does not state plainly: a record that is not an mRNA; one with
no CDS, which means the CDS was not annotated, not that the RNA is non-coding; a CDS the reader does not
interpret, a join (a programmed frameshift), one read from other than its first base, or one marked with
an exception; a translation that is not the CDS's; and a record of no human gene.
"""

from __future__ import annotations

import collections
import csv
import dataclasses
import gzip
import pathlib
import re
import sys
from collections.abc import Iterator

from weaver_data_provider import build
from weaver_data_provider.build import common, flatfile
from weaver_data_provider.v1 import bundle_pb2

_HUMAN = '9606'
_GENE2ACCESSION_COLUMNS = frozenset({'tax_id', 'GeneID', 'RNA_nucleotide_accession.version', 'Symbol'})
_VERSIONED = re.compile(r'^([A-Z]+_?\d+)\.(\d+)$')
_BASES = 'TCAG'
_AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
# The standard genetic code, NCBI's table 1; a codon with any other base translates as X.
_CODE = {
    a + b + c: _AMINO_ACIDS[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}


@dataclasses.dataclass(frozen=True)
class Release:
    """One cut of GenBank records, as local files."""

    release: str  # the caller's name for this cut, "GB_2026_11"
    records: tuple[pathlib.Path, ...]  # GenBank flat files, gzipped
    gene2accession: pathlib.Path  # NCBI's gene2accession, gzipped
    hgnc: pathlib.Path  # HGNC complete set, TSV

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [(f'records/{i}', path) for i, path in enumerate(self.records)]
        return [*roles, ('gene2accession', self.gene2accession), ('hgnc', self.hgnc)]


def _split(versioned: str) -> tuple[str, int]:
    match = _VERSIONED.match(versioned)
    if match is None:
        raise build.BuildError(f'{versioned!r} is not a versioned accession')
    return match.group(1), int(match.group(2))


def _read_records(paths: tuple[pathlib.Path, ...], left_out: collections.Counter[str]) -> dict[str, flatfile.Record]:
    """Every record of every file, by version; a record the reader takes no CDS from is counted in `left_out`.

    Raises:
        build.BuildError: If a version is in two files, or a file breaks what the reader relies on.
    """
    out: dict[str, flatfile.Record] = {}
    for path in paths:
        refused: dict[str, str] = {}
        records = flatfile.read_records(
            path, sequences=True, left_out=refused, leave_out=(flatfile.UninterpretableError,)
        )
        for versioned in (*records, *refused):
            if versioned in out:
                raise build.BuildError(f'{path}: {versioned} is in an earlier file too')
        out.update(records)
        if refused:
            left_out['a CDS the reader does not interpret'] += len(refused)
    return out


@dataclasses.dataclass(frozen=True)
class _Coding:
    """An mRNA record to bundle: its sequence, its CDS, and the protein the CDS names, with its translation."""

    sequence: bytes
    cds: flatfile.RecordCds
    protein_id: str
    translation: bytes


def _translates(coding: _Coding) -> bool:
    """Whether the record's translation is the conceptual translation of its CDS's bases.

    A closed start is the initiator, written M whatever its codon; a closed end is a stop codon, which the
    translation leaves out. A codon the code cannot read, with an ambiguous base, matches any residue.
    """
    cds = coding.cds
    bases = coding.sequence[cds.start : cds.end + 1].decode('ascii')
    residues = [_CODE.get(bases[k : k + 3], 'X') for k in range(0, len(bases) - 2, 3)]
    if not cds.end_open:
        if not residues or residues[-1] != '*':
            return False
        residues.pop()
    if not cds.start_open and residues:
        residues[0] = 'M'
    translation = coding.translation.decode('ascii')
    return len(translation) == len(residues) and all(r in ('X', t) for r, t in zip(residues, translation, strict=True))


def _coding(versioned: str, record: flatfile.Record) -> _Coding | str:
    """The record to bundle, or why it is left out.

    Raises:
        build.BuildError: If a record with a CDS lacks its sequence or protein, which the reader never returns.
    """
    if record.molecule != 'mRNA':
        return 'not an mRNA'
    if record.cds is None:
        return 'no CDS'
    if record.sequence is None or record.protein is None:
        raise build.BuildError(f'{versioned}: read with a CDS but without its sequence or protein')
    if record.joined:
        return 'a CDS that is a join'
    if record.codon_start != 1:
        return 'a CDS read from other than its first base'
    if record.exceptions:
        return 'a CDS marked with an exception'
    coding = _Coding(record.sequence, record.cds, *record.protein)
    if not _translates(coding):
        return "a translation that is not the CDS's"
    return coding


def _read_genes(path: pathlib.Path, wanted: set[str]) -> dict[str, dict[str, str]]:
    """The human genes `gene2accession` names for each wanted RNA version: GeneID to symbol, by version.

    Raises:
        build.BuildError: If the table lacks a column this reader relies on.
    """
    out: dict[str, dict[str, str]] = collections.defaultdict(dict)
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        header = fh.readline().lstrip('#').rstrip('\n').split('\t')
        if missing := _GENE2ACCESSION_COLUMNS - set(header):
            raise build.BuildError(f'{path}: columns {sorted(missing)} are not in the table; its schema has changed')
        for row in csv.DictReader(fh, fieldnames=header, delimiter='\t'):
            versioned = row['RNA_nucleotide_accession.version']
            if row['tax_id'] == _HUMAN and versioned in wanted:
                out[versioned][row['GeneID']] = row['Symbol']
    return dict(out)


def _bundle(
    gene_id: str,
    symbol: str,
    versions: list[str],
    records: dict[str, _Coding],
    hgnc: dict[str, list[dict[str, str]]],
) -> bundle_pb2.GeneBundle:
    row, others = common.hgnc_row(hgnc.get(gene_id), symbol)
    gene = common.gene_message(
        feature_id=f'GeneID:{gene_id}',
        symbol=symbol,
        hgnc_id='',
        ncbi_gene_id=gene_id,
        ensembl_gene_id='',
        description='',
        synonyms=others,
        hgnc=row,
    )
    bundle = bundle_pb2.GeneBundle(gene=gene)
    digests: dict[str, bundle_pb2.Sequence] = {}
    for versioned in versions:
        record = records[versioned]
        accession, version = _split(versioned)
        protein_accession, protein_version = _split(record.protein_id)
        cds = record.cds
        bundle.transcripts.add(
            accession=accession,
            version=version,
            publisher=bundle_pb2.PUBLISHER_GENBANK,
            status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
            biotype='mRNA',
            sequence_digest=common.sequence(digests, record.sequence, bundle_pb2.ALPHABET_NUCLEOTIDE),
            cds=bundle_pb2.Cds(
                start_index=cds.start, end_index_inclusive=cds.end, start_open=cds.start_open, end_open=cds.end_open
            ),
            protein_accession=protein_accession,
            protein_version=protein_version,
        )
        bundle.proteins.add(
            accession=protein_accession,
            version=protein_version,
            sequence_digest=common.sequence(digests, record.translation, bundle_pb2.ALPHABET_PROTEIN),
            transcript_accession=accession,
            transcript_version=version,
        )
    bundle.sequences.extend(digests.values())
    return bundle


def bundles(release: Release) -> Iterator[bundle_pb2.GeneBundle]:
    """The release's bundles, one per gene, by GeneID; the build's report goes to stderr.

    Raises:
        build.BuildError: If a version is in two files, a file breaks what its reader relies on, or a
            gene has no symbol in `gene2accession` or HGNC.
    """
    left_out: collections.Counter[str] = collections.Counter()
    records = _read_records(release.records, left_out)
    kept: dict[str, _Coding] = {}
    for versioned, record in records.items():
        coding = _coding(versioned, record)
        if isinstance(coding, str):
            left_out[coding] += 1
        else:
            kept[versioned] = coding
    genes = _read_genes(release.gene2accession, set(kept))
    by_gene: dict[str, tuple[str, list[str]]] = {}
    for versioned in kept:
        if versioned not in genes:
            left_out['no human gene in gene2accession'] += 1
            continue
        for gene_id, symbol in genes[versioned].items():
            by_gene.setdefault(gene_id, (symbol, []))[1].append(versioned)
    bundled = sum(1 for v in kept if v in genes)
    print(
        f'{release.release}: {bundled + left_out.total()} records, {bundled} bundled in {len(by_gene)} genes; '
        f'left out: {common.counted(left_out)}',
        file=sys.stderr,
    )
    hgnc = common.read_hgnc(release.hgnc, 'entrez_id')
    for gene_id in sorted(by_gene, key=int):
        symbol, versions = by_gene[gene_id]
        yield _bundle(gene_id, symbol, sorted(versions), kept, hgnc)
