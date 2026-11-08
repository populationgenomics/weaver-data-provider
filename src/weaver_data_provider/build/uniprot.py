"""Cut one UniProt release's human entries into gene bundles of proteins: each entry's canonical sequence and isoforms.

The inputs are a release's entries in UniProt's flat-file format — Swiss-Prot's reviewed human entries, and
TrEMBL's unreviewed ones if the caller gives them too — the release's isoform sequences (`varsplic`), and
HGNC's complete set. A UniProt protein has no transcript and no placement; the provider serves it by
name and by digest, and never lists it for a gene (`docs/design/genbank-and-uniprot.md`).

Each entry's canonical sequence is bundled under its accession and sequence version, `P38398.2`, the
spelling NCBI Protein uses. An entry with alternative isoforms has its canonical keyed again under the
displayed isoform's IDs, which are not always `-1`, and each other isoform under its ID alone,
`P38398-2`, with the sequence this release gives it: an isoform ID has no version. An entry goes into the
bundle of each gene it cross-references, by HGNC ID, or by GeneID where it names no gene HGNC holds.
Left out, and counted in the build's report: an entry of another species, an entry naming no gene or one
with no symbol, and an isoform whose sequence the release does not describe or the isoform file lacks.
An entry keys the isoform IDs written under its primary accession or under a secondary one it alone holds,
the IDs of an entry since merged into it; the displayed isoform also goes by another entry's ID where that
entry lists the ID as described elsewhere. An ID under another entry's primary accession is that entry's to
key; one under an accession several entries hold, an entry since split, names no one sequence and is left
out and counted. A name the release gives two sequences fails the build.
"""

from __future__ import annotations

import collections
import dataclasses
import gzip
import pathlib
import re
import sys
from collections.abc import Iterator

from weaver_data_provider import build
from weaver_data_provider.build import common
from weaver_data_provider.v1 import bundle_pb2

_HUMAN = '9606'
_TAXON = re.compile(r'NCBI_TaxID=(\d+)')
_SEQUENCE_VERSION = re.compile(r'sequence version (\d+)\.')
_GENE_NAME = re.compile(r'(?:^|; )Name=([^;{]+)')
_ORF_NAME = re.compile(r'ORFNames=([^;{,]+)')
_ISOFORM = re.compile(r'IsoId=([^;]+); Sequence=([^;]+);')
_ALTERNATIVE_PRODUCTS = 'ALTERNATIVE PRODUCTS:'


@dataclasses.dataclass(frozen=True)
class Release:
    """One UniProt release's entries and isoform sequences, as local files."""

    release: str  # UniProt's release identifier, "2026_04"
    entries: tuple[pathlib.Path, ...]  # flat-file entries, gzipped: Swiss-Prot's, and TrEMBL's if wanted
    isoforms: pathlib.Path  # the release's isoform sequences, uniprot_sprot_varsplic.fasta.gz
    hgnc: pathlib.Path  # HGNC complete set, TSV

    def inputs(self) -> list[tuple[str, pathlib.Path]]:
        """Every file read, by role, for the shard's record of what it was cut from."""
        roles = [(f'entries/{i}', path) for i, path in enumerate(self.entries)]
        return [*roles, ('isoforms', self.isoforms), ('hgnc', self.hgnc)]


@dataclasses.dataclass(frozen=True)
class _Isoform:
    ids: tuple[str, ...]
    sequence: str  # `Displayed`, `External`, `Not described`, or the variant spans that make it


@dataclasses.dataclass(frozen=True)
class _Entry:
    accession: str  # the primary accession
    secondaries: tuple[str, ...]  # accessions merged into this entry, under which old isoform IDs were written
    sequence_version: int
    taxon: str
    symbol: str  # the entry's gene name, else its ORF name, else its HGNC cross-reference's symbol; or empty
    hgnc_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    isoforms: tuple[_Isoform, ...]
    residues: bytes


def _isoforms(comments: list[str]) -> tuple[_Isoform, ...]:
    """The isoforms the entry's ALTERNATIVE PRODUCTS comment lists, from its CC lines."""
    text, inside = [], False
    for line in comments:
        if line.startswith('-!- '):
            inside = line.startswith(f'-!- {_ALTERNATIVE_PRODUCTS}')
        elif inside:
            text.append(line)
    return tuple(
        _Isoform(tuple(i.strip() for i in ids.split(',')), sequence.strip())
        for ids, sequence in _ISOFORM.findall(' '.join(text))
    )


def _entry(path: pathlib.Path, lines: list[str]) -> _Entry:
    """One entry from its lines, `//` excluded.

    Raises:
        build.BuildError: If the entry states no accession, sequence version, taxon or sequence.
    """
    by_code: dict[str, list[str]] = collections.defaultdict(list)
    residues: list[str] = []
    in_sequence = False
    for line in lines:
        if in_sequence:
            residues.append(''.join(line.split()))
        elif line.startswith('SQ   '):
            in_sequence = True
        else:
            by_code[line[:2]].append(line[5:].rstrip('\n'))
    accessions = ' '.join(by_code['AC']).replace(';', ' ').split()
    versions = [m.group(1) for line in by_code['DT'] if (m := _SEQUENCE_VERSION.search(line))]
    taxa = [m.group(1) for line in by_code['OX'] if (m := _TAXON.search(line))]
    if not accessions or len(versions) != 1 or len(taxa) != 1 or not residues:
        named = accessions[0] if accessions else 'an entry'
        raise build.BuildError(f'{path}: {named} states no single accession, sequence version, taxon and sequence')
    genes = ' '.join(by_code['GN'])
    names = [m.group(1).strip() for m in (_GENE_NAME.search(genes), _ORF_NAME.search(genes)) if m is not None]
    references = [line.split(';') for line in by_code['DR']]
    hgnc_symbols = [r[2].strip().rstrip('.') for r in references if r[0] == 'HGNC' and len(r) > 2]
    return _Entry(
        accession=accessions[0],
        secondaries=tuple(accessions[1:]),
        sequence_version=int(versions[0]),
        taxon=taxa[0],
        symbol=next(iter(names + hgnc_symbols), ''),
        hgnc_ids=tuple(r[1].strip() for r in references if r[0] == 'HGNC'),
        gene_ids=tuple(r[1].strip() for r in references if r[0] == 'GeneID'),
        isoforms=_isoforms([c.lstrip() for c in by_code['CC']]),
        residues=''.join(residues).upper().encode('ascii'),
    )


def _read_entries(paths: tuple[pathlib.Path, ...]) -> Iterator[_Entry]:
    """Every entry of every file.

    Raises:
        build.BuildError: If a file ends inside an entry, or an entry is malformed.
    """
    for path in paths:
        lines: list[str] = []
        with gzip.open(path, 'rt', encoding='utf-8') as fh:
            for line in fh:
                if line.startswith('//'):
                    yield _entry(path, lines)
                    lines = []
                else:
                    lines.append(line)
        if any(line.strip() for line in lines):
            raise build.BuildError(f'{path}: ends inside an entry, before its //')


def _read_isoforms(path: pathlib.Path) -> dict[str, bytes]:
    """The isoform file's sequences by isoform ID: `>sp|P38398-2|BRCA1_HUMAN ...` is `P38398-2`.

    Raises:
        build.BuildError: If a header is not `db|id|name`.
    """
    out: dict[str, bytes] = {}
    for header, residues in common.read_fasta(path).items():
        fields = header.split('|')
        if len(fields) != 3:
            raise build.BuildError(f'{path}: a header {header!r} that is not db|id|name')
        out[fields[1]] = residues
    return out


@dataclasses.dataclass
class _Gene:
    """The proteins bundled for one gene, and how it is identified."""

    hgnc: dict[str, str] | None
    symbol: str
    hgnc_id: str
    gene_id: str
    synonyms: list[str]
    proteins: dict[str, tuple[int, bytes]] = dataclasses.field(default_factory=dict)  # accession: (version, residues)


def _genes_of(
    entry: _Entry, by_hgnc_id: dict[str, list[dict[str, str]]], by_gene_id: dict[str, list[dict[str, str]]]
) -> list[tuple[str, _Gene]]:
    """The genes an entry cross-references, each keyed by its HGNC ID where HGNC knows it, else by GeneID.

    An HGNC ID HGNC does not hold is withdrawn, so it names no gene; the entry's GeneIDs are read instead.
    """
    out: list[tuple[str, _Gene]] = []
    for hgnc_id in entry.hgnc_ids:
        rows = by_hgnc_id.get(hgnc_id)
        if rows:
            row, others = common.hgnc_row(rows, entry.symbol)
            out.append(
                (hgnc_id, _Gene(row, entry.symbol, hgnc_id, row['entrez_id'] if row is not None else '', others))
            )
    if out:
        return out
    for gene_id in entry.gene_ids:
        row, others = common.hgnc_row(by_gene_id.get(gene_id), entry.symbol)
        hgnc_id = row['hgnc_id'] if row is not None else ''
        out.append((hgnc_id or f'GeneID:{gene_id}', _Gene(row, entry.symbol, hgnc_id, gene_id, others)))
    return out


def _bundle(gene: _Gene) -> bundle_pb2.GeneBundle:
    message = common.gene_message(
        feature_id=gene.hgnc_id or f'GeneID:{gene.gene_id}',
        symbol=gene.symbol,
        hgnc_id=gene.hgnc_id,
        ncbi_gene_id=gene.gene_id,
        ensembl_gene_id='',
        description='',
        synonyms=gene.synonyms,
        hgnc=gene.hgnc,
    )
    bundle = bundle_pb2.GeneBundle(gene=message)
    digests: dict[str, bundle_pb2.Sequence] = {}
    for accession, (version, residues) in gene.proteins.items():
        bundle.proteins.add(
            accession=accession,
            version=version,
            sequence_digest=common.sequence(digests, residues, bundle_pb2.ALPHABET_PROTEIN),
        )
    bundle.sequences.extend(digests.values())
    return bundle


@dataclasses.dataclass(frozen=True)
class _Owners:
    """Who may key an isoform ID, from every entry of the release."""

    primaries: frozenset[str]  # every entry's primary accession
    sole_secondaries: dict[str, str]  # a secondary accession one entry alone holds, to that entry's primary
    described_elsewhere: frozenset[str]  # IDs an entry lists under its own accession as another entry's canonical

    @classmethod
    def of(cls, entries: list[_Entry]) -> _Owners:
        primaries = frozenset(e.accession for e in entries)
        holders = collections.Counter(s for e in entries for s in e.secondaries)
        return cls(
            primaries=primaries,
            sole_secondaries={
                s: e.accession for e in entries for s in e.secondaries if holders[s] == 1 and s not in primaries
            },
            described_elsewhere=frozenset(
                i
                for e in entries
                for isoform in e.isoforms
                if isoform.sequence == 'External'
                for i in isoform.ids
                if i.startswith(f'{e.accession}-')
            ),
        )

    def owns(self, entry: _Entry, isoform_id: str) -> bool | None:
        """Whether the entry keys an isoform ID: True; False where another entry does; None where none can.

        An entry keys an ID written under its primary accession, or under a secondary accession it alone
        holds. An ID under another entry's primary is that entry's; one under an accession several entries
        hold, or none, names no one sequence.
        """
        accession = isoform_id.rpartition('-')[0]
        if accession == entry.accession or self.sole_secondaries.get(accession) == entry.accession:
            return True
        return False if accession in self.primaries else None


def _isoform_keys(
    entry: _Entry, owners: _Owners, isoforms: dict[str, bytes], left_out: collections.Counter[str]
) -> list[tuple[str, bytes]]:
    """The isoform IDs an entry keys and the residues each names; what it cannot key is counted in `left_out`.

    The canonical goes by every ID of the displayed isoform the entry owns, and by another entry's ID that
    entry lists as described here. Any other isoform goes by the IDs the entry owns, with the residues the
    isoform file gives the first of its IDs it holds. An isoform described in another entry is that entry's.
    """
    out: list[tuple[str, bytes]] = []
    for isoform in entry.isoforms:
        if isoform.sequence == 'External':
            continue
        owned = [i for i in isoform.ids if owners.owns(entry, i)]
        nobody = [i for i in isoform.ids if owners.owns(entry, i) is None]
        if isoform.sequence == 'Displayed':
            owned += [i for i in isoform.ids if i in owners.described_elsewhere and i not in owned]
            nobody = [i for i in nobody if i not in owners.described_elsewhere]
            out += [(i, entry.residues) for i in owned]
        elif isoform.sequence == 'Not described':
            left_out['isoforms the release does not describe'] += len(owned)
        else:
            residues = next((isoforms[i] for i in isoform.ids if i in isoforms), None)
            if residues is None:
                left_out['isoforms the isoform file lacks'] += len(owned)
            else:
                out += [(i, residues) for i in owned]
        left_out['isoform IDs of an accession no one entry holds'] += len(nobody)
    return out


def bundles(release: Release) -> Iterator[bundle_pb2.GeneBundle]:
    """The release's bundles, one per gene, in symbol order; the build's report goes to stderr.

    Raises:
        build.BuildError: If an entry or the isoform file is malformed, an accession is in two entries, or the
            release gives one name two sequences.
    """
    isoforms = _read_isoforms(release.isoforms)
    by_hgnc_id = common.read_hgnc(release.hgnc, 'hgnc_id')
    by_gene_id = common.read_hgnc(release.hgnc, 'entrez_id')
    genes: dict[str, _Gene] = {}
    keyed: dict[str, bytes] = {}  # every name keyed, to the residues it names
    left_out: collections.Counter[str] = collections.Counter()
    read = list(_read_entries(release.entries))
    if duplicated := [a for a, n in collections.Counter(e.accession for e in read).items() if n > 1]:
        raise build.BuildError(f'{duplicated[0]}: in two entries')
    owners = _Owners.of(read)
    isoform_ids = 0
    for entry in read:
        if entry.taxon != _HUMAN:
            left_out['entries of another species'] += 1
            continue
        named = _genes_of(entry, by_hgnc_id, by_gene_id)
        if not named:
            left_out['entries naming no gene'] += 1
            continue
        if any(not gene.symbol and gene.hgnc is None for _, gene in named):
            left_out['entries naming a gene with no symbol'] += 1
            continue
        keys = _isoform_keys(entry, owners, isoforms, left_out)
        proteins = [(entry.accession, entry.sequence_version, entry.residues)] + [(i, 0, r) for i, r in keys]
        isoform_ids += len(keys)
        for accession, _, residues in proteins:
            if keyed.setdefault(accession, residues) != residues:
                raise build.BuildError(f'{accession}: the release gives it two sequences')
        for key, gene in named:
            held = genes.setdefault(key, gene)
            held.proteins.update((accession, (version, residues)) for accession, version, residues in proteins)
    print(
        f'{release.release}: {len(read)} entries, bundled in {len(genes)} genes with {isoform_ids} isoform IDs; '
        f'left out: {common.counted(left_out)}',
        file=sys.stderr,
    )
    yield from sorted((_bundle(g) for g in genes.values()), key=lambda b: b.gene.symbol)
