"""The UniProt builder: synthetic flat-file entries, an isoform FASTA and HGNC in, bundles of proteins out."""

from __future__ import annotations

import gzip
import pathlib

import pytest
import weaver

from weaver_data_provider import build, refget, testing
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import cli, uniprot
from weaver_data_provider.v1 import bundle_pb2

CANONICAL = 'MKLVAAAA'
ISOFORM_2 = 'MKLVCCCC'
HGNC_COLUMNS = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']
HGNC = [
    ['HGNC:90001', 'FWD', 'forward gene', '', 'OLDFWD', '1001', ''],
    ['HGNC:90002', 'REV', 'reverse gene', '', '', '1002', ''],
]


def _entry(
    accession: str,
    residues: str,
    *,
    version: int = 1,
    taxon: str = '9606',
    gene: str = 'Name=FWD; Synonyms=OLDFWD;',
    references: tuple[str, ...] = ('GeneID; 1001; -.', 'HGNC; HGNC:90001; FWD.'),
    isoforms: tuple[tuple[str, str], ...] = (),
) -> str:
    """One entry as UniProt writes it, with a SEQUENCE CAUTION whose `Sequence=` names no isoform."""
    lines = [
        f'ID   {accession}_HUMAN               Reviewed;         {len(residues)} AA.',
        f'AC   {accession}; Q99999;',
        'DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.',
        f'DT   01-JAN-2001, sequence version {version}.',
        'DT   01-JAN-2026, entry version 50.',
        f'GN   {gene}',
        f'OX   NCBI_TaxID={taxon};',
    ]
    if isoforms:
        lines += [
            'CC   -!- ALTERNATIVE PRODUCTS:',
            f'CC       Event=Alternative splicing; Named isoforms={len(isoforms)};',
        ]
        for k, (ids, sequence) in enumerate(isoforms, 1):
            lines += [f'CC       Name={k};', f'CC         IsoId={ids}; Sequence={sequence};']
    lines += [
        'CC   -!- SEQUENCE CAUTION:',
        'CC       Sequence=AAA00001.1; Type=Erroneous initiation; Evidence={ECO:0000305};',
        *(f'DR   {r}' for r in references),
        f'SQ   SEQUENCE   {len(residues)} AA;  1000 MW;  0000000000000000 CRC64;',
        f'     {residues}',
        '//',
    ]
    return '\n'.join(lines) + '\n'


def _gz(path: pathlib.Path, text: str) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        fh.write(text)
    return path


def _release(tmp_path: pathlib.Path, entries: list[str], isoforms: dict[str, str]) -> uniprot.Release:
    fasta = ''.join(f'>sp|{i}|X_HUMAN Isoform of X OS=Homo sapiens OX=9606\n{r}\n' for i, r in isoforms.items())
    hgnc = tmp_path / 'hgnc.txt'
    hgnc.write_text('\n'.join('\t'.join(r) for r in [HGNC_COLUMNS, *HGNC]) + '\n', 'utf-8')
    return uniprot.Release(
        release='2026_01',
        entries=(_gz(tmp_path / 'entries.dat.gz', ''.join(entries)),),
        isoforms=_gz(tmp_path / 'varsplic.fasta.gz', fasta),
        hgnc=hgnc,
    )


def _names(bundle: bundle_pb2.GeneBundle) -> dict[str, str]:
    """Each protein's name, versioned or bare, and its digest."""
    return {store_mod.protein_name(p): p.sequence_digest for p in bundle.proteins}


_TWO_ISOFORMS = (('P00001-1', 'Displayed'), ('P00001-2', 'VSP_000001'))


def test_a_canonical_sequence_is_keyed_by_accession_and_sequence_version(tmp_path: pathlib.Path) -> None:
    (bundle,) = uniprot.bundles(_release(tmp_path, [_entry('P00001', CANONICAL, version=3)], {}))
    assert _names(bundle) == {'P00001.3': refget.digest(CANONICAL.encode())}


def test_the_displayed_isoform_id_is_a_second_key_for_the_canonical(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, version=3, isoforms=_TWO_ISOFORMS)
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2}))
    assert _names(bundle)['P00001-1'] == refget.digest(CANONICAL.encode())


def test_another_isoform_is_keyed_by_its_id_with_the_isoform_files_sequence(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, version=3, isoforms=_TWO_ISOFORMS)
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2}))
    assert _names(bundle)['P00001-2'] == refget.digest(ISOFORM_2.encode())


def test_a_canonical_displayed_as_another_isoform_leaves_1_its_own_sequence(tmp_path: pathlib.Path) -> None:
    # as release 2026_01 made P08195-2 the canonical in place of P08195-1
    isoforms = (('P00001-1', 'VSP_000001'), ('P00001-2', 'Displayed'))
    entry = _entry('P00001', CANONICAL, version=4, isoforms=isoforms)
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-1': ISOFORM_2}))
    assert _names(bundle) == {
        'P00001.4': refget.digest(CANONICAL.encode()),
        'P00001-2': refget.digest(CANONICAL.encode()),
        'P00001-1': refget.digest(ISOFORM_2.encode()),
    }


@pytest.mark.parametrize(
    'isoform',
    [
        ('P00001-2, P00009-1', 'VSP_000001'),  # identical to another entry's canonical: that entry keys its own
        ('P00009-2', 'External'),  # described in another entry
    ],
)
def test_an_isoform_another_entry_describes_is_not_keyed_here(tmp_path: pathlib.Path, isoform: tuple[str, str]) -> None:
    entry = _entry('P00001', CANONICAL, isoforms=(('P00001-1', 'Displayed'), isoform))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2, 'P00009-2': ISOFORM_2}))
    assert not {'P00009-1', 'P00009-2'} & set(_names(bundle))


@pytest.mark.parametrize(
    ('isoform', 'reason'),
    [
        (('P00001-2', 'Not described'), 'isoforms the release does not describe 1'),
        (('P00001-3', 'VSP_000003'), 'isoforms the isoform file lacks 1'),
    ],
)
def test_an_isoform_with_no_sequence_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], isoform: tuple[str, str], reason: str
) -> None:
    entry = _entry('P00001', CANONICAL, isoforms=(('P00001-1', 'Displayed'), isoform))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert set(_names(bundle)) == {'P00001.1', 'P00001-1'}
    assert reason in capsys.readouterr().err


def test_an_entry_naming_two_genes_is_in_each_bundle(tmp_path: pathlib.Path) -> None:
    references = ('HGNC; HGNC:90001; FWD.', 'HGNC; HGNC:90002; REV.')
    bundles = list(uniprot.bundles(_release(tmp_path, [_entry('P00001', CANONICAL, references=references)], {})))
    assert [(b.gene.symbol, list(_names(b))) for b in bundles] == [('FWD', ['P00001.1']), ('REV', ['P00001.1'])]


def test_an_entry_naming_its_gene_by_geneid_alone_is_bundled_under_hgncs_gene(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, references=('GeneID; 1002; -.',), gene='Name=SOMETHING;')
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert (bundle.gene.symbol, bundle.gene.hgnc_id) == ('REV', 'HGNC:90002')


def test_an_entry_of_a_gene_hgnc_does_not_know_keeps_its_own_gene_name(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, references=('GeneID; 7777; -.',), gene='Name=NOVEL1;')
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert (bundle.gene.symbol, bundle.gene.ncbi_gene_id) == ('NOVEL1', '7777')


@pytest.mark.parametrize(
    ('entry', 'reason'),
    [
        (_entry('P00002', CANONICAL, taxon='10090'), 'entries of another species 1'),
        (_entry('P00002', CANONICAL, references=()), 'entries naming no gene 1'),
        (
            _entry('P00002', CANONICAL, references=('GeneID; 7777; -.',), gene=''),
            'entries naming a gene with no symbol 1',
        ),
    ],
)
def test_an_entry_with_no_human_gene_to_bundle_it_under_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], entry: str, reason: str
) -> None:
    (bundle,) = uniprot.bundles(_release(tmp_path, [_entry('P00001', CANONICAL), entry], {}))
    assert list(_names(bundle)) == ['P00001.1']
    assert reason in capsys.readouterr().err


def test_an_accession_in_two_entries_is_refused(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, [_entry('P00001', CANONICAL), _entry('P00001', ISOFORM_2)], {})
    with pytest.raises(build.BuildError, match='P00001: in two entries'):
        list(uniprot.bundles(release))


@pytest.fixture
def served(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> provider_mod.BundleProvider:
    """A store the build commands made of the UniProt shard alone, beside the synthetic genome: FWD's P00001."""
    built = testing.build(tmp_path / 'reference')
    entry = _entry('P00001', CANONICAL, version=3, isoforms=_TWO_ISOFORMS)
    release = _release(tmp_path, [entry], {'P00001-2': ISOFORM_2})
    cli.main([
        'uniprot', '--release', release.release, '--entries', str(release.entries[0]),
        '--isoforms', str(release.isoforms), '--hgnc', str(release.hgnc), '--shards', str(tmp_path / 'shards'),
    ])  # fmt: skip
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    return provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'store')), built.genome)


# residue 5 is the canonical's Ala and isoform 2's Cys
@pytest.mark.parametrize('variant', ['P00001.3:p.Ala5Gly', 'P00001-1:p.Ala5Gly', 'P00001-2:p.Cys5Gly'])
def test_weaver_reads_a_name_on_each_key_the_entry_states(served: provider_mod.BundleProvider, variant: str) -> None:
    assert weaver.parse(variant).validate(served)


def test_weaver_refuses_a_residue_the_isoform_does_not_have(served: provider_mod.BundleProvider) -> None:
    assert not weaver.parse('P00001-2:p.Ala5Gly').validate(served)


def test_a_uniprot_protein_is_not_listed_for_its_gene(served: provider_mod.BundleProvider) -> None:
    assert served.get_symbol_accessions('FWD', 'gene', 'p') == []


def test_a_bare_accession_held_by_version_asks_for_the_version(served: provider_mod.BundleProvider) -> None:
    with pytest.raises(weaver.DataProviderError, match=r'P00001: .*give the versioned accession'):
        served.get_seq('P00001', 0, None, weaver.IdentifierType.ProteinAccession)
