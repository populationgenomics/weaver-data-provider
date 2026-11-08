"""The UniProt builder: synthetic flat-file entries, an isoform FASTA and HGNC in, bundles of proteins out."""

from __future__ import annotations

import gzip
import pathlib
import re

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
    secondaries: tuple[str, ...] = (),
) -> str:
    """One entry as UniProt writes it, with a SEQUENCE CAUTION whose `Sequence=` names no isoform."""
    lines = [
        f'ID   {accession}_HUMAN               Reviewed;         {len(residues)} AA.',
        'AC   ' + ' '.join(f'{a};' for a in (accession, *secondaries)),
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


def _counted(err: str, reason: str) -> int:
    """How many the build's report counts for a reason."""
    found = re.search(rf'\b{re.escape(reason)} (\d+)\b', err)
    return int(found.group(1)) if found else 0


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


def test_an_isoform_id_written_under_another_accession_is_keyed_where_it_is_described(tmp_path: pathlib.Path) -> None:
    # as P0DPK3 alone describes A0A0B4J2B3-2, an accession no entry holds
    entry = _entry('P00001', CANONICAL, isoforms=(('P00001-1', 'Displayed'), ('P00001-2, P00009-2', 'VSP_000001')))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2}))
    assert _names(bundle)['P00009-2'] == refget.digest(ISOFORM_2.encode())


@pytest.mark.parametrize(
    ('isoform', 'reason'),
    [
        (('P00001-2', 'Not described'), 'isoforms the release does not describe'),
        (('P00001-3', 'VSP_000003'), 'isoforms the isoform file lacks'),
    ],
)
def test_an_isoform_with_no_sequence_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], isoform: tuple[str, str], reason: str
) -> None:
    entry = _entry('P00001', CANONICAL, version=7, isoforms=(('P00001-1', 'Displayed'), isoform))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert set(_names(bundle)) == {'P00001.7', 'P00001-1'}
    assert _counted(capsys.readouterr().err, reason) == 1


def test_an_entry_naming_two_genes_is_in_each_bundle(tmp_path: pathlib.Path) -> None:
    references = ('HGNC; HGNC:90001; FWD.', 'HGNC; HGNC:90002; REV.')
    entry = _entry('P00001', CANONICAL, version=7, references=references)
    bundles = list(uniprot.bundles(_release(tmp_path, [entry], {})))
    assert [(b.gene.symbol, list(_names(b))) for b in bundles] == [('FWD', ['P00001.7']), ('REV', ['P00001.7'])]


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
        (_entry('P00002', CANONICAL, taxon='10090'), 'entries of another species'),
        (_entry('P00002', CANONICAL, references=()), 'entries naming no gene'),
        (
            _entry('P00002', CANONICAL, references=('GeneID; 7777; -.',), gene=''),
            'entries naming a gene with no symbol',
        ),
    ],
)
def test_an_entry_with_no_human_gene_to_bundle_it_under_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], entry: str, reason: str
) -> None:
    (bundle,) = uniprot.bundles(_release(tmp_path, [_entry('P00001', CANONICAL, version=7), entry], {}))
    assert list(_names(bundle)) == ['P00001.7']
    assert _counted(capsys.readouterr().err, reason) == 1


def test_a_withdrawn_hgnc_id_names_no_gene_so_the_geneid_is_read(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, references=('HGNC; HGNC:99999; OLDREV.', 'GeneID; 1002; -.'))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert (bundle.gene.symbol, bundle.gene.hgnc_id) == ('REV', 'HGNC:90002')


def test_an_entry_with_no_gene_name_takes_its_orf_name(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, references=('GeneID; 7777; -.',), gene='ORFNames=C1orf999;')
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert bundle.gene.symbol == 'C1orf999'


def test_an_entry_with_no_gene_or_orf_name_takes_its_hgnc_references_symbol(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, references=('HGNC; HGNC:99999; OLDSYM.', 'GeneID; 7777; -.'), gene='')
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert bundle.gene.symbol == 'OLDSYM'


# An entry as UniProt wraps it: a second AC line, a gene name with evidence and a second gene, an isoform ID list
# whose Sequence runs onto the next line, a variant list wrapped mid-list, and a `Sequence=` outside the isoforms.
_WRAPPED = """ID   FWD_HUMAN               Reviewed;         8 AA.
AC   P00001; Q00001; Q00002; Q00003; Q00004; Q00005; Q00006; Q00007; Q00008;
AC   Q00009;
DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.
DT   01-JAN-2001, sequence version 5.
DT   01-JAN-2026, entry version 50.
GN   Name=FWD {ECO:0000312|HGNC:HGNC:90001}; Synonyms=OLDFWD;
GN   and
GN   Name=REV {ECO:0000312|HGNC:HGNC:90002};
OX   NCBI_TaxID=9606 {ECO:0000312|EMBL:AAA00001.1};
CC   -!- ALTERNATIVE PRODUCTS:
CC       Event=Alternative splicing; Named isoforms=2;
CC       Name=1;
CC         IsoId=P00001-1;
CC         Sequence=Displayed;
CC       Name=2; Synonyms=Short {ECO:0000303|PubMed:1};
CC         IsoId=P00001-2; Sequence=VSP_000001,
CC         VSP_000002;
CC   -!- FUNCTION: Synthetic. Sequence=AAA00009.1 is not an isoform.
DR   GeneID; 1001; -.
DR   HGNC; HGNC:90001; FWD.
DR   HGNC; HGNC:90002; REV.
SQ   SEQUENCE   8 AA;  1000 MW;  0000000000000000 CRC64;
     MKLVAAAA
//
"""


def test_an_entry_wrapped_as_uniprot_writes_it_is_read_whole(tmp_path: pathlib.Path) -> None:
    bundles = list(uniprot.bundles(_release(tmp_path, [_WRAPPED], {'P00001-2': ISOFORM_2})))
    assert [(b.gene.symbol, sorted(_names(b))) for b in bundles] == [
        ('FWD', ['P00001-1', 'P00001-2', 'P00001.5']),
        ('REV', ['P00001-1', 'P00001-2', 'P00001.5']),
    ]


def test_an_isoform_id_its_own_entry_describes_elsewhere_is_keyed_where_it_is_displayed(tmp_path: pathlib.Path) -> None:
    # as P0DP91 displays IsoId=P0DP91-1, Q03468-2 and Q03468 lists Q03468-2 as External
    displaying = _entry('P00001', CANONICAL, isoforms=(('P00001-1, P00002-2', 'Displayed'),))
    describing = _entry(
        'P00002',
        ISOFORM_2,
        references=('HGNC; HGNC:90002; REV.',),
        isoforms=(('P00002-1', 'Displayed'), ('P00002-2, P00001-1', 'External')),
    )
    bundles = {b.gene.symbol: b for b in uniprot.bundles(_release(tmp_path, [displaying, describing], {}))}
    assert _names(bundles['FWD'])['P00002-2'] == refget.digest(CANONICAL.encode())


def test_an_isoform_id_under_a_secondary_accession_the_entry_alone_holds_is_keyed(tmp_path: pathlib.Path) -> None:
    # as P60953 (CDC42) displays IsoId=P60953-2, P21181-4, P21181 having been merged into it
    entry = _entry('P00001', CANONICAL, secondaries=('Q00007',), isoforms=(('P00001-1, Q00007-4', 'Displayed'),))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert _names(bundle)['Q00007-4'] == refget.digest(CANONICAL.encode())


def test_a_variant_isoforms_secondary_id_takes_the_sequence_the_file_gives_its_primary_id(
    tmp_path: pathlib.Path,
) -> None:
    isoforms = (('P00001-1', 'Displayed'), ('P00001-2, Q00007-2', 'VSP_000001'))
    entry = _entry('P00001', CANONICAL, secondaries=('Q00007',), isoforms=isoforms)
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2}))
    assert _names(bundle)['Q00007-2'] == refget.digest(ISOFORM_2.encode())


def test_isoform_ids_of_an_entry_split_among_several_are_each_keyed_where_described(tmp_path: pathlib.Path) -> None:
    # as O60344, a secondary accession of P0DPD6, P0DPD7 and P0DPD8, each of which displays one of its IDs
    first = _entry('P00001', CANONICAL, secondaries=('Q00007',), isoforms=(('P00001-1, Q00007-1', 'Displayed'),))
    second = _entry(
        'P00002',
        ISOFORM_2,
        secondaries=('Q00007',),
        references=('HGNC; HGNC:90002; REV.',),
        isoforms=(('P00002-1, Q00007-2', 'Displayed'), ('Q00007-1', 'External')),
    )
    bundles = {b.gene.symbol: _names(b) for b in uniprot.bundles(_release(tmp_path, [first, second], {}))}
    assert (bundles['FWD']['Q00007-1'], bundles['REV']['Q00007-2']) == (
        refget.digest(CANONICAL.encode()),
        refget.digest(ISOFORM_2.encode()),
    )


def test_an_isoform_id_entries_give_different_sequences_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # as P0DX52 and P0DX53 both display Q08AF8-1, the ID of the entry split between them
    first = _entry('P00001', CANONICAL, isoforms=(('P00001-1, Q00007-1', 'Displayed'),))
    second = _entry(
        'P00002', ISOFORM_2, references=('HGNC; HGNC:90002; REV.',), isoforms=(('P00002-1, Q00007-1', 'Displayed'),)
    )
    bundles = list(uniprot.bundles(_release(tmp_path, [first, second], {})))
    assert not any('Q00007-1' in _names(b) for b in bundles)
    assert _counted(capsys.readouterr().err, 'isoform IDs entries give different sequences') == 1


def test_an_isoform_id_entries_give_one_sequence_is_keyed_in_each(tmp_path: pathlib.Path) -> None:
    first = _entry('P00001', CANONICAL, isoforms=(('P00001-1, Q00007-1', 'Displayed'),))
    second = _entry(
        'P00002', CANONICAL, references=('HGNC; HGNC:90002; REV.',), isoforms=(('P00002-1, Q00007-1', 'Displayed'),)
    )
    bundles = list(uniprot.bundles(_release(tmp_path, [first, second], {})))
    assert [b.gene.symbol for b in bundles if 'Q00007-1' in _names(b)] == ['FWD', 'REV']


def test_an_isoform_the_entry_says_is_described_elsewhere_is_not_keyed_here(tmp_path: pathlib.Path) -> None:
    entry = _entry('P00001', CANONICAL, isoforms=(('P00001-1', 'Displayed'), ('P00001-2', 'External')))
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {'P00001-2': ISOFORM_2}))
    assert 'P00001-2' not in _names(bundle)


def test_an_entrys_gene_name_comes_before_its_hgnc_references_symbol(tmp_path: pathlib.Path) -> None:
    references = ('HGNC; HGNC:99999; OLDSYM.', 'GeneID; 7777; -.')
    entry = _entry('P00001', CANONICAL, references=references, gene='Name=NOVEL1;')
    (bundle,) = uniprot.bundles(_release(tmp_path, [entry], {}))
    assert bundle.gene.symbol == 'NOVEL1'


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
