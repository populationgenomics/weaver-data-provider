"""The GenBank builder: synthetic flat files, gene2accession and HGNC in, unplaced gene bundles out."""

from __future__ import annotations

import gzip
import pathlib

import pytest
import weaver

from weaver_data_provider import build, refget, testing
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import cli, genbank
from weaver_data_provider.v1 import bundle_pb2

# GG, then ATG AAA CTG GTG TAA (M K L V, stop) at bases 3..17, then CC
CODING = 'GGATGAAACTGGTGTAACC'
CODING_CDS = '3..17'
CODING_PROTEIN = 'MKLV'
GENE2ACCESSION_COLUMNS = [
    '#tax_id',
    'GeneID',
    'status',
    'RNA_nucleotide_accession.version',
    'RNA_nucleotide_gi',
    'protein_accession.version',
    'protein_gi',
    'Symbol',
]
HGNC_COLUMNS = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']
HGNC = [
    ['HGNC:90001', 'FWD', 'forward gene', '', 'OLDFWD', '1001', ''],
    ['HGNC:90002', 'REV', 'reverse gene', '', '', '1002', ''],
]


def _record(
    versioned: str,
    sequence: str,
    *,
    molecule: str = 'mRNA',
    cds: str | None = None,
    protein: tuple[str, str] | None = None,
    codon_start: int = 1,
    qualifiers: tuple[str, ...] = (),
    features: tuple[str, ...] = (),
) -> str:
    """One GenBank flat-file record as NCBI writes it: the CDS and its qualifiers, any further features, ORIGIN."""
    accession = versioned.split('.')[0]
    lines = [
        f'LOCUS       {accession}  {len(sequence)} bp    {molecule}    linear   PRI 01-JAN-2026',
        'DEFINITION  synthetic.',
        f'ACCESSION   {accession}',
        f'VERSION     {versioned}',
        'FEATURES             Location/Qualifiers',
        f'     source          1..{len(sequence)}',
        '                     /db_xref="taxon:9606"',
    ]
    if cds is not None:
        lines.append(f'     CDS             {cds}')
        lines += [f'                     {q}' for q in (f'/codon_start={codon_start}', *qualifiers)]
        if protein is not None:
            lines.append(f'                     /protein_id="{protein[0]}"')
            lines.append(f'                     /translation="{protein[1]}"')
    lines += features
    lines.append('ORIGIN')
    for k in range(0, len(sequence), 60):
        tens = [sequence[j : j + 10].lower() for j in range(k, min(k + 60, len(sequence)), 10)]
        lines.append(f'{k + 1:>9} {" ".join(tens)}')
    lines.append('//')
    return '\n'.join(lines) + '\n'


def _gz(path: pathlib.Path, text: str) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        fh.write(text)
    return path


def _gene2accession(path: pathlib.Path, rows: list[tuple[str, str, str, str]]) -> pathlib.Path:
    """`gene2accession` rows of (tax_id, GeneID, RNA accession.version, Symbol), the rest as NCBI leaves them."""
    lines = ['\t'.join(GENE2ACCESSION_COLUMNS)]
    lines += ['\t'.join([tax, gene, '-', rna, '-', '-', '-', symbol]) for tax, gene, rna, symbol in rows]
    return _gz(path, '\n'.join(lines) + '\n')


def _hgnc(path: pathlib.Path) -> pathlib.Path:
    path.write_text('\n'.join('\t'.join(r) for r in [HGNC_COLUMNS, *HGNC]) + '\n', 'utf-8')
    return path


def _release(tmp_path: pathlib.Path, records: list[str], genes: list[tuple[str, str, str, str]]) -> genbank.Release:
    return genbank.Release(
        release='GB_1',
        records=(_gz(tmp_path / 'records.gb.gz', ''.join(records)),),
        gene2accession=_gene2accession(tmp_path / 'gene2accession.gz', genes),
        hgnc=_hgnc(tmp_path / 'hgnc.txt'),
    )


def _coding(versioned: str = 'ZZ000001.1') -> str:
    return _record(versioned, CODING, cds=CODING_CDS, protein=('ZZZ00001.1', CODING_PROTEIN))


def test_a_coding_mrna_is_bundled_with_its_sequence_cds_and_protein(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, [_coding('ZZ000001.1')], [('9606', '1001', 'ZZ000001.1', 'FWD')])
    (bundle,) = genbank.bundles(release)
    (transcript,) = bundle.transcripts
    assert (transcript.accession, transcript.version) == ('ZZ000001', 1)
    assert transcript.publisher == bundle_pb2.PUBLISHER_GENBANK
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (2, 16)
    assert transcript.sequence_digest == refget.digest(CODING.encode())
    assert (transcript.protein_accession, transcript.protein_version) == ('ZZZ00001', 1)
    (protein,) = bundle.proteins
    assert protein.sequence_digest == refget.digest(CODING_PROTEIN.encode())
    assert (protein.transcript_accession, protein.transcript_version) == ('ZZ000001', 1)


def test_the_genes_identity_is_hgncs(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, [_coding('ZZ000001.1')], [('9606', '1001', 'ZZ000001.1', 'OLDFWD')])
    (bundle,) = genbank.bundles(release)
    assert (bundle.gene.symbol, bundle.gene.hgnc_id, bundle.gene.ncbi_gene_id) == ('FWD', 'HGNC:90001', '1001')


def test_a_gene_hgnc_does_not_know_keeps_the_symbol_gene2accession_gives(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, [_coding('ZZ000001.1')], [('9606', '7777', 'ZZ000001.1', 'NOVEL1')])
    (bundle,) = genbank.bundles(release)
    assert (bundle.gene.symbol, bundle.gene.ncbi_gene_id) == ('NOVEL1', '7777')


def test_a_record_ncbi_ties_to_two_genes_is_in_each_bundle(tmp_path: pathlib.Path) -> None:
    genes = [('9606', '1001', 'ZZ000001.1', 'FWD'), ('9606', '1002', 'ZZ000001.1', 'REV')]
    bundles = list(genbank.bundles(_release(tmp_path, [_coding('ZZ000001.1')], genes)))
    assert [(b.gene.symbol, [t.accession for t in b.transcripts]) for b in bundles] == [
        ('FWD', ['ZZ000001']),
        ('REV', ['ZZ000001']),
    ]


def test_an_open_start_is_bundled_open(tmp_path: pathlib.Path) -> None:
    # AAA CTG GTG TAA from the record's first base, its initiator not in the record
    record = _record('ZZ000002.1', 'AAACTGGTGTAA', cds='<1..12', protein=('ZZZ00002.1', 'KLV'))
    (bundle,) = genbank.bundles(_release(tmp_path, [record], [('9606', '1001', 'ZZ000002.1', 'FWD')]))
    cds = bundle.transcripts[0].cds
    assert (cds.start_index, cds.start_open, cds.end_open) == (0, True, False)


@pytest.mark.parametrize(
    ('sequence', 'cds', 'translation'),
    [
        ('CTGAAATAA', '1..9', 'MK'),  # CTG begins a CDS, written M
        ('ATGCTNTAA', '1..9', 'ML'),  # CTN is Leu whatever N is
        ('ATGAAATRA', '1..9', 'MK'),  # TRA is a stop whatever R is
        ('ATGAANTAA', '1..9', 'MN'),  # AAN is Lys or Asn, so a record stating either is consistent
        ('ATGAANTAA', '1..9', 'MX'),  # or writing it undecided
        ('ATGAAAGG', '1..>8', 'MKG'),  # at an open end, GG decides Gly whatever the third base
        ('ATGAAAAT', '1..>8', 'MK'),  # AT is Ile or Met, so it is not translated
    ],
)
def test_a_cds_translated_as_insdc_translates_it_is_bundled(
    tmp_path: pathlib.Path, sequence: str, cds: str, translation: str
) -> None:
    record = _record('ZZ000002.1', sequence, cds=cds, protein=('ZZZ00002.1', translation))
    (bundle,) = genbank.bundles(_release(tmp_path, [record], [('9606', '1001', 'ZZ000002.1', 'FWD')]))
    assert bundle.proteins[0].sequence_digest == refget.digest(translation.encode())


_SECOND_CDS = (
    '     CDS             3..8',
    '                     /protein_id="ZZZ00009.1"',
    '                     /translation="MK"',
)


@pytest.mark.parametrize(
    ('record', 'reason'),
    [
        (_record('NM_000003.1', CODING, cds=CODING_CDS, protein=('NP_000003.1', 'MKLV')), 'a RefSeq record'),
        (_record('ZZ000003.1', CODING), 'no CDS'),
        (_record('ZZ000003.1', CODING, molecule='DNA', cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV')), 'not an mRNA'),
        (_record('ZZ000003.1', CODING, cds=CODING_CDS), 'no single CDS naming its protein'),
        (
            _record('ZZ000003.1', CODING, cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV'), features=_SECOND_CDS),
            'no single CDS naming its protein',
        ),
        (_record('ZZ000003.1', CODING, cds='complement(3..17)'), 'a CDS the reader does not interpret'),
        (
            _record('ZZ000003.1', CODING, cds='join(3..8,9..17)', protein=('ZZZ00003.1', 'MKLV')),
            'a CDS that is a join',
        ),
        (
            # C, then AAA CTG GTG TAA: a partial CDS whose first codon starts at its second base
            _record('ZZ000003.1', 'CAAACTGGTGTAA', cds='<1..13', protein=('ZZZ00003.1', 'KLV'), codon_start=2),
            'a CDS read from other than its first base',
        ),
        (
            _record(
                'ZZ000003.1', CODING, cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV'), qualifiers=('/transl_table=2',)
            ),
            'a CDS in another genetic code',
        ),
        (
            _record(
                'ZZ000003.1', CODING, cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV'), qualifiers=('/exception="x"',)
            ),
            'a CDS marked with an exception',
        ),
        (
            _record(
                'ZZ000003.1',
                CODING,
                cds=CODING_CDS,
                protein=('ZZZ00003.1', 'MKLV'),
                qualifiers=('/transl_except=(pos:6..8,aa:Sec)',),
            ),
            'a CDS marked with an exception',
        ),
        (
            _record(
                'ZZ000003.1', CODING, cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV'), qualifiers=('/trans_splicing',)
            ),
            'a CDS marked with an exception',
        ),
        (
            _record('ZZ000003.1', CODING, cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLA')),
            "a translation that is not the CDS's",
        ),
        (
            # AAA AAA CTG GTG TAA: AAA does not begin a CDS, so its first residue is Lys, not Met
            _record('ZZ000003.1', 'GGAAAAAACTGGTGTAACC', cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLV')),
            "a translation that is not the CDS's",
        ),
        (
            # a closed CDS one base longer than its whole codons
            _record('ZZ000003.1', CODING, cds='3..18', protein=('ZZZ00003.1', 'MKLV')),
            "a translation that is not the CDS's",
        ),
        (
            # ATG AAA CTG GTG TAC: a closed end with no stop codon
            _record('ZZ000003.1', 'GGATGAAACTGGTGTACCC', cds=CODING_CDS, protein=('ZZZ00003.1', 'MKLVY')),
            "a translation that is not the CDS's",
        ),
        (
            # AAN is Lys or Asn, so a record translating it Glu states another protein
            _record('ZZ000003.1', 'ATGAANTAA', cds='1..9', protein=('ZZZ00003.1', 'ME')),
            "a translation that is not the CDS's",
        ),
        (
            # CTN is Leu, so a record translating it Trp states another protein
            _record('ZZ000003.1', 'ATGCTNTAA', cds='1..9', protein=('ZZZ00003.1', 'MW')),
            "a translation that is not the CDS's",
        ),
    ],
)
def test_a_record_whose_cds_or_protein_is_not_plain_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], record: str, reason: str
) -> None:
    genes = [
        ('9606', '1001', 'ZZ000001.1', 'FWD'),
        ('9606', '1001', 'ZZ000003.1', 'FWD'),
        ('9606', '1001', 'NM_000003.1', 'FWD'),
    ]
    (bundle,) = genbank.bundles(_release(tmp_path, [_coding('ZZ000001.1'), record], genes))
    assert [t.accession for t in bundle.transcripts] == ['ZZ000001']
    assert f'{reason} 1' in capsys.readouterr().err


def test_a_record_of_no_human_gene_is_left_out_and_counted(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    genes = [('9606', '1001', 'ZZ000001.1', 'FWD'), ('10090', '2001', 'ZZ000004.1', 'Fwd')]
    (bundle,) = genbank.bundles(_release(tmp_path, [_coding('ZZ000001.1'), _coding('ZZ000004.1')], genes))
    assert [t.accession for t in bundle.transcripts] == ['ZZ000001']
    assert 'no human gene in gene2accession 1' in capsys.readouterr().err


def test_a_version_in_two_files_is_refused(tmp_path: pathlib.Path) -> None:
    first = _gz(tmp_path / 'first.gb.gz', _coding('ZZ000001.1'))
    second = _gz(tmp_path / 'second.gb.gz', _coding('ZZ000001.1'))
    release = genbank.Release(
        release='GB_1',
        records=(first, second),
        gene2accession=_gene2accession(tmp_path / 'g2a.gz', [('9606', '1001', 'ZZ000001.1', 'FWD')]),
        hgnc=_hgnc(tmp_path / 'hgnc.txt'),
    )
    with pytest.raises(build.BuildError, match=r'ZZ000001\.1 is in an earlier file too'):
        list(genbank.bundles(release))


def test_a_version_one_file_holds_uninterpretable_and_another_holds_is_refused(tmp_path: pathlib.Path) -> None:
    first = _gz(tmp_path / 'first.gb.gz', _record('ZZ000001.1', CODING, cds=CODING_CDS))
    second = _gz(tmp_path / 'second.gb.gz', _coding('ZZ000001.1'))
    release = genbank.Release(
        release='GB_1',
        records=(first, second),
        gene2accession=_gene2accession(tmp_path / 'g2a.gz', [('9606', '1001', 'ZZ000001.1', 'FWD')]),
        hgnc=_hgnc(tmp_path / 'hgnc.txt'),
    )
    with pytest.raises(build.BuildError, match=r'ZZ000001\.1 is in an earlier file too'):
        list(genbank.bundles(release))


def test_a_gene2accession_without_its_symbol_column_is_refused(tmp_path: pathlib.Path) -> None:
    release = genbank.Release(
        release='GB_1',
        records=(_gz(tmp_path / 'records.gb.gz', _coding('ZZ000001.1')),),
        gene2accession=_gz(tmp_path / 'g2a.gz', '#tax_id\tGeneID\tRNA_nucleotide_accession.version\n'),
        hgnc=_hgnc(tmp_path / 'hgnc.txt'),
    )
    with pytest.raises(build.BuildError, match='Symbol'):
        list(genbank.bundles(release))


@pytest.fixture
def served(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> provider_mod.BundleProvider:
    """A store the build commands made of the GenBank shard alone, beside the synthetic genome: FWD's ZZ000001.1."""
    built = testing.build(tmp_path / 'reference')
    release = _release(tmp_path, [_coding('ZZ000001.1')], [('9606', '1001', 'ZZ000001.1', 'FWD')])
    cli.main([
        'genbank', '--release', 'GB_1', '--records', str(release.records[0]),
        '--gene2accession', str(release.gene2accession), '--hgnc', str(release.hgnc),
        '--shards', str(tmp_path / 'shards'),
    ])  # fmt: skip
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    return provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'store')), built.genome)


def test_weaver_predicts_a_bundled_records_protein(served: provider_mod.BundleProvider) -> None:
    # c.4 is the record's base 6, the first of codon AAA (Lys2); GAA is Glu
    mapper = weaver.VariantMapper(served)
    assert mapper.c_to_p(weaver.parse('ZZ000001.1:c.4A>G')).format() == 'ZZZ00001.1:p.(Lys2Glu)'


def test_weaver_reads_a_name_numbered_from_the_records_first_base(served: provider_mod.BundleProvider) -> None:
    # n.6 is the record's sixth base, as BIC numbers U14680.1
    assert weaver.parse('ZZ000001.1:n.6A>G').validate(served)
    assert not weaver.parse('ZZ000001.1:n.6C>G').validate(served)
