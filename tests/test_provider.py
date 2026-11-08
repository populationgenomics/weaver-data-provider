"""The weaver data provider over the synthetic store and genome, driven through weaver itself."""

from __future__ import annotations

import pathlib

import pytest
import weaver

from weaver_data_provider import genome as genome_mod
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import refget, testing
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2


@pytest.fixture
def built(tmp_path: pathlib.Path) -> testing.Reference:
    return testing.build(tmp_path)


@pytest.fixture
def provider(built: testing.Reference) -> provider_mod.BundleProvider:
    return built.provider()


def _base(genome: str, position_1based: int) -> str:
    return genome[position_1based - 1]


def test_transcript_model_is_weaver_shaped(provider: provider_mod.BundleProvider) -> None:
    model = provider.get_transcript('NM_000001.1', None)
    assert model['ac'] == 'NM_000001.1'
    assert model['gene'] == 'FWD'
    assert model.get('strand') == 1
    assert model.get('reference_accession') == testing.CHROM
    assert 'exons' in model, 'the model is unplaced'
    (a, b) = model['exons']
    assert (a['transcript_start'], a['transcript_end'], a['reference_start'], a['reference_end']) == (0, 10, 100, 109)
    assert (b['transcript_start'], b['transcript_end'], b['reference_start'], b['reference_end']) == (10, 22, 200, 211)
    assert model['cds_start_index'] == 3  # CDS starts at genome 104, four bases into exon 1
    assert model['cds_end_index'] == 14  # and ends at genome 205, five into exon 2, after exon 1's 10
    with pytest.raises(weaver.DataProviderError, match='versioned'):
        provider.get_transcript('NM_000001', None)
    with pytest.raises(weaver.DataProviderError, match='not in the reference data'):
        provider.get_transcript('NM_999999.9', None)


def test_sequences_come_from_the_record_and_the_genome(
    provider: provider_mod.BundleProvider, built: testing.Reference
) -> None:
    assert provider.get_seq('NM_000001.1', 0, 5, weaver.IdentifierType.TranscriptAccession) == built.fwd[:5]
    assert provider.get_seq('NM_000002.1', 0, None, 'transcript') == built.rev
    assert provider.get_seq('NP_000001.1', 0, 2, 'protein') == 'MK'
    genome = testing.genome_sequence()
    assert provider.get_seq(testing.CHROM, 100, 110, weaver.IdentifierType.GenomicAccession) == genome[100:110]


def test_weaver_projects_a_coding_variant_onto_the_genome(provider: provider_mod.BundleProvider) -> None:
    genome = testing.genome_sequence()
    mapper = weaver.VariantMapper(provider)
    ref = _base(genome, testing.FWD_CDS[0])  # c.1 is the CDS's first base, genome 104 on the plus strand
    alt = 'G' if ref != 'G' else 'C'
    variant = weaver.parse(f'NM_000001.1:c.1{ref}>{alt}')
    assert variant.validate(provider)
    assert mapper.c_to_g(variant, testing.CHROM).format() == f'{testing.CHROM}:g.{testing.FWD_CDS[0]}{ref}>{alt}'
    # the second exon: c.8 is the first base of exon 2 (exon 1 carries c.1..c.7)
    ref2 = _base(genome, testing.FWD_EXONS[1][0])
    alt2 = 'G' if ref2 != 'G' else 'C'
    assert mapper.c_to_g(weaver.parse(f'NM_000001.1:c.8{ref2}>{alt2}'), testing.CHROM).format() == (
        f'{testing.CHROM}:g.{testing.FWD_EXONS[1][0]}{ref2}>{alt2}'
    )


def test_weaver_projects_a_genomic_variant_onto_a_minus_strand_transcript(
    provider: provider_mod.BundleProvider,
) -> None:
    genome = testing.genome_sequence()
    mapper = weaver.VariantMapper(provider)
    # REV's start codon begins at genome 1106 (the CDS's highest coordinate on the minus strand)
    position = testing.REV_CDS[1]
    ref = _base(genome, position)
    alt = 'G' if ref != 'G' else 'C'
    projected = mapper.g_to_c(weaver.parse(f'{testing.CHROM}:g.{position}{ref}>{alt}'), 'NM_000002.1').format()
    comp = str.maketrans('ACGT', 'TGCA')
    assert projected == f'NM_000002.1:c.1{ref.translate(comp)}>{alt.translate(comp)}'


def test_the_reference_base_is_checked_against_the_record_not_the_genome(
    provider: provider_mod.BundleProvider, built: testing.Reference
) -> None:
    i = testing.REV_SUBSTITUTED_TRANSCRIPT_INDEX
    record_base, genome_base = built.rev[i], built.rev_genomic[i]
    # transcript index 2 lies in the 5' UTR: the CDS starts at index 4, so index 2 is c.-2
    assert weaver.parse(f'NM_000002.1:c.-2{record_base}>T').validate(provider)
    assert not weaver.parse(f'NM_000002.1:c.-2{genome_base}>T').validate(provider)


def test_symbol_and_accession_crosswalks(provider: provider_mod.BundleProvider) -> None:
    kinds = weaver.IdentifierType
    assert provider.get_symbol_accessions('FWD', 'gene', 'c') == [(kinds.TranscriptAccession, 'NM_000001.1')]
    assert provider.get_symbol_accessions('OLDFWD', 'gene', 'c') == [(kinds.TranscriptAccession, 'NM_000001.1')]
    assert provider.get_symbol_accessions('NM_000001.1', 'c', 'p') == [(kinds.ProteinAccession, 'NP_000001.1')]
    assert provider.get_symbol_accessions('REV', 'gene', 'p') == [(kinds.ProteinAccession, 'NP_000002.1')]
    assert provider.get_symbol_accessions('NOSUCH', 'gene', 'c') == []
    assert provider.get_symbol_accessions('NM_000001.1', 'c', 'g') == []


def test_identifier_kinds(provider: provider_mod.BundleProvider) -> None:
    kinds = weaver.IdentifierType
    assert provider.get_identifier_type(testing.CHROM) == kinds.GenomicAccession
    assert provider.get_identifier_type('NM_000001.1') == kinds.TranscriptAccession
    assert provider.get_identifier_type('NP_000001.1') == kinds.ProteinAccession
    assert provider.get_identifier_type('FWD') == kinds.GeneSymbol
    assert provider.get_identifier_type('nothing-here') == kinds.Unknown


def test_region_search_finds_a_transcript_over_its_intron(provider: provider_mod.BundleProvider) -> None:
    # FWD is placed over 100..211 (0-based), its intron 110..199
    assert provider.get_transcripts_for_region(testing.CHROM, 150, 150) == ['NM_000001.1']
    assert provider.get_transcripts_for_region(testing.CHROM, 211, 211) == ['NM_000001.1']
    assert provider.get_transcripts_for_region(testing.CHROM, 212, 999) == []
    assert provider.get_transcripts_for_region(testing.CHROM, 1000, 1110) == ['NM_000002.1']
    assert provider.get_transcripts_for_region(testing.CHROM, 0, 3000) == ['NM_000001.1', 'NM_000002.1']


def test_refget_identities_both_ways(provider: provider_mod.BundleProvider, built: testing.Reference) -> None:
    genome_digest = refget.digest(testing.genome_sequence().encode())
    assert provider.get_refget_accession(testing.CHROM) == genome_digest
    assert provider.get_accession_for_refget(genome_digest) == testing.CHROM
    assert provider.get_refget_accession('NM_000002.1') == refget.digest(built.rev.encode())
    assert provider.get_accession_for_refget(refget.digest(built.rev.encode())) == 'NM_000002.1'
    assert provider.get_refget_accession('NP_000001.1') == refget.digest(b'MKLV')
    assert provider.get_refget_accession('NM_999999.9') is None
    assert provider.get_accession_for_refget('SQ.' + 'B' * 32) is None


def test_a_store_and_a_genome_on_different_assemblies_are_refused(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    fasta = tmp_path / 'other.fa'
    fasta.write_text(f'>{testing.CHROM}\nACGT\n', 'ascii')
    genome_build.build_genome(fasta, tmp_path / 'grch37', assembly=bundle_pb2.ASSEMBLY_GRCH37)
    with pytest.raises(ValueError, match='GRCH38, genome on ASSEMBLY_GRCH37'):
        provider_mod.BundleProvider(built.store, genome_mod.Genome(str(tmp_path / 'grch37')))


def test_a_negative_start_on_a_record_reads_from_its_first_base(
    provider: provider_mod.BundleProvider, built: testing.Reference
) -> None:
    assert provider.get_seq('NM_000001.1', -3, 2, weaver.IdentifierType.TranscriptAccession) == built.fwd[:2]


def test_a_genomic_accession_the_genome_lacks_is_weavers_error(provider: provider_mod.BundleProvider) -> None:
    with pytest.raises(weaver.DataProviderError, match=r'NC_000001\.11'):
        provider.get_seq('NC_000001.11', 0, 10, weaver.IdentifierType.GenomicAccession)


def test_a_digest_with_the_ga4gh_prefix_resolves(provider: provider_mod.BundleProvider) -> None:
    genome_digest = refget.digest(testing.genome_sequence().encode())
    assert provider.get_accession_for_refget(f'ga4gh:{genome_digest}') == testing.CHROM


@pytest.mark.parametrize('malformed', ['SQ.abc', 'd41d8cd98f00b204e9800998ecf8427e', 'SQ.' + 'é' * 32])
def test_a_malformed_digest_names_nothing(provider: provider_mod.BundleProvider, malformed: str) -> None:
    assert provider.get_accession_for_refget(malformed) is None


def test_a_non_ascii_accession_is_not_in_the_reference_data(provider: provider_mod.BundleProvider) -> None:
    with pytest.raises(weaver.DataProviderError, match='not in the reference data'):
        provider.get_transcript('NM_000001.\u0661', None)  # an Arabic-Indic digit one


@pytest.fixture
def stacked(built: testing.Reference, tmp_path: pathlib.Path) -> provider_mod.BundleProvider:
    """Two releases holding the same versions; the second moves FWD's CDS start to index 5."""
    older, _, _, _ = testing.bundles()
    newer, _, _, _ = testing.bundles()
    newer[0].transcripts[0].cds.start_index = 5
    source = tmp_path / 'annotation.gff3'
    source.write_text('synthetic\n', 'utf-8')
    shards = [
        store_build.write_shard(older, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)]),
        store_build.write_shard(newer, tmp_path / 'shards', release='RS_2', inputs=[('annotation', source)]),
    ]
    store_build.write_index(tmp_path / 'stacked', shards, assembly=testing.ASSEMBLY)
    return provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'stacked')), built.genome)


def test_the_newest_release_holding_a_version_answers(stacked: provider_mod.BundleProvider) -> None:
    assert stacked.get_transcript('NM_000001.1', None)['cds_start_index'] == 5


def test_a_transcript_in_two_releases_is_found_once(stacked: provider_mod.BundleProvider) -> None:
    assert stacked.get_transcripts_for_region(testing.CHROM, 150, 150) == ['NM_000001.1']
    assert stacked.get_symbol_accessions('FWD', 'gene', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.1')
    ]
    assert stacked.get_symbol_accessions('FWD', 'gene', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.1')
    ]


@pytest.mark.parametrize(('protein', 'transcript'), [('NP_000001.1', 'NM_000001.1'), ('NP_000002.1', 'NM_000002.1')])
def test_a_protein_answers_the_transcript_that_encodes_it(
    provider: provider_mod.BundleProvider, protein: str, transcript: str
) -> None:
    assert provider.get_symbol_accessions(protein, 'p', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, transcript)
    ]


@pytest.mark.parametrize('protein', ['NP_999999.1', 'NP_000001.2'])
def test_a_protein_version_the_store_does_not_hold_has_no_transcript(
    provider: provider_mod.BundleProvider, protein: str
) -> None:
    assert provider.get_symbol_accessions(protein, 'p', 'c') == []


def test_a_bare_protein_accession_has_no_transcript(provider: provider_mod.BundleProvider) -> None:
    assert provider.get_symbol_accessions('NP_000001', 'p', 'c') == []


@pytest.fixture
def two_versions(built: testing.Reference, tmp_path: pathlib.Path) -> provider_mod.BundleProvider:
    """FWD's bundle holds NM_000001.1 encoding NP_000001.1 and NM_000001.2 encoding NP_000001.2."""
    made, _, _, _ = testing.bundles()
    fwd = made[0]
    transcript = fwd.transcripts.add()
    transcript.CopyFrom(fwd.transcripts[0])
    transcript.version = 2
    transcript.protein_version = 2
    protein = fwd.proteins.add()
    protein.CopyFrom(fwd.proteins[0])
    protein.version = 2
    protein.transcript_version = 2
    source = tmp_path / 'annotation.gff3'
    source.write_text('synthetic\n', 'utf-8')
    shard = store_build.write_shard(made, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)])
    store_build.write_index(tmp_path / 'two', [shard], assembly=testing.ASSEMBLY)
    return provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'two')), built.genome)


@pytest.mark.parametrize(('protein', 'transcript'), [('NP_000001.1', 'NM_000001.1'), ('NP_000001.2', 'NM_000001.2')])
def test_a_protein_version_answers_only_the_transcript_encoding_that_version(
    two_versions: provider_mod.BundleProvider, protein: str, transcript: str
) -> None:
    assert two_versions.get_symbol_accessions(protein, 'p', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, transcript)
    ]


def test_weaver_back_translates_a_protein_variant_without_being_given_a_transcript(
    provider: provider_mod.BundleProvider,
) -> None:
    mapper = weaver.VariantMapper(provider)
    coding, _ = mapper.p_to_c(weaver.parse('NP_000001.1:p.Lys2Arg'))
    assert coding.ac == 'NM_000001.1'


GENBANK_RECORD = 'ATGAAACTGGTGTAAGG'  # a synthetic mRNA: CDS ATG AAA CTG GTG TAA at indices 0..14, two bases beyond
GENBANK_PROTEIN = 'MKLV'  # NP_000001.1's residues, so the two share a digest


def _genbank_bundle(gene: bundle_pb2.Gene, *, record: str = GENBANK_RECORD) -> bundle_pb2.GeneBundle:
    """The gene's U00001.1, unplaced, encoding AAA00001.1."""
    residues, protein = record.encode(), GENBANK_PROTEIN.encode()
    return bundle_pb2.GeneBundle(
        gene=gene,
        transcripts=[
            bundle_pb2.Transcript(
                accession='U00001',
                version=1,
                publisher=bundle_pb2.PUBLISHER_GENBANK,
                status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
                biotype='mRNA',
                sequence_digest=refget.digest(residues),
                cds=bundle_pb2.Cds(start_index=0, end_index_inclusive=14),
                protein_accession='AAA00001',
                protein_version=1,
            )
        ],
        proteins=[
            bundle_pb2.Protein(
                accession='AAA00001',
                version=1,
                sequence_digest=refget.digest(protein),
                transcript_accession='U00001',
                transcript_version=1,
            )
        ],
        sequences=[
            bundle_pb2.Sequence(
                digest=refget.digest(residues), alphabet=bundle_pb2.ALPHABET_NUCLEOTIDE, residues=residues
            ),
            bundle_pb2.Sequence(digest=refget.digest(protein), alphabet=bundle_pb2.ALPHABET_PROTEIN, residues=protein),
        ],
    )


def _genbank_store(
    tmp_path: pathlib.Path, *, genbank_first: bool, record: str = GENBANK_RECORD
) -> store_mod.BundleStore:
    made, _, _, _ = testing.bundles()
    source = tmp_path / 'records.gb'
    source.write_text('synthetic\n', 'utf-8')
    refseq = store_build.write_shard(made, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)])
    genbank = store_build.write_shard(
        [_genbank_bundle(made[0].gene, record=record)],
        tmp_path / 'shards',
        release='GB_1',
        inputs=[('records', source)],
    )
    shards = [genbank, refseq] if genbank_first else [refseq, genbank]
    store_build.write_index(tmp_path / 'with-genbank', shards, assembly=testing.ASSEMBLY)
    return store_mod.BundleStore(str(tmp_path / 'with-genbank'))


@pytest.fixture
def genbank_store(tmp_path: pathlib.Path) -> store_mod.BundleStore:
    """The synthetic store with a GenBank shard stacked on it: FWD's U00001.1, unplaced, encoding AAA00001.1."""
    return _genbank_store(tmp_path, genbank_first=False)


@pytest.fixture
def with_genbank(built: testing.Reference, genbank_store: store_mod.BundleStore) -> provider_mod.BundleProvider:
    return provider_mod.BundleProvider(genbank_store, built.genome)


@pytest.mark.parametrize(
    ('identifier', 'kind'),
    [
        ('U00001.1', weaver.IdentifierType.TranscriptAccession),
        ('U00001', weaver.IdentifierType.TranscriptAccession),
        ('AAA00001.1', weaver.IdentifierType.ProteinAccession),
    ],
)
def test_an_accession_the_store_holds_is_the_kind_it_is_held_as(
    with_genbank: provider_mod.BundleProvider, identifier: str, kind: weaver.IdentifierType
) -> None:
    assert with_genbank.get_identifier_type(identifier) == kind


@pytest.mark.parametrize(
    ('identifier', 'kind'),
    [
        ('NM_999999.9', weaver.IdentifierType.TranscriptAccession),
        ('NP_999999.9', weaver.IdentifierType.ProteinAccession),
        ('U99999.1', weaver.IdentifierType.Unknown),
    ],
)
def test_an_accession_the_store_does_not_hold_is_the_kind_its_namespace_says(
    with_genbank: provider_mod.BundleProvider, identifier: str, kind: weaver.IdentifierType
) -> None:
    assert with_genbank.get_identifier_type(identifier) == kind


def test_a_protein_with_no_reference_prefix_serves_its_residues(with_genbank: provider_mod.BundleProvider) -> None:
    assert with_genbank.get_seq('AAA00001.1', 0, None, weaver.IdentifierType.ProteinAccession) == GENBANK_PROTEIN


def test_a_protein_with_no_reference_prefix_has_its_digest(with_genbank: provider_mod.BundleProvider) -> None:
    assert with_genbank.get_refget_accession('AAA00001.1') == refget.digest(GENBANK_PROTEIN.encode())


@pytest.mark.parametrize('protein', ['NP_999999', 'NP_999999.1'])
def test_a_protein_the_store_does_not_hold_is_not_in_the_reference_data(
    provider: provider_mod.BundleProvider, protein: str
) -> None:
    with pytest.raises(weaver.DataProviderError, match='not in the reference data'):
        provider.get_seq(protein, 0, None, weaver.IdentifierType.ProteinAccession)


def test_an_unplaced_transcript_is_given_without_a_placement_whatever_genome_is_asked_for(
    with_genbank: provider_mod.BundleProvider,
) -> None:
    model = with_genbank.get_transcript('U00001.1', testing.CHROM)
    assert {'strand', 'reference_accession', 'exons'}.isdisjoint(model)
    assert (model['cds_start_index'], model['cds_end_index']) == (0, 14)


def test_a_placed_transcript_asked_for_on_another_sequence_is_refused(provider: provider_mod.BundleProvider) -> None:
    with pytest.raises(weaver.DataProviderError, match=r'NM_000001\.1: no alignment on NC_000098\.1'):
        provider.get_transcript('NM_000001.1', 'NC_000098.1')


def test_weaver_checks_a_name_on_an_unplaced_record_against_its_own_bases(
    with_genbank: provider_mod.BundleProvider,
) -> None:
    # c.4 is the record's index 3, the first base of codon AAA
    assert weaver.parse('U00001.1:c.4A>G').validate(with_genbank)
    assert not weaver.parse('U00001.1:c.4C>G').validate(with_genbank)


def test_weaver_predicts_an_unplaced_records_protein_from_its_own_cds(
    with_genbank: provider_mod.BundleProvider,
) -> None:
    mapper = weaver.VariantMapper(with_genbank)
    # AAA (Lys2) to GAA (Glu)
    assert mapper.c_to_p(weaver.parse('U00001.1:c.4A>G')).format() == 'AAA00001.1:p.(Lys2Glu)'


def test_weaver_refuses_to_place_an_unplaced_record(with_genbank: provider_mod.BundleProvider) -> None:
    mapper = weaver.VariantMapper(with_genbank)
    with pytest.raises(weaver.UnsupportedOperationError, match=r'U00001\.1'):
        mapper.c_to_g(weaver.parse('U00001.1:c.4A>G'), testing.CHROM)


def test_a_refseq_transcript_bundled_without_an_alignment_is_given_unplaced(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    made, _, _, _ = testing.bundles()
    del made[1].transcripts[0].alignments[:]  # REV's NM_000002.1, as NCBI published no alignment for it
    source = tmp_path / 'annotation.gff3'
    source.write_text('synthetic\n', 'utf-8')
    shard = store_build.write_shard(made, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)])
    store_build.write_index(tmp_path / 'unaligned', [shard], assembly=testing.ASSEMBLY)
    provider = provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'unaligned')), built.genome)
    assert 'exons' not in provider.get_transcript('NM_000002.1', testing.CHROM)


def test_a_genes_transcripts_leave_out_a_genbank_record(with_genbank: provider_mod.BundleProvider) -> None:
    assert with_genbank.get_symbol_accessions('FWD', 'gene', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.1')
    ]


def test_a_genes_proteins_leave_out_a_genbank_records_protein(with_genbank: provider_mod.BundleProvider) -> None:
    assert with_genbank.get_symbol_accessions('FWD', 'gene', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.1')
    ]


def test_a_provider_made_to_list_genbank_lists_its_transcripts_for_a_gene(
    built: testing.Reference, genbank_store: store_mod.BundleStore
) -> None:
    publishers = {bundle_pb2.PUBLISHER_REFSEQ, bundle_pb2.PUBLISHER_GENBANK}
    provider = provider_mod.BundleProvider(genbank_store, built.genome, listed_publishers=publishers)
    assert provider.get_symbol_accessions('FWD', 'gene', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.1'),
        (weaver.IdentifierType.TranscriptAccession, 'U00001.1'),
    ]


def test_a_provider_made_to_list_genbank_lists_its_proteins_for_a_gene(
    built: testing.Reference, genbank_store: store_mod.BundleStore
) -> None:
    publishers = {bundle_pb2.PUBLISHER_REFSEQ, bundle_pb2.PUBLISHER_GENBANK}
    provider = provider_mod.BundleProvider(genbank_store, built.genome, listed_publishers=publishers)
    assert provider.get_symbol_accessions('FWD', 'gene', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.1'),
        (weaver.IdentifierType.ProteinAccession, 'AAA00001.1'),
    ]


def test_a_provider_made_to_list_only_ensembl_leaves_out_refseq(
    built: testing.Reference, genbank_store: store_mod.BundleStore
) -> None:
    provider = provider_mod.BundleProvider(
        genbank_store, built.genome, listed_publishers={bundle_pb2.PUBLISHER_ENSEMBL}
    )
    assert provider.get_symbol_accessions('FWD', 'gene', 'c') == []


def _with_retired(
    built: testing.Reference, tmp_path: pathlib.Path, *, include_retired: bool
) -> provider_mod.BundleProvider:
    """FWD's bundle holds NM_000001.2 encoding NP_000001.2, and NM_000001.1, superseded, encoding NP_000001.1."""
    made, _, _, _ = testing.bundles()
    fwd = made[0]
    retired = fwd.transcripts[0]
    retired.status = bundle_pb2.TRANSCRIPT_STATUS_SUPERSEDED
    del retired.tags[:]
    current = fwd.transcripts.add()
    current.CopyFrom(retired)
    current.version = 2
    current.status = bundle_pb2.TRANSCRIPT_STATUS_CURRENT
    current.protein_version = 2
    protein = fwd.proteins.add()
    protein.CopyFrom(fwd.proteins[0])
    protein.version = 2
    protein.transcript_version = 2
    source = tmp_path / 'annotation.gff3'
    source.write_text('synthetic\n', 'utf-8')
    shard = store_build.write_shard(made, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)])
    store_build.write_index(tmp_path / 'retired', [shard], assembly=testing.ASSEMBLY)
    store = store_mod.BundleStore(str(tmp_path / 'retired'))
    return provider_mod.BundleProvider(store, built.genome, include_retired=include_retired)


def test_a_genes_transcripts_leave_out_a_retired_version(built: testing.Reference, tmp_path: pathlib.Path) -> None:
    provider = _with_retired(built, tmp_path, include_retired=False)
    assert provider.get_symbol_accessions('FWD', 'gene', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.2')
    ]


def test_a_genes_proteins_leave_out_a_retired_versions_protein(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    provider = _with_retired(built, tmp_path, include_retired=False)
    assert provider.get_symbol_accessions('FWD', 'gene', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.2')
    ]


def test_a_retired_version_is_still_served_when_named(built: testing.Reference, tmp_path: pathlib.Path) -> None:
    provider = _with_retired(built, tmp_path, include_retired=False)
    assert provider.get_transcript('NM_000001.1', None)['ac'] == 'NM_000001.1'


def test_a_provider_made_to_include_retired_versions_lists_them_for_a_gene(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    provider = _with_retired(built, tmp_path, include_retired=True)
    assert provider.get_symbol_accessions('FWD', 'gene', 'c') == [
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.1'),
        (weaver.IdentifierType.TranscriptAccession, 'NM_000001.2'),
    ]


def test_a_provider_made_to_include_retired_versions_lists_them_over_a_region(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    provider = _with_retired(built, tmp_path, include_retired=True)
    assert provider.get_transcripts_for_region(testing.CHROM, 150, 150) == ['NM_000001.1', 'NM_000001.2']


@pytest.mark.parametrize('genbank_first', [False, True])
def test_a_digest_two_publishers_share_names_the_reference_publishers_accession(
    built: testing.Reference, tmp_path: pathlib.Path, genbank_first: bool
) -> None:
    store = _genbank_store(tmp_path, genbank_first=genbank_first)
    provider = provider_mod.BundleProvider(store, built.genome)
    assert provider.get_accession_for_refget(refget.digest(GENBANK_PROTEIN.encode())) == 'NP_000001.1'


def test_a_digest_only_genbank_holds_names_the_genbank_accession(with_genbank: provider_mod.BundleProvider) -> None:
    assert with_genbank.get_accession_for_refget(refget.digest(GENBANK_RECORD.encode())) == 'U00001.1'


def test_a_provider_made_to_include_retired_versions_lists_their_proteins_for_a_gene(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    provider = _with_retired(built, tmp_path, include_retired=True)
    assert provider.get_symbol_accessions('FWD', 'gene', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.1'),
        (weaver.IdentifierType.ProteinAccession, 'NP_000001.2'),
    ]


def test_a_digest_a_genbank_mrna_shares_with_a_refseq_transcript_names_the_refseq_accession(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    store = _genbank_store(tmp_path, genbank_first=True, record=built.fwd)  # U00001.1 with NM_000001.1's bases
    provider = provider_mod.BundleProvider(store, built.genome)
    assert provider.get_accession_for_refget(refget.digest(built.fwd.encode())) == 'NM_000001.1'


def test_a_provider_made_to_list_only_ensembl_lists_no_refseq_transcript_over_a_region(
    built: testing.Reference,
) -> None:
    provider = provider_mod.BundleProvider(built.store, built.genome, listed_publishers={bundle_pb2.PUBLISHER_ENSEMBL})
    assert provider.get_transcripts_for_region(testing.CHROM, 0, 3000) == []


@pytest.mark.parametrize('publishers', [set(), {bundle_pb2.PUBLISHER_UNSPECIFIED}])
def test_a_provider_listing_no_defined_publisher_is_refused(
    built: testing.Reference, publishers: set[bundle_pb2.Publisher]
) -> None:
    with pytest.raises(ValueError, match='listed_publishers'):
        provider_mod.BundleProvider(built.store, built.genome, listed_publishers=publishers)


def test_a_protein_with_no_version_serves_its_residues_by_its_accession(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    store_build.write_index(
        tmp_path / 'uniprot', [_isoform_release(tmp_path, 'UP_1', b'MKV')], assembly=testing.ASSEMBLY
    )
    provider = provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'uniprot')), built.genome)
    assert provider.get_seq('P00001-2', 0, None, weaver.IdentifierType.ProteinAccession) == 'MKV'


def _isoform_release(tmp_path: pathlib.Path, release: str, residues: bytes) -> pathlib.Path:
    """A UniProt-shaped shard holding FWD's isoform P00001-2 with these residues."""
    digest = refget.digest(residues)
    bundle = bundle_pb2.GeneBundle(
        gene=bundle_pb2.Gene(symbol='FWD'),
        proteins=[bundle_pb2.Protein(accession='P00001-2', version=0, sequence_digest=digest)],
        sequences=[bundle_pb2.Sequence(digest=digest, alphabet=bundle_pb2.ALPHABET_PROTEIN, residues=residues)],
    )
    source = tmp_path / f'{release}.dat'
    source.write_text('synthetic\n', 'utf-8')
    return store_build.write_shard([bundle], tmp_path / 'shards', release=release, inputs=[('entries', source)])


@pytest.fixture
def reused_isoform(built: testing.Reference, tmp_path: pathlib.Path) -> provider_mod.BundleProvider:
    """Two UniProt releases stacked, the later giving the isoform ID P00001-2 other residues, MKNEW for MKOLD."""
    shards = [_isoform_release(tmp_path, 'UP_1', b'MKOLD'), _isoform_release(tmp_path, 'UP_2', b'MKNEW')]
    store_build.write_index(tmp_path / 'reused', shards, assembly=testing.ASSEMBLY)
    return provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'reused')), built.genome)


def test_a_reused_isoform_id_serves_the_later_releases_residues(reused_isoform: provider_mod.BundleProvider) -> None:
    assert reused_isoform.get_seq('P00001-2', 0, None, weaver.IdentifierType.ProteinAccession) == 'MKNEW'


def test_a_reused_isoform_ids_digest_is_the_later_releases(reused_isoform: provider_mod.BundleProvider) -> None:
    assert reused_isoform.get_refget_accession('P00001-2') == refget.digest(b'MKNEW')


def test_a_digest_a_reused_name_no_longer_carries_names_nothing(reused_isoform: provider_mod.BundleProvider) -> None:
    assert reused_isoform.get_accession_for_refget(refget.digest(b'MKOLD')) is None


def test_a_uniprot_name_the_store_does_not_hold_is_a_protein_not_in_the_reference_data(
    reused_isoform: provider_mod.BundleProvider,
) -> None:
    with pytest.raises(weaver.DataProviderError, match=r'P21181-4: not in the reference data'):
        reused_isoform.get_seq('P21181-4', 0, None, weaver.IdentifierType.ProteinAccession)


def test_a_gene_symbol_in_the_shape_of_a_uniprot_accession_is_a_gene(
    built: testing.Reference, tmp_path: pathlib.Path
) -> None:
    made, _, _, _ = testing.bundles()
    made[0].gene.symbol = 'H2AC11'  # an HGNC symbol UniProt's accession format also fits
    source = tmp_path / 'annotation.gff3'
    source.write_text('synthetic\n', 'utf-8')
    shard = store_build.write_shard(made, tmp_path / 'shards', release='RS_1', inputs=[('annotation', source)])
    store_build.write_index(tmp_path / 'h2', [shard], assembly=testing.ASSEMBLY)
    provider = provider_mod.BundleProvider(store_mod.BundleStore(str(tmp_path / 'h2')), built.genome)
    assert provider.get_identifier_type('H2AC11') == weaver.IdentifierType.GeneSymbol
