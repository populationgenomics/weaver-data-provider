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
    assert model['strand'] == 1
    assert model['reference_accession'] == testing.CHROM
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
