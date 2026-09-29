"""The RefSeq builder, from files in NCBI's formats to bundles, and through the command to weaver.

The inputs are written here in the publisher's formats — GFF3, FASTA and TSV as text, the alignments
as a BAM that pysam writes — so no fixture comes from the builder itself. The synthetic release:

- PLUS, NM_000010.2 on NC_000099.1's plus strand: exons 101-120 and 201-230, CDS 105-210.
- MINUS, NM_000020.1 on the minus strand: exons 501-520 and 601-630, CDS 505-605; its record carries
  one base the genome lacks, after genome 610, so the alignment has an insertion.
- PAR, NM_000040.1 at 51-100 on both X and Y, as a pseudoautosomal transcript is.
- NONC, NR_000030.1, a non-coding transcript NCBI published no alignment for.
- MT-TF, a mitochondrial tRNA with no transcript accession; ALT, a gene on an alternate locus only.
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib
import random

import pysam.libcalignedsegment
import pysam.libcalignmentfile
import pysam.samtools
import pytest
import weaver

from weaver_data_provider import build, refget
from weaver_data_provider import genome as genome_mod
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import cli, refseq
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2

CHROM, X, Y = 'NC_000099.1', 'NC_000023.11', 'NC_000024.10'
_rng = random.Random(99)
GENOME = {CHROM: ''.join(_rng.choices('ACGT', k=1000)), X: ''.join(_rng.choices('ACGT', k=300))}
GENOME[Y] = GENOME[X]  # the pseudoautosomal region is the same sequence on both
INSERTED = 'A'  # the MINUS record's base the genome lacks, in genome orientation


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def _g(chrom: str, start: int, end: int) -> str:
    """Genome bases start..end, 1-based closed."""
    return GENOME[chrom][start - 1 : end]


_MINUS_GENOME_ORDER = _g(CHROM, 501, 520) + _g(CHROM, 601, 610) + INSERTED + _g(CHROM, 611, 630)
TRANSCRIPTS = {
    'NM_000010.2': _g(CHROM, 101, 120) + _g(CHROM, 201, 230),
    'NM_000020.1': _revcomp(_MINUS_GENOME_ORDER),
    'NR_000030.1': _g(CHROM, 801, 850),
    'NM_000040.1': _g(X, 51, 100),
    'NM_000050.1': 'ACGTACGTAC',
}
PROTEINS = {'NP_000010.1': 'MSEQPLUS', 'NP_000020.1': 'MSEQMINUS'}


def _row(seqid: str, kind: str, start: int, end: int, strand: str, **attrs: str) -> str:
    return '\t'.join(
        [seqid, 'RefSeq', kind, str(start), str(end), '.', strand, '.', ';'.join(f'{k}={v}' for k, v in attrs.items())]
    )


GFF = [
    '##gff-version 3',
    _row(
        CHROM,
        'gene',
        101,
        230,
        '+',
        ID='gene-PLUS',
        Dbxref='GeneID:10,HGNC:HGNC:10',
        Name='PLUS',
        description='plus%20gene',
        gene_synonym='PLS',
    ),
    _row(
        CHROM,
        'mRNA',
        101,
        230,
        '+',
        ID='rna-NM_000010.2',
        Parent='gene-PLUS',
        tag='RefSeq Select',
        transcript_id='NM_000010.2',
    ),
    _row(CHROM, 'exon', 101, 120, '+', ID='exon-NM_000010.2-1', Parent='rna-NM_000010.2'),
    _row(CHROM, 'exon', 201, 230, '+', ID='exon-NM_000010.2-2', Parent='rna-NM_000010.2'),
    _row(CHROM, 'CDS', 105, 120, '+', ID='cds-NP_000010.1', Parent='rna-NM_000010.2', protein_id='NP_000010.1'),
    _row(CHROM, 'CDS', 201, 210, '+', ID='cds-NP_000010.1', Parent='rna-NM_000010.2', protein_id='NP_000010.1'),
    _row(
        CHROM,
        'gene',
        501,
        630,
        '-',
        ID='gene-MINUS-OLD',
        Dbxref='GeneID:20',
        Name='MINUS-OLD',
        description='old%20name',
    ),
    _row(CHROM, 'mRNA', 501, 630, '-', ID='rna-NM_000020.1', Parent='gene-MINUS-OLD', transcript_id='NM_000020.1'),
    _row(CHROM, 'CDS', 601, 605, '-', ID='cds-NP_000020.1', Parent='rna-NM_000020.1', protein_id='NP_000020.1'),
    _row(CHROM, 'CDS', 505, 520, '-', ID='cds-NP_000020.1', Parent='rna-NM_000020.1', protein_id='NP_000020.1'),
    _row(
        CHROM,
        'gene',
        801,
        850,
        '+',
        ID='gene-NONC',
        Dbxref='GeneID:30',
        Name='NONC',
        description='a%20non-coding%20gene',
    ),
    _row(CHROM, 'lnc_RNA', 801, 850, '+', ID='rna-NR_000030.1', Parent='gene-NONC', transcript_id='NR_000030.1'),
    _row(X, 'gene', 51, 100, '+', ID='gene-PAR', Dbxref='GeneID:40', Name='PAR'),
    _row(X, 'lnc_RNA', 51, 100, '+', ID='rna-NM_000040.1', Parent='gene-PAR', transcript_id='NM_000040.1'),
    _row(Y, 'gene', 51, 100, '+', ID='gene-PAR-2', Dbxref='GeneID:40', Name='PAR'),
    _row(Y, 'lnc_RNA', 51, 100, '+', ID='rna-NM_000040.1-2', Parent='gene-PAR-2', transcript_id='NM_000040.1'),
    _row('NC_012920.1', 'gene', 1, 70, '+', ID='gene-MT-TF', Dbxref='GeneID:50', Name='MT-TF'),
    _row('NC_012920.1', 'tRNA', 1, 70, '+', ID='rna-MT-TF', Parent='gene-MT-TF'),
    _row('NT_187361.1', 'gene', 1, 10, '+', ID='gene-ALT', Dbxref='GeneID:60', Name='ALT'),
    _row('NT_187361.1', 'mRNA', 1, 10, '+', ID='rna-NM_000050.1', Parent='gene-ALT', transcript_id='NM_000050.1'),
]


@dataclasses.dataclass(frozen=True)
class Read:
    name: str
    chrom: str
    start: int  # 0-based leftmost genome position, as in a BAM
    cigar: str
    reverse: bool
    sequence: str  # genome orientation, as a BAM stores it
    secondary: bool = False


READS = [
    Read('NM_000010.2', CHROM, 100, '20=80N30=', False, TRANSCRIPTS['NM_000010.2']),
    Read('NM_000020.1', CHROM, 500, '20=80N10=1I20=', True, _MINUS_GENOME_ORDER),
    Read('NM_000040.1', X, 50, '50=', False, TRANSCRIPTS['NM_000040.1']),
    Read('NM_000040.1', Y, 50, '50=', False, TRANSCRIPTS['NM_000040.1']),
]
HGNC_COLUMNS = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']
HGNC = [
    ['HGNC:10', 'PLUS', 'plus gene, HGNC', 'PL1', 'OLDPLUS', '10', 'ENSG00000000010'],
    ['HGNC:20', 'MINUS', 'minus gene, HGNC', '', 'MINUS-OLD', '20', ''],
    ['HGNC:40', 'PAR', 'pseudoautosomal gene', '', '', '40', 'ENSG00000000040'],
]
MANE = [['NM_000010.2', 'ENST00000000010.1', 'MANE Select']]


def _fasta(path: pathlib.Path, records: dict[str, str]) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        for name, residues in records.items():
            fh.write(f'>{name} synthetic\n')
            for k in range(0, len(residues), 60):
                fh.write(residues[k : k + 60] + '\n')
    return path


def _bam(path: pathlib.Path, reads: list[Read]) -> pathlib.Path:
    order = [CHROM, X, Y]
    header = pysam.libcalignmentfile.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': len(GENOME[n])} for n in order]}
    )
    with pysam.libcalignmentfile.AlignmentFile(str(path), 'wb', header=header) as out:
        for read in sorted(reads, key=lambda r: (order.index(r.chrom), r.start)):
            segment = pysam.libcalignedsegment.AlignedSegment(header)
            segment.query_name = read.name
            segment.reference_name = read.chrom
            segment.reference_start = read.start
            segment.flag = (16 if read.reverse else 0) | (256 if read.secondary else 0)
            segment.mapping_quality = 60
            segment.cigarstring = read.cigar or None
            segment.query_sequence = read.sequence
            out.write(segment)
    pysam.samtools.index(str(path))
    return path


def _tsv(path: pathlib.Path, columns: list[str], rows: list[list[str]], *, gzipped: bool) -> pathlib.Path:
    text = '\n'.join('\t'.join(r) for r in [columns, *rows]) + '\n'
    if gzipped:
        with gzip.open(path, 'wt', encoding='utf-8') as fh:
            fh.write(text)
    else:
        path.write_text(text, 'utf-8')
    return path


def _release(
    tmp_path: pathlib.Path,
    *,
    gff: list[str] = GFF,
    reads: list[Read] = READS,
    transcripts: dict[str, str] = TRANSCRIPTS,
    proteins: dict[str, str] = PROTEINS,
    hgnc_columns: list[str] = HGNC_COLUMNS,
) -> refseq.Release:
    """The synthetic release written under `tmp_path`; each keyword replaces one input."""
    annotation = tmp_path / 'genomic.gff.gz'
    with gzip.open(annotation, 'wt', encoding='utf-8') as fh:
        fh.write('\n'.join(gff) + '\n')
    return refseq.Release(
        assembly=bundle_pb2.ASSEMBLY_GRCH38,
        release='RS_TEST',
        annotation=annotation,
        transcripts=_fasta(tmp_path / 'rna.fna.gz', transcripts),
        proteins=_fasta(tmp_path / 'protein.faa.gz', proteins),
        alignments=(_bam(tmp_path / 'alns.bam', reads),),
        hgnc=_tsv(tmp_path / 'hgnc.txt', hgnc_columns, [r[: len(hgnc_columns)] for r in HGNC], gzipped=False),
        mane=_tsv(tmp_path / 'mane.txt.gz', ['RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'], MANE, gzipped=True),
    )


def _by_symbol(release: refseq.Release) -> dict[str, bundle_pb2.GeneBundle]:
    return {b.gene.symbol: b for b in refseq.bundles(release)}


def _exons(alignment: bundle_pb2.Alignment) -> list[tuple[int, int, int, int, str]]:
    return [
        (e.transcript_start, e.transcript_end, e.genome_start, e.genome_end_inclusive, e.cigar) for e in alignment.exons
    ]


def test_a_plus_strand_transcript_takes_its_exons_from_the_alignment(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PLUS'].transcripts
    (alignment,) = transcript.alignments
    assert (alignment.chromosome, alignment.strand) == (CHROM, bundle_pb2.STRAND_PLUS)
    assert _exons(alignment) == [(0, 20, 100, 119, '20='), (20, 50, 200, 229, '30=')]


def test_a_plus_strand_cds_is_projected_onto_the_record(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PLUS'].transcripts
    # genome 105 is exon 1's fifth base; genome 210 is exon 2's tenth, after exon 1's 20
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 29)


def test_a_minus_strand_transcript_runs_5_to_3_with_its_insertion_in_transcript_order(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['MINUS'].transcripts
    (alignment,) = transcript.alignments
    assert alignment.strand == bundle_pb2.STRAND_MINUS
    assert _exons(alignment) == [(0, 31, 600, 629, '20=1I10='), (31, 51, 500, 519, '20=')]


def test_a_minus_strand_cds_is_projected_through_the_insertion(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['MINUS'].transcripts
    # genome 605 is 25 bases in from the exon's 5' end, genome 630; the inserted base adds one
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (26, 46)


def test_the_bundled_sequence_is_the_record_not_the_genome_splice(tmp_path: pathlib.Path) -> None:
    bundle = _by_symbol(_release(tmp_path))['MINUS']
    (transcript,) = bundle.transcripts
    residues = {s.digest: s.residues for s in bundle.sequences}
    assert residues[transcript.sequence_digest] == TRANSCRIPTS['NM_000020.1'].encode()
    assert transcript.sequence_digest == refget.digest(TRANSCRIPTS['NM_000020.1'].encode())


def test_a_protein_is_bundled_with_its_transcript(tmp_path: pathlib.Path) -> None:
    bundle = _by_symbol(_release(tmp_path))['PLUS']
    (protein,) = bundle.proteins
    assert (protein.accession, protein.version, protein.transcript_accession) == ('NP_000010', 1, 'NM_000010')
    assert {s.digest: s.residues for s in bundle.sequences}[protein.sequence_digest] == b'MSEQPLUS'


def test_a_pseudoautosomal_transcript_has_an_alignment_on_each_chromosome(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PAR'].transcripts
    assert sorted(a.chromosome for a in transcript.alignments) == [X, Y]


def test_a_noncoding_transcript_without_an_alignment_is_bundled_unplaced(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['NONC'].transcripts
    assert (transcript.accession, transcript.biotype) == ('NR_000030', 'lnc_RNA')
    assert not transcript.alignments
    assert not transcript.HasField('cds')


def test_only_transcripts_with_an_accession_on_a_chromosome_are_bundled(tmp_path: pathlib.Path) -> None:
    # MT-TF's tRNA names no accession; ALT is placed on an alternate locus only
    assert set(_by_symbol(_release(tmp_path))) == {'MINUS', 'NONC', 'PAR', 'PLUS'}


def test_bundles_come_in_symbol_order(tmp_path: pathlib.Path) -> None:
    symbols = [b.gene.symbol for b in refseq.bundles(_release(tmp_path))]
    assert symbols == sorted(symbols)


def test_hgnc_names_a_gene_it_has(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path))['MINUS'].gene
    assert (gene.symbol, gene.hgnc_id, gene.name) == ('MINUS', 'HGNC:20', 'minus gene, HGNC')
    assert list(gene.previous_symbols) == ['MINUS-OLD']


def test_the_annotation_names_a_gene_hgnc_lacks(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path))['NONC'].gene
    assert (gene.symbol, gene.hgnc_id, gene.ncbi_gene_id, gene.name) == ('NONC', '', '30', 'a non-coding gene')


def test_aliases_gather_hgnc_and_the_annotation(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path))['PLUS'].gene
    assert sorted(gene.alias_symbols) == ['PL1', 'PLS']
    assert list(gene.previous_symbols) == ['OLDPLUS']


def test_mane_select_is_tagged_with_its_ensembl_partner(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PLUS'].transcripts
    assert transcript.mane_partner == 'ENST00000000010.1'
    assert set(transcript.tags) == {bundle_pb2.TAG_REFSEQ_SELECT, bundle_pb2.TAG_MANE_SELECT}


def test_a_coding_transcript_without_an_alignment_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, reads=[r for r in READS if r.name != 'NM_000020.1'])
    with pytest.raises(build.BuildError, match=r'no alignment.*NM_000020\.1'):
        list(refseq.bundles(release))


def test_an_alignment_that_does_not_tile_the_record_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, transcripts={**TRANSCRIPTS, 'NM_000010.2': TRANSCRIPTS['NM_000010.2'] + 'A'})
    with pytest.raises(build.BuildError, match=r'NM_000010\.2.*covers transcript bases 0\.\.50 of a 51-base record'):
        list(refseq.bundles(release))


def test_an_alignment_outside_the_annotation_span_fails_the_build(tmp_path: pathlib.Path) -> None:
    shortened = _row(
        CHROM, 'mRNA', 101, 225, '+', ID='rna-NM_000010.2', Parent='gene-PLUS', transcript_id='NM_000010.2'
    )
    release = _release(tmp_path, gff=[shortened if 'ID=rna-NM_000010.2;' in line else line for line in GFF])
    with pytest.raises(build.BuildError, match=r'NM_000010\.2.*outside the annotation span 101\.\.225'):
        list(refseq.bundles(release))


def test_a_cds_bound_on_a_genome_only_base_fails_the_build(tmp_path: pathlib.Path) -> None:
    # the record lacks genome 105, the CDS's first base, so the alignment deletes it
    record = _g(CHROM, 101, 104) + _g(CHROM, 106, 120) + _g(CHROM, 201, 230)
    deleting = Read('NM_000010.2', CHROM, 100, '4=1D15=80N30=', False, record)
    release = _release(
        tmp_path,
        reads=[deleting if r.name == 'NM_000010.2' else r for r in READS],
        transcripts={**TRANSCRIPTS, 'NM_000010.2': record},
    )
    with pytest.raises(build.BuildError, match='genome position 104 is a genome-only base'):
        list(refseq.bundles(release))


def test_a_primary_alignment_with_no_cigar_fails_the_build(tmp_path: pathlib.Path) -> None:
    bare = Read('NM_000010.2', CHROM, 100, '', False, TRANSCRIPTS['NM_000010.2'])
    release = _release(tmp_path, reads=[bare if r.name == 'NM_000010.2' else r for r in READS])
    with pytest.raises(build.BuildError, match='no cigar'):
        list(refseq.bundles(release))


def test_a_transcript_missing_from_the_sequence_set_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, transcripts={k: v for k, v in TRANSCRIPTS.items() if k != 'NR_000030.1'})
    with pytest.raises(build.BuildError, match=r'no sequence.*NR_000030\.1'):
        list(refseq.bundles(release))


def test_a_protein_missing_from_the_protein_set_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, proteins={'NP_000020.1': PROTEINS['NP_000020.1']})
    with pytest.raises(build.BuildError, match=r'NP_000010\.1 is not in the protein set'):
        list(refseq.bundles(release))


def test_an_hgnc_table_without_a_needed_column_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, hgnc_columns=[c for c in HGNC_COLUMNS if c != 'ensembl_gene_id'])
    with pytest.raises(build.BuildError, match='ensembl_gene_id'):
        list(refseq.bundles(release))


def test_the_built_commands_feed_weaver(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> None:
    release = _release(tmp_path)
    fasta = _fasta(tmp_path / 'genome.fna.gz', GENOME)
    cli.main([
        'refseq', '--assembly', 'GRCh38', '--release', 'RS_TEST',
        '--annotation', str(release.annotation), '--transcripts', str(release.transcripts),
        '--proteins', str(release.proteins), '--alignments', str(release.alignments[0]),
        '--hgnc', str(release.hgnc), '--mane', str(release.mane), '--shards', str(tmp_path / 'shards'),
    ])  # fmt: skip
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    cli.main(['genome', '--assembly', 'GRCh38', '--fasta', str(fasta), '--out', str(tmp_path / 'genome')])
    provider = provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(tmp_path / 'genome'))
    )
    mapper = weaver.VariantMapper(provider)
    # MINUS's c.1 is genome 605 on the minus strand, so the bases are complemented
    genome_base = _g(CHROM, 605, 605)
    coding_base = _revcomp(genome_base)
    alt = 'G' if coding_base != 'G' else 'C'
    variant = weaver.parse(f'NM_000020.1:c.1{coding_base}>{alt}')
    assert variant.validate(provider)
    assert mapper.c_to_g(variant, CHROM).format() == f'{CHROM}:g.605{genome_base}>{_revcomp(alt)}'


def test_an_alignment_below_a_minus_strand_annotation_span_fails_the_build(tmp_path: pathlib.Path) -> None:
    # MINUS's 3' exon, 501-520, is its lowest on the genome and falls below a span narrowed to 505..630
    narrowed = _row(
        CHROM, 'mRNA', 505, 630, '-', ID='rna-NM_000020.1', Parent='gene-MINUS-OLD', transcript_id='NM_000020.1'
    )
    release = _release(tmp_path, gff=[narrowed if 'ID=rna-NM_000020.1;' in line else line for line in GFF])
    with pytest.raises(build.BuildError, match=r'NM_000020\.1.*outside the annotation span 505\.\.630'):
        list(refseq.bundles(release))


def test_a_secondary_alignment_is_not_the_placement(tmp_path: pathlib.Path) -> None:
    # it sorts ahead of the primary, so a builder taking the first read seen would place PLUS at 11..60
    secondary = Read('NM_000010.2', CHROM, 10, '50=', False, TRANSCRIPTS['NM_000010.2'], secondary=True)
    (transcript,) = _by_symbol(_release(tmp_path, reads=[secondary, *READS]))['PLUS'].transcripts
    (alignment,) = transcript.alignments
    assert _exons(alignment) == [(0, 20, 100, 119, '20='), (20, 50, 200, 229, '30=')]


def test_a_truncated_feature_line_fails_the_build(tmp_path: pathlib.Path) -> None:
    truncated = [line.rsplit('\t', 1)[0] if 'ID=gene-NONC;' in line else line for line in GFF]
    with pytest.raises(build.BuildError, match='8 fields, not 9'):
        list(refseq.bundles(_release(tmp_path, gff=truncated)))


def test_a_pseudoautosomal_transcript_is_modelled_on_the_chromosome_asked_for(tmp_path: pathlib.Path) -> None:
    shard = store_build.write_shard(
        refseq.bundles(_release(tmp_path)),
        tmp_path / 'shards',
        release='RS_TEST',
        inputs=[('annotation', tmp_path / 'genomic.gff.gz')],
    )
    store_build.write_index(tmp_path / 'store', [shard], assembly=bundle_pb2.ASSEMBLY_GRCH38)
    genome_build.build_genome(
        _fasta(tmp_path / 'genome.fna.gz', GENOME), tmp_path / 'genome', assembly=bundle_pb2.ASSEMBLY_GRCH38
    )
    provider = provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(tmp_path / 'genome'))
    )
    assert provider.get_transcript('NM_000040.1', Y)['reference_accession'] == Y
    assert provider.get_transcript('NM_000040.1', X)['reference_accession'] == X
