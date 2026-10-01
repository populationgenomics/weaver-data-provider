"""The Ensembl builder, from files in Ensembl's formats to bundles, and through the command to weaver.

The inputs are written here in the publisher's formats — Ensembl's GFF3 with its sequence-region rows and
`ID=transcript:` features, its GTF with the completeness tags, its cDNA, ncRNA and peptide FASTAs — over
a genome cut by the genome builder from a FASTA named by RefSeq accession, as NCBI's is. The synthetic
release, on sequence `22` (NC_000099.1) unless said otherwise:

- PLUSE, ENST00000000010.2 on the plus strand: exons 101-120 and 201-230, CDS 104-120 and 201-207,
  24 bases reading ATG, six codons and TAA.
- MINUSE, ENST00000000020.1 on the minus strand: exons 501-520 and 601-630, CDS 601-605 and 505-520,
  21 bases reading ATG, five codons and TAA in transcript orientation.
- NONCE, ENST00000000030.1, a lnc_RNA at 801-850 whose gene Ensembl gives no name.
- INCOMPLETE, ENST00000000040.1, an mRNA at 301-330 whose CDS the GTF tags cds_start_NF: it reads
  ATG, eight codons and TAA in frame, yet the annotation says the real start lies upstream.
- PATCHE, ENST00000000050.3 on the fix patch `HG1_PATCH` (NW_000001.1), MANE Select with the RefSeq
  partner NM_000070.1: the record's first five bases the chromosome lacks, so NCBI aligns the partner
  on the chromosome with a soft clip; on the patch, TT + record + TT, it sits whole at 3-37 with
  the CDS 3-26.
- NOSEQ, ENST00000000060.1, an unconfirmed_transcript Ensembl publishes no sequence for.
- MINUSP, ENST00000000080.1, a lnc_RNA on the patch's minus strand at 40-69, MANE-paired with
  NR_000080.1, which NCBI aligns on the chromosome at 950-980 with a mismatch and a base the record lacks.
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

from weaver_data_provider import build
from weaver_data_provider import genome as genome_mod
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import cli, ensembl
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2

CHROM, PATCH = 'NC_000099.1', 'NW_000001.1'
CHROM_NAME, PATCH_NAME = '22', 'HG1_PATCH'  # what the annotation calls them
_rng = random.Random(7)
_chrom = list(_rng.choices('ACGT', k=1000))


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


PLUSE_CDS = 'ATG' + 'GCTGAACAACCACTTTCT' + 'TAA'  # 24 bases: MAEQPLS
_chrom[103:120] = PLUSE_CDS[:17]
_chrom[200:207] = PLUSE_CDS[17:]
MINUSE_CDS = 'ATG' + 'GCTGAACAACCACTT' + 'TAA'  # 21 bases in transcript orientation: MAEQPL
_chrom[600:605] = _revcomp(MINUSE_CDS[:5])
_chrom[504:520] = _revcomp(MINUSE_CDS[5:])
INCOMPLETE_CDS = 'ATG' + 'GCTGAACAACCACTTTCTGAAGCT' + 'TAA'  # 30 bases, whole codons with a stop, at 301-330
_chrom[300:330] = INCOMPLETE_CDS
PATCHE_RECORD = 'ATG' + 'GCTGAACAACCACTTTCT' + 'TAA' + 'GGCCGGCCGGC'  # 35 bases; the chromosome lacks the first five
_chrom[900:930] = PATCHE_RECORD[5:]
MINUSP_RECORD = 'ACGTTGCAAGCTAGCTTACGGATCCATGCA'  # 30 bases, read off the patch's minus strand
_MINUSP_ON_GENOME = _revcomp(MINUSP_RECORD)  # the record in the genome's orientation
# the chromosome carries it with a base inserted after its tenth and its twentieth changed: 31 bases at 950-980
_chrom[949:980] = (
    _MINUSP_ON_GENOME[:10]
    + 'A'
    + _MINUSP_ON_GENOME[10:19]
    + ('C' if _MINUSP_ON_GENOME[19] != 'C' else 'G')
    + _MINUSP_ON_GENOME[20:]
)
GENOME = {CHROM: ''.join(_chrom)}
GENOME[PATCH] = 'TT' + PATCHE_RECORD + 'TT' + _MINUSP_ON_GENOME


def _g(chrom: str, start: int, end: int) -> str:
    """Genome bases start..end, 1-based closed."""
    return GENOME[chrom][start - 1 : end]


CDNA = {
    'ENST00000000010.2': _g(CHROM, 101, 120) + _g(CHROM, 201, 230),
    'ENST00000000020.1': _revcomp(_g(CHROM, 501, 520) + _g(CHROM, 601, 630)),
    'ENST00000000040.1': _g(CHROM, 301, 330),
    'ENST00000000050.3': PATCHE_RECORD,
}
NCRNA = {'ENST00000000030.1': _g(CHROM, 801, 850), 'ENST00000000080.1': MINUSP_RECORD}
PEPTIDES = {
    'ENSP00000000010.1': 'MAEQPLS',
    'ENSP00000000020.1': 'MAEQPL',
    'ENSP00000000040.1': 'MAEQPLSEA',
    'ENSP00000000050.2': 'MAEQPLS',
}


def _row(seqid: str, kind: str, start: int, end: int, strand: str, phase: str = '.', **attrs: str) -> str:
    columns = [seqid, 'ensembl', kind, str(start), str(end), '.', strand, phase]
    return '\t'.join([*columns, ';'.join(f'{k}={v}' for k, v in attrs.items())])


def _tx(
    seqid: str, kind: str, start: int, end: int, strand: str, enst: str, version: str, gene: str, **attrs: str
) -> str:
    ids = {'ID': f'transcript:{enst}', 'Parent': f'gene:{gene}', 'transcript_id': enst, 'version': version}
    return _row(seqid, kind, start, end, strand, **ids, **attrs)


def _cds(seqid: str, start: int, end: int, strand: str, phase: str, enst: str, ensp: str, version: str) -> str:
    ids = {'ID': f'CDS:{ensp}', 'Parent': f'transcript:{enst}', 'protein_id': ensp, 'version': version}
    return _row(seqid, 'CDS', start, end, strand, phase, **ids)


def _exon(seqid: str, start: int, end: int, strand: str, enst: str, ense: str, rank: str) -> str:
    return _row(seqid, 'exon', start, end, strand, Parent=f'transcript:{enst}', exon_id=ense, rank=rank)


def _gene(seqid: str, kind: str, start: int, end: int, strand: str, ensg: str, biotype: str, **attrs: str) -> str:
    return _row(seqid, kind, start, end, strand, ID=f'gene:{ensg}', biotype=biotype, gene_id=ensg, version='1', **attrs)


GFF = [
    '##gff-version 3',
    _row(CHROM_NAME, 'chromosome', 1, 1000, '.', ID=f'chromosome:{CHROM_NAME}', Alias=f'CM000099.1,chr22,{CHROM}'),
    _row(PATCH_NAME, 'scaffold', 1, len(GENOME[PATCH]), '.', ID=f'scaffold:{PATCH_NAME}', Alias=f'KN000001.1,{PATCH}'),
    _gene(
        CHROM_NAME, 'gene', 101, 230, '+', 'ENSG00000000010', 'protein_coding', Name='PLUSE', description='plus%20gene'
    ),
    _tx(
        CHROM_NAME,
        'mRNA',
        101,
        230,
        '+',
        'ENST00000000010',
        '2',
        'ENSG00000000010',
        biotype='protein_coding',
        tag='gencode_basic,Ensembl_canonical',
    ),
    _exon(CHROM_NAME, 101, 120, '+', 'ENST00000000010', 'ENSE00000000011', '1'),
    _exon(CHROM_NAME, 201, 230, '+', 'ENST00000000010', 'ENSE00000000012', '2'),
    _cds(CHROM_NAME, 104, 120, '+', '0', 'ENST00000000010', 'ENSP00000000010', '1'),
    _cds(CHROM_NAME, 201, 207, '+', '1', 'ENST00000000010', 'ENSP00000000010', '1'),
    _gene(CHROM_NAME, 'gene', 501, 630, '-', 'ENSG00000000020', 'protein_coding', Name='MINUSE'),
    _tx(
        CHROM_NAME,
        'mRNA',
        501,
        630,
        '-',
        'ENST00000000020',
        '1',
        'ENSG00000000020',
        biotype='protein_coding',
        tag='Ensembl_canonical',
    ),
    _exon(CHROM_NAME, 601, 630, '-', 'ENST00000000020', 'ENSE00000000021', '1'),
    _exon(CHROM_NAME, 501, 520, '-', 'ENST00000000020', 'ENSE00000000022', '2'),
    _cds(CHROM_NAME, 601, 605, '-', '0', 'ENST00000000020', 'ENSP00000000020', '1'),
    _cds(CHROM_NAME, 505, 520, '-', '1', 'ENST00000000020', 'ENSP00000000020', '1'),
    _gene(CHROM_NAME, 'ncRNA_gene', 801, 850, '+', 'ENSG00000000030', 'lncRNA'),
    _tx(CHROM_NAME, 'lnc_RNA', 801, 850, '+', 'ENST00000000030', '1', 'ENSG00000000030', biotype='lncRNA'),
    _exon(CHROM_NAME, 801, 850, '+', 'ENST00000000030', 'ENSE00000000031', '1'),
    _gene(CHROM_NAME, 'gene', 301, 330, '+', 'ENSG00000000040', 'protein_coding', Name='INCOMPLETE'),
    _tx(CHROM_NAME, 'mRNA', 301, 330, '+', 'ENST00000000040', '1', 'ENSG00000000040', biotype='protein_coding'),
    _exon(CHROM_NAME, 301, 330, '+', 'ENST00000000040', 'ENSE00000000041', '1'),
    _cds(CHROM_NAME, 301, 330, '+', '0', 'ENST00000000040', 'ENSP00000000040', '1'),
    _gene(PATCH_NAME, 'gene', 3, 37, '+', 'ENSG00000000050', 'protein_coding', Name='PATCHE'),
    _tx(
        PATCH_NAME,
        'mRNA',
        3,
        37,
        '+',
        'ENST00000000050',
        '3',
        'ENSG00000000050',
        biotype='protein_coding',
        tag='Ensembl_canonical,MANE_Select',
    ),
    _exon(PATCH_NAME, 3, 37, '+', 'ENST00000000050', 'ENSE00000000051', '1'),
    _cds(PATCH_NAME, 3, 26, '+', '0', 'ENST00000000050', 'ENSP00000000050', '2'),
    _gene(PATCH_NAME, 'ncRNA_gene', 40, 69, '-', 'ENSG00000000080', 'lncRNA', Name='MINUSP'),
    _tx(
        PATCH_NAME,
        'lnc_RNA',
        40,
        69,
        '-',
        'ENST00000000080',
        '1',
        'ENSG00000000080',
        biotype='lncRNA',
        tag='MANE_Select',
    ),
    _exon(PATCH_NAME, 40, 69, '-', 'ENST00000000080', 'ENSE00000000081', '1'),
    _gene(CHROM_NAME, 'gene', 951, 960, '+', 'ENSG00000000060', 'protein_coding', Name='NOSEQ'),
    _tx(
        CHROM_NAME,
        'unconfirmed_transcript',
        951,
        960,
        '+',
        'ENST00000000060',
        '1',
        'ENSG00000000060',
        biotype='protein_coding',
    ),
    _exon(CHROM_NAME, 951, 960, '+', 'ENST00000000060', 'ENSE00000000061', '1'),
]


def _gtf_row(enst: str, version: str, *tags: str) -> str:
    attrs = f'gene_id "ENSG"; gene_version "1"; transcript_id "{enst}"; transcript_version "{version}";'
    columns = ['22', 'ensembl', 'transcript', '1', '2', '.', '+', '.']
    return '\t'.join([*columns, attrs + ''.join(f' tag "{t}";' for t in tags)])


GTF = [
    '#!genome-build GRCh38',
    _gtf_row('ENST00000000010', '2', 'Ensembl_canonical'),
    _gtf_row('ENST00000000020', '1'),
    _gtf_row('ENST00000000030', '1'),
    _gtf_row('ENST00000000040', '1', 'cds_start_NF', 'basic'),
    _gtf_row('ENST00000000050', '3', 'MANE_Select'),
    _gtf_row('ENST00000000060', '1'),
    _gtf_row('ENST00000000080', '1'),
]


@dataclasses.dataclass(frozen=True)
class Read:
    name: str
    chrom: str
    start: int  # 0-based leftmost genome position, as in a BAM
    cigar: str
    reverse: bool
    sequence: str  # genome orientation, as a BAM stores it


READS = [
    Read('NM_000070.1', CHROM, 900, '5S30=', False, PATCHE_RECORD),
    # in genome orientation: ten matches, a genome-only base, nine matches, a mismatch, ten matches
    Read('NR_000080.1', CHROM, 949, '10=1D9=1X10=', True, _MINUSP_ON_GENOME),
]
HGNC_COLUMNS = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']
HGNC = [['HGNC:10', 'PLUSE', 'plus gene, HGNC', 'PL1', 'OLDPLUS', '10', 'ENSG00000000010']]
MANE = [['NM_000070.1', 'ENST00000000050.3', 'MANE Select'], ['NR_000080.1', 'ENST00000000080.1', 'MANE Select']]


def _fasta(path: pathlib.Path, records: dict[str, str], *, kind: str = 'cdna') -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        for name, residues in records.items():
            fh.write(f'>{name} {kind} chromosome:GRCh38:22:1:1000:1 gene:ENSG00000000000.1\n')
            for k in range(0, len(residues), 60):
                fh.write(residues[k : k + 60] + '\n')
    return path


def _bam(path: pathlib.Path, reads: list[Read], extra_references: dict[str, int]) -> pathlib.Path:
    lengths = {**{n: len(seq) for n, seq in GENOME.items()}, **extra_references}
    order = list(lengths)
    header = pysam.libcalignmentfile.AlignmentHeader.from_dict(
        {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': n, 'LN': lengths[n]} for n in order]}
    )
    with pysam.libcalignmentfile.AlignmentFile(str(path), 'wb', header=header) as out:
        for read in sorted(reads, key=lambda r: (order.index(r.chrom), r.start)):
            segment = pysam.libcalignedsegment.AlignedSegment(header)
            segment.query_name = read.name
            segment.reference_name = read.chrom
            segment.reference_start = read.start
            segment.flag = 16 if read.reverse else 0
            segment.mapping_quality = 60
            segment.cigarstring = read.cigar
            segment.query_sequence = read.sequence or None
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


def _genome(tmp_path: pathlib.Path, genome: dict[str, str]) -> pathlib.Path:
    fasta = tmp_path / 'genomic.fna.gz'
    with gzip.open(fasta, 'wt', encoding='ascii') as fh:
        for name, residues in genome.items():
            fh.write(f'>{name} synthetic\n')
            for k in range(0, len(residues), 70):
                fh.write(residues[k : k + 70] + '\n')
    genome_build.build_genome(fasta, tmp_path / 'genome', assembly=bundle_pb2.ASSEMBLY_GRCH38)
    return tmp_path / 'genome'


def _gzipped_lines(path: pathlib.Path, lines: list[str]) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    return path


def _release(
    tmp_path: pathlib.Path,
    *,
    gff: list[str] = GFF,
    gtf: list[str] = GTF,
    cdna: dict[str, str] = CDNA,
    ncrna: dict[str, str] = NCRNA,
    peptides: dict[str, str] = PEPTIDES,
    reads: list[Read] = READS,
    mane: list[list[str]] = MANE,
    extra_references: dict[str, int] | None = None,
    genome: dict[str, str] = GENOME,
    assembly: bundle_pb2.Assembly = bundle_pb2.ASSEMBLY_GRCH38,
    hgnc: list[list[str]] = HGNC,
) -> ensembl.Release:
    """The synthetic release written under `tmp_path`; each keyword replaces one input."""
    return ensembl.Release(
        assembly=assembly,
        release='116',
        annotation=_gzipped_lines(tmp_path / 'annotation.gff3.gz', gff),
        completeness=_gzipped_lines(tmp_path / 'annotation.gtf.gz', gtf),
        transcripts=(_fasta(tmp_path / 'cdna.fa.gz', cdna), _fasta(tmp_path / 'ncrna.fa.gz', ncrna, kind='ncrna')),
        proteins=_fasta(tmp_path / 'pep.fa.gz', peptides, kind='pep'),
        genome=_genome(tmp_path, genome),
        alignments=(_bam(tmp_path / 'alns.bam', reads, extra_references or {}),),
        hgnc=_tsv(tmp_path / 'hgnc.txt', HGNC_COLUMNS, hgnc, gzipped=False),
        mane=_tsv(tmp_path / 'mane.txt.gz', ['RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'], mane, gzipped=True),
    )


def _by_symbol(release: ensembl.Release) -> dict[str, bundle_pb2.GeneBundle]:
    return {b.gene.symbol: b for b in ensembl.bundles(release)}


def _exons(alignment: bundle_pb2.Alignment) -> list[tuple[int, int, int, int, str]]:
    return [
        (e.transcript_start, e.transcript_end, e.genome_start, e.genome_end_inclusive, e.cigar) for e in alignment.exons
    ]


def _residues(bundle: bundle_pb2.GeneBundle, digest: str) -> bytes:
    return {s.digest: s.residues for s in bundle.sequences}[digest]


def test_a_plus_strand_transcript_is_placed_at_its_exons_on_the_refseq_accession(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PLUSE'].transcripts
    (alignment,) = transcript.alignments
    assert (alignment.chromosome, alignment.strand, alignment.source) == (
        CHROM,
        bundle_pb2.STRAND_PLUS,
        bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION,
    )
    assert _exons(alignment) == [(0, 20, 100, 119, '20='), (20, 50, 200, 229, '30=')]


def test_a_plus_strand_cds_is_the_annotations_bounds_as_transcript_indices(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PLUSE'].transcripts
    # genome 104 is exon 1's fourth base; genome 207 is exon 2's seventh, after exon 1's 20
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (3, 26)
    assert not transcript.cds.start_open
    assert not transcript.cds.end_open


def test_a_minus_strand_transcript_runs_5_to_3_down_the_genome(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['MINUSE'].transcripts
    (alignment,) = transcript.alignments
    assert alignment.strand == bundle_pb2.STRAND_MINUS
    assert _exons(alignment) == [(0, 30, 600, 629, '30='), (30, 50, 500, 519, '20=')]


def test_a_minus_strand_cds_starts_at_its_highest_genome_base(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['MINUSE'].transcripts
    # genome 605 is 25 bases in from the exon's 5' end, 630; genome 505 is 15 in from 520, after exon 1's 30
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (25, 45)


def test_the_bundled_record_is_the_published_sequence(tmp_path: pathlib.Path) -> None:
    bundle = _by_symbol(_release(tmp_path))['MINUSE']
    (transcript,) = bundle.transcripts
    assert _residues(bundle, transcript.sequence_digest) == CDNA['ENST00000000020.1'].encode()


def test_an_ambiguity_code_in_the_genome_matches_the_n_the_record_writes(tmp_path: pathlib.Path) -> None:
    # the assembly says M, A or C, at NONCE's third base; Ensembl's sequence set writes N there
    genome = {**GENOME, CHROM: GENOME[CHROM][:802] + 'M' + GENOME[CHROM][803:]}
    record = NCRNA['ENST00000000030.1'][:2] + 'N' + NCRNA['ENST00000000030.1'][3:]
    bundle = _by_symbol(_release(tmp_path, genome=genome, ncrna={**NCRNA, 'ENST00000000030.1': record}))[
        'ENSG00000000030'
    ]
    (transcript,) = bundle.transcripts
    assert _residues(bundle, transcript.sequence_digest) == record.encode()


def test_a_record_that_is_not_the_genome_spliced_at_its_exons_fails_the_build(tmp_path: pathlib.Path) -> None:
    edited = CDNA['ENST00000000010.2']
    edited = edited[:10] + ('A' if edited[10] != 'A' else 'C') + edited[11:]
    with pytest.raises(build.BuildError, match=r'ENST00000000010\.2.*not the genome spliced at its exons'):
        list(ensembl.bundles(_release(tmp_path, cdna={**CDNA, 'ENST00000000010.2': edited})))


def test_a_sequence_with_no_refseq_accession_among_its_aliases_fails_the_build(tmp_path: pathlib.Path) -> None:
    unaliased = _row(CHROM_NAME, 'chromosome', 1, 1000, '.', ID=f'chromosome:{CHROM_NAME}', Alias='CM000099.1,chr22')
    gff = [unaliased if line.startswith(f'{CHROM_NAME}\tensembl\tchromosome') else line for line in GFF]
    with pytest.raises(build.BuildError, match=r'22: 0 RefSeq accessions among its aliases'):
        list(ensembl.bundles(_release(tmp_path, gff=gff)))


def test_a_noncoding_transcript_comes_from_the_ncrna_set(tmp_path: pathlib.Path) -> None:
    bundle = _by_symbol(_release(tmp_path))['ENSG00000000030']
    (transcript,) = bundle.transcripts
    assert _residues(bundle, transcript.sequence_digest) == NCRNA['ENST00000000030.1'].encode()


def test_a_noncoding_transcript_has_no_cds(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['ENSG00000000030'].transcripts
    assert not transcript.HasField('cds')


def test_the_biotype_is_ensembls_not_the_feature_type(tmp_path: pathlib.Path) -> None:
    bundles = _by_symbol(_release(tmp_path))
    (noncoding,) = bundles['ENSG00000000030'].transcripts
    (coding,) = bundles['PLUSE'].transcripts
    assert (noncoding.biotype, coding.biotype) == ('lncRNA', 'protein_coding')


def test_a_gene_with_no_name_is_named_by_its_stable_id(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path))['ENSG00000000030'].gene
    assert (gene.symbol, gene.ensembl_gene_id, gene.hgnc_id) == ('ENSG00000000030', 'ENSG00000000030', '')


def test_hgnc_names_a_gene_by_its_ensembl_id(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path))['PLUSE'].gene
    assert (gene.symbol, gene.hgnc_id, gene.ncbi_gene_id, gene.name) == ('PLUSE', 'HGNC:10', '10', 'plus gene, HGNC')
    assert list(gene.previous_symbols) == ['OLDPLUS']


@pytest.mark.parametrize(
    ('tags', 'open_ends'),
    [
        (('cds_start_NF',), (True, False)),
        (('cds_end_NF',), (False, True)),
        (('cds_start_NF', 'cds_end_NF'), (True, True)),
    ],
    ids=['start not found', 'end not found', 'neither found'],
)
def test_a_cds_end_the_gtf_tags_not_found_is_open(
    tmp_path: pathlib.Path, tags: tuple[str, ...], open_ends: tuple[bool, bool]
) -> None:
    gtf = [_gtf_row('ENST00000000040', '1', *tags) if 'ENST00000000040' in line else line for line in GTF]
    (transcript,) = _by_symbol(_release(tmp_path, gtf=gtf))['INCOMPLETE'].transcripts
    cds = transcript.cds
    assert (cds.start_index, cds.end_index_inclusive, cds.start_open, cds.end_open) == (0, 29, *open_ends)
    assert (transcript.protein_accession, transcript.protein_version) == ('ENSP00000000040', 1)


def test_completeness_is_the_gtfs_statement_whatever_the_cds_reads(tmp_path: pathlib.Path) -> None:
    # PLUSE's CDS cut one codon short no longer ends in a stop; the GTF still says it is complete
    short = _cds(CHROM_NAME, 201, 204, '+', '1', 'ENST00000000010', 'ENSP00000000010', '1')
    gff = [short if line.startswith(f'{CHROM_NAME}\tensembl\tCDS\t201\t') else line for line in GFF]
    (transcript,) = _by_symbol(_release(tmp_path, gff=gff))['PLUSE'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (3, 23)


def test_a_coding_transcript_the_gtf_does_not_name_fails_the_build(tmp_path: pathlib.Path) -> None:
    gtf = [line for line in GTF if 'transcript_id "ENST00000000020"' not in line]
    with pytest.raises(build.BuildError, match=r'ENST00000000020\.1: not in the GTF'):
        list(ensembl.bundles(_release(tmp_path, gtf=gtf)))


def test_a_coding_transcript_with_no_sequence_fails_the_build(tmp_path: pathlib.Path) -> None:
    cdna = {k: v for k, v in CDNA.items() if k != 'ENST00000000010.2'}
    with pytest.raises(build.BuildError, match=r'ENST00000000010\.2: a coding transcript with no sequence'):
        list(ensembl.bundles(_release(tmp_path, cdna=cdna)))


def test_a_transcript_ensembl_publishes_no_sequence_for_is_left_out(tmp_path: pathlib.Path) -> None:
    assert 'NOSEQ' not in _by_symbol(_release(tmp_path))


def test_a_transcript_on_a_sequence_the_assembly_no_longer_holds_is_left_out(tmp_path: pathlib.Path) -> None:
    # the annotation still carries a scaffold a later patch release retired; the genome has no such sequence
    retired = [
        _row('KI270000.1', 'scaffold', 1, 60, '.', ID='scaffold:KI270000.1', Alias='chrUn_KI270000v1,NT_187000.1'),
        _gene('KI270000.1', 'ncRNA_gene', 1, 50, '+', 'ENSG00000000070', 'lncRNA', Name='RETIRED'),
        _tx('KI270000.1', 'lnc_RNA', 1, 50, '+', 'ENST00000000070', '1', 'ENSG00000000070', biotype='lncRNA'),
        _exon('KI270000.1', 1, 50, '+', 'ENST00000000070', 'ENSE00000000071', '1'),
    ]
    ncrna = {**NCRNA, 'ENST00000000070.1': NCRNA['ENST00000000030.1']}
    assert 'RETIRED' not in _by_symbol(_release(tmp_path, gff=[*GFF, *retired], ncrna=ncrna))


def test_tags_and_the_mane_partner_come_from_the_annotation_and_the_summary(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PATCHE'].transcripts
    assert set(transcript.tags) == {bundle_pb2.TAG_ENSEMBL_CANONICAL, bundle_pb2.TAG_MANE_SELECT}
    assert transcript.mane_partner == 'NM_000070.1'


def test_a_mane_select_transcript_on_a_patch_gets_its_partners_chromosome_placement(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PATCHE'].transcripts
    assert [(a.chromosome, a.source, _exons(a)) for a in transcript.alignments] == [
        (CHROM, bundle_pb2.ALIGNMENT_SOURCE_MANE_PARTNER, [(0, 35, 900, 929, '5I30=')]),
        (PATCH, bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION, [(0, 35, 2, 36, '35=')]),
    ]


def test_a_mane_select_transcripts_cds_is_read_on_the_patch_that_carries_it_whole(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['PATCHE'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (0, 23)


def test_a_minus_strand_partner_alignment_with_a_mismatch_and_a_genome_only_base_is_copied(
    tmp_path: pathlib.Path,
) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['MINUSP'].transcripts
    assert [(a.chromosome, a.strand, _exons(a)) for a in transcript.alignments] == [
        (CHROM, bundle_pb2.STRAND_MINUS, [(0, 30, 949, 979, '10=1X9=1D10=')]),
        (PATCH, bundle_pb2.STRAND_MINUS, [(0, 30, 39, 68, '30=')]),
    ]


def test_a_transcript_on_a_chromosome_gets_no_partner_placements(tmp_path: pathlib.Path) -> None:
    # PLUSE paired with a partner NCBI aligns on the chromosome where PLUSE already is
    mane = [*MANE, ['NM_000010.2', 'ENST00000000010.2', 'MANE Select']]
    reads = [*READS, Read('NM_000010.2', CHROM, 100, '20=80N30=', False, CDNA['ENST00000000010.2'])]
    (transcript,) = _by_symbol(_release(tmp_path, reads=reads, mane=mane))['PLUSE'].transcripts
    assert [(a.chromosome, a.source) for a in transcript.alignments] == [
        (CHROM, bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION)
    ]


def test_a_partner_alignment_carrying_no_sequence_fails_the_build(tmp_path: pathlib.Path) -> None:
    reads = [Read('NM_000070.1', CHROM, 900, '5S30=', False, ''), READS[1]]
    with pytest.raises(build.BuildError, match=r'ENST00000000050\.3.*NM_000070\.1 on NC_000099\.1 carries no sequence'):
        list(ensembl.bundles(_release(tmp_path, reads=reads)))


def test_a_partner_ncbi_aligned_with_another_sequence_fails_the_build(tmp_path: pathlib.Path) -> None:
    other = PATCHE_RECORD[:-1] + ('A' if PATCHE_RECORD[-1] != 'A' else 'C')
    reads = [Read('NM_000070.1', CHROM, 900, '5S30=', False, other), READS[1]]
    with pytest.raises(build.BuildError, match=r'ENST00000000050\.3.*NM_000070\.1 on NC_000099\.1 is not this record'):
        list(ensembl.bundles(_release(tmp_path, reads=reads)))


def test_a_partner_alignment_that_does_not_describe_the_record_against_the_genome_fails_the_build(
    tmp_path: pathlib.Path,
) -> None:
    # the same record, aligned five bases upstream as if the chromosome carried its first five bases
    reads = [Read('NM_000070.1', CHROM, 895, '35=', False, PATCHE_RECORD), READS[1]]
    with pytest.raises(
        build.BuildError, match=r'ENST00000000050\.3.*NM_000070\.1 aligns on NC_000099\.1.*does not describe'
    ):
        list(ensembl.bundles(_release(tmp_path, reads=reads)))


def test_a_partners_alignment_on_a_chromosome_the_genome_lacks_fails_the_build(tmp_path: pathlib.Path) -> None:
    reads = [Read('NM_000070.1', 'NC_000098.1', 900, '5S30=', False, PATCHE_RECORD), READS[1]]
    with pytest.raises(build.BuildError, match=r'ENST00000000050\.3.*NC_000098\.1, which the genome lacks'):
        list(ensembl.bundles(_release(tmp_path, reads=reads, extra_references={'NC_000098.1': 1000})))


def test_a_protein_missing_from_the_peptide_set_fails_the_build(tmp_path: pathlib.Path) -> None:
    peptides = {k: v for k, v in PEPTIDES.items() if k != 'ENSP00000000020.1'}
    with pytest.raises(build.BuildError, match=r'ENST00000000020\.1.*ENSP00000000020\.1 is not in the peptide set'):
        list(ensembl.bundles(_release(tmp_path, peptides=peptides)))


def test_a_genome_on_another_assembly_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match='a genome on ASSEMBLY_GRCH38 for a ASSEMBLY_GRCH37 release'):
        list(ensembl.bundles(_release(tmp_path, assembly=bundle_pb2.ASSEMBLY_GRCH37)))


def test_bundles_come_in_symbol_order(tmp_path: pathlib.Path) -> None:
    symbols = [b.gene.symbol for b in ensembl.bundles(_release(tmp_path))]
    assert symbols == ['ENSG00000000030', 'INCOMPLETE', 'MINUSE', 'MINUSP', 'PATCHE', 'PLUSE']


def test_the_built_commands_feed_weaver(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> None:
    release = _release(tmp_path)
    cli.main([
        'ensembl', '--assembly', 'GRCh38', '--release', '116',
        '--annotation', str(release.annotation), '--completeness', str(release.completeness),
        '--transcripts', str(release.transcripts[0]), '--transcripts', str(release.transcripts[1]),
        '--proteins', str(release.proteins), '--genome', str(release.genome),
        '--alignments', str(release.alignments[0]),
        '--hgnc', str(release.hgnc), '--mane', str(release.mane), '--shards', str(tmp_path / 'shards'),
    ])  # fmt: skip
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    provider = provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(release.genome))
    )
    mapper = weaver.VariantMapper(provider)
    # PLUSE's c.1 is genome 104, the A of its ATG
    variant = weaver.parse('ENST00000000010.2:c.1A>G')
    assert variant.validate(provider)
    assert mapper.c_to_g(variant, CHROM).format() == f'{CHROM}:g.104A>G'


# two HGNC genes Ensembl annotates as one, as HGNC maps LINC00595 and LINC00856 to ENSG00000230417
_PLUSE_TWICE = [['HGNC:9', 'PLUSB', 'plus gene B', '', 'OLDPLUSB', '9', 'ENSG00000000010'], *HGNC]


def test_of_two_hgnc_rows_for_a_gene_the_one_named_as_the_annotation_names_it(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path, hgnc=_PLUSE_TWICE))['PLUSE'].gene
    assert (gene.symbol, gene.hgnc_id) == ('PLUSE', 'HGNC:10')  # not PLUSB, whose id is lower


def test_the_other_hgnc_rows_symbols_stay_findable_as_aliases(tmp_path: pathlib.Path) -> None:
    gene = _by_symbol(_release(tmp_path, hgnc=_PLUSE_TWICE))['PLUSE'].gene
    assert {'PLUSB', 'OLDPLUSB'} <= set(gene.alias_symbols)


def test_of_two_hgnc_rows_neither_named_as_the_annotation_the_lowest_id_names_the_gene(tmp_path: pathlib.Path) -> None:
    gff = [line.replace('Name=PLUSE;', 'Name=PLUSX;') for line in GFF]
    bundles = _by_symbol(_release(tmp_path, gff=gff, hgnc=_PLUSE_TWICE))
    assert (bundles['PLUSB'].gene.symbol, bundles['PLUSB'].gene.hgnc_id) == ('PLUSB', 'HGNC:9')


def test_a_description_loses_the_note_of_where_ensembl_took_it_from(tmp_path: pathlib.Path) -> None:
    # without HGNC's row the name is the annotation's description, which Ensembl ends with its source
    gff = [
        line.replace('description=plus%20gene', 'description=plus gene [Source:HGNC Symbol%3BAcc:HGNC:10]')
        for line in GFF
    ]
    assert _by_symbol(_release(tmp_path, gff=gff, hgnc=[]))['PLUSE'].gene.name == 'plus gene'


def _provider(tmp_path: pathlib.Path) -> provider_mod.BundleProvider:
    """The synthetic release built, indexed and opened as weaver's provider."""
    release = _release(tmp_path)
    shard = store_build.write_shard(
        ensembl.bundles(release), tmp_path / 'shards', release='116', inputs=release.inputs()
    )
    store_build.write_index(tmp_path / 'store', [shard], assembly=bundle_pb2.ASSEMBLY_GRCH38)
    return provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(release.genome))
    )


def test_ensembl_accessions_are_identified_as_transcripts_and_proteins(tmp_path: pathlib.Path) -> None:
    provider = _provider(tmp_path)
    kinds = (provider.get_identifier_type('ENST00000000010.2'), provider.get_identifier_type('ENSP00000000010.1'))
    assert kinds == (weaver.IdentifierType.TranscriptAccession, weaver.IdentifierType.ProteinAccession)


def test_an_ensembl_transcripts_protein_is_its_c_to_p_target(tmp_path: pathlib.Path) -> None:
    assert _provider(tmp_path).get_symbol_accessions('ENST00000000010.2', 'c', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'ENSP00000000010.1')
    ]


def test_an_ensembl_protein_is_read_by_its_accession(tmp_path: pathlib.Path) -> None:
    provider = _provider(tmp_path)
    assert provider.get_seq('ENSP00000000010.1', 0, None, weaver.IdentifierType.ProteinAccession) == 'MAEQPLS'
