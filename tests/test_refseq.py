"""The RefSeq builder, from files in NCBI's formats to bundles, and through the command to weaver.

The inputs are written here in the publisher's formats — GFF3, FASTA and TSV as text, the alignments
as a BAM that pysam writes — so no fixture comes from the builder itself. The synthetic release:

- PLUS, NM_000010.2 on NC_000099.1's plus strand: exons 101-120 and 201-230, CDS 105-210.
- MINUS, NM_000020.1 on the minus strand: exons 501-520 and 601-630, CDS 505-605; its record carries
  one base the genome lacks, after genome 610, so the alignment has an insertion.
- PAR, NM_000040.1 at 51-100 on both X and Y, as a pseudoautosomal transcript is.
- NONC, NR_000030.1, a non-coding transcript NCBI published no alignment for.
- ALT, NM_000050.1, a coding gene on an alternate locus only: exon 1-10 of NT_187361.1, CDS 3-8.
- PATCHED, NM_000060.1, placed on CHROM (exons 701-730 and 741-760, CDS 705-750, which reads ATG, ten
  codons and TAA) and again on the fix patch NW_000001.1, whose bases 1-150 are CHROM's 651-800, so the
  same exons sit at 51-80 and 91-110.
  Its second gene feature there carries a suffixed id and the same GeneID, as NCBI writes one, and a
  transcript NCBI places only on the patch, NR_000061.1 at 121-140. The patch rows come before the
  chromosome rows.
- CLIPPED, NM_000070.1: a record whose first five bases, GGGCC, the chromosome lacks; the rest is
  CHROM 901-930. On CHROM its alignment soft-clips them and its CDS, 901-920, is stated open at the
  start. On the fix patch NW_000002.1, whose bases are TT + the record + TT, it aligns whole at 3-37
  with the CDS 3-27 complete.
- MT-TF, a mitochondrial tRNA with no transcript accession.
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

CHROM, X, Y, ALT_LOCUS, PATCH, PATCH2 = (
    'NC_000099.1',
    'NC_000023.11',
    'NC_000024.10',
    'NT_187361.1',
    'NW_000001.1',
    'NW_000002.1',
)
_rng = random.Random(99)
GENOME = {
    CHROM: ''.join(_rng.choices('ACGT', k=1000)),
    X: ''.join(_rng.choices('ACGT', k=300)),
    ALT_LOCUS: 'ACGTACGTAC' + ''.join(_rng.choices('ACGT', k=90)),
}
PATCHED_CDS = (
    'ATG' + 'TCTGAACAACCAGCTACTTGTCATAAAGCT' + 'TAA'
)  # what NM_000060.1 reads over its CDS, 705-730 and 741-750
GENOME[CHROM] = GENOME[CHROM][:704] + PATCHED_CDS[:26] + GENOME[CHROM][730:740] + PATCHED_CDS[26:] + GENOME[CHROM][750:]
GENOME[Y] = GENOME[X]  # the pseudoautosomal region is the same sequence on both
GENOME[PATCH] = GENOME[CHROM][650:800]  # a fix patch carrying CHROM's 651-800
GENOME[PATCH2] = 'TT' + 'GGGCC' + GENOME[CHROM][900:930] + 'TT'  # a fix patch supplying CLIPPED's five 5' bases
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
    'NM_000050.1': _g(ALT_LOCUS, 1, 10),
    'NM_000060.1': _g(CHROM, 701, 730) + _g(CHROM, 741, 760),
    'NR_000061.1': _g(PATCH, 121, 140),
    'NM_000070.1': 'GGGCC' + _g(CHROM, 901, 930),
}
PROTEINS = {
    'NP_000010.1': 'MSEQPLUS',
    'NP_000020.1': 'MSEQMINUS',
    'NP_000050.1': 'MA',
    'NP_000060.1': 'MSEQPATCHKA',
    'NP_000070.1': 'MSEQCLIP',
}


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
    _row(ALT_LOCUS, 'gene', 1, 10, '+', ID='gene-ALT', Dbxref='GeneID:60', Name='ALT'),
    _row(ALT_LOCUS, 'mRNA', 1, 10, '+', ID='rna-NM_000050.1', Parent='gene-ALT', transcript_id='NM_000050.1'),
    _row(ALT_LOCUS, 'CDS', 3, 8, '+', ID='cds-NP_000050.1', Parent='rna-NM_000050.1', protein_id='NP_000050.1'),
    _row(PATCH, 'gene', 51, 110, '+', ID='gene-PATCHED-2', Dbxref='GeneID:70', Name='PATCHED'),
    _row(PATCH, 'mRNA', 51, 110, '+', ID='rna-NM_000060.1-2', Parent='gene-PATCHED-2', transcript_id='NM_000060.1'),
    _row(PATCH, 'CDS', 55, 80, '+', ID='cds-NP_000060.1-2', Parent='rna-NM_000060.1-2', protein_id='NP_000060.1'),
    _row(PATCH, 'CDS', 91, 100, '+', ID='cds-NP_000060.1-2', Parent='rna-NM_000060.1-2', protein_id='NP_000060.1'),
    _row(PATCH, 'lnc_RNA', 121, 140, '+', ID='rna-NR_000061.1', Parent='gene-PATCHED-2', transcript_id='NR_000061.1'),
    _row(CHROM, 'gene', 701, 760, '+', ID='gene-PATCHED', Dbxref='GeneID:70', Name='PATCHED'),
    _row(CHROM, 'mRNA', 701, 760, '+', ID='rna-NM_000060.1', Parent='gene-PATCHED', transcript_id='NM_000060.1'),
    _row(CHROM, 'CDS', 705, 730, '+', ID='cds-NP_000060.1', Parent='rna-NM_000060.1', protein_id='NP_000060.1'),
    _row(CHROM, 'CDS', 741, 750, '+', ID='cds-NP_000060.1', Parent='rna-NM_000060.1', protein_id='NP_000060.1'),
    _row(CHROM, 'gene', 901, 930, '+', ID='gene-CLIPPED', Dbxref='GeneID:80', Name='CLIPPED'),
    _row(CHROM, 'mRNA', 901, 930, '+', ID='rna-NM_000070.1', Parent='gene-CLIPPED', transcript_id='NM_000070.1'),
    _row(
        CHROM,
        'CDS',
        901,
        920,
        '+',
        ID='cds-NP_000070.1',
        Parent='rna-NM_000070.1',
        protein_id='NP_000070.1',
        partial='true',
        start_range='.,901',
    ),
    _row(PATCH2, 'gene', 3, 37, '+', ID='gene-CLIPPED-2', Dbxref='GeneID:80', Name='CLIPPED'),
    _row(PATCH2, 'mRNA', 3, 37, '+', ID='rna-NM_000070.1-2', Parent='gene-CLIPPED-2', transcript_id='NM_000070.1'),
    _row(PATCH2, 'CDS', 3, 27, '+', ID='cds-NP_000070.1-2', Parent='rna-NM_000070.1-2', protein_id='NP_000070.1'),
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
    Read('NM_000050.1', ALT_LOCUS, 0, '10=', False, TRANSCRIPTS['NM_000050.1']),
    Read('NM_000060.1', CHROM, 700, '30=10N20=', False, TRANSCRIPTS['NM_000060.1']),
    Read('NM_000060.1', PATCH, 50, '30=10N20=', False, TRANSCRIPTS['NM_000060.1']),
    Read('NR_000061.1', PATCH, 120, '20=', False, TRANSCRIPTS['NR_000061.1']),
    Read('NM_000070.1', CHROM, 900, '5S30=', False, TRANSCRIPTS['NM_000070.1']),
    Read('NM_000070.1', PATCH2, 2, '35=', False, TRANSCRIPTS['NM_000070.1']),
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
    order = list(GENOME)
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


def _transcript(bundle: bundle_pb2.GeneBundle, accession: str) -> bundle_pb2.Transcript:
    (found,) = [t for t in bundle.transcripts if t.accession == accession]
    return found


def _rows(gff: list[str], **attrs: str) -> list[str]:
    """The rows of a GFF carrying every one of these attribute values."""
    return [line for line in gff if all(f'{k}={v}' in line.split('\t')[-1].split(';') for k, v in attrs.items())]


def _replacing(gff: list[str], rows: list[str], with_rows: list[str]) -> list[str]:
    """`gff` with `rows` taken out and `with_rows` put where the first of them was."""
    first = min(gff.index(row) for row in rows)
    before = [line for line in gff[:first] if line not in rows]
    after = [line for line in gff[first:] if line not in rows]
    return [*before, *with_rows, *after]


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


def test_a_transcript_on_an_alternate_locus_only_is_bundled_with_its_placement_there(tmp_path: pathlib.Path) -> None:
    (transcript,) = _by_symbol(_release(tmp_path))['ALT'].transcripts
    (alignment,) = transcript.alignments
    assert alignment.chromosome == ALT_LOCUS
    assert _exons(alignment) == [(0, 10, 0, 9, '10=')]
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (2, 7)


def test_a_gene_placed_on_a_chromosome_and_a_patch_is_one_bundle(tmp_path: pathlib.Path) -> None:
    bundles = [b for b in refseq.bundles(_release(tmp_path)) if b.gene.symbol == 'PATCHED']
    (bundle,) = bundles
    assert sorted(f'{t.accession}.{t.version}' for t in bundle.transcripts) == ['NM_000060.1', 'NR_000061.1']


def test_a_transcript_placed_only_on_a_patch_keeps_that_placement(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_release(tmp_path))['PATCHED'], 'NR_000061')
    assert [(a.chromosome, _exons(a)) for a in transcript.alignments] == [(PATCH, [(0, 20, 120, 139, '20=')])]


def test_placements_are_kept_chromosome_first_whatever_order_the_annotation_lists_them(
    tmp_path: pathlib.Path,
) -> None:
    # the fixture lists PATCHED's patch rows before its chromosome rows
    transcript = _transcript(_by_symbol(_release(tmp_path))['PATCHED'], 'NM_000060')
    assert [a.chromosome for a in transcript.alignments] == [CHROM, PATCH]


def test_a_patch_placement_takes_its_exons_from_the_alignment_on_the_patch(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_release(tmp_path))['PATCHED'], 'NM_000060')
    (on_patch,) = [a for a in transcript.alignments if a.chromosome == PATCH]
    assert _exons(on_patch) == [(0, 30, 50, 79, '30='), (30, 50, 90, 109, '20=')]


def test_agreeing_cds_bounds_on_two_placements_project_to_one_cds(tmp_path: pathlib.Path) -> None:
    # genome 705 is exon 1's fifth base on CHROM, as 55 is on the patch; both project to index 4
    transcript = _transcript(_by_symbol(_release(tmp_path))['PATCHED'], 'NM_000060')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 39)


def _patch_cds_end(end: int) -> list[str]:
    """The fixture with PATCHED's last CDS row on the patch ending at `end` instead of 100."""
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000060.1-2') if '\tCDS\t91\t' in r]
    moved = _row(
        PATCH, 'CDS', 91, end, '+', ID='cds-NP_000060.1-2', Parent='rna-NM_000060.1-2', protein_id='NP_000060.1'
    )
    return _replacing(GFF, [row], [moved])


@pytest.mark.parametrize('patch_end', [99, 101, 97], ids=['one base short', 'one base long', 'one codon short'])
def test_placements_whose_cds_bounds_disagree_take_the_one_the_record_encodes_its_protein_over(
    tmp_path: pathlib.Path, patch_end: int
) -> None:
    # 4..39 through the chromosome reads ATG, ten codons and TAA; the patch's projection is off it by a base, which
    # breaks the frame, or by a codon, which drops the stop
    transcript = _transcript(_by_symbol(_release(tmp_path, gff=_patch_cds_end(patch_end)))['PATCHED'], 'NM_000060')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 39)


def test_the_record_can_confirm_a_patch_projection_over_the_chromosomes(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # the chromosome states the CDS one base short, 741-749; the patch states it whole
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000060.1') if '\tCDS\t741\t' in r]
    short = _row(CHROM, 'CDS', 741, 749, '+', ID='cds-NP_000060.1', Parent='rna-NM_000060.1', protein_id='NP_000060.1')
    transcript = _transcript(
        _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [short])))['PATCHED'], 'NM_000060'
    )
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 39)
    assert f'NM_000060.1 on {CHROM}' in capsys.readouterr().err


def test_a_placement_that_loses_the_cds_arbitration_keeps_its_alignment(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_release(tmp_path, gff=_patch_cds_end(99)))['PATCHED'], 'NM_000060')
    assert [a.chromosome for a in transcript.alignments] == [CHROM, PATCH]


def test_placements_whose_cds_bounds_disagree_and_the_protein_confirms_neither_fail_the_build(
    tmp_path: pathlib.Path,
) -> None:
    # a protein set naming a protein the record does not encode over either projection
    proteins = {**PROTEINS, 'NP_000060.1': 'MSEQPATCH'}
    with pytest.raises(build.BuildError, match=rf'NM_000060\.1.*4\.\.38 on {PATCH}, 4\.\.39 on {CHROM}.*over none'):
        list(refseq.bundles(_release(tmp_path, gff=_patch_cds_end(99), proteins=proteins)))


def test_placements_whose_cds_bounds_disagree_with_no_protein_in_the_set_fail_the_build(tmp_path: pathlib.Path) -> None:
    proteins = {k: v for k, v in PROTEINS.items() if k != 'NP_000060.1'}
    with pytest.raises(build.BuildError, match=r'NM_000060\.1: its protein NP_000060\.1 is not in the protein set'):
        list(refseq.bundles(_release(tmp_path, gff=_patch_cds_end(99), proteins=proteins)))


def test_a_determining_placement_with_no_alignment_leaves_the_cds_undetermined(tmp_path: pathlib.Path) -> None:
    # CLIPPED's whole CDS is stated on the patch only; without the patch alignment nothing projects it
    reads = [r for r in READS if not (r.name == 'NM_000070.1' and r.chrom == PATCH2)]
    transcript = _transcript(_by_symbol(_release(tmp_path, reads=reads))['CLIPPED'], 'NM_000070')
    assert transcript.cds_undetermined
    assert not transcript.HasField('cds')


def test_a_clipped_placement_whose_cds_start_lies_in_the_clip_does_not_determine_the_cds(
    tmp_path: pathlib.Path,
) -> None:
    # CLIPPED's CDS begins in the five bases the chromosome lacks; the patch states it whole: indices 0..24
    transcript = _transcript(_by_symbol(_release(tmp_path))['CLIPPED'], 'NM_000070')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (0, 24)
    assert not transcript.cds_undetermined


def test_a_clipped_placement_is_kept_with_the_clipped_bases_as_transcript_only(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_release(tmp_path))['CLIPPED'], 'NM_000070')
    assert [(a.chromosome, _exons(a)) for a in transcript.alignments] == [
        (CHROM, [(0, 35, 900, 929, '5I30=')]),
        (PATCH2, [(0, 35, 2, 36, '35=')]),
    ]


def test_a_cds_broken_internally_but_closed_at_both_ends_is_determined(tmp_path: pathlib.Path) -> None:
    # NCBI marks both rows either side of a frameshift partial, with the break's ends open; the outer ends stay closed.
    # PLUS has one placement, so nothing else could state the CDS.
    rows = [r for r in _rows(GFF, Parent='rna-NM_000010.2') if '\tCDS\t' in r]
    broken = [
        _row(
            CHROM,
            'CDS',
            105,
            120,
            '+',
            ID='cds-NP_000010.1',
            Parent='rna-NM_000010.2',
            protein_id='NP_000010.1',
            partial='true',
            end_range='120,.',
        ),
        _row(
            CHROM,
            'CDS',
            201,
            210,
            '+',
            ID='cds-NP_000010.1',
            Parent='rna-NM_000010.2',
            protein_id='NP_000010.1',
            partial='true',
            start_range='.,201',
        ),
    ]
    (transcript,) = _by_symbol(_release(tmp_path, gff=_replacing(GFF, rows, broken)))['PLUS'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 29)


def test_a_coding_transcript_no_placement_states_the_whole_cds_of_is_bundled_without_one(
    tmp_path: pathlib.Path,
) -> None:
    (row,) = _rows(GFF, Parent='rna-NM_000050.1')
    open_start = _row(
        ALT_LOCUS,
        'CDS',
        3,
        8,
        '+',
        ID='cds-NP_000050.1',
        Parent='rna-NM_000050.1',
        protein_id='NP_000050.1',
        partial='true',
        start_range='.,3',
    )
    (transcript,) = _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [open_start])))['ALT'].transcripts
    assert transcript.cds_undetermined
    assert not transcript.HasField('cds')
    assert (transcript.protein_accession, transcript.protein_version) == ('NP_000050', 1)


def test_a_minus_strand_cds_open_at_its_high_genome_end_does_not_determine_the_cds(tmp_path: pathlib.Path) -> None:
    # MINUS reads 5' to 3' downward on the genome, so an open start is `end_range` on its highest row, 601-605
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000020.1') if '\tCDS\t601\t' in r]
    open_5 = _row(
        CHROM,
        'CDS',
        601,
        605,
        '-',
        ID='cds-NP_000020.1',
        Parent='rna-NM_000020.1',
        protein_id='NP_000020.1',
        partial='true',
        end_range='605,.',
    )
    (transcript,) = _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [open_5])))['MINUS'].transcripts
    assert transcript.cds_undetermined


def test_a_minus_strand_cds_open_only_at_an_internal_break_is_determined(tmp_path: pathlib.Path) -> None:
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000020.1') if '\tCDS\t601\t' in r]
    inner_open = _row(
        CHROM,
        'CDS',
        601,
        605,
        '-',
        ID='cds-NP_000020.1',
        Parent='rna-NM_000020.1',
        protein_id='NP_000020.1',
        partial='true',
        start_range='.,601',
    )
    (transcript,) = _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [inner_open])))['MINUS'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (26, 46)


def test_a_transcript_on_a_scaffold_and_a_patch_only_is_kept_scaffold_first(tmp_path: pathlib.Path) -> None:
    # ALT's transcript placed again on PATCH, listed first, at 141-150
    patch_copy = [
        _row(PATCH, 'gene', 141, 150, '+', ID='gene-ALT-2', Dbxref='GeneID:60', Name='ALT'),
        _row(PATCH, 'mRNA', 141, 150, '+', ID='rna-NM_000050.1-2', Parent='gene-ALT-2', transcript_id='NM_000050.1'),
        _row(PATCH, 'CDS', 143, 148, '+', ID='cds-NP_000050.1-2', Parent='rna-NM_000050.1-2', protein_id='NP_000050.1'),
    ]
    read = Read('NM_000050.1', PATCH, 140, '10=', False, TRANSCRIPTS['NM_000050.1'])
    release = _release(tmp_path, gff=[*patch_copy, *GFF], reads=[*READS, read])
    (transcript,) = _by_symbol(release)['ALT'].transcripts
    assert [a.chromosome for a in transcript.alignments] == [ALT_LOCUS, PATCH]


def test_a_transcript_placed_twice_on_one_sequence_fails_the_build(tmp_path: pathlib.Path) -> None:
    again = _row(CHROM, 'mRNA', 801, 850, '+', ID='rna-NM_000010.2-2', Parent='gene-PLUS', transcript_id='NM_000010.2')
    with pytest.raises(build.BuildError, match=rf'NM_000010\.2: placed twice on {CHROM}'):
        list(refseq.bundles(_release(tmp_path, gff=[*GFF, again])))


def test_a_feature_with_no_transcript_accession_is_not_bundled(tmp_path: pathlib.Path) -> None:
    # MT-TF's tRNA names no accession
    assert 'MT-TF' not in _by_symbol(_release(tmp_path))


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


def _provider(tmp_path: pathlib.Path, *, gff: list[str] = GFF) -> provider_mod.BundleProvider:
    """The synthetic release built, indexed and cut, and opened as weaver's provider."""
    shard = store_build.write_shard(
        refseq.bundles(_release(tmp_path, gff=gff)),
        tmp_path / 'shards',
        release='RS_TEST',
        inputs=[('annotation', tmp_path / 'genomic.gff.gz')],
    )
    store_build.write_index(tmp_path / 'store', [shard], assembly=bundle_pb2.ASSEMBLY_GRCH38)
    genome_build.build_genome(
        _fasta(tmp_path / 'genome.fna.gz', GENOME), tmp_path / 'genome', assembly=bundle_pb2.ASSEMBLY_GRCH38
    )
    return provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(tmp_path / 'genome'))
    )


def test_a_pseudoautosomal_transcript_is_modelled_on_the_chromosome_asked_for(tmp_path: pathlib.Path) -> None:
    provider = _provider(tmp_path)
    assert provider.get_transcript('NM_000040.1', Y)['reference_accession'] == Y
    assert provider.get_transcript('NM_000040.1', X)['reference_accession'] == X


def test_a_transcript_also_on_a_patch_is_modelled_on_the_chromosome_when_none_is_named(
    tmp_path: pathlib.Path,
) -> None:
    assert _provider(tmp_path).get_transcript('NM_000060.1', None)['reference_accession'] == CHROM


def test_a_patch_a_transcript_is_placed_on_is_modelled_when_named(tmp_path: pathlib.Path) -> None:
    assert _provider(tmp_path).get_transcript('NM_000060.1', PATCH)['reference_accession'] == PATCH


def test_a_transcript_on_an_alternate_locus_only_projects_to_that_locus(tmp_path: pathlib.Path) -> None:
    # ALT's c.1 is NT_187361.1's third base, a G
    variant = weaver.parse('NM_000050.1:c.1G>C')
    assert weaver.VariantMapper(_provider(tmp_path)).c_to_g(variant, None).format() == f'{ALT_LOCUS}:g.3G>C'


def test_a_region_of_an_alternate_locus_finds_the_transcripts_placed_on_it(tmp_path: pathlib.Path) -> None:
    assert _provider(tmp_path).get_transcripts_for_region(ALT_LOCUS, 2, 2) == ['NM_000050.1']


def test_a_sequence_a_transcript_is_not_placed_on_is_refused_by_name(tmp_path: pathlib.Path) -> None:
    with pytest.raises(weaver.DataProviderError, match=f'no alignment on {CHROM}'):
        _provider(tmp_path).get_transcript('NM_000050.1', CHROM)


def _open_start_alt(tmp_path: pathlib.Path) -> provider_mod.BundleProvider:
    """The provider over a release whose ALT CDS is stated open at its start, so no placement determines it."""
    (row,) = _rows(GFF, Parent='rna-NM_000050.1')
    open_start = _row(
        ALT_LOCUS,
        'CDS',
        3,
        8,
        '+',
        ID='cds-NP_000050.1',
        Parent='rna-NM_000050.1',
        protein_id='NP_000050.1',
        partial='true',
        start_range='.,3',
    )
    return _provider(tmp_path, gff=_replacing(GFF, [row], [open_start]))


def test_a_coding_transcript_with_no_determined_cds_is_refused_rather_than_served_as_non_coding(
    tmp_path: pathlib.Path,
) -> None:
    with pytest.raises(weaver.DataProviderError, match=r'NM_000050\.1: coding, but no placement.*whole CDS'):
        _open_start_alt(tmp_path).get_transcript('NM_000050.1', None)


def test_a_coding_transcript_with_no_determined_cds_still_has_its_sequence(tmp_path: pathlib.Path) -> None:
    assert _open_start_alt(tmp_path).get_seq('NM_000050.1', 0, 10, 'c') == TRANSCRIPTS['NM_000050.1']
