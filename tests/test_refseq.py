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

Each coding record's CDS is written in a GenBank flat file, as NCBI's rna.gbff states it, in the
record's own coordinates. The same files, less the sequence sets and plus a status table, are NCBI's
historical set: the records carry the sequences and translations, and the table says which versions
Entrez holds to be replaced or suppressed.
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib
import random
import warnings
from collections.abc import Callable

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
# Each coding record's CDS as its GenBank record states it: 1-based, closed, in the record's coordinates.
RECORD_CDS = {
    'NM_000010.2': '5..30',
    'NM_000020.1': '27..47',
    'NM_000050.1': '3..8',
    'NM_000060.1': '5..40',
    'NM_000070.1': '1..25',
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


def _genbank(
    path: pathlib.Path,
    records: dict[str, str],
    cds: dict[str, str],
    *,
    proteins: dict[str, tuple[str, str]] | None = None,
    sequences: bool = False,
) -> pathlib.Path:
    """A GenBank flat file of these records, each coding one with its CDS location.

    With `proteins` (each coding record's protein id and residues) the CDS carries `/protein_id` and
    `/translation`, the latter wrapped onto further lines and closed on its last as NCBI writes it; a
    coding record then has a further feature with a qualifier of its own after the CDS, as NCBI's do.
    With `sequences` each record has its ORIGIN section, numbered and spaced as NCBI writes it.
    """
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        for versioned, residues in records.items():
            accession = versioned.split('.')[0]
            fh.write(f'LOCUS       {accession}  {len(residues)} bp    mRNA    linear   PRI 01-JAN-2026\n')
            fh.write(f'DEFINITION  synthetic.\nACCESSION   {accession}\nVERSION     {versioned}\n')
            fh.write('FEATURES             Location/Qualifiers\n')
            fh.write(f'     source          1..{len(residues)}\n')
            if versioned in cds:
                fh.write(f'     CDS             {cds[versioned]}\n                     /codon_start=1\n')
                if proteins is not None and versioned in proteins:
                    protein_id, translation = proteins[versioned]
                    fh.write(f'                     /protein_id="{protein_id}"\n')
                    pieces = [translation[k : k + 4] for k in range(0, len(translation), 4)]
                    fh.write(f'                     /translation="{pieces[0]}\n')
                    for piece in pieces[1:-1]:
                        fh.write(f'                     {piece}\n')
                    fh.write(f'                     {pieces[-1]}"\n' if len(pieces) > 1 else '                     "\n')
                fh.write(f'     polyA_site      {len(residues)}\n')
                fh.write('                     /note="synthetic, as a feature after the CDS"\n')
            if sequences:
                fh.write('ORIGIN\n')
                for k in range(0, len(residues), 60):
                    tens = [residues[j : j + 10].lower() for j in range(k, min(k + 60, len(residues)), 10)]
                    fh.write(f'{k + 1:>9} {" ".join(tens)}\n')
            fh.write('//\n')
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


# Each coding transcript's protein, as the GenBank record's CDS names and translates it.
RECORD_PROTEINS = {
    'NM_000010.2': ('NP_000010.1', PROTEINS['NP_000010.1']),
    'NM_000020.1': ('NP_000020.1', PROTEINS['NP_000020.1']),
    'NM_000050.1': ('NP_000050.1', PROTEINS['NP_000050.1']),
    'NM_000060.1': ('NP_000060.1', PROTEINS['NP_000060.1']),
    'NM_000070.1': ('NP_000070.1', PROTEINS['NP_000070.1']),
}
# Entrez's status of each version in the historical set: PLUS replaced, NONC suppressed, the rest live.
STATUS = [
    ['NM_000010.2', 'replaced', 'NM_000010.3'],
    ['NM_000020.1', 'live', ''],
    ['NR_000030.1', 'suppressed', ''],
    ['NM_000040.1', 'live', ''],
    ['NM_000050.1', 'live', ''],
    ['NM_000060.1', 'live', ''],
    ['NR_000061.1', 'live', ''],
    ['NM_000070.1', 'live', ''],
]
STATUS_COLUMNS = ['accession_version', 'status', 'replaced_by']


def _release(
    tmp_path: pathlib.Path,
    *,
    gff: list[str] = GFF,
    reads: list[Read] = READS,
    transcripts: dict[str, str] = TRANSCRIPTS,
    proteins: dict[str, str] = PROTEINS,
    hgnc_columns: list[str] = HGNC_COLUMNS,
    record_cds: dict[str, str] = RECORD_CDS,
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
        records=_genbank(tmp_path / 'rna.gbff.gz', transcripts, record_cds),
    )


def _historical(
    tmp_path: pathlib.Path,
    *,
    gff: list[str] = GFF,
    status: list[list[str]] = STATUS,
    status_columns: list[str] = STATUS_COLUMNS,
    record_proteins: dict[str, tuple[str, str]] = RECORD_PROTEINS,
    sequences: bool = True,
) -> refseq.Release:
    """The same set as NCBI's historical files: records with sequences and translations, a status table, no FASTA."""
    tmp_path.mkdir(exist_ok=True)
    annotation = tmp_path / 'genomic.gff.gz'
    with gzip.open(annotation, 'wt', encoding='utf-8') as fh:
        fh.write('\n'.join(gff) + '\n')
    return refseq.Release(
        assembly=bundle_pb2.ASSEMBLY_GRCH38,
        release='RS_TEST-historical',
        annotation=annotation,
        alignments=(_bam(tmp_path / 'alns.bam', READS),),
        hgnc=_tsv(tmp_path / 'hgnc.txt', HGNC_COLUMNS, HGNC, gzipped=False),
        mane=_tsv(tmp_path / 'mane.txt.gz', ['RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'], MANE, gzipped=True),
        records=_genbank(
            tmp_path / 'rna.gbff.gz', TRANSCRIPTS, RECORD_CDS, proteins=record_proteins, sequences=sequences
        ),
        status=_tsv(tmp_path / 'status.tsv', status_columns, [r[: len(status_columns)] for r in status], gzipped=False),
    )


def _by_symbol(release: refseq.Release) -> dict[str, bundle_pb2.GeneBundle]:
    return {b.gene.symbol: b for b in refseq.bundles(release)}


def _quietly(release: refseq.Release) -> dict[str, bundle_pb2.GeneBundle]:
    """The bundles by symbol, for a test about their contents rather than the warnings their build raises."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', build.BuildWarning)
        return _by_symbol(release)


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


def _alt_open(*, start: bool = False, end: bool = False) -> tuple[list[str], dict[str, str]]:
    """ALT's annotation and records with its CDS run off the record's start or end, as NCBI publishes such a model.

    The annotation's CDS row is open at the same end, so no placement states a bound there to check.
    """
    first, last = (1 if start else 3), (10 if end else 8)
    attrs = {'partial': 'true'} if start or end else {}
    if start:
        attrs['start_range'] = f'.,{first}'
    if end:
        attrs['end_range'] = f'{last},.'
    (row,) = _rows(GFF, ID='cds-NP_000050.1')
    opened = _row(ALT_LOCUS, 'CDS', first, last, '+', ID='cds-NP_000050.1', Parent='rna-NM_000050.1',
                  protein_id='NP_000050.1', **attrs)  # fmt: skip
    location = f'{"<" if start else ""}{first}..{">" if end else ""}{last}'
    return _replacing(GFF, [row], [opened]), {**RECORD_CDS, 'NM_000050.1': location}


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


def test_a_cds_both_placements_agree_with_raises_no_warning(tmp_path: pathlib.Path) -> None:
    # genome 705 is exon 1's fifth base on CHROM, as 55 is on the patch; both project to the record's 4..39
    with warnings.catch_warnings():
        warnings.simplefilter('error', build.BuildWarning)
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
def test_where_a_placement_projects_the_cds_elsewhere_the_records_is_kept(
    tmp_path: pathlib.Path, patch_end: int
) -> None:
    transcript = _transcript(_quietly(_release(tmp_path, gff=_patch_cds_end(patch_end)))['PATCHED'], 'NM_000060')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 39)


def test_a_placement_projecting_the_cds_elsewhere_warns_naming_it(tmp_path: pathlib.Path) -> None:
    expected = (
        rf'NM_000060\.1: CDS 4\.\.39 taken from its record; the annotation projects 4\.\.38 through .* on {PATCH}'
    )
    with pytest.warns(build.BuildWarning, match=expected):
        _by_symbol(_release(tmp_path, gff=_patch_cds_end(99)))


def test_a_chromosome_placement_projecting_the_cds_elsewhere_warns_too(tmp_path: pathlib.Path) -> None:
    # the chromosome states the CDS one base short, 741-749; the patch agrees with the record
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000060.1') if '\tCDS\t741\t' in r]
    short = _row(CHROM, 'CDS', 741, 749, '+', ID='cds-NP_000060.1', Parent='rna-NM_000060.1', protein_id='NP_000060.1')
    with pytest.warns(build.BuildWarning, match=rf'NM_000060\.1: .* through the alignment on {CHROM}'):
        _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [short])))


def test_a_placement_projecting_the_cds_elsewhere_keeps_its_alignment(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_quietly(_release(tmp_path, gff=_patch_cds_end(99)))['PATCHED'], 'NM_000060')
    assert [a.chromosome for a in transcript.alignments] == [CHROM, PATCH]


def test_the_records_cds_stands_where_no_aligned_placement_states_it_whole(tmp_path: pathlib.Path) -> None:
    # CLIPPED's whole CDS is stated on the patch only; without the patch alignment nothing can check it
    reads = [r for r in READS if not (r.name == 'NM_000070.1' and r.chrom == PATCH2)]
    transcript = _transcript(_by_symbol(_release(tmp_path, reads=reads))['CLIPPED'], 'NM_000070')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (0, 24)
    assert not transcript.cds.start_open
    assert not transcript.cds.end_open


def test_a_cds_starting_in_bases_a_placement_lacks_is_the_records(tmp_path: pathlib.Path) -> None:
    # CLIPPED's CDS begins in the five bases the chromosome lacks: indices 0..24, the record's
    transcript = _transcript(_by_symbol(_release(tmp_path))['CLIPPED'], 'NM_000070')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (0, 24)


def test_a_placement_open_at_its_cds_start_is_not_checked_against_the_record(tmp_path: pathlib.Path) -> None:
    # the chromosome's open row would project to index 5, the first base it has; it states no start to check
    with warnings.catch_warnings():
        warnings.simplefilter('error', build.BuildWarning)
        _by_symbol(_release(tmp_path))


def test_a_clipped_placement_is_kept_with_the_clipped_bases_as_transcript_only(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_release(tmp_path))['CLIPPED'], 'NM_000070')
    assert [(a.chromosome, _exons(a)) for a in transcript.alignments] == [
        (CHROM, [(0, 35, 900, 929, '5I30=')]),
        (PATCH2, [(0, 35, 2, 36, '35=')]),
    ]


def _broken_plus(end: int) -> list[str]:
    """PLUS's CDS rows marked partial either side of an internal break, outer ends closed, the last ending at `end`."""
    rows = [r for r in _rows(GFF, Parent='rna-NM_000010.2') if '\tCDS\t' in r]
    broken = [
        _row(CHROM, 'CDS', 105, 120, '+', ID='cds-NP_000010.1', Parent='rna-NM_000010.2', protein_id='NP_000010.1',
             partial='true', end_range='120,.'),
        _row(CHROM, 'CDS', 201, end, '+', ID='cds-NP_000010.1', Parent='rna-NM_000010.2', protein_id='NP_000010.1',
             partial='true', start_range='.,201'),
    ]  # fmt: skip
    return _replacing(GFF, rows, broken)


def test_a_cds_broken_internally_but_closed_at_both_ends_is_checked(tmp_path: pathlib.Path) -> None:
    # NCBI marks both rows either side of a frameshift partial; the outer ends stay closed, so the bounds are stated
    with pytest.warns(
        build.BuildWarning, match=r'NM_000010\.2: CDS 4\.\.29 taken from its record; the annotation projects 4\.\.28'
    ):
        _by_symbol(_release(tmp_path, gff=_broken_plus(209)))


def test_a_minus_strand_cds_open_at_its_high_genome_end_is_not_checked(tmp_path: pathlib.Path) -> None:
    # MINUS reads 5' to 3' downward on the genome, so an open start is `end_range` on its highest row, 601-605;
    # the row here ends at 603, which a check would find two bases short
    (row,) = [r for r in _rows(GFF, Parent='rna-NM_000020.1') if '\tCDS\t601\t' in r]
    open_5 = _row(CHROM, 'CDS', 601, 603, '-', ID='cds-NP_000020.1', Parent='rna-NM_000020.1',
                  protein_id='NP_000020.1', partial='true', end_range='603,.')  # fmt: skip
    with warnings.catch_warnings():
        warnings.simplefilter('error', build.BuildWarning)
        (transcript,) = _by_symbol(_release(tmp_path, gff=_replacing(GFF, [row], [open_5])))['MINUS'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (26, 46)


@pytest.mark.parametrize(
    ('start', 'end', 'bounds'),
    [(True, False, (0, 7)), (False, True, (2, 9)), (True, True, (0, 9))],
    ids=['open start', 'open end', 'open both'],
)
def test_a_record_cds_open_at_an_end_keeps_its_bounds_and_says_which_end(
    tmp_path: pathlib.Path, start: bool, end: bool, bounds: tuple[int, int]
) -> None:
    gff, record_cds = _alt_open(start=start, end=end)
    (transcript,) = _by_symbol(_release(tmp_path, gff=gff, record_cds=record_cds))['ALT'].transcripts
    cds = transcript.cds
    assert (cds.start_index, cds.end_index_inclusive, cds.start_open, cds.end_open) == (*bounds, start, end)


@pytest.mark.parametrize(
    ('location', 'error'),
    [('<3..8', r'open at its start, which is index 2, not 0'), ('3..>8', r'open at its end, which is index 7, not 9')],
    ids=['open start', 'open end'],
)
def test_an_open_cds_end_short_of_the_records_edge_fails_the_build(
    tmp_path: pathlib.Path, location: str, error: str
) -> None:
    # bases beyond an open end would be coding by the record yet numbered from nothing
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NM_000050.1': location})
    with pytest.raises(build.BuildError, match=rf'NM_000050\.1: CDS {error}'):
        list(refseq.bundles(release))


def test_a_transcript_with_an_open_cds_keeps_its_protein(tmp_path: pathlib.Path) -> None:
    gff, record_cds = _alt_open(start=True)
    (transcript,) = _by_symbol(_release(tmp_path, gff=gff, record_cds=record_cds))['ALT'].transcripts
    assert (transcript.protein_accession, transcript.protein_version) == ('NP_000050', 1)


def test_the_report_counts_the_records_that_leave_a_cds_end_open(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    gff, record_cds = _alt_open(start=True)
    list(refseq.bundles(_release(tmp_path, gff=gff, record_cds=record_cds)))
    assert "coding transcripts whose record leaves a CDS end open: 1 (5' 1);" in capsys.readouterr().err


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


def test_an_annotation_synonym_hgnc_lists_as_a_previous_symbol_is_not_an_alias_too(tmp_path: pathlib.Path) -> None:
    gff = [line.replace('gene_synonym=PLS', 'gene_synonym=PLS,OLDPLUS') for line in GFF]
    gene = _by_symbol(_release(tmp_path, gff=gff))['PLUS'].gene
    assert 'OLDPLUS' not in gene.alias_symbols


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


def _deleting_plus(tmp_path: pathlib.Path) -> refseq.Release:
    """PLUS with its record lacking genome 105, where the annotation starts the CDS, so the alignment deletes it."""
    record = _g(CHROM, 101, 104) + _g(CHROM, 106, 120) + _g(CHROM, 201, 230)
    deleting = Read('NM_000010.2', CHROM, 100, '4=1D15=80N30=', False, record)
    return _release(
        tmp_path,
        reads=[deleting if r.name == 'NM_000010.2' else r for r in READS],
        transcripts={**TRANSCRIPTS, 'NM_000010.2': record},
        record_cds={**RECORD_CDS, 'NM_000010.2': '5..29'},
    )


def test_a_cds_the_annotation_bounds_on_a_genome_only_base_is_the_records(tmp_path: pathlib.Path) -> None:
    (transcript,) = _quietly(_deleting_plus(tmp_path))['PLUS'].transcripts
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 28)


def test_an_annotation_cds_that_does_not_project_warns(tmp_path: pathlib.Path) -> None:
    with pytest.warns(build.BuildWarning, match='does not project: .*genome position 104 is a genome-only base'):
        _by_symbol(_deleting_plus(tmp_path))


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
        '--hgnc', str(release.hgnc), '--mane', str(release.mane), '--records', str(release.records),
        '--shards', str(tmp_path / 'shards'),
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


def _provider(
    tmp_path: pathlib.Path, *, gff: list[str] = GFF, record_cds: dict[str, str] = RECORD_CDS
) -> provider_mod.BundleProvider:
    """The synthetic release built, indexed and cut, and opened as weaver's provider."""
    shard = store_build.write_shard(
        refseq.bundles(_release(tmp_path, gff=gff, record_cds=record_cds)),
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


def _open_alt(tmp_path: pathlib.Path, *, start: bool = False, end: bool = False) -> provider_mod.BundleProvider:
    """The provider over a release whose ALT CDS runs off its record at the ends named."""
    gff, record_cds = _alt_open(start=start, end=end)
    return _provider(tmp_path, gff=gff, record_cds=record_cds)


@pytest.mark.parametrize(
    ('start', 'end', 'bounds'),
    [(False, False, (2, 7)), (True, False, (0, 7)), (False, True, (2, 9))],
    ids=['closed', 'open start', 'open end'],
)
def test_the_model_says_which_cds_ends_are_open(
    tmp_path: pathlib.Path, start: bool, end: bool, bounds: tuple[int, int]
) -> None:
    model = _open_alt(tmp_path, start=start, end=end).get_transcript('NM_000050.1', None)
    assert (model['cds_start_index'], model['cds_end_index']) == bounds
    assert (model.get('cds_start_open'), model.get('cds_end_open')) == (start, end)


def test_a_noncoding_model_has_no_open_cds_end(tmp_path: pathlib.Path) -> None:
    model = _provider(tmp_path).get_transcript('NR_000061.1', None)
    assert (model['cds_start_index'], model.get('cds_start_open'), model.get('cds_end_open')) == (None, False, False)


@pytest.mark.parametrize(
    ('open_end', 'variant', 'end'),
    [({'start': True}, 'NM_000050.1:c.1A>C', "5'"), ({'end': True}, 'NM_000050.1:c.*1A>C', "3'")],
    ids=['c. from an open start', 'c.* from an open end'],
)
def test_a_position_numbered_from_an_open_cds_end_is_refused(
    tmp_path: pathlib.Path, open_end: dict[str, bool], variant: str, end: str
) -> None:
    mapper = weaver.VariantMapper(_open_alt(tmp_path, **open_end))
    with pytest.raises(weaver.ValidationError, match=f'open at the {end} end'):
        mapper.c_to_g(weaver.parse(variant), None)


@pytest.mark.parametrize(
    ('open_end', 'variant', 'projected'),
    [({'end': True}, 'NM_000050.1:c.1G>C', 'g.3G>C'), ({'start': True}, 'NM_000050.1:c.*1A>C', 'g.9A>C')],
    ids=['c. on an open end', 'c.* on an open start'],
)
def test_a_position_numbered_from_the_closed_end_of_an_open_cds_projects(
    tmp_path: pathlib.Path, open_end: dict[str, bool], variant: str, projected: str
) -> None:
    mapped = weaver.VariantMapper(_open_alt(tmp_path, **open_end)).c_to_g(weaver.parse(variant), None)
    assert mapped.format() == f'{ALT_LOCUS}:{projected}'


# ---- the GenBank records -------------------------------------------------------------------------------


def test_a_transcript_coding_by_the_annotation_but_not_its_record_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, record_cds={k: v for k, v in RECORD_CDS.items() if k != 'NM_000010.2'})
    with pytest.raises(build.BuildError, match=r'state none in their record.*NM_000010\.2'):
        list(refseq.bundles(release))


def test_a_transcript_coding_only_by_its_record_must_still_be_aligned(tmp_path: pathlib.Path) -> None:
    # NONC has no alignment; a CDS in its record makes it coding, and a coding transcript needs a placement
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NR_000030.1': '5..40'})
    with pytest.raises(build.BuildError, match=r'no alignment.*NR_000030\.1'):
        list(refseq.bundles(release))


def test_a_record_cds_where_the_annotation_gives_none_warns(tmp_path: pathlib.Path) -> None:
    # PAR is aligned and the annotation calls it non-coding, but its record states a CDS
    with pytest.warns(
        build.BuildWarning, match=r'NM_000040\.1: CDS 4\.\.39 taken from its record; the annotation gives it no CDS'
    ):
        _by_symbol(_release(tmp_path, record_cds={**RECORD_CDS, 'NM_000040.1': '5..40'}))


def test_a_record_cds_outside_its_sequence_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': '5..300'})
    with pytest.raises(build.BuildError, match=r'NM_000010\.2: CDS 4\.\.299 outside its 50-base sequence'):
        list(refseq.bundles(release))


def test_a_closed_cds_read_from_other_than_its_first_base_fails_the_build(tmp_path: pathlib.Path) -> None:
    # a whole CDS starts with its start codon, so it reads from its first base
    shifted = '5..30\n                     /codon_start=2'
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': shifted})
    with pytest.raises(build.BuildError, match=r'NM_000010\.2: a CDS with a closed start but /codon_start=2'):
        list(refseq.bundles(release))


def test_an_open_cds_may_read_from_its_second_base(tmp_path: pathlib.Path) -> None:
    gff, record_cds = _alt_open(start=True)
    record_cds['NM_000050.1'] += '\n                     /codon_start=2'
    (transcript,) = _by_symbol(_release(tmp_path, gff=gff, record_cds=record_cds))['ALT'].transcripts
    assert transcript.cds.start_open


def test_a_frameshifted_cds_is_taken_by_its_outer_bounds(tmp_path: pathlib.Path) -> None:
    # a programmed frameshift is written as a join that skips a base; c. numbering needs the outer bounds
    cds = _quietly(_release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': 'join(5..17,19..30)'}))['PLUS']
    assert (cds.transcripts[0].cds.start_index, cds.transcripts[0].cds.end_index_inclusive) == (4, 29)


def test_a_frameshifted_cds_wrapped_onto_a_second_line_is_read_whole(tmp_path: pathlib.Path) -> None:
    wrapped = 'join(5..17,\n                     19..30)'
    transcript = _quietly(_release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': wrapped}))['PLUS'].transcripts[0]
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (4, 29)


def test_a_record_stating_two_cds_features_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': '5..30\n     CDS             7..30'})
    with pytest.raises(build.BuildError, match=r'NM_000010\.2 states 2 CDS features, 0 of them naming a protein'):
        list(refseq.bundles(release))


def test_a_cds_in_a_record_with_no_version_line_fails_the_build(tmp_path: pathlib.Path) -> None:
    # otherwise the CDS would be credited to the record before it
    release = _release(tmp_path)
    text = gzip.decompress(release.records.read_bytes()).decode('ascii')
    release.records.write_bytes(gzip.compress(text.replace('VERSION     NM_000020.1\n', '', 1).encode('ascii')))
    with pytest.raises(build.BuildError, match='a CDS in a record with no VERSION line'):
        list(refseq.bundles(release))


@pytest.mark.parametrize(
    'location',
    [
        'complement(5..30)',  # a minus-strand CDS on a transcript record
        'order(5..17,19..30)',  # not one contiguous reading
        'join(NM_000001.2:5..17,19..30)',  # a range on another record, whose digits are no bound of this one
        '5^6',  # a site between two bases
        '5',  # a single base
        '30..5',  # backwards
        'join(19..30,5..17)',  # a join whose ranges run backwards
        'join(5..>17,19..30)',  # open at an inner bound
    ],
)
def test_a_cds_location_of_another_shape_fails_the_build(tmp_path: pathlib.Path, location: str) -> None:
    release = _release(tmp_path, record_cds={**RECORD_CDS, 'NM_000010.2': location})
    with pytest.raises(build.BuildError, match=r'NM_000010\.2: CDS location'):
        list(refseq.bundles(release))


def test_a_file_ending_inside_a_cds_location_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path)
    text = gzip.decompress(release.records.read_bytes()).decode('ascii')
    cut = text.index('     CDS             5..30')
    release.records.write_bytes(gzip.compress((text[:cut] + '     CDS             join(5..17,\n').encode('ascii')))
    with pytest.raises(build.BuildError, match=r'ends inside the CDS location of NM_000010\.2'):
        list(refseq.bundles(release))


def test_a_version_line_with_no_accession_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path)
    text = gzip.decompress(release.records.read_bytes()).decode('ascii')
    blanked = text.replace('VERSION     NM_000020.1\n', 'VERSION\n')
    release.records.write_bytes(gzip.compress(blanked.encode('ascii')))
    with pytest.raises(build.BuildError, match='a VERSION line with no accession'):
        list(refseq.bundles(release))


def test_the_shard_names_the_records_it_was_cut_from(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path)
    shard = store_build.write_shard(
        refseq.bundles(release), tmp_path / 'shards', release='RS_TEST', inputs=release.inputs()
    )
    assert 'records' in {i.role for i in store_build.shard_record(shard).inputs}


# ---- the command's warning report ------------------------------------------------------------------------


def _refseq_command(release: refseq.Release, shards: pathlib.Path) -> list[str]:
    return [
        'refseq', '--assembly', 'GRCh38', '--release', 'RS_TEST',
        '--annotation', str(release.annotation), '--transcripts', str(release.transcripts),
        '--proteins', str(release.proteins), '--alignments', str(release.alignments[0]),
        '--hgnc', str(release.hgnc), '--mane', str(release.mane), '--records', str(release.records),
        '--shards', str(shards),
    ]  # fmt: skip


def test_the_command_reports_each_build_warning_on_a_line(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # the command reports build warnings whatever filters the caller runs it under
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        cli.main(_refseq_command(_release(tmp_path, gff=_patch_cds_end(99)), tmp_path / 'shards'))
    report = capsys.readouterr().err.splitlines()
    assert report[-2].startswith('warning: NM_000060.1: CDS 4..39 taken from its record; the annotation projects 4..38')
    assert report[-1] == '1 build warning'


def test_a_crashed_build_still_reports_its_warnings(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    # stands in for a defect or an I/O error during the build, which is not a BuildError
    write_shard = store_build.write_shard

    def crash(*args: object, **kwargs: object) -> pathlib.Path:
        write_shard(*args, **kwargs)  # type: ignore[arg-type]
        raise OSError('disk full')

    monkeypatch.setattr(store_build, 'write_shard', crash)
    with pytest.raises(OSError, match='disk full'):
        cli.main(_refseq_command(_release(tmp_path, gff=_patch_cds_end(99)), tmp_path / 'shards'))
    report = capsys.readouterr().err.splitlines()
    assert report[-2].startswith('warning: NM_000060.1: CDS 4..39 taken from its record')
    assert report[-1] == '1 build warning'


def test_a_failed_build_reports_its_warnings_before_the_failure(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # genes are bundled in the annotation's order: PATCHED's warning is raised before CLIPPED's missing protein
    proteins = {k: v for k, v in PROTEINS.items() if k != 'NP_000070.1'}
    with pytest.raises(SystemExit):
        cli.main(_refseq_command(_release(tmp_path, gff=_patch_cds_end(99), proteins=proteins), tmp_path / 'shards'))
    report = capsys.readouterr().err.splitlines()
    assert report[-3].startswith('warning: NM_000060.1: CDS 4..39 taken from its record')
    assert report[-2] == '1 build warning'
    assert report[-1] == 'FAILED: NM_000070.1: its protein NP_000070.1 is not in the protein set'


def test_a_warning_of_another_kind_is_shown_and_the_command_finishes(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    # stands in for a library warning raised during the build; showing it while recording hung the command
    write_shard = store_build.write_shard

    def warning_first(*args: object, **kwargs: object) -> pathlib.Path:
        warnings.warn('a library warning', RuntimeWarning, stacklevel=1)
        return write_shard(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(store_build, 'write_shard', warning_first)
    with pytest.warns(RuntimeWarning, match='a library warning'):  # re-shown once the build is over
        cli.main(_refseq_command(_release(tmp_path), tmp_path / 'shards'))
    assert capsys.readouterr().out.strip().endswith('.bagz')  # the command finished and printed its shard


# ---- the historical set ---------------------------------------------------------------------------


def test_the_historical_set_takes_sequences_and_proteins_from_the_records(tmp_path: pathlib.Path) -> None:
    bundles = _by_symbol(_historical(tmp_path))
    transcript = _transcript(bundles['MINUS'], 'NM_000020')
    sequences = {s.digest: bytes(s.residues) for s in bundles['MINUS'].sequences}
    assert sequences[transcript.sequence_digest] == TRANSCRIPTS['NM_000020.1'].encode()
    (protein,) = bundles['MINUS'].proteins
    assert (protein.accession, protein.version) == ('NP_000020', 1)
    assert sequences[protein.sequence_digest] == PROTEINS['NP_000020.1'].encode()
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (26, 46)


def test_each_version_carries_the_status_entrez_gives_it(tmp_path: pathlib.Path) -> None:
    rows = [
        ['NM_000010.2', 'replaced', 'NM_000010.3'],
        ['NR_000030.1', 'suppressed', ''],
        *[r for r in STATUS if r[0] not in ('NM_000010.2', 'NR_000030.1')],
    ]
    bundles = _by_symbol(_historical(tmp_path, status=rows))
    assert _transcript(bundles['PLUS'], 'NM_000010').status == bundle_pb2.TRANSCRIPT_STATUS_SUPERSEDED
    assert _transcript(bundles['NONC'], 'NR_000030').status == bundle_pb2.TRANSCRIPT_STATUS_SUPPRESSED
    assert _transcript(bundles['MINUS'], 'NM_000020').status == bundle_pb2.TRANSCRIPT_STATUS_CURRENT


def test_a_retired_version_carries_no_tags_and_no_mane_partner(tmp_path: pathlib.Path) -> None:
    # PLUS is tagged RefSeq Select in the annotation and MANE Select in the summary; replaced, it keeps neither
    transcript = _transcript(_by_symbol(_historical(tmp_path))['PLUS'], 'NM_000010')
    assert list(transcript.tags) == []
    assert transcript.mane_partner == ''


def test_a_version_with_no_status_row_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _historical(tmp_path, status=[r for r in STATUS if r[0] != 'NM_000040.1'])
    with pytest.raises(build.BuildError, match=r'no row in the status table.*NM_000040\.1'):
        list(refseq.bundles(release))


def test_a_status_word_that_is_not_one_of_entrez_s_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _historical(tmp_path, status=[['NM_000010.2', 'retired', ''], *STATUS[1:]])
    with pytest.raises(build.BuildError, match=r"NM_000010\.2 has status 'retired'"):
        list(refseq.bundles(release))


def test_a_status_table_without_its_columns_fails_the_build(tmp_path: pathlib.Path) -> None:
    release = _historical(tmp_path, status_columns=['accession_version'])
    with pytest.raises(build.BuildError, match="'status'"):
        list(refseq.bundles(release))


def test_a_record_without_a_translation_fails_a_build_reading_sequences_from_records(tmp_path: pathlib.Path) -> None:
    proteins = {k: v for k, v in RECORD_PROTEINS.items() if k != 'NM_000050.1'}
    with pytest.raises(build.BuildError, match=r'NM_000050\.1 names no protein_id or carries no translation'):
        list(refseq.bundles(_historical(tmp_path, record_proteins=proteins)))


def test_a_record_without_an_origin_fails_a_build_reading_sequences_from_records(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match='no ORIGIN section'):
        list(refseq.bundles(_historical(tmp_path, sequences=False)))


def test_a_release_with_one_sequence_set_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='sequence sets'):
        dataclasses.replace(_release(tmp_path), proteins=None)


def test_a_release_with_a_status_table_beside_its_sequence_sets_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='status table'):
        dataclasses.replace(_release(tmp_path), status=tmp_path / 'status.tsv')


def test_a_historical_shard_under_the_current_one_serves_a_retired_version(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Built by the command and stacked under the current shard, a retired version is served when named only."""
    historical = _historical(tmp_path / 'historical')
    cli.main([
        'historical', '--assembly', 'GRCh38', '--release', 'RS_TEST-historical',
        '--annotation', str(historical.annotation), '--records', str(historical.records),
        '--alignments', str(historical.alignments[0]), '--status', str(historical.status),
        '--hgnc', str(historical.hgnc), '--mane', str(historical.mane), '--shards', str(tmp_path / 'shards'),
    ])  # fmt: skip
    historical_shard = capsys.readouterr().out.strip()
    # the current release no longer carries PLUS's NM_000010.2 or NONC's NR_000030.1
    retired = ('NM_000010.2', 'NR_000030.1')
    (tmp_path / 'current').mkdir()
    current = _release(
        tmp_path / 'current',
        gff=[row for row in GFF if not any(r in row for r in retired)],
        reads=[r for r in READS if r.name not in retired],
        transcripts={k: v for k, v in TRANSCRIPTS.items() if k not in retired},
        record_cds={k: v for k, v in RECORD_CDS.items() if k not in retired},
    )
    current_shard = store_build.write_shard(
        refseq.bundles(current), tmp_path / 'shards', release='RS_TEST', inputs=current.inputs()
    )
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), historical_shard, str(current_shard)])
    fasta = _fasta(tmp_path / 'genome.fna.gz', GENOME)
    cli.main(['genome', '--assembly', 'GRCh38', '--fasta', str(fasta), '--out', str(tmp_path / 'genome')])
    provider = provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(tmp_path / 'genome'))
    )
    assert provider.get_transcript('NM_000010.2', None)['cds_start_index'] == 4
    ref = _g(CHROM, 105, 105)  # PLUS's c.1
    assert weaver.parse(f'NM_000010.2:c.1{ref}>{"G" if ref != "G" else "C"}').validate(provider)
    assert provider.get_transcripts_for_region(CHROM, 100, 230) == []  # PLUS's only transcript is retired
    assert provider.get_transcripts_for_region(CHROM, 800, 850) == []  # NONC's is suppressed
    assert provider.get_transcripts_for_region(CHROM, 500, 630) == ['NM_000020.1']


def _with_second_cds(tmp_path: pathlib.Path, *, named: bool) -> refseq.Release:
    """The historical set with a second CDS feature on ALT's record, naming a protein or not."""
    release = _historical(tmp_path)
    extra = '     CDS             3..5\n                     /codon_start=1\n'
    if named:
        extra += '                     /protein_id="NP_000051.1"\n                     /translation="M"\n'
    with gzip.open(release.records, 'rt', encoding='ascii') as fh:
        text = fh.read()
    marker = 'ORIGIN'
    head, _, tail = text.partition('VERSION     NM_000050.1\n')
    tail = tail.replace(marker, extra + marker, 1)
    with gzip.open(release.records, 'wt', encoding='ascii') as fh:
        fh.write(head + 'VERSION     NM_000050.1\n' + tail)
    return release


def test_a_record_with_a_second_cds_naming_no_protein_takes_the_one_that_does(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_symbol(_with_second_cds(tmp_path, named=False))['ALT'], 'NM_000050')
    assert (transcript.cds.start_index, transcript.cds.end_index_inclusive) == (2, 7)
    assert (transcript.protein_accession, transcript.protein_version) == ('NP_000050', 1)


def test_a_record_with_two_cds_features_naming_proteins_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match=r'NM_000050\.1 states 2 CDS features, 2 of them naming a protein'):
        list(refseq.bundles(_with_second_cds(tmp_path, named=True)))


def _complemented(tmp_path: pathlib.Path, historical: bool) -> refseq.Release:
    """ALT's record stating its CDS as a `complement` location, in the historical set or in a release."""
    cds = {**RECORD_CDS, 'NM_000050.1': 'complement(3..8)'}
    if not historical:
        return _release(tmp_path, record_cds=cds)
    release = _historical(tmp_path)
    _genbank(release.records, TRANSCRIPTS, cds, proteins=RECORD_PROTEINS, sequences=True)
    return release


def test_the_historical_set_leaves_out_a_record_whose_cds_statement_the_reader_refuses(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    bundles = _by_symbol(_complemented(tmp_path, historical=True))
    assert 'ALT' not in bundles  # its only transcript is left out, so the gene has no bundle
    assert 'NM_000040' in [t.accession for t in bundles['PAR'].transcripts]  # the rest are bundled
    report = capsys.readouterr().err
    assert 'left out for a CDS statement the reader refuses: 1' in report
    assert 'NM_000050.1: CDS location' in report


def test_a_release_record_whose_cds_statement_the_reader_refuses_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match=r'NM_000050\.1: CDS location .* is not a range'):
        list(refseq.bundles(_complemented(tmp_path, historical=False)))


def test_two_records_translating_one_protein_version_differently_fail_the_build(tmp_path: pathlib.Path) -> None:
    # ALT's record names PLUS's protein and translates it to something else
    proteins = {**RECORD_PROTEINS, 'NM_000050.1': ('NP_000010.1', 'MA')}
    with pytest.raises(build.BuildError, match=r'NM_000050\.1: translates NP_000010\.1 to other residues'):
        list(refseq.bundles(_historical(tmp_path, record_proteins=proteins)))


def _rewritten_records(release: refseq.Release, edit: Callable[[str], str]) -> refseq.Release:
    """The release with its GenBank file's text put through `edit`."""
    with gzip.open(release.records, 'rt', encoding='ascii') as fh:
        text = fh.read()
    with gzip.open(release.records, 'wt', encoding='ascii') as fh:
        fh.write(edit(text))
    return release


def test_a_record_whose_origin_has_only_blank_lines_fails_a_build_reading_sequences(tmp_path: pathlib.Path) -> None:
    def blank_origin(text: str) -> str:
        head, _, tail = text.partition('VERSION     NR_000030.1\n')
        body, _, rest = tail.partition('//\n')
        return head + 'VERSION     NR_000030.1\n' + body[: body.index('ORIGIN')] + 'ORIGIN\n\n//\n' + rest

    with pytest.raises(build.BuildError, match=r'NR_000030\.1 has no ORIGIN section, or an empty one'):
        list(refseq.bundles(_rewritten_records(_historical(tmp_path), blank_origin)))


def test_a_record_not_closed_before_the_next_fails_the_build(tmp_path: pathlib.Path) -> None:
    # without its //, NM_000010.2's ORIGIN would swallow the record after it
    def unclosed(text: str) -> str:
        head, _, tail = text.partition('VERSION     NM_000010.2\n')
        return head + 'VERSION     NM_000010.2\n' + tail.replace('//\n', '', 1)

    with pytest.raises(build.BuildError, match=r'NM_000010\.2 is not closed by // before the next record'):
        list(refseq.bundles(_rewritten_records(_historical(tmp_path), unclosed)))


def test_a_file_ending_inside_a_record_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match=r'ends inside the record of NM_000070\.1'):
        list(refseq.bundles(_rewritten_records(_release(tmp_path), lambda text: text[: text.rindex('//\n')])))


def test_the_historical_set_leaves_out_a_record_whose_open_cds_end_is_not_its_last_base(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ALT's record, 10 bases, marks its CDS open at an end two bases short of the record's
    release = _historical(tmp_path)
    cds = {**RECORD_CDS, 'NM_000050.1': '3..>8'}
    _genbank(release.records, TRANSCRIPTS, cds, proteins=RECORD_PROTEINS, sequences=True)
    assert 'ALT' not in _by_symbol(release)
    assert 'NM_000050.1: CDS open at its end, which is index 7, not 9' in capsys.readouterr().err


def test_the_protein_is_the_records_where_the_annotation_names_another(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # NCBI's historical GFF3 names another version's protein on a quarter of its CDS rows; the record is the statement.
    # ALT's rows name PLUS's protein, which the set holds with other residues, so nothing would fail by itself.
    gff = [row.replace('protein_id=NP_000050.1', 'protein_id=NP_000010.1') for row in GFF]
    proteins = {**RECORD_PROTEINS, 'NM_000050.1': ('NP_000050.1', 'MA')}
    bundles = _by_symbol(_historical(tmp_path, gff=gff, record_proteins=proteins))
    transcript = _transcript(bundles['ALT'], 'NM_000050')
    assert (transcript.protein_accession, transcript.protein_version) == ('NP_000050', 1)
    (protein,) = bundles['ALT'].proteins
    residues = {s.digest: bytes(s.residues) for s in bundles['ALT'].sequences}
    assert (protein.accession, protein.version, residues[protein.sequence_digest]) == ('NP_000050', 1, b'MA')
    assert "names another protein than their record, the record's bundled: 1" in capsys.readouterr().err


def test_a_record_protein_id_that_is_not_versioned_fails_the_build_naming_the_record(tmp_path: pathlib.Path) -> None:
    proteins = {**RECORD_PROTEINS, 'NM_000050.1': ('NP_000050', 'MA')}
    with pytest.raises(build.BuildError, match=r"NM_000050\.1: protein_id 'NP_000050' is not a versioned accession"):
        list(refseq.bundles(_historical(tmp_path, record_proteins=proteins)))
