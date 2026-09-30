"""The Ensembl builder, from files in Ensembl's formats to bundles, and through the command to weaver.

The inputs are written here as Ensembl and NCBI publish them — the GFF3, FASTA and TSV as text, the
assembly report naming Ensembl's chromosome `99` as NC_000099.1 — so no fixture comes from the builder
itself. The synthetic release:

- FWD, ENSG00000000001: ENST00000000001.2 on the plus strand, exons 101-120 and 201-230, CDS 105-210,
  MANE Select paired with NM_000010.2; and ENST00000000008.1, a retained intron at 101-120, listed
  first in the annotation.
- REV, ENSG00000000002: ENST00000000002.1 on the minus strand, exons 501-520 and 601-630, CDS 505-605,
  protein ENSP00000000002.3; its record differs from the genome at one base, its sixth, as Ensembl's
  LYPD8-207 does on GRCh38. HGNC has two rows for the gene, REV and REVB.
- ENSG00000000003, a lncRNA gene the annotation names nothing, with one transcript at 801-850, whose
  record is in the ncRNA set, and a TEC model at 861-870 that has no record.
- ENST00000000099.1, a MANE record whose transcript Ensembl places off the primary chromosomes.
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib
import random
import warnings

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

CHROM = 'NC_000099.1'
GENOME = ''.join(random.Random(7).choices('ACGT', k=1000))
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def _g(start: int, end: int) -> str:
    """Genome bases start..end, 1-based closed."""
    return GENOME[start - 1 : end]


def _revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


_REV_SPLICED = _revcomp(_g(501, 520) + _g(601, 630))
RECORDS = {
    'ENST00000000001.2': _g(101, 120) + _g(201, 230),
    'ENST00000000008.1': _g(101, 120),
    'ENST00000000002.1': _REV_SPLICED[:5] + ('A' if _REV_SPLICED[5] != 'A' else 'C') + _REV_SPLICED[6:],
    'ENST00000000003.1': _g(801, 850),
    'ENST00000000099.1': 'ACGTACGTACGT',
}
NCRNA = frozenset({'ENST00000000003.1'})  # written to the ncRNA set, the rest to the cDNA set
PROTEINS = {'ENSP00000000001.2': 'MSEQFWD', 'ENSP00000000002.3': 'MSEQREV'}


def _row(kind: str, start: int, end: int, strand: str, phase: str = '.', **attrs: str) -> str:
    """A GFF3 line; a CDS states its phase, 0 unless given, and every other feature `.`."""
    if kind == 'CDS' and phase == '.':
        phase = '0'
    fields = ['99', 'ensembl', kind, str(start), str(end), '.', strand, phase]
    return '\t'.join([*fields, ';'.join(f'{k}={v}' for k, v in attrs.items())])


GFF = [
    '##gff-version 3',
    _row('gene', 101, 230, '+', ID='gene:ENSG00000000001', Name='FWD', biotype='protein_coding',
         description='forward gene [Source:HGNC Symbol%3BAcc:HGNC:1]', gene_id='ENSG00000000001', version='3'),
    _row('transcript', 101, 120, '+', ID='transcript:ENST00000000008', Parent='gene:ENSG00000000001',
         biotype='retained_intron', transcript_id='ENST00000000008', version='1'),
    _row('exon', 101, 120, '+', Parent='transcript:ENST00000000008', exon_id='ENSE00000000009', rank='1', version='1'),
    _row('mRNA', 101, 230, '+', ID='transcript:ENST00000000001', Parent='gene:ENSG00000000001', Name='FWD-201',
         biotype='protein_coding', tag='gencode_basic,Ensembl_canonical,MANE_Select',
         transcript_id='ENST00000000001', version='2'),
    _row('exon', 101, 120, '+', Parent='transcript:ENST00000000001', exon_id='ENSE00000000001', rank='1', version='1'),
    _row('CDS', 105, 120, '+', ID='CDS:ENSP00000000001', Parent='transcript:ENST00000000001',
         protein_id='ENSP00000000001', version='2'),
    _row('exon', 201, 230, '+', Parent='transcript:ENST00000000001', exon_id='ENSE00000000002', rank='2', version='1'),
    _row('CDS', 201, 210, '+', ID='CDS:ENSP00000000001', Parent='transcript:ENST00000000001',
         protein_id='ENSP00000000001', version='2'),
    _row('gene', 501, 630, '-', ID='gene:ENSG00000000002', Name='REV', biotype='protein_coding',
         description='reverse gene', gene_id='ENSG00000000002', version='1'),
    _row('mRNA', 501, 630, '-', ID='transcript:ENST00000000002', Parent='gene:ENSG00000000002', Name='REV-201',
         biotype='protein_coding', transcript_id='ENST00000000002', version='1'),
    _row('exon', 601, 630, '-', Parent='transcript:ENST00000000002', exon_id='ENSE00000000003', rank='1', version='1'),
    _row('CDS', 601, 605, '-', ID='CDS:ENSP00000000002', Parent='transcript:ENST00000000002',
         protein_id='ENSP00000000002', version='3'),
    _row('exon', 501, 520, '-', Parent='transcript:ENST00000000002', exon_id='ENSE00000000004', rank='2', version='1'),
    _row('CDS', 505, 520, '-', ID='CDS:ENSP00000000002', Parent='transcript:ENST00000000002',
         protein_id='ENSP00000000002', version='3'),
    _row('ncRNA_gene', 801, 870, '+', ID='gene:ENSG00000000003', biotype='lncRNA',
         description='long%2C non-coding', gene_id='ENSG00000000003', version='1'),
    _row('lnc_RNA', 801, 850, '+', ID='transcript:ENST00000000003', Parent='gene:ENSG00000000003', biotype='lncRNA',
         transcript_id='ENST00000000003', version='1'),
    _row('exon', 801, 850, '+', Parent='transcript:ENST00000000003', exon_id='ENSE00000000005', rank='1', version='1'),
    _row('unconfirmed_transcript', 861, 870, '+', ID='transcript:ENST00000000004', Parent='gene:ENSG00000000003',
         biotype='TEC', transcript_id='ENST00000000004', version='1'),
    _row('exon', 861, 870, '+', Parent='transcript:ENST00000000004', exon_id='ENSE00000000006', rank='1', version='1'),
]  # fmt: skip
HGNC_COLUMNS = ['hgnc_id', 'symbol', 'name', 'alias_symbol', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']
HGNC = [
    ['HGNC:1', 'FWD', 'forward gene, HGNC', 'FW1', 'OLDFWD', '11', 'ENSG00000000001'],
    ['HGNC:3', 'REVB', 'reverse gene B', '', 'OLDREVB', '13', 'ENSG00000000002'],
    ['HGNC:4', 'REV', 'reverse gene, HGNC', '', '', '12', 'ENSG00000000002'],
]
MANE = [
    ['NM_000010.2', 'ENST00000000001.2', 'MANE Select'],
    ['NM_000099.1', 'ENST00000000099.1', 'MANE Select'],
]
REPORT = [
    '# Assembly name:  SYNTHETIC',
    '# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\t'
    'Relationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name',
    f'99\tassembled-molecule\t99\tChromosome\tCM000099.1\t=\t{CHROM}\tPrimary Assembly\t1000\tchr99',
]


@dataclasses.dataclass(frozen=True)
class Inputs:
    """Each input a keyword can replace, as text before it is written."""

    gff: tuple[str, ...] = tuple(GFF)
    records: tuple[tuple[str, str], ...] = tuple(RECORDS.items())
    proteins: tuple[tuple[str, str], ...] = tuple(PROTEINS.items())
    report: tuple[str, ...] = tuple(REPORT)
    genome_assembly: bundle_pb2.Assembly = bundle_pb2.ASSEMBLY_GRCH38


def _fasta(path: pathlib.Path, records: dict[str, str]) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        for name, residues in records.items():
            fh.write(f'>{name} synthetic\n')
            for k in range(0, len(residues), 60):
                fh.write(residues[k : k + 60] + '\n')
    return path


def _tsv(path: pathlib.Path, rows: list[list[str]], *, gzipped: bool) -> pathlib.Path:
    text = '\n'.join('\t'.join(r) for r in rows) + '\n'
    if gzipped:
        with gzip.open(path, 'wt', encoding='utf-8') as fh:
            fh.write(text)
    else:
        path.write_text(text, 'utf-8')
    return path


def _release(tmp_path: pathlib.Path, inputs: Inputs = Inputs()) -> ensembl.Release:  # noqa: B008
    """The synthetic release written under `tmp_path`, its genome built from the synthetic chromosome."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    annotation = tmp_path / 'annotation.gff3.gz'
    with gzip.open(annotation, 'wt', encoding='utf-8') as fh:
        fh.write('\n'.join(inputs.gff) + '\n')
    (tmp_path / 'report.txt').write_text('\n'.join(inputs.report) + '\n', 'utf-8')
    genome_fasta = _fasta(tmp_path / 'genome.fna.gz', {CHROM: GENOME})
    genome_build.build_genome(genome_fasta, tmp_path / 'genome', assembly=inputs.genome_assembly)
    records = dict(inputs.records)
    return ensembl.Release(
        assembly=bundle_pb2.ASSEMBLY_GRCH38,
        release='ensembl_TEST',
        annotation=annotation,
        transcripts=(
            _fasta(tmp_path / 'cdna.fa.gz', {k: v for k, v in records.items() if k not in NCRNA}),
            _fasta(tmp_path / 'ncrna.fa.gz', {k: v for k, v in records.items() if k in NCRNA}),
        ),
        proteins=_fasta(tmp_path / 'pep.fa.gz', dict(inputs.proteins)),
        hgnc=_tsv(tmp_path / 'hgnc.txt', [HGNC_COLUMNS, *HGNC], gzipped=False),
        mane=_tsv(tmp_path / 'mane.txt.gz', [['RefSeq_nuc', 'Ensembl_nuc', 'MANE_status'], *MANE], gzipped=True),
        assembly_report=tmp_path / 'report.txt',
        genome=tmp_path / 'genome',
    )


def _by_gene_id(release: ensembl.Release) -> dict[str, bundle_pb2.GeneBundle]:
    return {b.gene.ensembl_gene_id: b for b in ensembl.bundles(release)}


def _transcript(bundles: dict[str, bundle_pb2.GeneBundle], accession: str) -> bundle_pb2.Transcript:
    (found,) = [t for b in bundles.values() for t in b.transcripts if t.accession == accession]
    return found


def _line(kind: str, containing: str) -> str:
    """The one annotation line of this kind containing this text."""
    (found,) = [line for line in GFF if line.split('\t')[2:3] == [kind] and containing in line]
    return found


def _replacing(old: str, new: str) -> tuple[str, ...]:
    """The annotation with one exact line replaced."""
    if old not in GFF:
        raise ValueError(f'no such line: {old!r}')
    return tuple(new if line == old else line for line in GFF)


def _exons(alignment: bundle_pb2.Alignment) -> list[tuple[int, int, int, int, str]]:
    return [
        (e.transcript_start, e.transcript_end, e.genome_start, e.genome_end_inclusive, e.cigar) for e in alignment.exons
    ]


# ---- placements --------------------------------------------------------------------------------------


def test_a_plus_strand_transcript_is_its_annotated_exons(tmp_path: pathlib.Path) -> None:
    (alignment,) = _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000001').alignments
    assert (alignment.chromosome, alignment.strand) == (CHROM, bundle_pb2.STRAND_PLUS)
    assert _exons(alignment) == [(0, 20, 100, 119, '20='), (20, 50, 200, 229, '30=')]


def test_an_annotated_placement_says_where_it_came_from(tmp_path: pathlib.Path) -> None:
    transcript = _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000001')
    assert transcript.alignments[0].source == bundle_pb2.ALIGNMENT_SOURCE_ANNOTATION


def test_a_transcript_names_ensembl_as_its_publisher(tmp_path: pathlib.Path) -> None:
    assert _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000001').publisher == bundle_pb2.PUBLISHER_ENSEMBL


def test_a_record_base_the_genome_does_not_have_is_a_mismatch_in_the_cigar(tmp_path: pathlib.Path) -> None:
    (alignment,) = _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000002').alignments
    assert alignment.strand == bundle_pb2.STRAND_MINUS
    assert _exons(alignment) == [(0, 30, 600, 629, '5=1X24='), (30, 50, 500, 519, '20=')]


def test_records_are_read_from_every_sequence_set(tmp_path: pathlib.Path) -> None:
    # ENST00000000003.1's record is in the ncRNA set, the second --transcripts file
    assert _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000003').biotype == 'lnc_RNA'


def test_a_transcript_keeps_the_annotations_feature_type_as_its_biotype(tmp_path: pathlib.Path) -> None:
    bundles = _by_gene_id(_release(tmp_path))
    biotypes = {a: _transcript(bundles, a).biotype for a in ('ENST00000000001', 'ENST00000000008')}
    assert biotypes == {'ENST00000000001': 'mRNA', 'ENST00000000008': 'transcript'}


# ---- the CDS ----------------------------------------------------------------------------------------


def test_the_cds_is_the_annotations_projected_through_the_exons(tmp_path: pathlib.Path) -> None:
    bundles = _by_gene_id(_release(tmp_path))
    fwd, rev = _transcript(bundles, 'ENST00000000001').cds, _transcript(bundles, 'ENST00000000002').cds
    # FWD: genome 105 is exon 1's fifth base, 210 exon 2's tenth; REV: 605 is 25 in from 630, 505 is 15 in from 520
    assert (fwd.start_index, fwd.end_index_inclusive) == (4, 29)
    assert (rev.start_index, rev.end_index_inclusive) == (25, 45)


_FWD_FIRST_CDS = _line('CDS', '\t105\t')


def test_a_cds_beginning_mid_codon_is_left_off_its_transcript(tmp_path: pathlib.Path) -> None:
    mid_codon = _FWD_FIRST_CDS.replace('\t0\tID=', '\t1\tID=')
    transcript = _transcript(
        _by_gene_id(_release(tmp_path, Inputs(gff=_replacing(_FWD_FIRST_CDS, mid_codon)))), 'ENST00000000001'
    )
    assert not transcript.HasField('cds')


def test_a_transcript_whose_cds_begins_mid_codon_keeps_its_protein(tmp_path: pathlib.Path) -> None:
    mid_codon = _FWD_FIRST_CDS.replace('\t0\tID=', '\t1\tID=')
    bundle = _by_gene_id(_release(tmp_path, Inputs(gff=_replacing(_FWD_FIRST_CDS, mid_codon))))['ENSG00000000001']
    assert [p.accession for p in bundle.proteins] == ['ENSP00000000001']


@pytest.mark.parametrize(
    ('span', 'keeps_cds'),
    [
        ('\t601\t', False),  # REV's 5'-most CDS span on the minus strand, the higher one
        ('\t505\t', True),  # its 3'-most: a phase there is a codon split across the intron
    ],
)
def test_on_the_minus_strand_the_phase_that_matters_is_the_highest_spans(
    tmp_path: pathlib.Path, span: str, keeps_cds: bool
) -> None:
    line = _line('CDS', span)
    phased = Inputs(gff=_replacing(line, line.replace('\t0\tID=', '\t2\tID=')))
    assert _transcript(_by_gene_id(_release(tmp_path, phased)), 'ENST00000000002').HasField('cds') == keeps_cds


def test_the_report_counts_the_cdss_left_off_for_beginning_mid_codon(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    mid_codon = _FWD_FIRST_CDS.replace('\t0\tID=', '\t1\tID=')
    list(ensembl.bundles(_release(tmp_path, Inputs(gff=_replacing(_FWD_FIRST_CDS, mid_codon)))))
    assert 'bundled without a CDS for beginning mid-codon: 1;' in capsys.readouterr().err


def test_a_cds_feature_without_a_phase_fails_the_build(tmp_path: pathlib.Path) -> None:
    unphased = _FWD_FIRST_CDS.replace('\t0\tID=', '\t.\tID=')
    with pytest.raises(build.BuildError, match=r"a CDS feature of transcript:ENST00000000001 with phase '\.'"):
        list(ensembl.bundles(_release(tmp_path, Inputs(gff=_replacing(_FWD_FIRST_CDS, unphased)))))


# ---- proteins ------------------------------------------------------------------------------------------


def test_a_protein_is_bundled_with_its_own_version_and_its_transcripts(tmp_path: pathlib.Path) -> None:
    (protein,) = _by_gene_id(_release(tmp_path))['ENSG00000000002'].proteins
    assert (protein.accession, protein.version, protein.transcript_accession, protein.transcript_version) == (
        'ENSP00000000002',
        3,
        'ENST00000000002',
        1,
    )


def test_a_protein_carries_its_residues(tmp_path: pathlib.Path) -> None:
    bundle = _by_gene_id(_release(tmp_path))['ENSG00000000002']
    (protein,) = bundle.proteins
    residues = {s.digest: (s.alphabet, s.residues) for s in bundle.sequences}
    assert residues[protein.sequence_digest] == (bundle_pb2.ALPHABET_PROTEIN, b'MSEQREV')


# ---- genes ----------------------------------------------------------------------------------------------


def test_hgnc_names_a_gene_it_has(tmp_path: pathlib.Path) -> None:
    gene = _by_gene_id(_release(tmp_path))['ENSG00000000001'].gene
    assert (gene.symbol, gene.hgnc_id, gene.ncbi_gene_id, gene.name) == ('FWD', 'HGNC:1', '11', 'forward gene, HGNC')
    assert (list(gene.alias_symbols), list(gene.previous_symbols)) == (['FW1'], ['OLDFWD'])


def test_of_two_hgnc_rows_for_a_gene_the_one_named_as_the_annotation_names_it(tmp_path: pathlib.Path) -> None:
    gene = _by_gene_id(_release(tmp_path))['ENSG00000000002'].gene
    assert (gene.symbol, gene.hgnc_id) == ('REV', 'HGNC:4')  # not REVB, whose id is lower


def test_the_other_hgnc_rows_symbols_stay_findable_as_aliases(tmp_path: pathlib.Path) -> None:
    gene = _by_gene_id(_release(tmp_path))['ENSG00000000002'].gene
    assert sorted(gene.alias_symbols) == ['OLDREVB', 'REVB']


def test_of_two_hgnc_rows_neither_named_as_the_annotation_the_lowest_id_names_the_gene(tmp_path: pathlib.Path) -> None:
    rev = _line('gene', 'gene:ENSG00000000002;')
    gene = _by_gene_id(_release(tmp_path, Inputs(gff=_replacing(rev, rev.replace('Name=REV;', 'Name=REV1;')))))[
        'ENSG00000000002'
    ].gene
    assert (gene.symbol, gene.hgnc_id) == ('REVB', 'HGNC:3')


def test_a_gene_the_annotation_names_nothing_is_named_by_its_stable_id(tmp_path: pathlib.Path) -> None:
    assert _by_gene_id(_release(tmp_path))['ENSG00000000003'].gene.symbol == 'ENSG00000000003'


def test_the_annotations_description_loses_its_source_note(tmp_path: pathlib.Path) -> None:
    # without HGNC's row the name is the annotation's description
    release = _release(tmp_path, Inputs(gff=tuple(line.replace('ENSG00000000001', 'ENSG00000000077') for line in GFF)))
    assert _by_gene_id(release)['ENSG00000000077'].gene.name == 'forward gene'


def test_the_annotations_description_is_percent_decoded(tmp_path: pathlib.Path) -> None:
    assert _by_gene_id(_release(tmp_path))['ENSG00000000003'].gene.name == 'long, non-coding'


def test_bundles_come_in_symbol_order(tmp_path: pathlib.Path) -> None:
    assert [b.gene.symbol for b in ensembl.bundles(_release(tmp_path))] == ['ENSG00000000003', 'FWD', 'REV']


def test_a_genes_transcripts_come_in_accession_order(tmp_path: pathlib.Path) -> None:
    # ENST00000000008 comes first in the annotation
    bundle = _by_gene_id(_release(tmp_path))['ENSG00000000001']
    assert [t.accession for t in bundle.transcripts] == ['ENST00000000001', 'ENST00000000008']


# ---- MANE ------------------------------------------------------------------------------------------------


def test_a_mane_transcript_carries_its_refseq_partner(tmp_path: pathlib.Path) -> None:
    assert _transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000001').mane_partner == 'NM_000010.2'


def test_the_annotations_tags_are_the_transcripts(tmp_path: pathlib.Path) -> None:
    tags = set(_transcript(_by_gene_id(_release(tmp_path)), 'ENST00000000001').tags)
    assert tags == {bundle_pb2.TAG_MANE_SELECT, bundle_pb2.TAG_ENSEMBL_CANONICAL}


def test_the_mane_summary_tags_a_transcript_the_annotation_does_not(tmp_path: pathlib.Path) -> None:
    tagged = _line('mRNA', 'transcript:ENST00000000001;')
    untagged = tagged.replace('tag=gencode_basic,Ensembl_canonical,MANE_Select;', 'tag=gencode_basic;')
    transcript = _transcript(
        _by_gene_id(_release(tmp_path, Inputs(gff=_replacing(tagged, untagged)))), 'ENST00000000001'
    )
    assert set(transcript.tags) == {bundle_pb2.TAG_MANE_SELECT}


# ---- what is left out, and what fails -------------------------------------------------------------------


def test_a_non_coding_transcript_with_no_record_is_left_out(tmp_path: pathlib.Path) -> None:
    bundle = _by_gene_id(_release(tmp_path))['ENSG00000000003']
    assert [t.accession for t in bundle.transcripts] == ['ENST00000000003']  # not the TEC model ENST00000000004


def test_the_report_counts_what_it_leaves_out_and_names_the_mane_among_it(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    list(ensembl.bundles(_release(tmp_path)))
    report = capsys.readouterr().err
    assert 'left out for having no record: 1 (unconfirmed_transcript 1)' in report
    assert 'placed off these chromosomes, left out: 1, of them MANE 1: ENST00000000099.1' in report


def test_a_coding_transcript_with_no_record_fails_the_build(tmp_path: pathlib.Path) -> None:
    records = tuple((k, v) for k, v in RECORDS.items() if k != 'ENST00000000002.1')
    with pytest.raises(build.BuildError, match=r'coding transcripts have no record.*ENST00000000002\.1'):
        list(ensembl.bundles(_release(tmp_path, Inputs(records=records))))


def test_a_record_whose_length_is_not_its_exons_fails_the_build(tmp_path: pathlib.Path) -> None:
    records = tuple((k, v + 'A' if k == 'ENST00000000001.2' else v) for k, v in RECORDS.items())
    with pytest.raises(build.BuildError, match=r'ENST00000000001\.2: a 51-base record for 50 bases of annotated exons'):
        list(ensembl.bundles(_release(tmp_path, Inputs(records=records))))


def test_a_chromosome_the_assembly_report_does_not_name_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match=r"sequences \['99'\] have no accession"):
        list(ensembl.bundles(_release(tmp_path, Inputs(report=tuple(REPORT[:2])))))


def test_a_chromosome_named_for_a_sequence_the_genome_lacks_fails_the_build(tmp_path: pathlib.Path) -> None:
    report = (*REPORT[:2], REPORT[2].replace(CHROM, 'NC_000098.1'))
    with pytest.raises(build.BuildError, match=r"sequences \['99'\] have no accession"):
        list(ensembl.bundles(_release(tmp_path, Inputs(report=report))))


def test_a_genome_on_another_assembly_fails_the_build(tmp_path: pathlib.Path) -> None:
    with pytest.raises(build.BuildError, match='the genome is on ASSEMBLY_GRCH37, the release on ASSEMBLY_GRCH38'):
        list(ensembl.bundles(_release(tmp_path, Inputs(genome_assembly=bundle_pb2.ASSEMBLY_GRCH37))))


def test_a_protein_missing_from_the_peptide_set_fails_the_build(tmp_path: pathlib.Path) -> None:
    proteins = (('ENSP00000000001.2', 'MSEQFWD'),)
    with pytest.raises(build.BuildError, match=r'ENSP00000000002\.3 is not in the protein set'):
        list(ensembl.bundles(_release(tmp_path, Inputs(proteins=proteins))))


def test_a_truncated_feature_line_fails_the_build(tmp_path: pathlib.Path) -> None:
    rev = _line('gene', 'gene:ENSG00000000002;')
    truncated = Inputs(gff=_replacing(rev, rev.rsplit('\t', 1)[0]))
    with pytest.raises(build.BuildError, match='8 fields, not 9'):
        list(ensembl.bundles(_release(tmp_path, truncated)))


# ---- the shard and the command ----------------------------------------------------------------------------


def test_the_shard_names_every_file_it_was_cut_from(tmp_path: pathlib.Path) -> None:
    release = _release(tmp_path)
    shard = store_build.write_shard(
        ensembl.bundles(release), tmp_path / 'shards', release='ensembl_TEST', inputs=release.inputs()
    )
    roles = {i.role for i in store_build.shard_record(shard).inputs}
    assert roles == {
        'annotation',
        'transcripts/0',
        'transcripts/1',
        'proteins',
        'hgnc',
        'mane',
        'assembly_report',
        'genome',
    }


def _command(release: ensembl.Release, shards: pathlib.Path, *extra: str) -> list[str]:
    return [
        'ensembl', '--assembly', 'GRCh38', '--release', 'ensembl_TEST',
        '--annotation', str(release.annotation),
        '--transcripts', str(release.transcripts[0]), '--transcripts', str(release.transcripts[1]),
        '--proteins', str(release.proteins), '--hgnc', str(release.hgnc), '--mane', str(release.mane),
        '--assembly-report', str(release.assembly_report), '--genome', str(release.genome), '--shards', str(shards),
        *extra,
    ]  # fmt: skip


def test_the_command_names_the_shard_for_its_release_and_prefix(
    tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cli.main(_command(_release(tmp_path), tmp_path / 'shards', '--prefix', '02-'))
    shard = pathlib.Path(capsys.readouterr().out.strip())
    assert shard.name.startswith('02-ensembl_TEST-')
    assert store_build.shard_record(shard).release == 'ensembl_TEST'


def test_the_command_reads_both_sequence_sets(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> None:
    cli.main(_command(_release(tmp_path), tmp_path / 'shards'))
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    assert store_mod.BundleStore(str(tmp_path / 'store')).by_accession('ENST00000000003.1')


@pytest.fixture
def provider(tmp_path: pathlib.Path, capsys: pytest.CaptureFixture[str]) -> provider_mod.BundleProvider:
    release = _release(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        cli.main(_command(release, tmp_path / 'shards'))
    shard = capsys.readouterr().out.strip()
    cli.main(['index', '--assembly', 'GRCh38', '--out', str(tmp_path / 'store'), shard])
    return provider_mod.BundleProvider(
        store_mod.BundleStore(str(tmp_path / 'store')), genome_mod.Genome(str(release.genome))
    )


def test_weaver_projects_an_ensembl_coding_variant_onto_the_genome(provider: provider_mod.BundleProvider) -> None:
    # REV's c.1 is genome 605 on the minus strand, so the bases are complemented
    genome_base = _g(605, 605)
    coding_base = _revcomp(genome_base)
    alt = 'G' if coding_base != 'G' else 'C'
    variant = weaver.parse(f'ENST00000000002.1:c.1{coding_base}>{alt}')
    assert variant.validate(provider)
    assert (
        weaver.VariantMapper(provider).c_to_g(variant, CHROM).format() == f'{CHROM}:g.605{genome_base}>{_revcomp(alt)}'
    )


def test_ensembl_accessions_are_identified_as_transcripts_and_proteins(provider: provider_mod.BundleProvider) -> None:
    assert provider.get_identifier_type('ENST00000000001.2') == weaver.IdentifierType.TranscriptAccession
    assert provider.get_identifier_type('ENSP00000000001.2') == weaver.IdentifierType.ProteinAccession


def test_an_ensembl_transcripts_protein_is_its_c_to_p_target(provider: provider_mod.BundleProvider) -> None:
    assert provider.get_symbol_accessions('ENST00000000002.1', 'c', 'p') == [
        (weaver.IdentifierType.ProteinAccession, 'ENSP00000000002.3')
    ]


def test_an_ensembl_protein_is_read_by_its_accession(provider: provider_mod.BundleProvider) -> None:
    assert provider.get_seq('ENSP00000000002.3', 0, None, weaver.IdentifierType.ProteinAccession) == 'MSEQREV'
