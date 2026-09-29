"""A small, internally consistent reference set for tests, here and in weaver's own suite.

The genome is a 3,000-base repeat of ACGT. FWD is a plus-strand coding gene of two exons whose record
is the spliced genome. REV is a minus-strand coding gene of two exons whose record differs from the
genome by one substitution, so a variant at that base validates against the record and not the
genome. Every coordinate is derived from those definitions, so a test asserts arithmetic, not
remembered numbers.

    reference = testing.build(tmp_path)
    provider = reference.provider()
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib

from weaver_data_provider import genome as genome_mod
from weaver_data_provider import provider as provider_mod
from weaver_data_provider import refget
from weaver_data_provider import store as store_mod
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2

CHROM = 'NC_000099.1'
GENOME_LENGTH = 3000
FWD_EXONS = ((101, 110), (201, 212))  # 1-based closed
FWD_CDS = (104, 205)
REV_EXONS = ((1001, 1010), (1101, 1110))
REV_CDS = (1003, 1106)
REV_SUBSTITUTED_TRANSCRIPT_INDEX = 2  # transcript base 3 differs from the genome
ASSEMBLY = bundle_pb2.ASSEMBLY_GRCH38


def genome_sequence() -> str:
    return ('ACGT' * (GENOME_LENGTH // 4))[:GENOME_LENGTH]


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


@dataclasses.dataclass(frozen=True)
class Reference:
    """The built store and genome, with the sequences a test compares against."""

    store: store_mod.BundleStore
    genome: genome_mod.Genome
    fwd: str  # FWD's record, plus strand
    rev: str  # REV's record, minus strand, one base off the genome
    rev_genomic: str  # what REV would be if it matched the genome

    def provider(self) -> provider_mod.BundleProvider:
        return provider_mod.BundleProvider(self.store, self.genome)


def _exons(spans: tuple[tuple[int, int], ...], *, minus: bool) -> list[bundle_pb2.Exon]:
    ordered = sorted(spans, reverse=minus)
    out, t = [], 0
    for start, end in ordered:
        n = end - start + 1
        out.append(
            bundle_pb2.Exon(
                transcript_start=t,
                transcript_end=t + n,
                genome_start=start - 1,
                genome_end_inclusive=end - 1,
                cigar=f'{n}=',
            )
        )
        t += n
    return out


def _cds_indices(exons: list[bundle_pb2.Exon], cds: tuple[int, int], *, minus: bool) -> tuple[int, int]:
    low, high = cds[0] - 1, cds[1] - 1
    first, last = (high, low) if minus else (low, high)

    def index(g: int) -> int:
        for e in exons:
            if e.genome_start <= g <= e.genome_end_inclusive:
                return e.transcript_start + (e.genome_end_inclusive - g if minus else g - e.genome_start)
        raise ValueError(g)

    return index(first), index(last)


def _bundle(
    symbol: str,
    accession: str,
    protein: str,
    residues: str,
    exons: list[bundle_pb2.Exon],
    cds: tuple[int, int],
    *,
    minus: bool,
    aliases: tuple[str, ...] = (),
) -> bundle_pb2.GeneBundle:
    digest = refget.digest(residues.encode())
    protein_residues = b'MKLV'
    start, end = _cds_indices(exons, cds, minus=minus)
    return bundle_pb2.GeneBundle(
        gene=bundle_pb2.Gene(symbol=symbol, alias_symbols=list(aliases), ncbi_gene_id='1'),
        transcripts=[
            bundle_pb2.Transcript(
                accession=accession,
                version=1,
                publisher=bundle_pb2.PUBLISHER_REFSEQ,
                status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
                biotype='mRNA',
                tags=[bundle_pb2.TAG_MANE_SELECT],
                sequence_digest=digest,
                cds=bundle_pb2.Cds(start_index=start, end_index_inclusive=end),
                protein_accession=protein,
                protein_version=1,
                alignments=[
                    bundle_pb2.Alignment(
                        assembly=ASSEMBLY,
                        chromosome=CHROM,
                        strand=bundle_pb2.STRAND_MINUS if minus else bundle_pb2.STRAND_PLUS,
                        source=bundle_pb2.ALIGNMENT_SOURCE_NCBI_BAM,
                        exons=exons,
                    )
                ],
            )
        ],
        proteins=[
            bundle_pb2.Protein(
                accession=protein,
                version=1,
                sequence_digest=refget.digest(protein_residues),
                transcript_accession=accession,
                transcript_version=1,
            )
        ],
        sequences=[
            bundle_pb2.Sequence(digest=digest, alphabet=bundle_pb2.ALPHABET_NUCLEOTIDE, residues=residues.encode()),
            bundle_pb2.Sequence(
                digest=refget.digest(protein_residues), alphabet=bundle_pb2.ALPHABET_PROTEIN, residues=protein_residues
            ),
        ],
    )


def bundles() -> tuple[list[bundle_pb2.GeneBundle], str, str, str]:
    """The two genes' bundles, and the FWD record, the REV record and REV's genome-spliced form."""
    genome = genome_sequence()
    fwd = ''.join(genome[s - 1 : e] for s, e in FWD_EXONS)
    rev_genomic = revcomp(''.join(genome[s - 1 : e] for s, e in REV_EXONS))
    i = REV_SUBSTITUTED_TRANSCRIPT_INDEX
    rev = rev_genomic[:i] + ('A' if rev_genomic[i] != 'A' else 'C') + rev_genomic[i + 1 :]
    made = [
        _bundle(
            'FWD',
            'NM_000001',
            'NP_000001',
            fwd,
            _exons(FWD_EXONS, minus=False),
            FWD_CDS,
            minus=False,
            aliases=('OLDFWD',),
        ),
        _bundle('REV', 'NM_000002', 'NP_000002', rev, _exons(REV_EXONS, minus=True), REV_CDS, minus=True),
    ]
    return made, fwd, rev, rev_genomic


def build(root: pathlib.Path) -> Reference:
    """Write the store under `root/store` and the genome under `root/genome`, and open both."""
    made, fwd, rev, rev_genomic = bundles()
    inputs_file = root / 'inputs.txt'
    root.mkdir(parents=True, exist_ok=True)
    inputs_file.write_text('synthetic\n', 'utf-8')
    shard_path = store_build.write_shard(
        made, root / 'shards', release='SYNTHETIC_1', inputs=[('annotation', inputs_file)]
    )
    store_build.write_index(root / 'store', [shard_path], assembly=ASSEMBLY)

    genome = genome_sequence()
    fasta = root / 'genome.fa.gz'
    with gzip.open(fasta, 'wt', encoding='ascii') as fh:
        fh.write(f'>{CHROM} a synthetic chromosome\n')
        for k in range(0, GENOME_LENGTH, 70):
            fh.write(genome[k : k + 70] + '\n')
    genome_build.build_genome(fasta, root / 'genome', assembly=ASSEMBLY)
    return Reference(
        store=store_mod.BundleStore(str(root / 'store')),
        genome=genome_mod.Genome(str(root / 'genome')),
        fwd=fwd,
        rev=rev,
        rev_genomic=rev_genomic,
    )
