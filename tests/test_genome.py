"""The derived genome: the blocks file, the catalogue, and the reader over them."""

from __future__ import annotations

import concurrent.futures
import gzip
import pathlib
import random

import pytest

from weaver_data_provider import _files, build, refget, testing
from weaver_data_provider import genome as genome_mod
from weaver_data_provider.build import genome as genome_build
from weaver_data_provider.v1 import bundle_pb2, genome_pb2

_ASSEMBLY = bundle_pb2.ASSEMBLY_GRCH38


@pytest.fixture
def built(tmp_path: pathlib.Path) -> testing.Reference:
    return testing.build(tmp_path)


def _random_genome(tmp_path: pathlib.Path, lengths: dict[str, int]) -> dict[str, str]:
    """A FASTA of random sequences, soft masking included, at `tmp_path/g.fa.gz`; returns them as published."""
    rng = random.Random(3)
    sequences = {name: ''.join(rng.choices('ACGTacgtN', k=n)) for name, n in lengths.items()}
    with gzip.open(tmp_path / 'g.fa.gz', 'wt', encoding='ascii') as fh:
        for name, sequence in sequences.items():
            fh.write(f'>{name}\n')
            for k in range(0, len(sequence), 61):
                fh.write(sequence[k : k + 61] + '\n')
    return sequences


def _rewrite_catalogue(root: pathlib.Path, catalogue: genome_pb2.SequenceCatalogue) -> None:
    (root / genome_mod.CATALOGUE).unlink()
    _files.write_record(root / genome_mod.CATALOGUE, catalogue.SerializeToString())


def test_slices_come_back_as_the_genome_has_them(built: testing.Reference) -> None:
    genome = testing.genome_sequence()
    assert built.genome.fetch(testing.CHROM, 0, 8) == genome[:8]
    assert built.genome.fetch(testing.CHROM, 1103, 1110) == genome[1103:1110]
    assert built.genome.fetch(testing.CHROM, 2995) == genome[2995:]  # to the end
    assert built.genome.fetch(testing.CHROM, 2995, 5000) == genome[2995:]  # clipped, not an error
    assert built.genome.fetch(testing.CHROM, 10, 10) == ''


def test_slices_across_block_boundaries_match_the_source(tmp_path: pathlib.Path) -> None:
    # Blocks of 7 bases put a boundary inside nearly every slice, and sequence ends mid-block.
    sequences = _random_genome(tmp_path, {'a': 500, 'b': 7, 'c': 64, 'd': 1})
    genome_build.build_genome(tmp_path / 'g.fa.gz', tmp_path / 'out', assembly=_ASSEMBLY, block_size=7)
    genome = genome_mod.Genome(str(tmp_path / 'out'))
    rng = random.Random(5)
    for name, sequence in sequences.items():
        for _ in range(200):
            start = rng.randrange(len(sequence) + 3)
            end = start + rng.randrange(40)
            assert genome.fetch(name, start, end) == sequence[start:end].upper()
        assert genome.fetch(name, 0) == sequence.upper()


def test_a_gather_answers_each_range_as_a_fetch_would(tmp_path: pathlib.Path) -> None:
    sequences = _random_genome(tmp_path, {'a': 900})
    genome_build.build_genome(tmp_path / 'g.fa.gz', tmp_path / 'out', assembly=_ASSEMBLY, block_size=16)
    genome = genome_mod.Genome(str(tmp_path / 'out'))
    ranges: list[tuple[int, int | None]] = [(850, 910), (3, 40), (0, 0), (100, 101), (90, 400), (880, None), (5, -1)]
    expected = sequences['a'].upper()
    assert (
        genome.fetch_many('a', ranges)
        == [
            expected[850:900],  # clipped at the sequence's end
            expected[3:40],
            '',
            expected[100:101],
            expected[90:400],
            expected[880:],
            expected[5:],  # a negative end runs to the end, as None does
        ]
    )


def test_one_genome_serves_many_threads_through_a_small_cache(tmp_path: pathlib.Path) -> None:
    # A two-block cache over 64 blocks evicts on nearly every read, so threads race on the eviction.
    sequences = _random_genome(tmp_path, {'a': 1024})
    genome_build.build_genome(tmp_path / 'g.fa.gz', tmp_path / 'out', assembly=_ASSEMBLY, block_size=16)
    genome = genome_mod.Genome(str(tmp_path / 'out'), cache_blocks=2)
    expected = sequences['a'].upper()

    def reads(seed: int) -> bool:
        rng = random.Random(seed)
        ranges = [(s, s + rng.randrange(1, 50)) for s in rng.choices(range(1000), k=300)]
        return all(genome.fetch('a', s, e) == expected[s:e] for s, e in ranges)

    with concurrent.futures.ThreadPoolExecutor(12) as pool:
        assert all(pool.map(reads, range(12)))


def test_the_catalogue_names_lengths_and_refget_digests(built: testing.Reference) -> None:
    genome = testing.genome_sequence()
    assert built.genome.names == [testing.CHROM]
    assert testing.CHROM in built.genome
    assert built.genome.length(testing.CHROM) == testing.GENOME_LENGTH
    assert built.genome.digest(testing.CHROM) == refget.digest(genome.encode())
    assert built.genome.name_for_digest(refget.digest(genome.encode())) == testing.CHROM
    assert built.genome.name_for_digest('SQ.' + 'A' * 32) is None
    assert built.genome.assembly == _ASSEMBLY


def test_soft_masked_residues_digest_as_their_upper_case(tmp_path: pathlib.Path) -> None:
    fasta = tmp_path / 'g.fa.gz'
    with gzip.open(fasta, 'wt') as fh:
        fh.write('>chrT masked\nacgtACGTacgt\n')
    genome_build.build_genome(fasta, tmp_path / 'out', assembly=_ASSEMBLY)
    genome = genome_mod.Genome(str(tmp_path / 'out'))
    assert genome.digest('chrT') == refget.digest(b'ACGTACGTACGT')
    assert genome.fetch('chrT', 0, 12) == 'ACGTACGTACGT'  # slices are upper-cased on the way out too


def test_a_plain_fasta_and_a_multi_record_one_are_cut(tmp_path: pathlib.Path) -> None:
    fasta = tmp_path / 'g.fa'
    fasta.write_text('>one first\nACGT\n>two\nGGGGCCCC\nAA\n', 'ascii')
    catalogue = genome_build.build_genome(fasta, tmp_path / 'out', assembly=_ASSEMBLY)
    assert [(s.name, s.length, s.description) for s in catalogue.sequences] == [('one', 4, 'first'), ('two', 10, '')]
    assert genome_mod.Genome(str(tmp_path / 'out')).fetch('two', 6) == 'CCAA'


@pytest.mark.parametrize(
    ('text', 'match'),
    [
        ('>empty\n>ok\nACGT\n', 'empty'),
        ('>one\nACGT\n>one\nACGT\n', 'twice'),
        ('ACGT\n>one\nACGT\n', 'before the first header'),
    ],
)
def test_a_malformed_fasta_fails_the_build(tmp_path: pathlib.Path, text: str, match: str) -> None:
    fasta = tmp_path / 'g.fa'
    fasta.write_text(text, 'ascii')
    with pytest.raises(build.BuildError, match=match):
        genome_build.build_genome(fasta, tmp_path / 'out', assembly=_ASSEMBLY)
    assert not (tmp_path / 'out' / genome_mod.CATALOGUE).exists()  # no catalogue, so no genome


def test_a_tar_archive_genome_is_refused(tmp_path: pathlib.Path) -> None:
    archive = tmp_path / 'genome.tar.gz'
    archive.write_bytes(b'')
    with pytest.raises(build.BuildError, match='tar'):
        genome_build.build_genome(archive, tmp_path / 'out', assembly=_ASSEMBLY)


def test_a_catalogue_implying_other_blocks_is_refused(built: testing.Reference) -> None:
    root = pathlib.Path(built.genome.root)
    catalogue = genome_mod.read_catalogue(str(root))
    catalogue.sequences.add(name='ghost', length=1, digest='SQ.' + 'A' * 32)
    _rewrite_catalogue(root, catalogue)
    with pytest.raises(genome_mod.GenomeError, match='catalogue implies'):
        genome_mod.Genome(str(root))


def test_a_rewritten_blocks_file_is_refused(tmp_path: pathlib.Path) -> None:
    # Same sequence count and lengths, different bases: only the digest tells them apart.
    for text, out in (('>a\nACGT\n', 'first'), ('>a\nTTTT\n', 'second')):
        (tmp_path / 'g.fa').write_text(text, 'ascii')
        genome_build.build_genome(tmp_path / 'g.fa', tmp_path / out, assembly=_ASSEMBLY)
    (tmp_path / 'first' / genome_mod.BLOCKS).write_bytes((tmp_path / 'second' / genome_mod.BLOCKS).read_bytes())
    with pytest.raises(genome_mod.GenomeError, match='crc32c'):
        genome_mod.Genome(str(tmp_path / 'first'))


def test_an_unknown_format_version_is_refused(built: testing.Reference) -> None:
    root = pathlib.Path(built.genome.root)
    catalogue = genome_mod.read_catalogue(str(root))
    catalogue.format_version += 1
    _rewrite_catalogue(root, catalogue)
    with pytest.raises(genome_mod.GenomeError, match='format version'):
        genome_mod.Genome(str(root))


def test_a_negative_start_reads_from_the_sequences_first_base(tmp_path: pathlib.Path) -> None:
    # blocks of 4, so a start of -2 on `b` would otherwise reach into `a`'s last block
    sequences = _random_genome(tmp_path, {'a': 8, 'b': 10})
    genome_build.build_genome(tmp_path / 'g.fa.gz', tmp_path / 'out', assembly=_ASSEMBLY, block_size=4)
    genome = genome_mod.Genome(str(tmp_path / 'out'))
    assert genome.fetch('b', -2, 3) == sequences['b'][:3].upper()
    assert genome.fetch('a', -2, 3) == sequences['a'][:3].upper()


def test_a_catalogue_naming_a_digest_this_reader_cannot_check_is_refused(built: testing.Reference) -> None:
    root = pathlib.Path(built.genome.root)
    catalogue = genome_mod.read_catalogue(str(root))
    catalogue.digest, catalogue.digest_algorithm = '0' * 32, 'md5'
    _rewrite_catalogue(root, catalogue)
    with pytest.raises(genome_mod.GenomeError, match='cannot check a md5 digest'):
        genome_mod.Genome(str(root))
