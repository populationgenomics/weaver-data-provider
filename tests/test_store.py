"""The key tables, the store's writer and reader, and the checks the reader makes against its manifest."""

from __future__ import annotations

import concurrent.futures
import pathlib

import bagz
import pytest

from weaver_data_provider import _files, build, keytable, refget, store
from weaver_data_provider.build import store as store_build
from weaver_data_provider.v1 import bundle_pb2, index_pb2, store_pb2

_ASSEMBLY = bundle_pb2.ASSEMBLY_GRCH38


def _bundle(
    symbol: str,
    accession: str,
    *,
    chromosome: str,
    start: int,
    residues: bytes,
    aliases: tuple[str, ...] = (),
    version: int = 2,
    assembly: bundle_pb2.Assembly = _ASSEMBLY,
) -> bundle_pb2.GeneBundle:
    digest = refget.digest(residues)
    return bundle_pb2.GeneBundle(
        gene=bundle_pb2.Gene(symbol=symbol, alias_symbols=list(aliases)),
        transcripts=[
            bundle_pb2.Transcript(
                accession=accession,
                version=version,
                publisher=bundle_pb2.PUBLISHER_REFSEQ,
                status=bundle_pb2.TRANSCRIPT_STATUS_CURRENT,
                sequence_digest=digest,
                protein_accession='NP_1',
                protein_version=1,
                alignments=[
                    bundle_pb2.Alignment(
                        assembly=assembly,
                        chromosome=chromosome,
                        strand=bundle_pb2.STRAND_PLUS,
                        source=bundle_pb2.ALIGNMENT_SOURCE_NCBI_BAM,
                        exons=[
                            bundle_pb2.Exon(
                                transcript_start=0,
                                transcript_end=4,
                                genome_start=start,
                                genome_end_inclusive=start + 3,
                                cigar='4=',
                            ),
                            bundle_pb2.Exon(
                                transcript_start=4,
                                transcript_end=len(residues),
                                genome_start=start + 100,
                                genome_end_inclusive=start + 100 + len(residues) - 5,
                                cigar=f'{len(residues) - 4}=',
                            ),
                        ],
                    )
                ],
            )
        ],
        proteins=[
            bundle_pb2.Protein(
                accession='NP_1',
                version=1,
                sequence_digest=refget.digest(b'MK'),
                transcript_accession=accession,
                transcript_version=version,
            )
        ],
        sequences=[
            bundle_pb2.Sequence(digest=digest, alphabet=bundle_pb2.ALPHABET_NUCLEOTIDE, residues=residues),
            bundle_pb2.Sequence(digest=refget.digest(b'MK'), alphabet=bundle_pb2.ALPHABET_PROTEIN, residues=b'MK'),
        ],
    )


_BUNDLES = [
    _bundle('ARSA', 'NM_000487', chromosome='NC_000022.11', start=5000, residues=b'GGGGCCCCAA', aliases=('ASA',)),
    _bundle('BRCA2', 'NM_000059', chromosome='NC_000013.11', start=1000, residues=b'ACGTACGTAC'),
    # an alias ending in a Greek alpha, as HGNC records some: it reduces to its ASCII stem
    _bundle(
        'HNF1A',
        'NM_000545',
        chromosome='NC_000012.12',
        start=300,
        residues=b'CCCCGGGGAA',
        aliases=('HNF1' + chr(0x3B1),),
    ),
    # a long gene starting first on chromosome 13, so a position past the short genes' ends still finds it
    _bundle('LONG', 'NM_000001', chromosome='NC_000013.11', start=100, residues=b'A' * 2000),
    _bundle('XPC', 'NM_004628', chromosome='NC_000003.12', start=200, residues=b'TTTTTTTTAA'),
]


def _inputs(tmp_path: pathlib.Path) -> list[tuple[str, pathlib.Path]]:
    path = tmp_path / 'annotation.gff3'
    path.write_text('an input\n', 'utf-8')
    return [('annotation', path)]


def _write(
    tmp_path: pathlib.Path, shards: list[tuple[str, list[bundle_pb2.GeneBundle]]], *, name: str = 'store'
) -> pathlib.Path:
    written = [
        store_build.write_shard(bundles, tmp_path / 'shards', release=release, inputs=_inputs(tmp_path))
        for release, bundles in shards
    ]
    store_build.write_index(tmp_path / name, written, assembly=_ASSEMBLY)
    return tmp_path / name


@pytest.fixture
def written(tmp_path: pathlib.Path) -> store.BundleStore:
    return store.BundleStore(str(_write(tmp_path, [('RS_1', _BUNDLES)])))


# ---- key tables ----------------------------------------------------------------------------------


def test_key_table_round_trips_every_key_and_answers_absent_ones_empty() -> None:
    writer = keytable.KeyTableWriter('test')
    keys = {f'key-{i:04d}'.encode(): [i, i + 1000] for i in range(500)}
    for key, records in keys.items():
        writer.add(key, records)
        writer.add(key, records)  # a repeated add is idempotent
    table = keytable.KeyTable(writer.message())
    assert len(table) == 500
    assert all(table.lookup(key) == tuple(records) for key, records in keys.items())
    assert table.lookup(b'key-0500') == ()
    assert table.lookup(b'key-') == ()  # a prefix of a key is not the key


def test_key_table_prefix_query_returns_every_version_once() -> None:
    writer = keytable.KeyTableWriter('accession')
    writer.add(b'NM_000059.3', [7])
    writer.add(b'NM_000059.4', [7])
    writer.add(b'NM_0000590.1', [8])  # shares the digits but not the `NM_000059.` prefix
    writer.add(b'NM_000060.1', [9])
    table = keytable.KeyTable(writer.message())
    assert table.prefixed(b'NM_000059.') == (7,)
    assert table.prefixed(b'NM_00005') == (7, 8)
    assert table.prefixed(b'NM_1') == ()


def test_an_empty_key_table_is_valid_and_answers_nothing() -> None:
    table = keytable.KeyTable(keytable.KeyTableWriter('empty').message())
    assert len(table) == 0
    assert table.lookup(b'anything') == ()
    assert table.prefixed(b'') == ()


def test_an_empty_key_is_refused() -> None:
    with pytest.raises(ValueError, match='empty key'):
        keytable.KeyTableWriter('test').add(b'', [1])


def test_refget_raw_form_round_trips() -> None:
    digest = refget.digest(b'ACGT')
    assert digest == 'SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'  # the refget spec's own example
    assert len(refget.raw(digest)) == 24
    assert refget.identifier(refget.raw(digest)) == digest
    with pytest.raises(ValueError, match='not a refget'):
        refget.raw('ga4gh:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2')


# ---- lookups -----------------------------------------------------------------------------------


def test_symbol_lookup_normalises_and_reaches_aliases(written: store.BundleStore) -> None:
    assert [b.gene.symbol for b in written.by_symbol('brca2')] == ['BRCA2']
    assert [b.gene.symbol for b in written.by_symbol('XP-C')] == ['XPC']
    assert [b.gene.symbol for b in written.by_symbol('ASA')] == ['ARSA']  # withdrawn symbol, as papers print it
    assert [b.gene.symbol for b in written.by_symbol('HNF1')] == ['HNF1A']  # the Greek-lettered alias, ASCII part
    assert written.by_symbol('NOSUCHGENE') == []
    assert written.by_symbol('') == []


def test_accession_lookup_accepts_versioned_and_unversioned(written: store.BundleStore) -> None:
    assert [b.gene.symbol for b in written.by_accession('NM_000059.2')] == ['BRCA2']
    assert [b.gene.symbol for b in written.by_accession('NM_000059')] == ['BRCA2']
    assert written.by_accession('NM_000059.9') == []
    assert written.by_accession('NM_00005') == []  # a bare accession is a whole accession, not a prefix


def test_protein_and_digest_lookups(written: store.BundleStore) -> None:
    assert len(written.by_protein('NP_1')) == 5  # the fixture gives every gene the same protein
    assert len(written.by_protein('NP_1.1')) == 5
    found = written.by_digest(refget.digest(b'acgtacgtac'))  # case-normalised
    assert [b.gene.symbol for b in found] == ['BRCA2']


def test_position_lookup_is_by_chromosome(written: store.BundleStore) -> None:
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 1050)] == ['BRCA2', 'LONG']
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 2000)] == ['LONG']
    assert written.overlapping('NC_000013.11', 99) == []
    assert written.overlapping('NC_000099.1', 1050) == []  # a chromosome with no placements


def test_interval_lookup_finds_every_placement_it_touches(written: store.BundleStore) -> None:
    # BRCA2's placement is 1000..1105 (0-based), LONG's 100..2105
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 0, 99)] == []
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 0, 100)] == ['LONG']
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 1105, 5000)] == ['BRCA2', 'LONG']
    assert [b.gene.symbol for b in written.overlapping('NC_000013.11', 1106, 5000)] == ['LONG']
    with pytest.raises(ValueError, match='precedes'):
        written.overlapping('NC_000013.11', 10, 5)


# ---- the cache -----------------------------------------------------------------------------------


class _CountingReads:
    """Stands in for a store's bundle reader, recording the record numbers each batch read asks for."""

    def __init__(self, reader: bagz.Reader) -> None:
        self._reader = reader
        self.batches: list[list[int]] = []

    def __len__(self) -> int:
        return len(self._reader)

    def read_indices(self, records: list[int]) -> list[bytes]:
        self.batches.append(list(records))
        return self._reader.read_indices(records)


def _count_reads(opened: store.BundleStore, monkeypatch: pytest.MonkeyPatch) -> _CountingReads:
    # The read is the behaviour a cache exists to save, so it is counted where the store makes it.
    counting = _CountingReads(opened._bundles)
    monkeypatch.setattr(opened, '_bundles', counting)
    return counting


def test_a_cached_bundle_is_not_read_again(written: store.BundleStore, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = _count_reads(written, monkeypatch)
    written.by_symbol('BRCA2')
    written.by_accession('NM_000059')
    assert reads.batches == [[1]]


def test_a_returned_bundle_is_the_callers_own(written: store.BundleStore) -> None:
    written.by_symbol('BRCA2')[0].gene.symbol = 'CHANGED'
    assert [b.gene.symbol for b in written.by_accession('NM_000059')] == ['BRCA2']


def test_a_batch_larger_than_the_cache_is_still_returned_whole(tmp_path: pathlib.Path) -> None:
    small = store.BundleStore(str(_write(tmp_path, [('RS_1', _BUNDLES)])), cache_size=2)
    found = small.bundles([0, 1, 2, 3, 0])
    assert [b.gene.symbol for b in found] == ['ARSA', 'BRCA2', 'HNF1A', 'LONG', 'ARSA']


def test_the_cache_keeps_the_most_recently_used(tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    small = store.BundleStore(str(_write(tmp_path, [('RS_1', _BUNDLES)])), cache_size=2)
    small.bundles([0])
    small.bundles([1])
    small.bundles([0])  # a hit renews record 0
    small.bundles([2])  # evicts record 1, the least recently used
    reads = _count_reads(small, monkeypatch)
    small.bundles([0])
    small.bundles([1])
    assert reads.batches == [[1]]


def test_one_store_serves_many_threads_with_a_cache_smaller_than_the_working_set(tmp_path: pathlib.Path) -> None:
    shared = store.BundleStore(str(_write(tmp_path, [('RS_1', _BUNDLES)])), cache_size=2)
    symbols = [b.gene.symbol for b in _BUNDLES]

    def lookups(seed: int) -> list[str]:
        return [shared.by_symbol(symbols[(seed + i) % len(symbols)])[0].gene.symbol for i in range(200)]

    with concurrent.futures.ThreadPoolExecutor(16) as pool:
        results = list(pool.map(lookups, range(16)))
    for seed, got in enumerate(results):
        assert got == [symbols[(seed + i) % len(symbols)] for i in range(200)]


# ---- shards and the manifest ---------------------------------------------------------------------


def test_record_numbers_run_across_shards_and_each_record_knows_its_release(tmp_path: pathlib.Path) -> None:
    newer = _bundle('BRCA2', 'NM_000059', chromosome='NC_000013.11', start=1000, residues=b'ACGTACGTAC', version=3)
    stacked = store.BundleStore(str(_write(tmp_path, [('RS_1', _BUNDLES), ('RS_2', [newer])])))
    assert len(stacked) == len(_BUNDLES) + 1
    assert stacked.release_of(0) == 'RS_1'
    assert stacked.release_of(len(_BUNDLES)) == 'RS_2'
    with pytest.raises(IndexError):
        stacked.release_of(len(stacked))
    # a symbol held in both releases answers both records, oldest release first
    held = stacked.by_symbol('BRCA2')
    assert [t.version for b in held for t in b.transcripts] == [2, 3]
    assert [b.transcripts[0].version for b in stacked.by_accession('NM_000059')] == [2, 3]


def test_a_shard_is_named_for_its_contents_and_written_once(tmp_path: pathlib.Path) -> None:
    path = store_build.write_shard(_BUNDLES, tmp_path, release='RS_1', inputs=_inputs(tmp_path), prefix='02-')
    assert path.name == f'02-RS_1-{_files.crc32c_of_local(path)}.bagz'
    # the same contents name the same path, and a shard there is not written over
    with pytest.raises(FileExistsError, match='already written'):
        store_build.write_shard(_BUNDLES, tmp_path, release='RS_1', inputs=_inputs(tmp_path), prefix='02-')
    # a release cut again with different contents lands beside it, not over it
    recut = store_build.write_shard(_BUNDLES[:2], tmp_path, release='RS_1', inputs=_inputs(tmp_path), prefix='02-')
    assert recut != path
    assert path.exists()


def test_a_release_name_that_is_not_a_plain_file_name_is_refused(tmp_path: pathlib.Path) -> None:
    with pytest.raises(ValueError, match='file name'):
        store_build.write_shard(_BUNDLES, tmp_path, release='../RS_1', inputs=_inputs(tmp_path))


def test_a_bundle_breaking_the_schema_is_refused_on_write(tmp_path: pathlib.Path) -> None:
    broken = bundle_pb2.GeneBundle()
    broken.CopyFrom(_BUNDLES[0])
    broken.transcripts[0].sequence_digest = 'not-a-digest'
    with pytest.raises(build.BuildError, match='sequence_digest'):
        store_build.write_shard([broken], tmp_path / 'shards', release='RS_1', inputs=_inputs(tmp_path))
    assert list((tmp_path / 'shards').iterdir()) == []  # nothing half-written is left behind


def test_a_protein_citing_a_sequence_its_bundle_does_not_carry_is_refused_on_write(tmp_path: pathlib.Path) -> None:
    dangling = bundle_pb2.GeneBundle()
    dangling.CopyFrom(_BUNDLES[0])
    dangling.proteins[0].sequence_digest = refget.digest(b'NOTCARRIED')
    with pytest.raises(build.BuildError, match='cites sequences it does not carry'):
        store_build.write_shard([dangling], tmp_path / 'shards', release='RS_1', inputs=_inputs(tmp_path))


def test_a_bundle_citing_a_sequence_it_does_not_carry_is_refused_on_write(tmp_path: pathlib.Path) -> None:
    dangling = bundle_pb2.GeneBundle()
    dangling.CopyFrom(_BUNDLES[0])
    dangling.transcripts[0].sequence_digest = refget.digest(b'NOTCARRIED')  # well formed, but not in the bundle
    with pytest.raises(build.BuildError, match='cites sequences it does not carry'):
        store_build.write_shard([dangling], tmp_path / 'shards', release='RS_1', inputs=_inputs(tmp_path))


def test_an_alignment_on_another_assembly_is_refused_by_the_index(tmp_path: pathlib.Path) -> None:
    grch37 = _bundle(
        'X', 'NM_9', chromosome='NC_000001.10', start=1, residues=b'ACGTACGT', assembly=bundle_pb2.ASSEMBLY_GRCH37
    )
    path = store_build.write_shard([grch37], tmp_path, release='RS_1', inputs=_inputs(tmp_path))
    with pytest.raises(build.BuildError, match='GRCH37 in a store on ASSEMBLY_GRCH38'):
        store_build.write_index(tmp_path / 'store', [path], assembly=_ASSEMBLY)


def test_a_store_without_its_manifest_is_not_a_store(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    (root / store.MANIFEST).unlink()
    with pytest.raises(FileNotFoundError):
        store.BundleStore(str(root))


def test_an_unknown_format_version_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    manifest = store.read_manifest(str(root))
    manifest.format_version += 1
    (root / store.MANIFEST).unlink()
    _files.write_record(root / store.MANIFEST, manifest.SerializeToString())
    with pytest.raises(store.StoreError, match='format version'):
        store.BundleStore(str(root))


def test_a_shard_rewritten_under_the_index_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    shard = _files.join(str(root), store.read_manifest(str(root)).shards[0].path)
    # the same number of records at the same path, in another order: an index would misread every one
    pathlib.Path(shard).unlink()
    with bagz.Writer(shard, bagz.Writer.Options(compression=bagz.CompressionZstd())) as writer:
        for bundle in reversed(_BUNDLES):
            writer.write(bundle.SerializeToString())
    with pytest.raises(store.StoreError, match='crc32c'):
        store.BundleStore(str(root))


def test_an_index_over_a_changed_shard_is_refused_at_build(tmp_path: pathlib.Path) -> None:
    path = store_build.write_shard(_BUNDLES, tmp_path, release='RS_1', inputs=_inputs(tmp_path))
    # the record beside the shard no longer matches the shard's bytes
    changed = store_build.shard_record(path)
    changed.digest = '0' * 8
    record = store_build.shard_record_path(path)
    record.unlink()
    _files.write_record(record, changed.SerializeToString())
    with pytest.raises(build.BuildError, match='crc32c'):
        store_build.write_index(tmp_path / 'store', [path], assembly=_ASSEMBLY)


def test_a_shard_leaves_the_record_an_index_needs_beside_it(tmp_path: pathlib.Path) -> None:
    path = store_build.write_shard(_BUNDLES, tmp_path, release='RS_1', inputs=_inputs(tmp_path))
    record = store_build.shard_record(path)
    assert (record.path, record.release, record.records) == (path.name, 'RS_1', len(_BUNDLES))
    assert record.digest == _files.crc32c_of_local(path)
    assert [i.role for i in record.inputs] == ['annotation']


def _only_symbols(root: pathlib.Path) -> None:
    """Rewrite the store's key file to hold the symbol table alone."""
    (root / store.KEYS).unlink()
    _files.write_record(root / store.KEYS, keytable.KeyTableWriter(store.SYMBOL).message().SerializeToString())


def _rewrite_manifest(root: pathlib.Path, manifest: store_pb2.StoreManifest) -> None:
    (root / store.MANIFEST).unlink()
    _files.write_record(root / store.MANIFEST, manifest.SerializeToString())


def test_a_key_file_rewritten_under_its_manifest_is_refused(tmp_path: pathlib.Path) -> None:
    # what an index rebuilt in place leaves when it is interrupted before its manifest is written
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    _only_symbols(root)
    with pytest.raises(store.StoreError, match=r'keys\.bagz: crc32c'):
        store.BundleStore(str(root))


def test_a_manifest_recording_no_index_digests_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    manifest = store.read_manifest(str(root))
    del manifest.index_files[:]
    _rewrite_manifest(root, manifest)
    with pytest.raises(store.StoreError, match='records no digest'):
        store.BundleStore(str(root))


def test_a_key_file_missing_a_table_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    _only_symbols(root)
    manifest = store.read_manifest(str(root))
    keys = next(f for f in manifest.index_files if f.path == store.KEYS)
    keys.digest = _files.crc32c_of_local(root / store.KEYS)
    _rewrite_manifest(root, manifest)
    with pytest.raises(store.StoreError, match='lacks'):
        store.BundleStore(str(root))


# ---- paths ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('root', 'part', 'expected'),
    [
        ('gs://b/ref/v/build1', 'manifest.bagz', 'gs://b/ref/v/build1/manifest.bagz'),
        ('gs://b/ref/v/build1', '../shards/s.bagz', 'gs://b/ref/v/shards/s.bagz'),
        ('/tmp/x/idx', '../shards/s.bagz', '/tmp/x/shards/s.bagz'),
    ],
)
def test_paths_join_relative_to_a_root_on_disk_or_in_a_bucket(root: str, part: str, expected: str) -> None:
    assert _files.join(root, part) == expected


@pytest.mark.parametrize(('root', 'part'), [('gs://b/a', '../../x'), ('/tmp', '/etc/passwd'), ('/tmp', 'gs://b/x')])
def test_paths_that_escape_or_are_absolute_are_refused(root: str, part: str) -> None:
    with pytest.raises(ValueError, match=r'climbs|absolute'):
        _files.join(root, part)


def test_crc32c_is_the_castagnoli_checksum(tmp_path: pathlib.Path) -> None:
    # the check value every CRC-32C implementation publishes for these nine bytes
    path = tmp_path / 'check'
    path.write_bytes(b'123456789')
    assert _files.crc32c_of_local(path) == 'e3069283'


def test_a_manifest_recording_more_records_than_its_shards_hold_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    manifest = store.read_manifest(str(root))
    manifest.shards[0].records += 1
    _rewrite_manifest(root, manifest)
    with pytest.raises(store.StoreError, match='shards hold 5 records, the manifest 6'):
        store.BundleStore(str(root))


def test_an_interval_table_on_another_assembly_is_refused(tmp_path: pathlib.Path) -> None:
    root = _write(tmp_path, [('RS_1', _BUNDLES)])
    table = index_pb2.ChromosomeIntervals(
        assembly=bundle_pb2.ASSEMBLY_GRCH37, chromosome='NC_000013.10', start=[0], end_inclusive=[9], record=[0]
    )
    (root / store.INTERVALS).unlink()
    _files.write_record(root / store.INTERVALS, table.SerializeToString())
    manifest = store.read_manifest(str(root))
    intervals = next(f for f in manifest.index_files if f.path == store.INTERVALS)
    intervals.digest = _files.crc32c_of_local(root / store.INTERVALS)
    _rewrite_manifest(root, manifest)
    with pytest.raises(store.StoreError, match='an interval table on'):
        store.BundleStore(str(root))


def test_an_index_over_a_shard_holding_other_than_its_recorded_count_is_refused(tmp_path: pathlib.Path) -> None:
    path = store_build.write_shard(_BUNDLES, tmp_path, release='RS_1', inputs=_inputs(tmp_path))
    miscounted = store_build.shard_record(path)
    miscounted.records += 1
    record = store_build.shard_record_path(path)
    record.unlink()
    _files.write_record(record, miscounted.SerializeToString())
    with pytest.raises(build.BuildError, match='different number of records'):
        store_build.write_index(tmp_path / 'store', [path], assembly=_ASSEMBLY)
