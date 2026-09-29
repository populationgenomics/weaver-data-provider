"""Read a store: gene bundles by symbol, accession, protein, digest or position.

A store is a directory or `gs://` prefix holding a manifest, the bundle shards it names, and the key
and interval tables over them, written by `weaver_data_provider.build`. Open reads the manifest,
refuses a format version it does not know, checks every shard and both tables against the digests the
manifest records and the shards against its record count, and loads the tables whole. A lookup is
then an in-memory search and one ranged read.
"""

from __future__ import annotations

import bisect
import collections
import itertools
import re
import threading
from collections.abc import Sequence
from typing import TYPE_CHECKING

import bagz

import weaver_data_provider
from weaver_data_provider import _files, keytable, refget
from weaver_data_provider.v1 import bundle_pb2, index_pb2, store_pb2

if TYPE_CHECKING:  # the client doubles the package's import time, and only gs:// paths use it
    from google.cloud import storage

MANIFEST = 'manifest.bagz'
KEYS = 'keys.bagz'
INTERVALS = 'intervals.bagz'
SYMBOL, ACCESSION, PROTEIN, DIGEST = 'symbol', 'accession', 'protein', 'digest'
KINDS = (SYMBOL, ACCESSION, PROTEIN, DIGEST)
_NOT_SYMBOL = re.compile(r'[^A-Z0-9]')
# Record limits in memory: a bundle read is then one request, not a limits read and a data read.
_BUNDLE_OPTIONS = bagz.Reader.Options(limits_storage=bagz.LimitsStorage.IN_MEMORY)


class StoreError(ValueError):
    """A store that does not match its own manifest, or a manifest this version cannot read."""


def symbol_key(symbol: str) -> str:
    """A gene symbol in the shape the index stores it: upper-case ASCII alphanumerics only.

    Papers print `DNM1l` and `XP-C` for symbols stored as `DNM1L` and `XPC`, and an HGNC alias that
    ends in a Greek letter reduces to its ASCII stem; a caller normalises what a paper prints the same way.
    """
    return _NOT_SYMBOL.sub('', symbol.upper())


def read_manifest(root: str) -> store_pb2.StoreManifest:
    """A store's manifest, refused if its format version is not the one this package reads."""
    manifest = store_pb2.StoreManifest()
    manifest.ParseFromString(_files.read_record(_files.join(root, MANIFEST)))
    if manifest.format_version != weaver_data_provider.FORMAT_VERSION:
        raise StoreError(
            f'{root}: format version {manifest.format_version}; '
            f'this package reads version {weaver_data_provider.FORMAT_VERSION}'
        )
    return manifest


def _check_digest(path: str, digest: str, algorithm: str, storage_client: storage.Client | None) -> None:
    if algorithm != 'crc32c':
        raise StoreError(f'{path}: cannot check a {algorithm} digest')
    if (found := _files.crc32c_of(path, storage_client)) != digest:
        raise StoreError(f'{path}: crc32c {found}, the manifest records {digest}')


class _Intervals:
    """One chromosome's placements sorted by start, with the running maximum of their ends.

    A position query bisects twice: `starts` bounds the placements that begin at or before the
    position, and the running maximum bounds those that could still reach it, so no scan runs over
    the placements that ended long before.
    """

    def __init__(self, table: index_pb2.ChromosomeIntervals) -> None:
        self.starts = list(table.start)
        self.ends = list(table.end_inclusive)
        self.records = list(table.record)
        self.max_end_so_far = list(itertools.accumulate(self.ends, max))

    def overlapping(self, start: int, end: int) -> list[int]:
        """Records of placements overlapping the 0-based closed interval [start, end]."""
        hi = bisect.bisect_right(self.starts, end)
        lo = bisect.bisect_left(self.max_end_so_far, start, 0, hi)
        return sorted({self.records[k] for k in range(lo, hi) if self.ends[k] >= start})


class BundleStore:
    """One assembly's gene bundles and the tables over them. Safe to share between threads."""

    def __init__(self, root: str, *, cache_size: int = 4096, storage_client: storage.Client | None = None) -> None:
        """Open a store and check it against its manifest.

        Args:
            root: The store's directory or `gs://` prefix.
            cache_size: Bundles kept decoded, by record number.
            storage_client: A `google.cloud.storage.Client`, used on `gs://` only, to read each shard's
                recorded crc32c. One is made when omitted.

        Raises:
            StoreError: If the manifest's format version is unknown, or a shard's or index file's
                digest, or the record count, disagrees with the manifest.
        """
        if cache_size < 1:
            raise ValueError('cache_size must be at least 1')
        self.root = root.rstrip('/')
        self.manifest = read_manifest(self.root)
        self.assembly = self.manifest.assembly

        paths = [_files.join(self.root, shard.path) for shard in self.manifest.shards]
        for shard, path in zip(self.manifest.shards, paths, strict=True):
            _check_digest(path, shard.digest, shard.digest_algorithm, storage_client)
        index_files = {f.path: f for f in self.manifest.index_files}
        if missing_files := {KEYS, INTERVALS} - index_files.keys():
            raise StoreError(f'{self.root}: the manifest records no digest for {sorted(missing_files)}')
        for name in (KEYS, INTERVALS):
            entry = index_files[name]
            _check_digest(_files.join(self.root, name), entry.digest, entry.digest_algorithm, storage_client)
        # One reader over every shard, in manifest order; record numbers run across them in that order.
        self._bundles = bagz.Reader(','.join(paths), _BUNDLE_OPTIONS)
        expected = sum(shard.records for shard in self.manifest.shards)
        if len(self._bundles) != expected:
            raise StoreError(f'{self.root}: shards hold {len(self._bundles)} records, the manifest {expected}')
        self._shard_ends = list(itertools.accumulate(shard.records for shard in self.manifest.shards))

        self._tables: dict[str, keytable.KeyTable] = {}
        for raw in bagz.Reader(_files.join(self.root, KEYS)):
            message = index_pb2.KeyTable()
            message.ParseFromString(raw)
            self._tables[message.kind] = keytable.KeyTable(message)
        if missing := set(KINDS) - self._tables.keys():
            raise StoreError(f'{self.root}: {KEYS} lacks the {sorted(missing)} tables')

        self._intervals: dict[str, _Intervals] = {}
        for raw in bagz.Reader(_files.join(self.root, INTERVALS)):
            table = index_pb2.ChromosomeIntervals()
            table.ParseFromString(raw)
            if table.assembly != self.assembly:
                raise StoreError(f'{self.root}: an interval table on {table.assembly} in a store on {self.assembly}')
            self._intervals[table.chromosome] = _Intervals(table)

        # Serialized, so no caller holds a message another caller is also given.
        self._cache: collections.OrderedDict[int, bytes] = collections.OrderedDict()
        self._cache_size = cache_size
        # Held only around cache bookkeeping, never across a read.
        self._lock = threading.Lock()

    def __len__(self) -> int:
        return len(self._bundles)

    def release_of(self, record: int) -> str:
        """The annotation release a record number came from, via the shard it falls in."""
        if not 0 <= record < len(self):
            raise IndexError(f'record {record} is outside a store of {len(self)}')
        return self.manifest.shards[bisect.bisect_right(self._shard_ends, record)].release

    def bundles(self, records: Sequence[int]) -> list[bundle_pb2.GeneBundle]:
        """The bundles at these record numbers, each a message of the caller's own to keep or change.

        The ones not cached are read in one batch.
        """
        wanted = list(dict.fromkeys(records))
        found: dict[int, bytes] = {}
        with self._lock:
            for record in wanted:
                cached = self._cache.get(record)
                if cached is not None:
                    self._cache.move_to_end(record)
                    found[record] = cached
        missing = [r for r in wanted if r not in found]
        if missing:
            read = {r: bytes(raw) for r, raw in zip(missing, self._bundles.read_indices(missing), strict=True)}
            found.update(read)
            with self._lock:
                for record, raw in read.items():
                    self._cache[record] = raw
                    self._cache.move_to_end(record)
                while len(self._cache) > self._cache_size:
                    self._cache.popitem(last=False)
        return [bundle_pb2.GeneBundle.FromString(found[r]) for r in records]

    def _versioned_or_all(self, kind: str, accession: str) -> tuple[int, ...]:
        """Records under an exact `NM_000059.4`, or under every version of a bare `NM_000059`."""
        accession = accession.strip()
        if not accession.isascii():
            return ()  # every key is ASCII
        if '.' in accession:
            return self._tables[kind].lookup(accession.encode('ascii'))
        return self._tables[kind].prefixed(f'{accession}.'.encode('ascii'))

    def by_symbol(self, symbol: str) -> list[bundle_pb2.GeneBundle]:
        """Bundles whose approved, previous or alias symbol matches, symbol normalised.

        More than one is a real answer, not an error: a spelling can be one gene's symbol and another's alias.
        """
        key = symbol_key(symbol)
        return self.bundles(self._tables[SYMBOL].lookup(key.encode('ascii'))) if key else []

    def by_accession(self, accession: str) -> list[bundle_pb2.GeneBundle]:
        """Bundles holding a transcript with this accession, versioned or not."""
        return self.bundles(self._versioned_or_all(ACCESSION, accession))

    def by_protein(self, accession: str) -> list[bundle_pb2.GeneBundle]:
        """Bundles holding a protein with this accession, versioned or not."""
        return self.bundles(self._versioned_or_all(PROTEIN, accession))

    def by_digest(self, digest: str) -> list[bundle_pb2.GeneBundle]:
        """Bundles holding a sequence with this refget `SQ.` digest."""
        return self.bundles(self._tables[DIGEST].lookup(refget.raw(digest)))

    def overlapping(self, chromosome: str, start: int, end: int | None = None) -> list[bundle_pb2.GeneBundle]:
        """Bundles with a transcript placed over a 0-based position, or the closed interval [start, end].

        A chromosome with no placements answers empty; the store is on one assembly, so there is no
        other assembly for a caller to have meant.
        """
        if end is None:
            end = start
        if end < start:
            raise ValueError(f'interval end {end} precedes start {start}')
        table = self._intervals.get(chromosome)
        return self.bundles(table.overlapping(start, end)) if table is not None else []
