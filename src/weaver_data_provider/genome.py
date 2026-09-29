"""Read a derived genome: any slice of any sequence, and the refget identity of each.

A genome is a directory or `gs://` prefix holding a blocks file — every sequence cut into fixed-size
blocks of bases, one zstd-compressed bagz record each — and a catalogue naming every sequence with its
length and refget digest, written by `weaver_data_provider.build`. A slice is the blocks it touches,
read in one batch; `fetch_many` gathers several slices of one sequence, a transcript's exons say, in
one batch too. Decoded blocks are kept in a bounded cache, since the reads one variant needs cluster
on a few blocks.
"""

from __future__ import annotations

import collections
import threading
from collections.abc import Iterable, Sequence
from typing import TYPE_CHECKING

import bagz

import weaver_data_provider
from weaver_data_provider import _files
from weaver_data_provider.v1 import genome_pb2

if TYPE_CHECKING:  # the client doubles the package's import time, and only gs:// paths use it
    from google.cloud import storage

BLOCKS = 'genome.bagz'
CATALOGUE = 'sequences.bagz'


class GenomeError(ValueError):
    """A genome this version cannot read, or one that disagrees with its own catalogue."""


def read_catalogue(root: str) -> genome_pb2.SequenceCatalogue:
    catalogue = genome_pb2.SequenceCatalogue()
    catalogue.ParseFromString(_files.read_record(_files.join(root, CATALOGUE)))
    if catalogue.format_version != weaver_data_provider.FORMAT_VERSION:
        raise GenomeError(
            f'{root}: format version {catalogue.format_version}; '
            f'this package reads version {weaver_data_provider.FORMAT_VERSION}'
        )
    return catalogue


class Genome:
    """One assembly's sequences, sliceable by name, with their lengths and refget digests.

    Safe to share between threads.
    """

    def __init__(self, root: str, *, cache_blocks: int = 256, storage_client: storage.Client | None = None) -> None:
        """Open a genome and check its blocks file against its catalogue.

        Args:
            root: The genome's directory or `gs://` prefix.
            cache_blocks: Decoded blocks kept, least recently used evicted first.
            storage_client: A `google.cloud.storage.Client`, used on `gs://` only, to read the blocks
                file's recorded digest. One is made when omitted.

        Raises:
            GenomeError: If the catalogue's format version is unknown, or the blocks file's digest or
                block count disagrees with the catalogue.
        """
        if cache_blocks < 1:
            raise ValueError('cache_blocks must be at least 1')
        self.root = root.rstrip('/')
        catalogue = read_catalogue(self.root)
        self.assembly = catalogue.assembly
        self.block_size = catalogue.block_size
        path = _files.join(self.root, BLOCKS)
        if catalogue.digest_algorithm != 'crc32c':
            raise GenomeError(f'{path}: cannot check a {catalogue.digest_algorithm} digest')
        if (found := _files.crc32c_of(path, storage_client)) != catalogue.digest:
            raise GenomeError(f'{path}: crc32c {found}, the catalogue records {catalogue.digest}')
        self._blocks = bagz.Reader(path)

        self._entries = {entry.name: entry for entry in catalogue.sequences}
        self._by_digest = {entry.digest: entry.name for entry in catalogue.sequences}
        self._first_block: dict[str, int] = {}
        total = 0
        for entry in catalogue.sequences:
            self._first_block[entry.name] = total
            total += -(-entry.length // self.block_size)
        if len(self._blocks) != total:
            raise GenomeError(f'{path}: holds {len(self._blocks)} blocks, the catalogue implies {total}')

        self._cache: collections.OrderedDict[int, str] = collections.OrderedDict()
        self._cache_blocks = cache_blocks
        # Held only around cache bookkeeping, never across a read.
        self._lock = threading.Lock()

    @property
    def names(self) -> list[str]:
        return list(self._entries)

    def __contains__(self, name: str) -> bool:
        return name in self._entries

    def length(self, name: str) -> int:
        return self._entries[name].length

    def digest(self, name: str) -> str:
        """The `SQ.` refget identifier of a sequence."""
        return self._entries[name].digest

    def name_for_digest(self, digest: str) -> str | None:
        return self._by_digest.get(digest)

    def _read(self, records: Iterable[int]) -> dict[int, str]:
        """These blocks, upper-cased; the ones not cached are read in one batch."""
        wanted = list(dict.fromkeys(records))
        found: dict[int, str] = {}
        with self._lock:
            for record in wanted:
                cached = self._cache.get(record)
                if cached is not None:
                    self._cache.move_to_end(record)
                    found[record] = cached
        missing = [r for r in wanted if r not in found]
        if missing:
            read = {
                record: bytes(raw).decode('ascii').upper()
                for record, raw in zip(missing, self._blocks.read_indices(missing), strict=True)
            }
            found.update(read)
            with self._lock:
                self._cache.update(read)
                while len(self._cache) > self._cache_blocks:
                    self._cache.popitem(last=False)
        return found

    def fetch_many(self, name: str, ranges: Sequence[tuple[int, int | None]]) -> list[str]:
        """Residues of `name` over each range, with every block the ranges touch read in one batch.

        Each range is 0-based and half-open; an end of None, or a negative one, runs to the end of the
        sequence. Residues are upper-cased and clipped to the sequence at both ends, so a negative start
        reads from its first base.

        Raises:
            KeyError: If the catalogue has no such sequence.
        """
        length = self._entries[name].length
        first, size = self._first_block[name], self.block_size
        clipped = [(max(start, 0), length if end is None or end < 0 else min(end, length)) for start, end in ranges]
        blocks = self._read(
            first + b for start, stop in clipped if start < stop for b in range(start // size, (stop - 1) // size + 1)
        )
        out = []
        for start, stop in clipped:
            if start >= stop:
                out.append('')
                continue
            lo, hi = start // size, (stop - 1) // size
            joined = ''.join(blocks[first + b] for b in range(lo, hi + 1))
            out.append(joined[start - lo * size : stop - lo * size])
        return out

    def fetch(self, name: str, start: int, end: int | None = None) -> str:
        """Residues of `name` over the 0-based half-open range, upper-cased; clipped to the sequence.

        Raises:
            KeyError: If the catalogue has no such sequence.
        """
        return self.fetch_many(name, [(start, end)])[0]
