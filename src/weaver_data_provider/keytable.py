"""Sorted key tables: a store's indexes, written as one record each and searched in memory.

A table is every key in bytewise order, concatenated, with offset arrays into the keys and into
the bundle record numbers they map to. Loaded, it is three byte strings viewed as arrays; a lookup
is a binary search that slices one key out per comparison, and a prefix query is the run of keys
between two such searches.
"""

from __future__ import annotations

import array
import bisect
import collections
import sys
from collections.abc import Iterable

from weaver_data_provider.v1 import index_pb2


def _uint32s(values: Iterable[int]) -> bytes:
    if sys.byteorder != 'little':
        raise RuntimeError('key tables are written little-endian; this host is big-endian')
    return array.array('I', values).tobytes()


class KeyTableWriter:
    """Collects key -> record numbers for one kind of key, then emits the sorted table."""

    def __init__(self, kind: str) -> None:
        self.kind = kind
        self._records: dict[bytes, set[int]] = collections.defaultdict(set)

    def add(self, key: bytes, records: Iterable[int]) -> None:
        if not key:
            raise ValueError(f'{self.kind}: an empty key indexes nothing')
        self._records[key].update(records)

    def __len__(self) -> int:
        return len(self._records)

    def message(self) -> index_pb2.KeyTable:
        keys = sorted(self._records)
        key_offsets, record_offsets, records = [0], [0], []
        for key in keys:
            key_offsets.append(key_offsets[-1] + len(key))
            records.extend(sorted(self._records[key]))
            record_offsets.append(len(records))
        return index_pb2.KeyTable(
            kind=self.kind,
            keys=b''.join(keys),
            key_offsets=_uint32s(key_offsets),
            records=_uint32s(records),
            record_offsets=_uint32s(record_offsets),
        )


class KeyTable:
    """A loaded table; lookups are binary searches over the concatenated keys."""

    def __init__(self, message: index_pb2.KeyTable) -> None:
        if sys.byteorder != 'little':
            raise RuntimeError('key tables are little-endian; this host is big-endian')
        self.kind = message.kind
        self._keys = message.keys
        self._key_offsets = memoryview(message.key_offsets).cast('I')
        self._records = memoryview(message.records).cast('I')
        self._record_offsets = memoryview(message.record_offsets).cast('I')
        if len(self._key_offsets) != len(self._record_offsets) or len(self._key_offsets) == 0:
            raise ValueError(f'{self.kind}: offset arrays are inconsistent')

    def __len__(self) -> int:
        return len(self._key_offsets) - 1

    def key(self, i: int) -> bytes:
        return self._keys[self._key_offsets[i] : self._key_offsets[i + 1]]

    def records(self, i: int) -> tuple[int, ...]:
        return tuple(self._records[self._record_offsets[i] : self._record_offsets[i + 1]])

    def _first_at_or_after(self, key: bytes) -> int:
        return bisect.bisect_left(range(len(self)), key, key=self.key)

    def lookup(self, key: bytes) -> tuple[int, ...]:
        """The record numbers under exactly `key`, ascending; empty when absent."""
        i = self._first_at_or_after(key)
        if i < len(self) and self.key(i) == key:
            return self.records(i)
        return ()

    def prefixed(self, prefix: bytes) -> tuple[int, ...]:
        """The record numbers under every key starting with `prefix`, ascending and distinct."""
        found: set[int] = set()
        i = self._first_at_or_after(prefix)
        while i < len(self) and self.key(i).startswith(prefix):
            found.update(self.records(i))
            i += 1
        return tuple(sorted(found))
