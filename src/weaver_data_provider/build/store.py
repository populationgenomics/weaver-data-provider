"""Write a store: bundle shards, and the index and manifest over an ordered list of them.

A shard is one release's bundles, written once and never rewritten. An index is derived: it is
rebuilt whole by reading its shards back in the order given, so stacking a new release is writing
one shard and rebuilding the index over the longer list. The manifest is written last and records
the digest of every shard and index file, so a reader refuses a store whose files are not the ones
its manifest was written over — an index rebuilt in place and interrupted before its manifest.
"""

from __future__ import annotations

import collections
import itertools
import os
import pathlib
import re
from collections.abc import Iterable, Sequence

import bagz

import weaver_data_provider
from weaver_data_provider import _files, build, keytable, refget
from weaver_data_provider import store as store_mod
from weaver_data_provider.v1 import bundle_pb2, index_pb2, store_pb2

_ZSTD = bagz.Writer.Options(compression=bagz.CompressionZstd())


def _input(role: str, path: pathlib.Path) -> store_pb2.Input:
    return store_pb2.Input(role=role, digest=_files.md5_of_local(path), digest_algorithm='md5')


def shard_record_path(path: pathlib.Path) -> pathlib.Path:
    """Where a shard's record sits: `RS_1.bagz` has `RS_1.shard.bagz` beside it.

    bagz decides whether to decompress from the file name when it reads, and only decompresses a name
    ending in `.bagz`; a record written as `RS_1.bagz.shard` would read back as compressed bytes.
    """
    return path.with_name(f'{path.stem}.shard.bagz')


def shard_record(path: pathlib.Path) -> store_pb2.Shard:
    """The record `write_shard` left beside a shard: its release, record count, digest and inputs."""
    shard = store_pb2.Shard()
    shard.ParseFromString(_files.read_record(str(shard_record_path(path))))
    return shard


_RELEASE_NAME = re.compile(r'^[A-Za-z0-9._-]+$')


def _check_digests_resolve(bundle: bundle_pb2.GeneBundle, what: str) -> None:
    """Every digest a transcript or protein cites is one of the bundle's sequences.

    Raises:
        build.BuildError: Naming the accessions whose digests the bundle does not carry.
    """
    carried = {sequence.digest for sequence in bundle.sequences}
    dangling = [
        f'{item.accession}.{item.version}'
        for item in (*bundle.transcripts, *bundle.proteins)
        if item.sequence_digest not in carried
    ]
    if dangling:
        raise build.BuildError(f'{what}: cites sequences it does not carry, for {dangling}')


_KIND_ORDER = list(build.SEQUENCE_KINDS)


def _check_alignment_order(bundle: bundle_pb2.GeneBundle, what: str) -> None:
    """Each transcript's alignments are one per sequence: chromosomes first, then scaffolds, then patches, by accession.

    Raises:
        build.BuildError: Naming a transcript whose alignments are on a sequence of no known kind, on one sequence
            twice, or not in that order.
    """
    for transcript in bundle.transcripts:
        names = [a.chromosome for a in transcript.alignments]
        if any(name[:3] not in build.SEQUENCE_KINDS for name in names):
            raise build.BuildError(f'{what}: {transcript.accession}.{transcript.version} is aligned on {names}')
        keys = [(_KIND_ORDER.index(name[:3]), name) for name in names]
        if any(a >= b for a, b in itertools.pairwise(keys)):
            raise build.BuildError(
                f'{what}: {transcript.accession}.{transcript.version} alignments are ordered {names}, not one per '
                'sequence, chromosomes first, then scaffolds, then patches, each by accession'
            )


def write_shard(
    bundles: Iterable[bundle_pb2.GeneBundle],
    directory: pathlib.Path,
    *,
    release: str,
    inputs: Sequence[tuple[str, pathlib.Path]],
    prefix: str = '',
) -> pathlib.Path:
    """Write one release's bundles as a shard in `directory`, each checked against the bundle schema.

    The shard is named for its contents, `<prefix><release>-<crc32c>.bagz`, so a release cut again —
    a new builder, a fixed encoding — lands at a new path rather than over the one an index already
    names. Its record (release, record count, digest, inputs) is written beside it, so an index built
    later, over any list of shards, needs only their paths.

    Args:
        bundles: The release's bundles, in the order they are to be numbered.
        directory: Where the shard is written, created if absent.
        release: The publisher's release identifier, which begins the shard's name.
        inputs: Each file the release was cut from, by role, recorded with its md5.
        prefix: Prepended to the name, for a caller whose layout orders shards by name (`02-`).

    Returns:
        The shard's path.

    Raises:
        build.BuildError: If a bundle violates the schema, or the release is empty.
        FileExistsError: If a shard with these exact contents is already there.
        ValueError: If `release` or `prefix` would not make a plain file name.
    """
    if not _RELEASE_NAME.match(release) or (prefix and not _RELEASE_NAME.match(prefix)):
        raise ValueError(f'{prefix}{release!r} is not a plain file name: use letters, digits, ".", "_" and "-" only')
    directory.mkdir(parents=True, exist_ok=True)
    # Written under a name no final shard can have, then renamed once its digest is known.
    pending = directory / f'.{prefix}{release}.pending.bagz'
    records = 0
    try:
        with bagz.Writer(str(pending), _ZSTD) as writer:
            for bundle in bundles:
                what = f'{release} bundle {records} ({bundle.gene.symbol})'
                build.validated(bundle, what)
                _check_digests_resolve(bundle, what)
                _check_alignment_order(bundle, what)
                writer.write(bundle.SerializeToString())
                records += 1
        if records == 0:
            raise build.BuildError(f'{release}: no bundles to write')
        digest = _files.crc32c_of_local(pending)
        path = directory / f'{prefix}{release}-{digest}.bagz'
        if path.exists():
            raise FileExistsError(f'{path}: a shard with these contents is already written')
        pending.rename(path)
    finally:
        pending.unlink(missing_ok=True)
    shard = store_pb2.Shard(
        path=path.name,
        release=release,
        records=records,
        digest=digest,
        digest_algorithm='crc32c',
        inputs=[_input(role, p) for role, p in inputs],
    )
    build.validated(shard, f'{path}: shard record')
    _files.write_record(shard_record_path(path), shard.SerializeToString())
    return path


def _index_bundle(
    bundle: bundle_pb2.GeneBundle,
    record: int,
    tables: dict[str, keytable.KeyTableWriter],
    intervals: dict[str, list[tuple[int, int, int]]],
    *,
    assembly: bundle_pb2.Assembly,
) -> None:
    """Key one bundle under its symbols, accessions, proteins and digests, and place its alignments."""
    gene = bundle.gene
    for name in (gene.symbol, *gene.previous_symbols, *gene.alias_symbols):
        if key := store_mod.symbol_key(name):
            tables[store_mod.SYMBOL].add(key.encode('ascii'), [record])
    for transcript in bundle.transcripts:
        tables[store_mod.ACCESSION].add(f'{transcript.accession}.{transcript.version}'.encode('ascii'), [record])
        for alignment in transcript.alignments:
            if alignment.assembly != assembly:
                name = bundle_pb2.Assembly.Name
                raise build.BuildError(
                    f'{transcript.accession}.{transcript.version}: aligned on '
                    f'{name(alignment.assembly)} in a store on {name(assembly)}'
                )
            span_start = min(e.genome_start for e in alignment.exons)
            span_end = max(e.genome_end_inclusive for e in alignment.exons)
            intervals[alignment.chromosome].append((span_start, span_end, record))
    for protein in bundle.proteins:
        tables[store_mod.PROTEIN].add(f'{protein.accession}.{protein.version}'.encode('ascii'), [record])
    for sequence in bundle.sequences:
        tables[store_mod.DIGEST].add(refget.raw(sequence.digest), [record])


def _write_tables(
    out: pathlib.Path,
    tables: dict[str, keytable.KeyTableWriter],
    intervals: dict[str, list[tuple[int, int, int]]],
    *,
    assembly: bundle_pb2.Assembly,
) -> None:
    with bagz.Writer(str(out / store_mod.KEYS), _ZSTD) as writer:
        for kind in store_mod.KINDS:
            message = tables[kind].message()
            build.validated(message, f'the {kind} key table')
            writer.write(message.SerializeToString())
    with bagz.Writer(str(out / store_mod.INTERVALS), _ZSTD) as writer:
        for chromosome, rows in sorted(intervals.items()):
            rows.sort()
            table = index_pb2.ChromosomeIntervals(
                assembly=assembly,
                chromosome=chromosome,
                start=[r[0] for r in rows],
                end_inclusive=[r[1] for r in rows],
                record=[r[2] for r in rows],
            )
            build.validated(table, f'intervals on {chromosome}')
            writer.write(table.SerializeToString())


def write_index(
    out: pathlib.Path, shards: Sequence[pathlib.Path], *, assembly: bundle_pb2.Assembly
) -> store_pb2.StoreManifest:
    """Write the key tables, interval tables and manifest in `out` over these shards, in this order.

    Each shard is read with the record `write_shard` left beside it. Record numbers run across the
    shards in the order given; the manifest names each by its path relative to `out`. Accessions and
    proteins are keyed by their versioned form only; the reader answers a bare accession with a prefix
    search over the versions.

    Raises:
        build.BuildError: If a shard's bytes no longer match its digest, or an alignment is on another assembly.
    """
    if not shards:
        raise build.BuildError('an index spans at least one shard')
    out.mkdir(parents=True, exist_ok=True)
    tables = {kind: keytable.KeyTableWriter(kind) for kind in store_mod.KINDS}
    intervals: dict[str, list[tuple[int, int, int]]] = collections.defaultdict(list)
    manifest = store_pb2.StoreManifest(
        format_version=weaver_data_provider.FORMAT_VERSION, assembly=assembly, builder=build.BUILDER
    )
    record = 0
    for path in shards:
        shard = shard_record(path)
        if (found := _files.crc32c_of_local(path)) != shard.digest:
            raise build.BuildError(f'{path}: crc32c {found}, but the shard was recorded as {shard.digest}')
        entry = manifest.shards.add()
        entry.CopyFrom(shard)
        entry.path = os.path.relpath(path, out)
        for raw in bagz.Reader(str(path)):
            bundle = bundle_pb2.GeneBundle()
            bundle.ParseFromString(raw)
            _index_bundle(bundle, record, tables, intervals, assembly=assembly)
            record += 1
        if record != sum(s.records for s in manifest.shards):
            raise build.BuildError(f'{path}: holds a different number of records than recorded ({shard.records})')
    _write_tables(out, tables, intervals, assembly=assembly)
    for name in (store_mod.KEYS, store_mod.INTERVALS):
        manifest.index_files.add(path=name, digest=_files.crc32c_of_local(out / name), digest_algorithm='crc32c')
    build.validated(manifest, 'manifest')
    _files.write_record(out / store_mod.MANIFEST, manifest.SerializeToString())
    return manifest
