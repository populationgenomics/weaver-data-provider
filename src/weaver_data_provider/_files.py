"""Paths, single-record files and object digests, on local disk or `gs://`, through bagz where possible.

Metadata files (a store's manifest, a genome's catalogue) are single-record bagz files so they are read
through the same GCS client as the bundles themselves. Only a data file's digest (a shard, a genome's
blocks) needs object metadata, and only on `gs://`; that is the one place a Google Cloud Storage
client is used.
"""

from __future__ import annotations

import base64
import hashlib
import pathlib
import posixpath
from typing import TYPE_CHECKING

import bagz
import google_crc32c

if TYPE_CHECKING:  # the client doubles the package's import time, and only gs:// paths use it
    from google.cloud import storage

_ZSTD = bagz.Writer.Options(compression=bagz.CompressionZstd())


def join(root: str, *parts: str) -> str:
    """`root` joined with relative `parts`, `..` resolved, on a local path or a `gs://` URL alike.

    A manifest names its shards relative to itself, and a layout may keep shards beside the directory
    the manifest is in; nothing resolves `..` inside a URL, so it is resolved here.

    Raises:
        ValueError: If a part is absolute, or the result climbs above the bucket.
    """
    for part in parts:
        if part.startswith('/') or '://' in part:
            raise ValueError(f'{part!r} is absolute; a store names its files relative to itself')
    scheme, sep, rest = root.partition('://') if '://' in root else ('', '', root)
    joined = posixpath.normpath(posixpath.join(rest, *parts))
    # Under a URL the first component is the bucket; if normalising changed it, a `..` climbed out.
    if sep and joined.split('/')[0] != rest.split('/')[0]:
        raise ValueError(f'{root} joined with {parts} climbs out of the bucket')
    return f'{scheme}{sep}{joined}'


def is_remote(path: str) -> bool:
    return path.startswith('gs://')


def check_readable(path: str) -> None:
    """A `gs://` path needs the `gcs` extra; said before any read, so the message names the extra, not a backend.

    Raises:
        ImportError: If `path` is remote and the extra is not installed.
    """
    if not is_remote(path):
        return
    try:
        from google.cloud import storage  # noqa: F401 — installed by bagz[gcs], with bagz's GCS backend
    except ImportError as error:
        raise ImportError(f"{path}: reading gs:// needs the gcs extra: pip install 'hgvs-weaver-data[gcs]'") from error


def read_record(path: str) -> bytes:
    """The only record of a single-record bagz file.

    Raises:
        ValueError: If the file holds other than exactly one record.
    """
    check_readable(path)
    reader = bagz.Reader(path)
    if len(reader) != 1:
        raise ValueError(f'{path}: expected exactly one record, found {len(reader)}')
    return bytes(reader[0])


def write_record(path: pathlib.Path, data: bytes) -> None:
    with bagz.Writer(str(path), _ZSTD) as writer:
        writer.write(data)


def md5_of_local(path: pathlib.Path) -> str:
    h = hashlib.md5()  # noqa: S324 — identifies a build input, not a security digest
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def crc32c_of_local(path: pathlib.Path) -> str:
    """A file's crc32c as 8 hex characters: the checksum GCS records on every object."""
    h = google_crc32c.Checksum()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.digest().hex()


def crc32c_of(path: str, storage_client: storage.Client | None = None) -> str:
    """A file's crc32c: computed for a local file, read from object metadata on `gs://`.

    crc32c rather than md5, because GCS records no md5 for a composite object, which is what
    `gcloud storage cp` makes of a large upload.

    Raises:
        FileNotFoundError: If a `gs://` object does not exist.
        ValueError: If a `gs://` object records no crc32c.
    """
    if not is_remote(path):
        return crc32c_of_local(pathlib.Path(path))
    check_readable(path)
    from google.cloud import storage  # imported on first use, for the reason the TYPE_CHECKING import gives

    client = storage_client if storage_client is not None else storage.Client()
    bucket, _, name = path.removeprefix('gs://').partition('/')
    blob = client.bucket(bucket).get_blob(name)
    if blob is None:
        raise FileNotFoundError(path)
    if blob.crc32c is None:
        raise ValueError(f'{path}: the object records no crc32c')
    return base64.b64decode(blob.crc32c).hex()
