"""Cut a genome into blocks of bases, and catalogue its sequences.

A publisher's gzip FASTA cannot be sliced, so this writes every sequence as fixed-size blocks of
bases, one zstd-compressed bagz record each, and a catalogue naming every sequence with its length and
refget digest. Residues are written as read, soft masking kept. The catalogue is written last, which
makes its presence the mark of a complete genome.
"""

from __future__ import annotations

import gzip
import hashlib
import pathlib
from collections.abc import Iterator

import bagz

import weaver_data_provider
from weaver_data_provider import _files, build, refget
from weaver_data_provider import genome as genome_mod
from weaver_data_provider.v1 import bundle_pb2, genome_pb2, store_pb2

BLOCK_SIZE = 64 << 10
# Decompression costs the same at any level, so the slow compression is paid once, at build.
_OPTIONS = bagz.Writer.Options(compression=bagz.CompressionZstd(level=19))


def _parse_header(header: bytes) -> tuple[str, str]:
    name, _, description = header[1:].decode('ascii').strip().partition(' ')
    return name, description


def _records(path: pathlib.Path) -> Iterator[tuple[str, str, list[bytes]]]:
    """Each FASTA record as (name, description, residue lines).

    Raises:
        build.BuildError: If residues come before the first header.
    """
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as fh:
        header: bytes | None = None
        lines: list[bytes] = []
        for raw in fh:
            if raw.startswith(b'>'):
                if header is not None:
                    yield (*_parse_header(header), lines)
                header, lines = raw, []
            elif header is None:
                raise build.BuildError(f'{path}: residues before the first header')
            else:
                lines.append(raw.strip())
        if header is not None:
            yield (*_parse_header(header), lines)


def build_genome(
    source: pathlib.Path, out: pathlib.Path, *, assembly: bundle_pb2.Assembly, block_size: int = BLOCK_SIZE
) -> genome_pb2.SequenceCatalogue:
    """Write `out/genome.bagz` from a plain or gzip FASTA, then its catalogue.

    Args:
        source: The publisher's FASTA; gzipped when its name ends in `.gz`.
        out: The genome's directory, created if absent.
        assembly: The assembly the FASTA is, recorded in the catalogue.
        block_size: Bases per block. A slice decompresses every block it touches, so larger blocks
            compress better and slice slower.

    Raises:
        build.BuildError: If the source is a tar archive, or holds an empty or duplicate sequence.
    """
    if source.name.endswith(('.tar.gz', '.tgz', '.tar')):
        raise build.BuildError(f'{source}: a tar archive; extract the FASTA first')
    out.mkdir(parents=True, exist_ok=True)
    catalogue = genome_pb2.SequenceCatalogue(
        format_version=weaver_data_provider.FORMAT_VERSION,
        assembly=assembly,
        source=store_pb2.Input(role='genome', digest=_files.md5_of_local(source), digest_algorithm='md5'),
        builder=build.BUILDER,
        block_size=block_size,
        digest_algorithm='crc32c',
    )
    names: set[str] = set()
    blocks = out / genome_mod.BLOCKS
    with bagz.Writer(str(blocks), _OPTIONS) as writer:
        for name, description, lines in _records(source):
            if name in names:
                raise build.BuildError(f'{source}: {name} appears twice')
            names.add(name)
            hasher = hashlib.sha512()
            length = 0
            carry = bytearray()
            for line in lines:
                hasher.update(line.upper())
                length += len(line)
                carry += line
                while len(carry) >= block_size:
                    writer.write(bytes(carry[:block_size]))
                    del carry[:block_size]
            if length == 0:
                raise build.BuildError(f'{name}: an empty sequence')
            if carry:
                writer.write(bytes(carry))
            catalogue.sequences.add(
                name=name,
                length=length,
                digest=refget.identifier(hasher.digest()[: refget.RAW_LENGTH]),
                description=description,
            )
    catalogue.digest = _files.crc32c_of_local(blocks)
    build.validated(catalogue, f'{source}: catalogue')
    _files.write_record(out / genome_mod.CATALOGUE, catalogue.SerializeToString())
    return catalogue
