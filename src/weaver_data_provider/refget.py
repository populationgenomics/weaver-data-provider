"""The refget sequence identifier: sha512t24u over the normalised residues, and its raw form."""

from __future__ import annotations

import base64
import hashlib

PREFIX = 'SQ.'
RAW_LENGTH = 24


def digest(residues: bytes) -> str:
    """The `SQ.` identifier refget and VRS give a sequence.

    Args:
        residues: The sequence as ASCII residues; case is normalised, so a soft-masked FASTA digests
            as its uppercase form, which is the form every refget server has indexed.
    """
    return PREFIX + base64.urlsafe_b64encode(hashlib.sha512(residues.upper()).digest()[:RAW_LENGTH]).decode('ascii')


def raw(identifier: str) -> bytes:
    """The 24 digest bytes behind an `SQ.` identifier, the fixed-width form an index keys on.

    Raises:
        ValueError: If the identifier is not an `SQ.` digest.
    """
    if not identifier.startswith(PREFIX):
        raise ValueError(f'{identifier!r} is not a refget SQ identifier')
    decoded = base64.urlsafe_b64decode(identifier[len(PREFIX) :] + '=')
    if len(decoded) != RAW_LENGTH:
        raise ValueError(f'{identifier!r} does not decode to {RAW_LENGTH} bytes')
    return decoded


def identifier(raw_digest: bytes) -> str:
    """The `SQ.` identifier for 24 digest bytes; the inverse of `raw`."""
    if len(raw_digest) != RAW_LENGTH:
        raise ValueError(f'a refget digest is {RAW_LENGTH} bytes, not {len(raw_digest)}')
    return PREFIX + base64.urlsafe_b64encode(raw_digest).decode('ascii').rstrip('=')
