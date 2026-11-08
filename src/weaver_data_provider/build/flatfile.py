"""Read INSDC flat files, GenBank's record format: each record's CDS in its own coordinates, sequence and protein.

RefSeq's transcript records and GenBank's submitted ones are the same format, so one reader serves the RefSeq builder,
which takes each transcript's CDS and, for the historical set, its sequence and protein from the record, and the
GenBank builder, which takes everything from it.
"""

from __future__ import annotations

import dataclasses
import gzip
import pathlib
import re

from weaver_data_provider import build
from weaver_data_provider.build import common
from weaver_data_provider.v1 import bundle_pb2

_ACCESSION = re.compile(r'^([A-Z]{2}_\d+)\.(\d+)$')
_CDS_RANGE = re.compile(r'<?(\d+)\.\.>?(\d+)')
_CDS_LOCATION = re.compile(rf'{_CDS_RANGE.pattern}|join\({_CDS_RANGE.pattern}(?:,{_CDS_RANGE.pattern})+\)')


@dataclasses.dataclass(frozen=True)
class RecordCds:
    """A GenBank record's CDS, 0-based and inclusive, and whether each end runs off the record (`<`, `>`)."""

    start: int
    end: int
    start_open: bool
    end_open: bool


class CdsStatementError(build.BuildError):
    """A record states its CDS in a form the reader does not interpret; the historical set leaves such a record out."""


@dataclasses.dataclass
class Record:
    """What a transcript's GenBank record states: its CDS, and when read for them, its sequence and protein."""

    cds: RecordCds | None = None
    protein_id: str | None = None  # the CDS's protein_id
    sequence: bytes | None = None
    protein: tuple[str, bytes] | None = None  # the CDS's protein_id and translation, when read for sequences


def parse_cds(version: str, location: str) -> RecordCds:
    """One CDS location, `a..b` or a forward `join` of such ranges, either end possibly open.

    Raises:
        CdsStatementError: If the location is any other shape — `complement`, `order`, a remote
            accession, a single base — or its ranges do not run forward.
    """
    if _CDS_LOCATION.fullmatch(location) is None:
        raise CdsStatementError(f'{version}: CDS location {location!r} is not a range or a forward join of ranges')
    ranges = [(int(a), int(b)) for a, b in _CDS_RANGE.findall(location)]
    bounds = [n for r in ranges for n in r]
    if bounds != sorted(bounds):
        raise CdsStatementError(f'{version}: CDS location {location!r} does not run forward')
    if '<' in location.partition(',')[2] or '>' in location.rpartition(',')[0]:
        raise CdsStatementError(f'{version}: CDS location {location!r} is open at an inner bound')
    return RecordCds(ranges[0][0] - 1, ranges[-1][1] - 1, start_open='<' in location, end_open='>' in location)


def _quoted(version: str, qualifier: str, text: str) -> str:
    """The value of a `/name="..."` qualifier, its continuation lines already joined.

    Raises:
        build.BuildError: If the value is not quoted.
    """
    value = text.partition('=')[2]
    if len(value) < 2 or value[0] != '"' or value[-1] != '"':
        raise build.BuildError(f'{version}: CDS qualifier {qualifier} is not a quoted value: {text[:60]!r}')
    return value[1:-1]


class _CdsQualifiers:
    """The qualifiers of one record's CDS feature, read a line at a time.

    Without `protein`, the translation is not kept: a release build takes its proteins from the protein
    set, and buffering every translation would hold them all for nothing.
    """

    def __init__(self, version: str, cds: RecordCds, *, protein: bool) -> None:
        self.version = version
        self.cds = cds
        self._protein = protein
        self._lines: list[str] = []
        self.protein_id: str | None = None
        self.translation: str | None = None

    def add(self, line: str) -> None:
        text = line.strip()
        if text.startswith('/'):
            self._close()
            self._lines = [] if text.startswith('/translation=') and not self._protein else [text]
        elif self._lines:
            self._lines.append(text)

    def _close(self) -> None:
        if not self._lines:
            return
        text = ''.join(self._lines)
        name = text.partition('=')[0]
        if name == '/codon_start' and text != '/codon_start=1' and not self.cds.start_open:
            raise CdsStatementError(f'{self.version}: a CDS with a closed start but {text}')
        if name == '/protein_id':
            protein_id = _quoted(self.version, name, text)
            if _ACCESSION.match(protein_id) is None:
                raise build.BuildError(f'{self.version}: protein_id {protein_id!r} is not a versioned accession')
            self.protein_id = protein_id
        elif name == '/translation':
            self.translation = _quoted(self.version, name, text)
        self._lines = []

    def finish(self) -> None:
        """Close the last qualifier; called once no further qualifier line can follow."""
        self._close()

    def protein(self) -> tuple[str, bytes] | None:
        """The protein id and translation, once finished; None unless both were stated."""
        if self.protein_id is None or self.translation is None:
            return None
        return self.protein_id, self.translation.encode('ascii')


def read_records(path: pathlib.Path, *, sequences: bool, left_out: dict[str, str] | None = None) -> dict[str, Record]:
    """Each GenBank record's CDS in its own coordinates and, when `sequences`, its sequence and protein.

    With `left_out`, a record whose CDS statement the reader refuses — a location of another shape, a
    closed start read from other than its first base, an open end that is not the record's — is left
    out of the result and noted there by version with the reason, instead of failing the read.

    A CDS written as a join — a programmed frameshift, `join(114..317,319..801)` — is taken by its outer
    bounds, the first and last coding bases, which is what c. numbering needs. The sequence is the
    ORIGIN section upper-cased; the protein is the CDS's `/protein_id` with its `/translation`.

    A record with several CDS features — an old record annotating a pseudogene's reading frame beside
    the protein's — has as its CDS the one that names a protein; several naming one is refused.

    Raises:
        build.BuildError: If a record states several CDS features naming a protein, or several and
            none naming one, a version appears twice, a CDS or an ORIGIN comes before any VERSION
            line, a VERSION line has no accession, a location is of a shape `parse_cds`
            refuses, a record is not closed by `//` before the next or before the file ends, a CDS
            with a closed start reads from other than its first base (`/codon_start` 2 or 3), which
            no start codon can do, or, when read for sequences, a record has no ORIGIN section or an
            empty one, or a CDS no protein id or translation.
    """
    out: dict[str, Record] = {}
    version: str | None = None
    features: list[_CdsQualifiers] = []  # the record's CDS features, each with its qualifiers
    qualifiers: _CdsQualifiers | None = None
    origin: list[str] | None = None

    def refuse(error: CdsStatementError) -> None:
        if left_out is None or version is None:
            raise error
        left_out[version] = str(error)

    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if qualifiers is not None and not line.startswith(' ' * 21):
                qualifiers.finish()
                qualifiers = None
            if qualifiers is not None:
                try:
                    qualifiers.add(line)
                except CdsStatementError as error:
                    refuse(error)
                    qualifiers = None
            elif line.startswith('LOCUS'):
                if version is not None:
                    raise build.BuildError(f'{path}: {version} is not closed by // before the next record')
                features = []
            elif origin is not None and not line.startswith('//'):
                origin.append(''.join(line.split()[1:]))
            elif line.startswith('VERSION'):
                fields = line.split()
                if len(fields) < 2:
                    raise build.BuildError(f'{path}: a VERSION line with no accession')
                version = fields[1]
                if version in out:
                    raise build.BuildError(f'{path}: {version} appears twice')
            elif line.startswith('     CDS '):
                if version is None:
                    raise build.BuildError(f'{path}: a CDS in a record with no VERSION line: {line.strip()[:60]!r}')
                location = line[21:].strip()
                while location.count('(') > location.count(')'):  # a join may wrap onto the next lines
                    continued = next(fh, None)
                    if continued is None:
                        raise build.BuildError(f'{path}: ends inside the CDS location of {version}')
                    location += continued.strip()
                try:
                    qualifiers = _CdsQualifiers(version, parse_cds(version, location), protein=sequences)
                except CdsStatementError as error:
                    refuse(error)
                else:
                    features.append(qualifiers)
            elif line.startswith('ORIGIN') and sequences:
                if version is None:
                    raise build.BuildError(f'{path}: an ORIGIN in a record with no VERSION line')
                origin = []
            elif line.startswith('//') and version is not None:
                if qualifiers is not None:
                    qualifiers.finish()
                    qualifiers = None
                if left_out is None or version not in left_out:
                    try:
                        out[version] = _record(path, version, features, origin, sequences=sequences)
                    except CdsStatementError as error:
                        refuse(error)
                version, features, origin = None, [], None
    if version is not None:
        raise build.BuildError(f'{path}: ends inside the record of {version}, before its //')
    return out


def _record(
    path: pathlib.Path, version: str, features: list[_CdsQualifiers], origin: list[str] | None, *, sequences: bool
) -> Record:
    """One record from its CDS features and ORIGIN, as `read_records` documents.

    Raises:
        CdsStatementError: If, read for sequences, the CDS lies outside the sequence or is open at an
            end that is not the sequence's.
        build.BuildError: If the CDS features do not single one out, or, when read for sequences, the
            record has no ORIGIN section or an empty one, or its CDS no protein id or translation.
    """
    named = [f for f in features if f.protein_id is not None]
    if len(named) > 1 or (not named and len(features) > 1):
        raise build.BuildError(
            f'{path}: {version} states {len(features)} CDS features, {len(named)} of them naming a protein'
        )
    chosen = named[0] if named else features[0] if features else None
    record = Record(
        cds=chosen.cds if chosen is not None else None,
        protein_id=chosen.protein_id if chosen is not None else None,
        protein=chosen.protein() if chosen is not None else None,
    )
    if not sequences:
        return record
    residues = ''.join(origin or ()).upper().encode('ascii')
    if not residues:
        raise build.BuildError(f'{path}: {version} has no ORIGIN section, or an empty one, to read its sequence from')
    record.sequence = residues
    if record.cds is not None:
        if record.protein is None:
            raise build.BuildError(f'{path}: the CDS of {version} names no protein_id or carries no translation')
        cds = record.cds
        message = bundle_pb2.Cds(
            start_index=cds.start, end_index_inclusive=cds.end, start_open=cds.start_open, end_open=cds.end_open
        )
        try:
            common.check_cds(version, message, len(residues))
        except build.BuildError as error:
            raise CdsStatementError(str(error)) from error
    return record
