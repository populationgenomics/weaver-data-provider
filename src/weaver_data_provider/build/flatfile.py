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

# A versioned INSDC or RefSeq accession: `AAA73985.1`, `NP_009225.1`.
_ACCESSION = re.compile(r'^[A-Z]+_?\d+\.\d+$')
# Qualifiers saying a CDS's translation is not the conceptual translation of its bases.
_EXCEPTIONS = frozenset(
    {'/exception', '/transl_except', '/ribosomal_slippage', '/trans_splicing', '/artificial_location'}
)
# A LOCUS line's length and molecule type, wherever a long locus name has pushed them: `5711 bp    mRNA`.
_LOCUS_MOLECULE = re.compile(r'\s\d+ bp\s+(\S+)')
_CDS_RANGE = re.compile(r'<?(\d+)\.\.>?(\d+)')
_CDS_LOCATION = re.compile(rf'{_CDS_RANGE.pattern}|join\({_CDS_RANGE.pattern}(?:,{_CDS_RANGE.pattern})+\)')


@dataclasses.dataclass(frozen=True)
class RecordCds:
    """A record's CDS as it states it: its bounds, and what its qualifiers say about reading its bases.

    The bounds are 0-based and inclusive; an open end (`<`, `>`) is one the coding sequence runs past.
    """

    start: int
    end: int
    start_open: bool
    end_open: bool
    joined: bool  # a join of ranges, such as a programmed frameshift, taken by its outer bounds
    codon_start: int  # which of the CDS's first three bases begins its first whole codon, from 1
    translation_table: int  # NCBI's genetic code number; 1, the standard code, where the record states none
    exceptions: tuple[str, ...]  # qualifiers saying its translation is not its bases' conceptual one


@dataclasses.dataclass(frozen=True)
class _Location:
    start: int
    end: int
    start_open: bool
    end_open: bool
    joined: bool


class UninterpretableError(build.BuildError):
    """A record the reader takes no CDS from; a reader of an archive may leave it out rather than fail."""


class CdsStatementError(UninterpretableError):
    """A record states its CDS in a form the reader does not interpret; the historical set leaves such a record out."""


class CdsChoiceError(UninterpretableError):
    """A record states no single CDS with its protein: several naming one, several naming none, or no translation."""


@dataclasses.dataclass(frozen=True)
class Record:
    """What a GenBank record states: its molecule, its CDS, and when read for them, its sequence and protein."""

    molecule: str  # the LOCUS line's molecule type: `mRNA`, `DNA`, `RNA`, ...; empty where it states none
    cds: RecordCds | None
    protein_id: str | None  # the CDS's protein_id
    sequence: bytes | None  # None unless read for sequences
    protein: tuple[str, bytes] | None  # the CDS's protein_id and translation, when read for sequences


def _parse_location(version: str, location: str) -> _Location:
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
    return _Location(
        ranges[0][0] - 1,
        ranges[-1][1] - 1,
        start_open='<' in location,
        end_open='>' in location,
        joined=location.startswith('join('),
    )


def _quoted(version: str, qualifier: str, text: str) -> str:
    """The value of a `/name="..."` qualifier, its continuation lines already joined.

    Raises:
        build.BuildError: If the value is not quoted.
    """
    value = text.partition('=')[2]
    if len(value) < 2 or value[0] != '"' or value[-1] != '"':
        raise build.BuildError(f'{version}: CDS qualifier {qualifier} is not a quoted value: {text[:60]!r}')
    return value[1:-1]


def _number(version: str, text: str) -> int:
    """The value of a `/name=N` qualifier.

    Raises:
        build.BuildError: If the value is not a whole number.
    """
    value = text.partition('=')[2]
    if not value.isdigit():
        raise build.BuildError(f'{version}: CDS qualifier {text[:60]!r} is not a whole number')
    return int(value)


class _CdsQualifiers:
    """The qualifiers of one record's CDS feature, read a line at a time.

    Without `protein`, the translation is not kept: a release build takes its proteins from the protein
    set, and buffering every translation would hold them all for nothing. A codon start or genetic code
    the feature does not state is INSDC's default, 1.
    """

    def __init__(self, version: str, location: _Location, *, protein: bool) -> None:
        self.version = version
        self.location = location
        self._protein = protein
        self._lines: list[str] = []
        self.protein_id: str | None = None
        self.translation: str | None = None
        self.codon_start = 1
        self.translation_table = 1
        self.exceptions: list[str] = []

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
        if name == '/codon_start' and text != '/codon_start=1' and not self.location.start_open:
            raise CdsStatementError(f'{self.version}: a CDS with a closed start but {text}')
        if name in _EXCEPTIONS:
            self.exceptions.append(name)
        if name == '/codon_start':
            self.codon_start = _number(self.version, text)
        elif name == '/transl_table':
            self.translation_table = _number(self.version, text)
        elif name == '/protein_id':
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

    def cds(self) -> RecordCds:
        """The CDS, once finished."""
        loc = self.location
        return RecordCds(
            loc.start,
            loc.end,
            start_open=loc.start_open,
            end_open=loc.end_open,
            joined=loc.joined,
            codon_start=self.codon_start,
            translation_table=self.translation_table,
            exceptions=tuple(self.exceptions),
        )

    def protein(self) -> tuple[str, bytes] | None:
        """The protein id and translation, once finished; None unless both were stated."""
        if self.protein_id is None or self.translation is None:
            return None
        return self.protein_id, self.translation.encode('ascii')


def read_records(
    path: pathlib.Path,
    *,
    sequences: bool,
    left_out: dict[str, UninterpretableError] | None = None,
    leave_out: tuple[type[UninterpretableError], ...] = (CdsStatementError,),
) -> dict[str, Record]:
    """Each GenBank record's CDS in its own coordinates and, when `sequences`, its sequence and protein.

    With `left_out`, a record refused with one of the `leave_out` errors is left out of the result and
    noted there by version with the error, instead of failing the read. A `CdsStatementError` is a CDS
    statement the reader refuses — a location of another shape, a closed start read from other than its
    first base, an open end that is not the record's; a `CdsChoiceError`, a record with no single CDS
    naming a protein and carrying its translation.

    A CDS written as a join — a programmed frameshift, `join(114..317,319..801)` — is taken by its outer
    bounds, the first and last coding bases, which is what c. numbering needs. The sequence is the
    ORIGIN section upper-cased; the protein is the CDS's `/protein_id` with its `/translation`.

    A record with several CDS features — an old record annotating a pseudogene's reading frame beside
    the protein's — has as its CDS the one that names a protein; several naming one is refused.

    Raises:
        build.BuildError: If a record states several CDS features naming a protein, or several and
            none naming one, a version appears twice, a CDS or an ORIGIN comes before any VERSION
            line, a VERSION line has no accession, a location is of a shape `_parse_location`
            refuses, a record is not closed by `//` before the next or before the file ends, a CDS
            with a closed start reads from other than its first base (`/codon_start` 2 or 3), which
            no start codon can do, a codon start or genetic code is not a number, or, when read for
            sequences, a record has no ORIGIN section or an empty one, or a CDS no protein id or
            translation.
    """
    out: dict[str, Record] = {}
    version: str | None = None
    molecule = ''
    features: list[_CdsQualifiers] = []  # the record's CDS features, each with its qualifiers
    qualifiers: _CdsQualifiers | None = None
    origin: list[str] | None = None

    def refuse(error: UninterpretableError) -> None:
        if left_out is None or version is None or not isinstance(error, leave_out):
            raise error
        left_out[version] = error

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
                found = _LOCUS_MOLECULE.search(line)
                molecule = found.group(1) if found is not None else ''
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
                    qualifiers = _CdsQualifiers(version, _parse_location(version, location), protein=sequences)
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
                        out[version] = _record(path, version, molecule, features, origin, sequences=sequences)
                    except UninterpretableError as error:
                        refuse(error)
                version, features, origin = None, [], None
    if version is not None:
        raise build.BuildError(f'{path}: ends inside the record of {version}, before its //')
    return out


def _record(
    path: pathlib.Path,
    version: str,
    molecule: str,
    features: list[_CdsQualifiers],
    origin: list[str] | None,
    *,
    sequences: bool,
) -> Record:
    """One record from its CDS features and ORIGIN, as `read_records` documents.

    Raises:
        CdsStatementError: If, read for sequences, the CDS lies outside the sequence or is open at an
            end that is not the sequence's.
        CdsChoiceError: If the CDS features do not single one out, or, when read for sequences, the CDS
            names no protein id or carries no translation.
        build.BuildError: If, read for sequences, the record has no ORIGIN section or an empty one.
    """
    named = [f for f in features if f.protein_id is not None]
    if len(named) > 1 or (not named and len(features) > 1):
        raise CdsChoiceError(
            f'{path}: {version} states {len(features)} CDS features, {len(named)} of them naming a protein'
        )
    chosen = named[0] if named else features[0] if features else None
    cds = chosen.cds() if chosen is not None else None
    protein_id = chosen.protein_id if chosen is not None else None
    protein = chosen.protein() if chosen is not None else None
    if not sequences:
        return Record(molecule=molecule, cds=cds, protein_id=protein_id, sequence=None, protein=protein)
    residues = ''.join(origin or ()).upper().encode('ascii')
    if not residues:
        raise build.BuildError(f'{path}: {version} has no ORIGIN section, or an empty one, to read its sequence from')
    if cds is not None:
        if protein is None:
            raise CdsChoiceError(f'{path}: the CDS of {version} names no protein_id or carries no translation')
        message = bundle_pb2.Cds(
            start_index=cds.start, end_index_inclusive=cds.end, start_open=cds.start_open, end_open=cds.end_open
        )
        try:
            common.check_cds(version, message, len(residues))
        except build.BuildError as error:
            raise CdsStatementError(str(error)) from error
    return Record(molecule=molecule, cds=cds, protein_id=protein_id, sequence=residues, protein=protein)
