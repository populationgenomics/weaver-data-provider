"""NCBI's current status of each RefSeq transcript version, fetched from Entrez as a table for the builder.

A GenBank record says when it was replaced, but not when it was suppressed, and a record written
before its retirement says nothing at all; what NCBI holds a version to be today is in Entrez. This
asks E-utilities' `esummary` for every version named in a GenBank file and writes one row per
version: `accession_version`, `status` (`live`, `replaced` or `suppressed`) and `replaced_by` (the
successor NCBI names, possibly unversioned or another accession; empty otherwise). The builder's
`historical` command takes the table as `--status`. This is the one builder step that reads from the
network rather than a local file, since no file NCBI publishes states a version's status; it is a
command of its own, so a build proper still reads only what it was given.

Entrez's summary carries a `status` field only for a retired record; a live one has none, and `live`
here is read from that absence. A sentinel version known to be replaced is checked first, so a change
in how Entrez states status fails the run rather than writing every version as live.

Entrez takes about ten seconds to answer for 500 versions, so requests run three at a time, which
stays well inside NCBI's limit of three requests a second without an API key (ten with one, read from
`NCBI_API_KEY`). The fetch of 190,000 versions takes about twenty minutes. Rows are appended to the
table as each batch is answered, and a rerun with the same table fetches only the versions it lacks,
so an interrupted fetch resumes. A version Entrez answers nothing for is reported at the end and the
run fails, with every other row kept; such a version has no status and the build refuses the set
until it is resolved by hand.
"""

from __future__ import annotations

import concurrent.futures
import gzip
import json
import pathlib
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable

BATCH = 500
ATTEMPTS = 6
WORKERS = 3
HEADER = ['accession_version', 'status', 'replaced_by']
SENTINEL = ('NM_000059.3', 'replaced')  # BRCA2's previous version, replaced by NM_000059.4 in 2019
# One esummary request: the versions asked for, the API key, the pacer; the summaries by uid.
Post = Callable[[list[str], 'str | None', '_Pace'], dict[str, dict[str, str]]]


class _Pace:
    """The least interval between two requests, kept across threads."""

    def __init__(self, interval: float) -> None:
        self._interval = interval
        self._next = 0.0
        self._lock = threading.Lock()

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            start = max(now, self._next)
            self._next = start + self._interval
        time.sleep(start - now)


def versions_in(records: pathlib.Path) -> list[str]:
    """Every VERSION line's accession.version in a gzipped GenBank file, in file order, each once."""
    seen: dict[str, None] = {}
    with gzip.open(records, 'rt', encoding='ascii', errors='replace') as fh:
        for line in fh:
            if line.startswith('VERSION'):
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(f'{records}: a VERSION line with no accession')
                seen[fields[1]] = None
    return list(seen)


def _fetch(body: bytes, pace: _Pace) -> dict[str, object]:
    """One esummary request, paced, retried on a rate-limit or server error."""
    for attempt in range(1, ATTEMPTS + 1):
        pace.wait()
        try:
            request = urllib.request.Request(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi', data=body, method='POST'
            )
            with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310 — a fixed https URL
                return json.load(response)
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt == ATTEMPTS:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == ATTEMPTS:
                raise
        time.sleep(2**attempt)
    raise AssertionError('unreachable: every attempt returns or raises')


def _post(ids: list[str], api_key: str | None, pace: _Pace) -> dict[str, dict[str, str]]:
    """The per-record summaries of one batch, by uid."""
    form = {'db': 'nuccore', 'id': ','.join(ids), 'retmode': 'json'}
    if api_key:
        form['api_key'] = api_key
    payload = _fetch(urllib.parse.urlencode(form).encode('ascii'), pace)
    result = payload.get('result')
    if not isinstance(result, dict) or 'uids' not in result:
        raise RuntimeError(f'esummary answered without a result: {str(payload)[:200]}')
    return {uid: result[uid] for uid in result['uids']}


def _row(summary: dict[str, str]) -> tuple[str, tuple[str, str]]:
    """A summary as `(version, (status, replaced_by))`; live is the absence of a status."""
    version = summary.get('accessionversion')
    if not version:
        raise RuntimeError(f'an esummary record with no accessionversion: {str(summary)[:200]}')
    status, replaced_by = summary.get('status'), summary.get('replacedby') or ''
    if status is None and replaced_by:
        raise RuntimeError(f'{version}: Entrez names a successor {replaced_by} but no status')
    return version, (status or 'live', replaced_by)


def _batch(
    ids: list[str], api_key: str | None, pace: _Pace, post: Post
) -> tuple[dict[str, tuple[str, str]], list[str]]:
    """One batch's `(status, replaced_by)` by version, and the versions Entrez answered nothing for.

    A version missing from the batch's answer is asked for once more on its own: esummary now and
    then drops one from a large batch it answers alone.
    """
    out: dict[str, tuple[str, str]] = {}
    for summary in post(ids, api_key, pace).values():
        if 'error' not in summary:
            out.update([_row(summary)])
    missing = [v for v in ids if v not in out]
    if missing:
        for summary in post(missing, api_key, pace).values():
            if 'error' not in summary:
                out.update([_row(summary)])
        missing = [v for v in ids if v not in out]
    return out, missing


def check_sentinel(api_key: str | None, pace: _Pace, post: Post) -> None:
    """Entrez still states a retired record's status as this module reads it.

    Raises:
        RuntimeError: If the sentinel's status is not the one known.
    """
    found, _ = _batch([SENTINEL[0]], api_key, pace, post)
    if found.get(SENTINEL[0], ('', ''))[0] != SENTINEL[1]:
        raise RuntimeError(
            f'{SENTINEL[0]} is known to be {SENTINEL[1]}, but Entrez answered {found.get(SENTINEL[0])}; '
            'the summary no longer states status as this module reads it'
        )


def read_table(path: pathlib.Path) -> dict[str, tuple[str, str]]:
    """The rows a table already holds, by version; empty if there is no table."""
    if not path.exists():
        return {}
    out: dict[str, tuple[str, str]] = {}
    with path.open(encoding='ascii') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header != HEADER:
            raise ValueError(f'{path}: header {header} is not {HEADER}')
        for line in fh:
            version, status, replaced_by = line.rstrip('\n').split('\t')
            out[version] = (status, replaced_by)
    return out


def fetch_into(
    path: pathlib.Path, versions: list[str], api_key: str | None, *, log: bool = True, post: Post | None = None
) -> list[str]:
    """Append each version's row to the table as it is answered, skipping versions already there.

    Args:
        path: The table, created if absent, resumed if present.
        versions: Every version wanted in it.
        api_key: An NCBI API key, for the faster rate; None for the anonymous one.
        log: Whether to report progress on stderr.
        post: How one esummary request is made; Entrez itself unless a test passes its own.

    Returns:
        The versions Entrez answered nothing for.

    Raises:
        RuntimeError: If the sentinel's status is not the one known, or an answer is not of the shape read.
        ValueError: If the table's header is not this module's.
    """
    post = post or _post
    done = read_table(path)
    wanted = [v for v in versions if v not in done]
    batches = [wanted[k : k + BATCH] for k in range(0, len(wanted), BATCH)]
    pace = _Pace(0.1 if api_key else 0.34)  # the least time between two requests NCBI allows
    check_sentinel(api_key, pace, post)
    unanswered: list[str] = []

    def fetch(index: int) -> tuple[dict[str, tuple[str, str]], list[str]]:
        return _batch(batches[index], api_key, pace, post)

    with path.open('a', encoding='ascii') as fh:
        if not done:
            fh.write('\t'.join(HEADER) + '\n')
        with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
            for n, (found, missing) in enumerate(pool.map(fetch, range(len(batches))), start=1):
                for version in sorted(found):
                    status, replaced_by = found[version]
                    fh.write(f'{version}\t{status}\t{replaced_by}\n')
                fh.flush()
                unanswered.extend(missing)
                if log:
                    print(f'{len(done) + min(n * BATCH, len(wanted))}/{len(versions)}', file=sys.stderr, end='\r')
    return unanswered


def summary(path: pathlib.Path) -> str:
    """`189348 versions, live 86975, replaced 98226, suppressed 4147`, from the table as it stands."""
    counts: dict[str, int] = {}
    for status, _ in read_table(path).values():
        counts[status] = counts.get(status, 0) + 1
    return f'{sum(counts.values())} versions, ' + ', '.join(f'{status} {n}' for status, n in sorted(counts.items()))
