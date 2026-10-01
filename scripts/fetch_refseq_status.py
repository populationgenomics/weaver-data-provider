"""Fetch NCBI's current status of each RefSeq transcript version from Entrez, as a table for the builder.

A GenBank record says when it was replaced, but not when it was suppressed, and a record written
before its retirement says nothing at all; what NCBI holds a version to be today is in Entrez. This
asks E-utilities' `esummary` for every version named in a GenBank file and writes one row per
version: `accession_version`, `status` (`live`, `replaced` or `suppressed`) and `replaced_by` (the
successor NCBI names, possibly unversioned or another accession; empty otherwise). The builder's
`historical` command takes the table as `--status`.

    uv run python scripts/fetch_refseq_status.py --records knownrefseq_rna.gbff.gz --out status.tsv

Entrez takes about ten seconds to answer for 500 versions, so requests run three at a time, which
stays well inside NCBI's limit of three requests a second without an API key (ten with one, read from
`NCBI_API_KEY`). The fetch of 190,000 versions takes about twenty minutes.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import json
import os
import pathlib
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ESUMMARY = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
BATCH = 500
ATTEMPTS = 6
WORKERS = 3


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


def _fetch(body: bytes) -> dict[str, object]:
    """One esummary request, retried on a rate-limit or server error."""
    for attempt in range(1, ATTEMPTS + 1):
        try:
            request = urllib.request.Request(ESUMMARY, data=body, method='POST')  # noqa: S310 — a fixed https URL
            with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310
                return json.load(response)
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt == ATTEMPTS:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == ATTEMPTS:
                raise
        time.sleep(2**attempt)
    raise AssertionError('unreachable: every attempt returns or raises')


def _post(ids: list[str], api_key: str | None) -> dict[str, dict[str, str]]:
    """The per-record summaries of one batch, by uid."""
    form = {'db': 'nuccore', 'id': ','.join(ids), 'retmode': 'json'}
    if api_key:
        form['api_key'] = api_key
    payload = _fetch(urllib.parse.urlencode(form).encode('ascii'))
    result = payload.get('result')
    if not isinstance(result, dict) or 'uids' not in result:
        raise RuntimeError(f'esummary answered without a result: {str(payload)[:200]}')
    return {uid: result[uid] for uid in result['uids']}


def _batch(ids: list[str], api_key: str | None) -> dict[str, tuple[str, str]]:
    """One batch's `(status, replaced_by)` by version; a version Entrez does not answer for is an error."""
    out: dict[str, tuple[str, str]] = {}
    for summary in _post(ids, api_key).values():
        if 'error' in summary:
            continue
        out[summary['accessionversion']] = (summary.get('status') or 'live', summary.get('replacedby') or '')
    missing = [v for v in ids if v not in out]
    if missing:
        raise RuntimeError(f'esummary answered nothing for {len(missing)} versions, e.g. {missing[:5]}')
    return out


def statuses(versions: list[str], api_key: str | None, *, log: bool = True) -> dict[str, tuple[str, str]]:
    """Each version's `(status, replaced_by)`, fetched in batches a few at a time."""
    out: dict[str, tuple[str, str]] = {}
    batches = [versions[k : k + BATCH] for k in range(0, len(versions), BATCH)]
    interval = 0.1 if api_key else 0.34  # the least time between two requests NCBI allows

    def fetch(index: int) -> dict[str, tuple[str, str]]:
        time.sleep(index % WORKERS * interval)  # stagger the first requests; later ones are paced by the answers
        return _batch(batches[index], api_key)

    with concurrent.futures.ThreadPoolExecutor(WORKERS) as pool:
        for done, found in enumerate(pool.map(fetch, range(len(batches))), start=1):
            out.update(found)
            if log:
                print(f'{min(done * BATCH, len(versions))}/{len(versions)}', file=sys.stderr, end='\r')
    return out


def write(path: pathlib.Path, found: dict[str, tuple[str, str]]) -> None:
    with path.open('w', encoding='ascii') as fh:
        fh.write('accession_version\tstatus\treplaced_by\n')
        for version in sorted(found):
            status, replaced_by = found[version]
            fh.write(f'{version}\t{status}\t{replaced_by}\n')


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--records', type=pathlib.Path, required=True, help='a gzipped GenBank file naming the versions'
    )
    parser.add_argument('--out', type=pathlib.Path, required=True, help='the TSV to write')
    args = parser.parse_args(argv)
    versions = versions_in(args.records)
    found = statuses(versions, os.environ.get('NCBI_API_KEY'))
    write(args.out, found)
    counts: dict[str, int] = {}
    for status, _ in found.values():
        counts[status] = counts.get(status, 0) + 1
    print(
        f'{args.out}: {len(found)} versions, ' + ', '.join(f'{s} {n}' for s, n in sorted(counts.items())),
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
