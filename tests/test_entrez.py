"""The status fetch: Entrez's answers as the table the historical build reads, through a fake esummary."""

from __future__ import annotations

import gzip
import pathlib

import pytest

from weaver_data_provider.build import cli, entrez

# What esummary answers for a version: a live record has no status; a retired one names its status.
ANSWERS: dict[str, dict[str, str]] = {
    'NM_000059.3': {'accessionversion': 'NM_000059.3', 'status': 'replaced', 'replacedby': 'NM_000059.4'},
    'NM_000010.1': {'accessionversion': 'NM_000010.1'},
    'NM_000010.2': {'accessionversion': 'NM_000010.2', 'status': 'replaced', 'replacedby': 'NM_000010.3'},
    'NR_000030.1': {'accessionversion': 'NR_000030.1', 'status': 'suppressed'},
}


def _fake(
    answers: dict[str, dict[str, str]], *, drop_once: frozenset[str] = frozenset()
) -> tuple[entrez.Post, list[list[str]]]:
    """An esummary that answers from `answers`, dropping the versions in `drop_once` from their first batch."""
    requests: list[list[str]] = []
    dropped: set[str] = set()

    def post(ids: list[str], api_key: str | None, pace: entrez._Pace) -> dict[str, dict[str, str]]:
        requests.append(list(ids))
        out = {}
        for uid, version in enumerate(ids):
            if version in drop_once and version not in dropped:
                dropped.add(version)
                continue
            if version in answers:
                out[str(uid)] = answers[version]
        return out

    return post, requests


def _records(path: pathlib.Path, versions: list[str]) -> pathlib.Path:
    with gzip.open(path, 'wt', encoding='ascii') as fh:
        for version in versions:
            fh.write(f'LOCUS       {version.split(".")[0]}\nVERSION     {version}\n//\n')
    return path


def test_the_table_holds_each_versions_status_with_live_read_from_the_absence_of_one(tmp_path: pathlib.Path) -> None:
    post, _ = _fake(ANSWERS)
    unanswered = entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1', 'NM_000010.2', 'NR_000030.1'], None, post=post)
    assert unanswered == []
    assert entrez.read_table(tmp_path / 's.tsv') == {
        'NM_000010.1': ('live', ''),
        'NM_000010.2': ('replaced', 'NM_000010.3'),
        'NR_000030.1': ('suppressed', ''),
    }


def test_a_rerun_fetches_only_the_versions_the_table_lacks(tmp_path: pathlib.Path) -> None:
    post, requests = _fake(ANSWERS)
    entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1'], None, post=post)
    entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1', 'NM_000010.2'], None, post=post)
    asked = [ids for ids in requests if ids != [entrez.SENTINEL[0]]]
    assert asked == [['NM_000010.1'], ['NM_000010.2']]
    assert set(entrez.read_table(tmp_path / 's.tsv')) == {'NM_000010.1', 'NM_000010.2'}


def test_a_version_dropped_from_a_batch_is_asked_for_again_and_one_never_answered_is_reported(
    tmp_path: pathlib.Path,
) -> None:
    post, _ = _fake(ANSWERS, drop_once=frozenset({'NM_000010.2'}))
    unanswered = entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1', 'NM_000010.2', 'NM_999999.9'], None, post=post)
    assert unanswered == ['NM_999999.9']
    assert set(entrez.read_table(tmp_path / 's.tsv')) == {'NM_000010.1', 'NM_000010.2'}


def test_a_sentinel_entrez_no_longer_calls_replaced_fails_the_fetch_before_any_row(tmp_path: pathlib.Path) -> None:
    post, _ = _fake({**ANSWERS, 'NM_000059.3': {'accessionversion': 'NM_000059.3'}})
    with pytest.raises(RuntimeError, match='no longer states status'):
        entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1'], None, post=post)
    assert not (tmp_path / 's.tsv').exists()


def test_a_summary_naming_a_successor_but_no_status_is_refused(tmp_path: pathlib.Path) -> None:
    post, _ = _fake({**ANSWERS, 'NM_000010.1': {'accessionversion': 'NM_000010.1', 'replacedby': 'NM_000010.2'}})
    with pytest.raises(RuntimeError, match=r'names a successor NM_000010\.2 but no status'):
        entrez.fetch_into(tmp_path / 's.tsv', ['NM_000010.1'], None, post=post)


def test_a_table_with_another_header_is_refused(tmp_path: pathlib.Path) -> None:
    (tmp_path / 's.tsv').write_text('version\tstate\n', 'ascii')
    with pytest.raises(ValueError, match='header'):
        entrez.read_table(tmp_path / 's.tsv')


def test_the_status_command_writes_the_table_for_the_records_versions(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    post, _ = _fake(ANSWERS)
    monkeypatch.setattr(entrez, '_post', post)
    records = _records(tmp_path / 'rna.gbff.gz', ['NM_000010.1', 'NR_000030.1'])
    cli.main(['status', '--records', str(records), '--out', str(tmp_path / 'status.tsv')])
    expected = {'NM_000010.1': ('live', ''), 'NR_000030.1': ('suppressed', '')}
    assert entrez.read_table(tmp_path / 'status.tsv') == expected
    assert '2 versions, live 1, suppressed 1 of 2 named' in capsys.readouterr().err


def test_the_status_command_fails_when_a_version_goes_unanswered(
    tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    post, _ = _fake(ANSWERS)
    monkeypatch.setattr(entrez, '_post', post)
    records = _records(tmp_path / 'rna.gbff.gz', ['NM_000010.1', 'NM_999999.9'])
    with pytest.raises(SystemExit):
        cli.main(['status', '--records', str(records), '--out', str(tmp_path / 'status.tsv')])
    assert set(entrez.read_table(tmp_path / 'status.tsv')) == {'NM_000010.1'}
