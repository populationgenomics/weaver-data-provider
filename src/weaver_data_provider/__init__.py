"""Reference data for hgvs-weaver: gene bundles built from a publisher's release, and a DataProvider over them."""

from __future__ import annotations

import importlib.metadata

__version__ = importlib.metadata.version('hgvs-weaver-data')

# The at-rest formats this version reads: any from the first to the second. A writer stamps its files with the lowest
# that reads them correctly, so an older reader refuses only what it would misread; the floor rises only in a major
# release (docs/design/store.md).
MIN_FORMAT_VERSION = 1
FORMAT_VERSION = 2


def format_versions() -> str:
    """The versions this package reads, for an error refusing another."""
    return f'this package reads versions {MIN_FORMAT_VERSION} to {FORMAT_VERSION}'
