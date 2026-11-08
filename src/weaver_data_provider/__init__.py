"""Reference data for hgvs-weaver: gene bundles built from a publisher's release, and a DataProvider over them."""

from __future__ import annotations

import importlib.metadata

__version__ = importlib.metadata.version('hgvs-weaver-data')

# The at-rest formats this version reads, oldest to newest; see docs/design/store.md.
MIN_FORMAT_VERSION = 1
FORMAT_VERSION = 2
