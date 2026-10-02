"""Reference data for hgvs-weaver: gene bundles built from a publisher's release, and a DataProvider over them."""

from __future__ import annotations

import importlib.metadata

__version__ = importlib.metadata.version('hgvs-weaver-data')

# The at-rest format this version writes and reads. A reader refuses any other, rather than misreading it.
FORMAT_VERSION = 1
