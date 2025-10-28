"""Test configuration for SpectrumPY-Flight.

This file ensures that the repository's root directory is available on
``sys.path`` before any test modules are imported.  In developer
environments where the project has not been installed as a package,
pytest may otherwise fail to resolve imports such as ``idex_packet`` or
``olivine_metrics``.  By explicitly appending the root directory to the
module search path we keep the tests importable without requiring an
editable install first.
"""

from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
