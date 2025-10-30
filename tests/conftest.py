"""Test configuration for SpectrumPY-Flight.

The tests expect the package to be importable either from an installed
location or from the repository's ``src`` tree.  When running the suite
locally without installing the package first we augment ``sys.path`` to
point at the ``src`` directory so Python can resolve
``spectrumpy_flight`` just as it would after ``pip install``.
"""

from __future__ import annotations

import sys
from pathlib import Path


SRC_ROOT = Path(__file__).resolve().parent.parent / "src"

if SRC_ROOT.exists():
    sys.path.insert(0, str(SRC_ROOT))
