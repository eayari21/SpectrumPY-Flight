"""Allow ``python -m imfpy.gui`` to behave like ``python -m imfpy.gui.main``."""

from __future__ import annotations

import sys

from .main import main

if __name__ == "__main__":
    sys.exit(main())
