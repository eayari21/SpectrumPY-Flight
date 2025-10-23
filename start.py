"""Convenience entry point for launching the SpectrumPY UI."""
import os
os.environ.pop("QT_DEBUG_PLUGINS", None)
from spectrum_launcher import main


if __name__ == "__main__":
    raise SystemExit(main())
