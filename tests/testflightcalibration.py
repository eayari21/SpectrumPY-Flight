"""Manual driver for the flight calibration analysis.

The module lives inside ``tests/`` because it exercises the full
end-to-end pipeline used by :mod:`pytest`.  It is intentionally marked as a
skipped test so that the heavy-weight decoding logic does not run during the
regular unit-test suite.  Invoke it directly from the repository root:

.. code-block:: console

    $ python tests/testflightcalibration.py --output-dir reports

The script ensures that every ``ois_output_*`` capture under ``Data/`` has an
associated ``HDF5/`` analysis file before generating
``flight_calibration_report.pdf``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable, Sequence

# Ensure ``spectrumpy_flight`` can be imported when executing the script directly
# via ``python tests/testflightcalibration.py``.  Python sets ``sys.path[0]`` to
# the directory containing the script (``tests/``), so the repository root is
# not automatically available on ``sys.path``.  Inject it manually before
# importing project modules so the helper works without installing the package.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

try:  # pragma: no cover - makes pytest ignore the helper module
    import pytest
except Exception:  # pragma: no cover - pytest is not required to run the script
    pytest = None  # type: ignore[assignment]
else:  # pragma: no cover
    pytestmark = pytest.mark.skip("Manual helper script — invoke via __main__.")

from spectrumpy_flight.flight_calibration import generate_flight_calibration_report
from spectrumpy_flight.idex_packet import IDEXEvent, _resolve_output_path


def _iter_data_products(data_root: Path) -> Iterable[Path]:
    """Yield every capture (decoded or raw) inside ``data_root``."""

    return sorted(data_root.glob("ois_output_*"))


def ensure_hdf_products(data_root: Path) -> list[Path]:
    """Decode every raw capture under ``data_root`` into ``HDF5/`` outputs."""

    generated: list[Path] = []
    for path in _iter_data_products(data_root):
        if path.suffix.lower() == ".h5" and path.exists():
            generated.append(path)
            continue

        target = _resolve_output_path(str(path))
        if target.exists():
            generated.append(target)
            continue

        event = IDEXEvent(str(path))
        target.parent.mkdir(parents=True, exist_ok=True)
        event.write_to_hdf5(event.data, str(target))
        if not target.exists():
            raise RuntimeError(f"Failed to write decoded analysis file: {target}")
        generated.append(target)

    return generated


def main(argv: Sequence[str] | None = None) -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_data = repo_root / "Data"
    default_output = repo_root / "reports"

    parser = argparse.ArgumentParser(description="Rebuild the flight calibration report using repository data.")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=default_data,
        help="Directory containing raw Data/ois_output captures.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help="Directory where the report PDF and JSON summary will be written.",
    )
    parser.add_argument(
        "--report-name",
        default="flight_calibration_report.pdf",
        help="Optional custom filename for the generated PDF report.",
    )
    args = parser.parse_args(argv)

    data_root = args.data_root.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not data_root.exists():
        raise SystemExit(f"Data directory not found: {data_root}")

    decoded = ensure_hdf_products(data_root)
    print(f"Decoded or verified {len(decoded)} analysis products under {data_root}.")

    results = generate_flight_calibration_report(
        data_root,
        output_dir,
        report_name=args.report_name,
    )
    print(f"Report saved to {results['pdf']}")
    print(f"Summary saved to {results['summary']}")
    return 0


if __name__ == "__main__":  # pragma: no cover - manual invocation helper
    raise SystemExit(main())
