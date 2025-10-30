"""Command line entry point for SpectrumPY graphical interfaces."""

from __future__ import annotations

import argparse
import importlib.util
from importlib import resources
import sys
from types import ModuleType
from typing import Optional


def _load_quicklook_module() -> ModuleType:
    """Return the dynamically imported ``IDEX_quicklook`` module."""

    module_name = "spectrumpy_flight.IDEX_quicklook"
    resource = resources.files("spectrumpy_flight") / "IDEX-quicklook.py"
    if not resource.is_file():
        msg = "Unable to locate 'IDEX-quicklook.py' inside the installed package."
        raise FileNotFoundError(msg)

    with resources.as_file(resource) as module_path:
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            msg = "Failed to create a module spec for 'IDEX-quicklook.py'."
            raise ImportError(msg)

        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        sys.modules.setdefault("IDEX_quicklook", module)
        spec.loader.exec_module(module)
        return module


def _create_application() -> "QApplication":
    """Initialise a :class:`~PyQt6.QtWidgets.QApplication` with SpectrumPY defaults."""

    from PyQt6.QtGui import QFont
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    try:
        app.setStyle("Fusion")
    except Exception:
        # Some Qt builds may not expose the Fusion style; ignore the error.
        pass

    base_font = QFont(app.font())
    point_size = base_font.pointSize()
    if point_size <= 0:
        point_size = 12
    else:
        point_size = max(point_size, 12)
    base_font.setPointSize(point_size)
    app.setFont(base_font)
    return app


def _launch_quicklook(filename: Optional[str], eventnumber: Optional[int]) -> int:
    """Launch the IDEX Quicklook window."""

    module = _load_quicklook_module()
    MainWindow = getattr(module, "MainWindow")

    app = _create_application()
    window = MainWindow(filename=filename, eventnumber=eventnumber)
    window.show()
    return app.exec()


def _launch_welcome() -> int:
    """Launch the welcome screen defined in :mod:`spectrum_launcher`."""

    from .spectrum_launcher import main as launcher_main

    return launcher_main()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments for the launcher."""

    parser = argparse.ArgumentParser(
        description="Launch the SpectrumPY quicklook tools from a single entry point.",
    )
    subparsers = parser.add_subparsers(dest="command")

    launcher_parser = subparsers.add_parser(
        "launcher",
        help="Open the welcome screen (default).",
    )
    launcher_parser.set_defaults(func=lambda args: _launch_welcome())

    quicklook_parser = subparsers.add_parser(
        "quicklook",
        help="Open the IDEX Quicklook window directly.",
    )
    quicklook_parser.add_argument(
        "--filename",
        metavar="PATH",
        default=None,
        help="Optional path to a science data product (HDF5/CDF/trace).",
    )
    quicklook_parser.add_argument(
        "--eventnumber",
        type=int,
        default=None,
        help="1-based event index to select after loading a file.",
    )
    quicklook_parser.set_defaults(func=lambda args: _launch_quicklook(args.filename, args.eventnumber))

    parser.set_defaults(func=lambda args: _launch_welcome())

    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Entry point compatible with ``python -m imfpy.gui.main``."""

    args = parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
