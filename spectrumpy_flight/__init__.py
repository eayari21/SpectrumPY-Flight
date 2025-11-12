"""Runtime shim to allow importing the package directly from the repository root."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

__all__ = []


def _load_package() -> ModuleType:
    """Load ``src/spectrumpy_flight`` and return it as the active package."""

    repo_root = Path(__file__).resolve().parent.parent
    package_dir = repo_root / "src" / "spectrumpy_flight"
    init_path = package_dir / "__init__.py"

    if not init_path.is_file():
        raise ModuleNotFoundError(
            "Cannot locate 'spectrumpy_flight' inside the repository."
        )

    spec = importlib.util.spec_from_file_location(
        __name__,
        init_path,
        submodule_search_locations=[str(package_dir)],
    )
    if spec is None or spec.loader is None:
        raise ModuleNotFoundError("Unable to create a module spec for 'spectrumpy_flight'.")

    module = importlib.util.module_from_spec(spec)
    sys.modules[__name__] = module
    spec.loader.exec_module(module)
    return module


_module = _load_package()
globals().update(_module.__dict__)

