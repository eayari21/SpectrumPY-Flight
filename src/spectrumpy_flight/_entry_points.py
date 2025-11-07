"""Console-script adapters for legacy top-level scripts."""

from __future__ import annotations

import importlib.util
from importlib import resources
import sys
from types import ModuleType

_SCRIPT_ALIASES = {
    "IDEX-quicklook.py": "spectrumpy_flight.IDEX_quicklook",
    "Scope-IDEX-quicklook.py": "spectrumpy_flight.Scope_IDEX_quicklook",
}


def _load_script(script: str) -> ModuleType:
    """Load a script stored alongside the package and return the module."""

    resource = resources.files("spectrumpy_flight") / script
    if not resource.is_file():
        raise FileNotFoundError(f"Unable to locate '{script}' inside the installed package")

    alias = _SCRIPT_ALIASES.get(script, f"spectrumpy_flight.{script.replace('-', '_').replace('.py', '')}")

    with resources.as_file(resource) as module_path:
        spec = importlib.util.spec_from_file_location(alias, module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Failed to create an import spec for '{script}'")

        module = importlib.util.module_from_spec(spec)
        sys.modules.setdefault(alias, module)
        spec.loader.exec_module(module)
        return module


def _call_main(module: ModuleType) -> int:
    """Invoke ``main`` from *module* and normalise the return value."""

    main = getattr(module, "main", None)
    if main is None:
        raise AttributeError(f"Module '{module.__name__}' does not define a main() function")

    result = main()
    return int(result) if result is not None else 0


def quicklook() -> int:
    """Entry point for the desktop quicklook GUI."""

    module = _load_script("IDEX-quicklook.py")
    return _call_main(module)


def scope_quicklook() -> int:
    """Entry point for the oscilloscope-tuned quicklook variant."""

    module = _load_script("Scope-IDEX-quicklook.py")
    return _call_main(module)


__all__ = ["quicklook", "scope_quicklook"]
