"""Compatibility wrapper that exposes :mod:`spectrumpy_flight` as ``spectrumpy``.

The PyPI distribution is named ``spectrumpy`` so users can install the toolkit
with ``pip install spectrumpy``.  Historically the code lived under the
``spectrumpy_flight`` package.  Importing this module re-exports the public API
and registers aliases for the high-traffic submodules so existing notebooks and
scripts keep working.
"""

from __future__ import annotations

import sys
from importlib import import_module
from typing import Any

_spectrumpy_flight = import_module("spectrumpy_flight")

__version__ = getattr(_spectrumpy_flight, "__version__", "0.0.0")
__all__ = list(getattr(_spectrumpy_flight, "__all__", []))


def __getattr__(name: str) -> Any:
    """Proxy attribute lookups to :mod:`spectrumpy_flight`.

    When ``name`` matches a submodule we also register a ``spectrumpy`` alias in
    :data:`sys.modules` so ``import spectrumpy.<module>`` returns the original
    implementation.
    """

    try:
        value = getattr(_spectrumpy_flight, name)
    except AttributeError:
        module_name = f"spectrumpy_flight.{name}"
        module = import_module(module_name)
        sys.modules.setdefault(f"spectrumpy.{name}", module)
        return module
    else:
        return value


def __dir__() -> list[str]:
    return sorted(set(list(globals()) + dir(_spectrumpy_flight)))


for _name in __all__:
    globals()[_name] = getattr(_spectrumpy_flight, _name)


for _alias in ("idex_packet", "drive_idex_packet", "start"):
    module = import_module(f"spectrumpy_flight.{_alias}")
    sys.modules.setdefault(f"spectrumpy.{_alias}", module)
