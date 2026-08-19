"""Auto-discovering converter registry.

Scans this directory for BaseConverter subclasses at import time.
Add a new converter by dropping a gseXXXXX.py file here — no registration needed.
"""

import importlib
import pkgutil
from pathlib import Path

from .base import BaseConverter

_registry = None  # dict[str, type[BaseConverter]] | None


def _discover():
    """Import all converter modules and collect BaseConverter subclasses."""
    global _registry
    if _registry is not None:
        return _registry

    _registry = {}

    package_dir = Path(__file__).parent
    for module_info in pkgutil.iter_modules([str(package_dir)]):
        if module_info.name in ("base", "__init__"):
            continue
        try:
            importlib.import_module(f".{module_info.name}", package=__package__)
        except Exception:
            import sys
            print(
                f"  Warning: skipping converter {module_info.name} (import error)",
                file=sys.stderr,
            )

    # Collect all concrete subclasses (skip intermediate abstract classes)
    for cls in BaseConverter.__subclasses__():
        # Name: GSE139129Converter → gse139129
        name = cls.__name__.replace("Converter", "").lower()
        _registry[name] = cls

    return _registry


def get_converter(name: str):
    registry = _discover()
    if name not in registry:
        available = sorted(registry.keys())
        raise ValueError(
            f"Unknown converter: {name}. Available: {available}"
        )
    return registry[name]


def detect_format(raw_dir: Path) -> str:
    """Try each converter's detect() and return the first match."""
    registry = _discover()
    for name, cls in registry.items():
        try:
            if cls.detect(raw_dir):
                return name
        except Exception:
            continue
    available = sorted(registry.keys())
    raise ValueError(
        f"No converter could detect data format in {raw_dir}. "
        f"Available converters: {available}"
    )


def list_converters():
    return dict(_discover())
