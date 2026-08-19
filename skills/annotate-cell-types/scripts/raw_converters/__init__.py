"""Raw converter auto-discovery and registry.

Scans this directory for BaseRawConverter subclasses.
Priority: dataset-specific converters before generic ones.
"""

import importlib
import inspect
import logging
import pkgutil
from pathlib import Path

logger = logging.getLogger(__name__)

from .base import BaseRawConverter
from .common import (
    TenXMEXConverter, TenXH5Converter,
    SingleCSVConverter, MultiCSVConverter, H5adPassthroughConverter,
)
from .emtab9154 import EMTAB9154RawConverter
try:
    from .gse236030 import GSE236030RawConverter
except ImportError:
    GSE236030RawConverter = None
try:
    from .zenodo_nested_mex import ZenodoNestedMEXConverter
except ImportError:
    ZenodoNestedMEXConverter = None
try:
    from .gse190260 import GSE190260RawConverter
except ImportError:
    GSE190260RawConverter = None
try:
    from .gse245577 import GSE245577Converter
except ImportError:
    GSE245577Converter = None
try:
    from .gse260997 import GSE260997Converter
except ImportError:
    GSE260997Converter = None
try:
    from .gse190447 import GSE190447Converter
except ImportError:
    GSE190447Converter = None
from .gse157526 import GSE157526RawConverter
from .gse172138 import GSE172138Converter
from .gse186814 import GSE186814Converter

# Registry: converter_name -> converter_class
# Built-in converters are registered explicitly.
# Additional converters in this directory are auto-discovered.

_REGISTRY = {}

# Explicit registration order (dataset-specific first, then generic):
_BUILTIN = [
    EMTAB9154RawConverter,
    GSE236030RawConverter,
    ZenodoNestedMEXConverter,
    GSE190260RawConverter,
    GSE186814Converter,
    GSE157526RawConverter,
    GSE172138Converter,
    GSE260997Converter,
    GSE245577Converter,
    GSE190447Converter,
    TenXMEXConverter,
    TenXH5Converter,
    SingleCSVConverter,
    MultiCSVConverter,
    H5adPassthroughConverter,
]
# Filter out None entries from failed imports
_BUILTIN = [c for c in _BUILTIN if c is not None]


def _register(converter_cls):
    name = converter_cls.__name__
    _REGISTRY[name] = converter_cls


def _auto_discover():
    """Scan raw_converters directory for BaseRawConverter subclasses."""
    _REGISTRY.clear()
    for cls in _BUILTIN:
        _register(cls)
    # Also scan for additional modules
    try:
        import raw_converters as pkg
        pkg_path = Path(pkg.__path__[0])
        for _, mod_name, _ in pkgutil.iter_modules([str(pkg_path)]):
            if mod_name in ("base", "common", "__init__", "emtab9154", "gse245577", "gse190447", "gse157526", "gse172138"):
                continue
            try:
                mod = importlib.import_module(f".{mod_name}", package="raw_converters")
                for _, obj in inspect.getmembers(mod, inspect.isclass):
                    if (issubclass(obj, BaseRawConverter)
                            and obj is not BaseRawConverter
                            and obj.__name__ not in _REGISTRY):
                        _register(obj)
                        logger.info("Auto-discovered converter: %s", obj.__name__)
            except Exception as e:
                logger.debug("Skipping %s: %s", mod_name, e)
    except Exception:
        pass


_auto_discover()


def detect_format(raw_dir: str) -> str:
    """Detect which converter can handle a raw directory.

    Returns converter name, or raises ValueError.
    Dataset-specific converters are tried before generic ones.
    """
    path = Path(raw_dir)
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {raw_dir}")

    for name in [c.__name__ for c in _BUILTIN]:
        if name not in _REGISTRY:
            continue
        cls = _REGISTRY[name]
        try:
            if cls.detect(path):
                return name
        except Exception:
            continue

    # Try any remaining auto-discovered converters
    for name, cls in list(_REGISTRY.items()):
        if name not in [c.__name__ for c in _BUILTIN]:
            try:
                if cls.detect(path):
                    return name
            except Exception:
                continue

    raise ValueError(
        f"No raw converter found for: {raw_dir}. "
        f"Available converters: {list(_REGISTRY.keys())}"
    )


def get_converter(name: str) -> type:
    """Get converter class by name."""
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown converter: {name}. Available: {list(_REGISTRY.keys())}"
        )
    return _REGISTRY[name]


def list_converters() -> dict:
    """Return all registered converters."""
    return dict(_REGISTRY)
