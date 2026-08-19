"""Preflight dependency check for annotate-cell-types.

Checks that all required Python packages are installed before the pipeline
runs. Reports missing packages clearly instead of cryptic import errors.
"""

import importlib
import sys


REQUIRED_PACKAGES = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("scipy", "scipy"),
    ("anndata", "anndata"),
    ("scanpy", "scanpy"),
    ("celltypist", "celltypist"),
    ("openpyxl", "openpyxl"),
]

# Optional packages — warn if missing but don't fail
OPTIONAL_PACKAGES = [
    ("sklearn", "scikit-learn"),
]


def check_dependencies() -> dict:
    """Check required and optional packages.

    Returns:
        dict with keys:
            all_ok: bool
            missing: list[str]
            optional_missing: list[str]
            details: list[dict]
    """
    missing = []
    optional_missing = []
    details = []

    for import_name, pkg_name in REQUIRED_PACKAGES:
        try:
            mod = importlib.import_module(import_name)
            ver = getattr(mod, "__version__", "?")
            details.append({
                "package": pkg_name,
                "status": "ok",
                "version": str(ver),
            })
        except ImportError:
            missing.append(pkg_name)
            details.append({
                "package": pkg_name,
                "status": "missing",
                "version": None,
            })

    for import_name, pkg_name in OPTIONAL_PACKAGES:
        try:
            importlib.import_module(import_name)
            details.append({
                "package": pkg_name,
                "status": "ok_optional",
                "version": None,
            })
        except ImportError:
            optional_missing.append(pkg_name)
            details.append({
                "package": pkg_name,
                "status": "missing_optional",
                "version": None,
            })

    return {
        "all_ok": len(missing) == 0,
        "missing": missing,
        "optional_missing": optional_missing,
        "details": details,
        "python_version": sys.version,
        "python_executable": sys.executable,
    }


def check_and_report() -> int:
    """Check dependencies and print a report. Returns exit code (0 = OK)."""
    result = check_dependencies()

    print("=== annotate-cell-types dependency check ===")
    print(f"  Python: {result['python_version']}")
    print(f"  Path:   {result['python_executable']}")
    print()

    for d in result["details"]:
        status_icon = {
            "ok": "[OK]",
            "missing": "[MISSING]",
            "ok_optional": "[OPT]",
            "missing_optional": "[OPT-MISS]",
        }.get(d["status"], "[?]")
        ver_str = f" ({d['version']})" if d["version"] else ""
        print(f"  {status_icon} {d['package']}{ver_str}")

    print()

    if result["missing"]:
        print("MISSING REQUIRED PACKAGES:")
        for pkg in result["missing"]:
            print(f"  - {pkg}")
        print()
        print("Install with:")
        print(f"  pip install {' '.join(result['missing'])}")
        print()
        print("Or from requirements file:")
        print("  pip install -r skills/annotate-cell-types/requirements.txt")
        return 1

    if result["optional_missing"]:
        print("Optional packages not installed:")
        for pkg in result["optional_missing"]:
            print(f"  - {pkg}")
        print("  (not required for default operation)")

    print("All required dependencies OK.")
    return 0


if __name__ == "__main__":
    sys.exit(check_and_report())
