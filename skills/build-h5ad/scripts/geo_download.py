#!/usr/bin/env python3
"""Backward-compatible wrapper — delegates to download.py.

Usage (still works):
    python3 geo_download.py GSE139129
    python3 geo_download.py GSE139129 --output-dir raw_data/GSE139129/
    python3 geo_download.py GSE139129 --skip-download

Now also works for non-GEO accessions:
    python3 geo_download.py E-MTAB-13502
"""

import sys
from pathlib import Path

# Ensure sibling scripts are importable
_scripts_dir = Path(__file__).resolve().parent
if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))

from download import main

if __name__ == "__main__":
    main()
