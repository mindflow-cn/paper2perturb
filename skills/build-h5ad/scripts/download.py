#!/usr/bin/env python3
"""Download scRNA-seq data from public repositories (GEO, ArrayExpress, etc.).

Usage:
    python3 download.py GSE139129
    python3 download.py E-MTAB-13502
    python3 download.py GSE139129 --output-dir raw_data/GSE139129/
    python3 download.py GSE139129 --skip-download  # extract only
"""

import argparse
import ftplib
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from urllib.parse import urljoin, urlparse
from urllib.request import urlopen, Request

# ---- Provider detection ----

# Provider registry: (pattern, display_name, download_function)
# Ordered by priority — first match wins.
_PROVIDERS = []


def _provider(pattern, name):
    """Decorator to register a download function for an accession pattern."""
    def decorator(fn):
        _PROVIDERS.append((re.compile(pattern, re.IGNORECASE), name, fn))
        return fn
    return decorator


def detect_provider(accession: str) -> tuple[str, callable]:
    """Return (provider_name, download_fn) for an accession, or raise ValueError."""
    for pattern, name, fn in _PROVIDERS:
        if pattern.match(accession):
            return name, fn
    raise ValueError(
        f"Unknown accession prefix: {accession}. "
        f"Supported: GEO (GSE*), ArrayExpress/BioStudies (E-MTAB-*, E-GEOD-*, E-MEXP-*), "
        f"Zenodo (Zenodo <record_id>), CELLxGENE (cellxgene:<collection_uuid>), "
        f"GitHub (GitHub: owner/repo), EGA (EGAS*, EGAD*), Figshare (figshare: <article_id>)"
    )


# ---- Utility functions ----

def _download_file(url: str, dest: Path) -> bool:
    """Download a file via wget or curl. Skips if dest exists."""
    if dest.exists():
        print(f"    Skipping (exists): {dest.name}")
        return True

    print(f"    Downloading: {dest.name}")
    try:
        subprocess.run(
            ["wget", "-q", "--show-progress", "-O", str(dest), url],
            check=True,
        )
        return True
    except subprocess.CalledProcessError:
        try:
            subprocess.run(
                ["curl", "-L", "-o", str(dest), url],
                check=True,
            )
            return True
        except subprocess.CalledProcessError:
            print(f"    ERROR: Failed to download {url}", file=sys.stderr)
            return False


def _extract_archive(archive_path: Path, output_dir: Path) -> list[Path]:
    """Extract .tar or .tar.gz archive. Returns list of extracted file paths."""
    extracted = []
    try:
        with tarfile.open(archive_path, "r:*") as tar:
            members = tar.getmembers()
            print(f"  Extracting {len(members)} file(s) from {archive_path.name}...")
            tar.extractall(path=output_dir)
            extracted = [output_dir / m.name for m in members if m.isfile()]
    except tarfile.TarError as e:
        print(f"  ERROR extracting {archive_path.name}: {e}", file=sys.stderr)
    return extracted


def _extract_gz_files(directory: Path) -> list[Path]:
    """Gunzip all .gz files in directory (skip .tar.gz). Returns extracted paths."""
    extracted = []
    for gz_path in sorted(directory.glob("*.gz")):
        if gz_path.name.endswith(".tar.gz"):
            continue
        out_path = gz_path.with_suffix("")
        if out_path.exists():
            print(f"    Skipping (exists): {out_path.name}")
            extracted.append(out_path)
            continue
        print(f"    Gunzipping: {gz_path.name}")
        with gzip.open(gz_path, "rb") as f_in:
            with open(out_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        extracted.append(out_path)
    return extracted


def _post_download_extract(output_dir: Path):
    """Extract archives and gunzip raw .gz files (when converters can't read them)."""
    for archive in sorted(output_dir.glob("*.tar*")):
        _extract_archive(archive, output_dir)

    _gz_files = list(output_dir.glob("*.gz"))
    _has_data_files = bool(
        list(output_dir.glob("*.txt.gz"))
        or list(output_dir.glob("*.csv.gz"))
        or list(output_dir.glob("*.tsv.gz"))
    )
    if not _has_data_files:
        _extract_gz_files(output_dir)


# ---- GEO provider ----

GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"
GEO_HTTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _gse_nnn(gse: str) -> str:
    """Extract first 3 digits of GSE ID for path construction."""
    digits = re.search(r"\d+", gse)
    if not digits:
        raise ValueError(f"Cannot parse GSE ID: {gse}")
    return digits.group()[:3]


def _geo_suppl_path(gse: str) -> str:
    return f"/geo/series/GSE{_gse_nnn(gse)}nnn/{gse}/suppl"


def _geo_suppl_url(gse: str) -> str:
    return f"{GEO_HTTP_BASE}/GSE{_gse_nnn(gse)}nnn/{gse}/suppl"


def _geo_ftp_list(gse: str) -> list[str]:
    """List supplementary files via FTP."""
    path = _geo_suppl_path(gse)
    print(f"  Listing FTP: {GEO_FTP_HOST}{path}")
    try:
        ftp = ftplib.FTP(GEO_FTP_HOST, timeout=30)
        ftp.login()
        ftp.cwd(path)
        files = []
        ftp.retrlines("LIST", files.append)
        ftp.quit()

        result = []
        for line in files:
            parts = line.split()
            if len(parts) >= 9:
                fname = " ".join(parts[8:])
                if not fname.startswith("."):
                    result.append(fname)
        return result
    except Exception as e:
        print(f"  FTP listing failed: {e}")
        return []


@_provider(r"^GSE\d+", "NCBI GEO")
def download_geo(accession: str, output_dir: Path):
    """Download supplementary files for a GEO accession.

    Strategy: List via FTP, then download each file via HTTP.
    Falls back to downloading the RAW tar via the GEO download endpoint.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    files = _geo_ftp_list(accession)
    if files:
        print(f"  Found {len(files)} file(s) via FTP")
        base_url = _geo_suppl_url(accession)
        for fname in files:
            url = f"{base_url}/{fname}"
            dest = output_dir / fname
            if _download_file(url, dest):
                downloaded.append(dest)
        if downloaded:
            return downloaded

    # Fallback: direct HTTP RAW tar download
    print("  Trying GEO download endpoint...")
    raw_tar = f"{accession}_RAW.tar"
    url = f"{_geo_suppl_url(accession)}/{raw_tar}"
    dest = output_dir / raw_tar
    if _download_file(url, dest):
        downloaded.append(dest)

    return downloaded


# ---- ArrayExpress / BioStudies provider ----

BIOSTUDIES_BASE = "https://www.ebi.ac.uk/biostudies/files"


def _biostudies_list_files(accession: str) -> list[str]:
    """Fetch file listing from BioStudies JSON metadata."""
    url = f"{BIOSTUDIES_BASE}/{accession}/{accession}.json"
    print(f"  Fetching metadata: {url}")
    try:
        req = Request(url, headers={"User-Agent": "build-h5ad/1.0"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())

        def find_files(obj):
            if isinstance(obj, dict):
                if obj.get("type") == "file" and "path" in obj:
                    yield obj["path"]
                for v in obj.values():
                    yield from find_files(v)
            elif isinstance(obj, list):
                for v in obj:
                    yield from find_files(v)

        return list(find_files(data))
    except Exception as e:
        print(f"  BioStudies API error: {e}")
        return []


def _is_processed_matrix(filename: str) -> bool:
    """Check if a file is a processed count matrix (MTX, features, barcodes)."""
    return bool(re.search(r"_(raw|filtered)_(matrix\.mtx|barcodes|features)\.(tsv|csv|txt)\.gz$", filename))


def _is_raw_data_tar(filename: str) -> bool:
    """Check if a file is a raw data archive."""
    return bool(
        re.search(r"\.(tar|tar\.gz)$", filename)
        or re.search(r"\.(fastq|fq|bam)\.gz$", filename)
    )


@_provider(r"^E-(MTAB|GEOD|MEXP)-\d+", "EBI ArrayExpress / BioStudies")
def download_arrayexpress(accession: str, output_dir: Path):
    """Download processed data files from ArrayExpress via BioStudies API.

    Prioritises processed count matrices (.mtx.gz + barcodes + features).
    Falls back to all files if no processed data found.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    all_files = _biostudies_list_files(accession)
    if not all_files:
        print(f"  No files found for {accession} in BioStudies")
        return downloaded

    print(f"  Found {len(all_files)} file(s) in BioStudies")

    # Prioritise processed matrix files; fall back to all non-raw files
    processed = [f for f in all_files if _is_processed_matrix(f)]
    candidates = processed if processed else [
        f for f in all_files if not _is_raw_data_tar(f)
    ]
    # Always include SDRF/IDF metadata files
    metadata_files = [f for f in all_files if f.endswith(('.sdrf.txt', '.idf.txt'))]
    for mf in metadata_files:
        if mf not in candidates:
            candidates.append(mf)

    if not candidates:
        # Last resort: download everything
        candidates = all_files

    skipped_raw = len(all_files) - len(candidates)
    if skipped_raw > 0:
        print(f"  Skipping {skipped_raw} raw-data file(s) (FASTQ/BAM/tar)")

    print(f"  Downloading {len(candidates)} file(s)...")
    for fname in sorted(candidates):
        url = f"{BIOSTUDIES_BASE}/{accession}/{fname}"
        dest = output_dir / fname
        if _download_file(url, dest):
            downloaded.append(dest)

    return downloaded


# ---- Zenodo provider ----

ZENODO_API_BASE = "https://zenodo.org/api/records"


def _parse_zenodo_id(accession: str) -> str:
    """Extract Zenodo record ID from various accession formats.

    Supported formats:
      - "Zenodo 7942968"             → "7942968"
      - "10.5281/zenodo.7942968"     → "7942968" (from DOI)
      - "Zenodo: 10.5281/zenodo.19364848" → "19364848"
    """
    # "Zenodo 7942968" — space-separated
    m = re.search(r"(?i)zenodo\s+(\d+)\s*$", accession)
    if m:
        return m.group(1)
    # "Zenodo: 7942968" or "Zenodo:7942968" — colon-separated
    m = re.search(r"(?i)zenodo:\s*(\d+)", accession)
    if m:
        return m.group(1)
    # DOI format — "zenodo.<id>" or "/zenodo.<id>"
    m = re.search(r"(?i)zenodo[./](\d+)", accession)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot extract Zenodo record ID from: {accession}")


@_provider(r"(?i)zenodo", "Zenodo")
def download_zenodo(accession: str, output_dir: Path):
    """Download files from a Zenodo record via the REST API.

    Uses https://zenodo.org/api/records/{record_id} to discover files,
    then downloads each one.
    """
    record_id = _parse_zenodo_id(accession)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    api_url = f"{ZENODO_API_BASE}/{record_id}"
    print(f"  Fetching metadata: {api_url}")
    try:
        req = Request(api_url, headers={"User-Agent": "build-h5ad/1.0"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  Zenodo API error: {e}", file=sys.stderr)
        return downloaded

    files = data.get("files", [])
    if not files:
        print(f"  No files found in Zenodo record {record_id}")
        return downloaded

    print(f"  Found {len(files)} file(s) in Zenodo record {record_id}")
    for f in files:
        fname = f["key"]
        download_url = f["links"]["self"]
        size_mb = f.get("size", 0) / (1024 * 1024)
        print(f"    {fname} ({size_mb:.1f} MB)")
        dest = output_dir / fname
        if _download_file(download_url, dest):
            downloaded.append(dest)

    return downloaded


# ---- CELLxGENE (CZ CELLxGENE Discover) provider ----

CELLXGENE_COLLECTION_API = "https://api.cellxgene.cziscience.com/curation/v1/collections"


def _parse_cellxgene_id(accession: str) -> str:
    """Extract CELLxGENE collection UUID from accession string.

    Supported formats:
      - "cellxgene:c2879de0-affc-496b-8e2b-f57ed9ec3c34"
      - "cellxgene c2879de0-affc-496b-8e2b-f57ed9ec3c34"
    """
    m = re.search(r"(?i)cellxgene[:\s]+([a-f0-9-]{36})", accession)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot extract CELLxGENE collection ID from: {accession}")


@_provider(r"(?i)cellxgene", "CZ CELLxGENE Discover")
def download_cellxgene(accession: str, output_dir: Path):
    """Download h5ad files from a CZ CELLxGENE Discover collection.

    Queries the CELLxGENE curation API for the collection, then downloads
    all h5ad assets from each dataset in the collection.
    """
    collection_id = _parse_cellxgene_id(accession)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    api_url = f"{CELLXGENE_COLLECTION_API}/{collection_id}"
    print(f"  Fetching metadata: {api_url}")
    try:
        req = Request(api_url, headers={"User-Agent": "build-h5ad/1.0"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  CELLxGENE API error: {e}", file=sys.stderr)
        return downloaded

    datasets = data.get("datasets", [])
    if not datasets:
        print(f"  No datasets found in collection {collection_id}")
        return downloaded

    print(f"  Collection: {data.get('name', 'N/A')[:80]}")
    print(f"  Datasets: {len(datasets)}")

    for ds in datasets:
        ds_id = ds.get("dataset_id", "unknown")
        title = ds.get("title", "")[:80]
        cell_count = ds.get("cell_count", "?")
        print(f"    Dataset {ds_id}: {title} ({cell_count} cells)")

        for asset in ds.get("assets", []):
            fname_base = f"{ds_id}"
            filetype = asset.get("filetype", "").lower()
            ext = ".h5ad" if filetype == "h5ad" else f".{filetype}"

            fname = f"{fname_base}{ext}"
            download_url = asset["url"]
            size_mb = asset.get("filesize", 0) / (1024 * 1024)
            print(f"      {filetype.upper()}: {download_url} ({size_mb:.1f} MB)")
            dest = output_dir / fname
            if _download_file(download_url, dest):
                downloaded.append(dest)

    return downloaded


# ---- Google Drive provider ----


def _parse_gdrive_id(accession: str) -> str:
    """Extract Google Drive file/folder ID from accession string.

    Supported formats:
      - "custom_gdrive_18-KInmm43wKdBX95Gq9zbuzAQwtLjgE9"
      - "custom_gdrive_1abc123xyz"
    """
    m = re.search(r"(?i)custom_gdrive[_:]?([a-zA-Z0-9_-]{20,})", accession)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot extract Google Drive ID from: {accession}")


@_provider(r"(?i)custom_gdrive", "Google Drive")
def download_gdrive(accession: str, output_dir: Path):
    """Download files from a Google Drive file or folder.

    Uses gdown to handle Google Drive downloads, including large files
    that require confirmation tokens.
    """
    import gdown

    file_id = _parse_gdrive_id(accession)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    url = f"https://drive.google.com/uc?id={file_id}"
    print(f"  Google Drive ID: {file_id}")

    # Try as a file first; fall back to folder
    try:
        output_path = str(output_dir / f"{file_id}.download")
        result = gdown.download(url, output=output_path, quiet=False)
        if result:
            downloaded.append(Path(result))
            print(f"  Downloaded: {Path(result).name}")
            return downloaded
    except Exception as e:
        print(f"  File download failed: {e}")

    # Try as a folder
    try:
        folder_url = f"https://drive.google.com/drive/folders/{file_id}"
        print(f"  Trying folder download: {folder_url}")
        result = gdown.download_folder(folder_url, output=str(output_dir),
                                       quiet=False, use_cookies=False)
        if result:
            for path in result:
                downloaded.append(Path(path))
            return downloaded
    except Exception as e:
        print(f"  Folder download failed: {e}")

    return downloaded


# ---- Figshare provider ----

FIGSHARE_API_BASE = "https://api.figshare.com/v2/articles"


def _parse_figshare_id(accession: str) -> str:
    """Extract Figshare article ID from accession string.

    Supported formats:
      - "figshare: 28777181"    → "28777181"
      - "figshare:28777181"     → "28777181"
      - "figshare 28777181"     → "28777181"
    """
    m = re.search(r"(?i)figshare[:\s]+(\d+)", accession)
    if m:
        return m.group(1)
    raise ValueError(f"Cannot extract Figshare article ID from: {accession}")


@_provider(r"(?i)figshare", "Figshare")
def download_figshare(accession: str, output_dir: Path):
    """Download files from a Figshare article via the REST API.

    Uses https://api.figshare.com/v2/articles/{article_id} to discover files,
    then downloads each one.
    """
    article_id = _parse_figshare_id(accession)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    api_url = f"{FIGSHARE_API_BASE}/{article_id}"
    print(f"  Fetching metadata: {api_url}")
    try:
        req = Request(api_url, headers={"User-Agent": "build-h5ad/1.0"})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  Figshare API error: {e}", file=sys.stderr)
        return downloaded

    files = data.get("files", [])
    if not files:
        print(f"  No files found in Figshare article {article_id}")
        return downloaded

    print(f"  Found {len(files)} file(s) in Figshare article {article_id}")
    for f in files:
        fname = f["name"]
        download_url = f["download_url"]
        size_mb = f.get("size", 0) / (1024 * 1024)
        print(f"    {fname} ({size_mb:.1f} MB)")
        dest = output_dir / fname
        if _download_file(download_url, dest):
            downloaded.append(dest)

    return downloaded


# ---- GitHub provider ----

GITHUB_API_BASE = "https://api.github.com"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com"


def _parse_github_repo(accession: str) -> str:
    """Extract owner/repo from various GitHub accession formats.

    Supported formats:
      - "GitHub: federicogiorgi/panobinostat"  → "federicogiorgi/panobinostat"
      - "GitHub federicogiorgi/panobinostat"    → "federicogiorgi/panobinostat"
      - "github.com/federicogiorgi/panobinostat" → "federicogiorgi/panobinostat"
    """
    # "GitHub: owner/repo" or "GitHub owner/repo"
    m = re.search(r"(?i)github[:\s]+([a-zA-Z0-9_.-]+/[a-zA-Z0-9_.-]+)", accession)
    if m:
        return m.group(1)
    # "github.com/owner/repo" (with optional trailing .git or /)
    m = re.search(r"(?i)github\.com/([a-zA-Z0-9_.-]+/[a-zA-Z0-9_.-]+)", accession)
    if m:
        return m.group(1).rstrip("/").removesuffix(".git")
    raise ValueError(f"Cannot extract GitHub owner/repo from: {accession}")


def _github_get_default_branch(owner: str, repo: str) -> str:
    """Get the default branch of a GitHub repo."""
    url = f"{GITHUB_API_BASE}/repos/{owner}/{repo}"
    try:
        req = Request(url, headers={
            "User-Agent": "build-h5ad/1.0",
            "Accept": "application/vnd.github+json",
        })
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
        return data.get("default_branch", "main")
    except Exception:
        return "main"


def _github_list_files(owner: str, repo: str, branch: str | None = None) -> list[dict]:
    """List all files in a GitHub repository via the Git Trees API."""
    if branch is None:
        branch = _github_get_default_branch(owner, repo)

    url = f"{GITHUB_API_BASE}/repos/{owner}/{repo}/git/trees/{branch}?recursive=1"
    print(f"  Listing files: {url}")
    try:
        req = Request(url, headers={
            "User-Agent": "build-h5ad/1.0",
            "Accept": "application/vnd.github+json",
        })
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
    except Exception as e:
        print(f"  GitHub API error: {e}", file=sys.stderr)
        return []

    files = []
    for item in data.get("tree", []):
        if item.get("type") == "blob":
            files.append({
                "path": item["path"],
                "size": item.get("size", 0),
                "url": f"{GITHUB_RAW_BASE}/{owner}/{repo}/{branch}/{item['path']}",
            })
    return files


def _github_is_lfs_pointer(raw_url: str) -> tuple[bool, str | None, int | None]:
    """Check if a remote file is a Git LFS pointer. Returns (is_lfs, oid, size)."""
    try:
        req = Request(raw_url, headers={"User-Agent": "build-h5ad/1.0"})
        with urlopen(req, timeout=10) as resp:
            content = resp.read(1024)
        text = content.decode("utf-8", errors="ignore")
        if text.startswith("version https://git-lfs.github.com/spec/v1"):
            m_oid = re.search(r"oid sha256:([a-f0-9]+)", text)
            m_size = re.search(r"size (\d+)", text)
            oid = m_oid.group(1) if m_oid else None
            size = int(m_size.group(1)) if m_size else None
            return True, oid, size
    except Exception:
        pass
    return False, None, None


def _github_clone_lfs_files(owner: str, repo: str, file_paths: list[str],
                             output_dir: Path) -> list[Path]:
    """Shallow-clone a GitHub repo with git-lfs and copy specified files out.

    Handles repos whose LFS objects are only served to authenticated git-lfs
    clients (batch API alone may return 404 for unauthenticated requests).
    """
    import tempfile

    clone_url = f"https://github.com/{owner}/{repo}.git"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmppath = Path(tmpdir)

        print(f"    Cloning {owner}/{repo} (shallow) with git-lfs...")
        result = subprocess.run(
            ["git", "clone", "--depth", "1", clone_url, str(tmppath)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"    Clone failed: {result.stderr.strip()}", file=sys.stderr)
            return []

        # Pull LFS objects (git clone will have already smudged them;
        # this is a safety net if smudge was skipped)
        subprocess.run(
            ["git", "-C", str(tmppath), "lfs", "pull"],
            capture_output=True, text=True,
        )

        downloaded = []
        for rel_path in file_paths:
            src = tmppath / rel_path
            dest = output_dir / Path(rel_path).name
            if src.exists() and src.stat().st_size > 200:
                shutil.copy2(src, dest)
                downloaded.append(dest)
            else:
                print(f"    Skipping {rel_path}: missing or still an LFS pointer",
                      file=sys.stderr)

        return downloaded


@_provider(r"(?i)github", "GitHub")
def download_github(accession: str, output_dir: Path):
    """Download files from a GitHub repository.

    Non-LFS files are downloaded directly via raw.githubusercontent.com.
    LFS-backed files are obtained via shallow git clone with git-lfs.
    """
    owner_repo = _parse_github_repo(accession)
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    owner, repo = owner_repo.split("/", 1)
    print(f"  Repository: {owner}/{repo}")

    files = _github_list_files(owner, repo)
    if not files:
        print(f"  No files found in {owner}/{repo}")
        return downloaded

    # Prioritise data files, skip code/scripts
    _data_exts = {".rda", ".rds", ".rdata", ".h5ad", ".h5",
                  ".csv", ".tsv", ".txt", ".mtx", ".gz", ".zip", ".tar.gz"}
    data_files = [f for f in files
                  if any(f["path"].lower().endswith(ext) for ext in _data_exts)]

    if not data_files:
        # Fall back to all non-script files
        _code_exts = {".R", ".Rmd", ".md", ".Rproj", ".py", ".sh"}
        data_files = [f for f in files
                      if not any(f["path"].endswith(ext) for ext in _code_exts)]

    print(f"  Found {len(data_files)} data file(s) among {len(files)} total files")

    # Separate LFS and non-LFS files
    lfs_paths = []
    for f in data_files:
        is_lfs, oid, lfs_size = _github_is_lfs_pointer(f["url"])
        if is_lfs and oid:
            size_mb = (lfs_size or 0) / (1024 * 1024)
            print(f"    {f['path']} [LFS] ({size_mb:.1f} MB)")
            lfs_paths.append(f["path"])
        else:
            size_mb = f.get("size", 0) / (1024 * 1024)
            print(f"    {f['path']} ({size_mb:.1f} MB)")
            dest = output_dir / Path(f["path"]).name
            if _download_file(f["url"], dest):
                downloaded.append(dest)

    # Download LFS files in one shallow clone
    if lfs_paths:
        lfs_downloaded = _github_clone_lfs_files(owner, repo, lfs_paths, output_dir)
        downloaded.extend(lfs_downloaded)

    return downloaded


# ---- EGA (European Genome-Phenome Archive) provider ----


def _ega_parse_study_datasets(accession: str) -> list[str]:
    """Parse EGA study page to find scRNA-seq datasets.

    Returns list of dataset IDs (EGAD...) from the study page HTML.
    Filters for datasets whose description mentions scRNA/RNA/expression/single cell.
    """
    import urllib.request

    study_url = f"https://ega-archive.org/studies/{accession}"
    print(f"  Fetching study page: {study_url}")
    try:
        req = urllib.request.Request(
            study_url,
            headers={"User-Agent": "build-h5ad/1.0"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            html = resp.read().decode("utf-8", errors="ignore")
    except Exception as e:
        print(f"  ERROR fetching study page: {e}", file=sys.stderr)
        return []

    # Extract dataset IDs from links: /datasets/EGAD00001011155
    dataset_ids = list(set(re.findall(r'/datasets/(EGAD\d+)', html)))
    if not dataset_ids:
        print(f"  No datasets found on study page")
        return []

    # Filter for scRNA-seq datasets by looking at surrounding text for each link
    sc_datasets = []
    for ds_id in dataset_ids:
        # Find the label/description near this dataset link
        pattern = rf'<a[^>]*href="[^"]*{ds_id}[^"]*"[^>]*>(.*?)</a>.*?<label[^>]*for="{ds_id}"[^>]*>(.*?)</label>'
        m = re.search(pattern, html, re.DOTALL)
        description = ""
        if m:
            description = m.group(1) + " " + m.group(2)

        rna_keywords = ["scRNA", "single cell RNA", "gene expression", "cDNA",
                       "RNA sequencing", "expression", "transcript"]
        is_rna_seq = any(kw.lower() in description.lower() for kw in rna_keywords)
        is_not_dna = "DNA" not in description or "cDNA" in description
        is_not_bulk = "bulk" not in description.lower()

        if is_rna_seq and is_not_dna and is_not_bulk:
            sc_datasets.append(ds_id)
            print(f"    Found scRNA-seq dataset: {ds_id}")
        else:
            print(f"    Skipping non-scRNA dataset: {ds_id}")

    return sc_datasets


def _find_ega_credentials() -> dict | None:
    """Find EGA credentials from config files or environment variables.

    Checks (in order):
      1. ~/.ega/credential.json
      2. ~/.egac.json
      3. EGA_USERNAME + EGA_PASSWORD env vars
    """
    cred_file = Path.home() / ".ega" / "credential.json"
    if cred_file.exists():
        try:
            with open(cred_file) as f:
                return json.load(f)
        except Exception:
            pass

    egac_file = Path.home() / ".egac.json"
    if egac_file.exists():
        try:
            with open(egac_file) as f:
                return json.load(f)
        except Exception:
            pass

    username = os.environ.get("EGA_USERNAME")
    password = os.environ.get("EGA_PASSWORD")
    if username and password:
        return {"username": username, "password": password}

    return None


@_provider(r"^EGAS\d+", "EGA Study")
def download_ega_study(accession: str, output_dir: Path):
    """Download scRNA-seq data from an EGA study.

    Parses the study page for relevant datasets, then downloads via pyega3.
    Requires EGA credentials in ~/.ega/credential.json or env vars.
    """
    import tempfile

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    # Find credentials
    creds = _find_ega_credentials()
    if not creds:
        print(
            "  ERROR: EGA credentials not found.\n"
            "  Create ~/.ega/credential.json with:\n"
            '    {"username": "your_email@example.com", "password": "your_password"}\n'
            "  Or set EGA_USERNAME and EGA_PASSWORD environment variables.\n"
            "  Register at https://ega-archive.org/ and request access to the dataset.",
            file=sys.stderr,
        )
        return downloaded

    # Write credentials to temp file for pyega3
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as cf:
        json.dump(creds, cf)
        cred_path = cf.name

    try:
        # Find scRNA-seq datasets from study page
        datasets = _ega_parse_study_datasets(accession)
        if not datasets:
            print(f"  WARNING: No scRNA-seq datasets found. Trying to list authorized datasets...")
            result = subprocess.run(
                ["pyega3", "-cf", cred_path, "datasets"],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                for line in result.stdout.strip().split("\n"):
                    line = line.strip()
                    if re.match(r"EGAD\d+", line):
                        datasets.append(line.split()[0])
                print(f"  Found {len(datasets)} authorized dataset(s)")
            else:
                print(f"  pyega3 error: {result.stderr.strip()}", file=sys.stderr)

        if not datasets:
            print("  No datasets found or authorized", file=sys.stderr)
            return downloaded

        # For each dataset, list and download files
        for ds_id in datasets:
            print(f"\n  Dataset: {ds_id}")

            # List files in dataset
            result = subprocess.run(
                ["pyega3", "-cf", cred_path, "files", ds_id],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                print(f"    pyega3 error: {result.stderr.strip()}", file=sys.stderr)
                continue

            files = []
            for line in result.stdout.strip().split("\n"):
                line = line.strip()
                if line and re.match(r"EGAF\d+", line):
                    parts = line.split()
                    if len(parts) >= 1:
                        files.append(parts[0])

            if not files:
                print(f"    No files found in dataset")
                continue

            print(f"    Found {len(files)} file(s)")

            # Download files with pyega3 fetch
            for file_id in files:
                print(f"    Fetching dataset {ds_id} to {output_dir}...")
                result = subprocess.run(
                    ["pyega3", "-cf", cred_path, "fetch", ds_id,
                     "--output-dir", str(output_dir)],
                    capture_output=True, text=True,
                )
                if result.returncode == 0:
                    for f in output_dir.iterdir():
                        if f.is_file() and f not in downloaded:
                            downloaded.append(f)
                else:
                    err = result.stderr.strip()
                    if "401" in err or "403" in err or "unauthorized" in err.lower():
                        print(
                            f"    Authorization failed. You may need to request access to {ds_id}\n"
                            f"    at https://ega-archive.org/datasets/{ds_id}",
                            file=sys.stderr,
                        )
                    else:
                        print(f"    pyega3 error: {err}", file=sys.stderr)
                    break  # Don't keep trying if auth fails

    finally:
        Path(cred_path).unlink(missing_ok=True)

    return downloaded


@_provider(r"^EGAD\d+", "EGA Dataset")
def download_ega_dataset(accession: str, output_dir: Path):
    """Download files from a specific EGA dataset via pyega3.

    Requires EGA credentials in ~/.ega/credential.json or env vars.
    """
    import tempfile

    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []

    creds = _find_ega_credentials()
    if not creds:
        print(
            "  ERROR: EGA credentials not found.\n"
            "  Create ~/.ega/credential.json with:\n"
            '    {"username": "your_email@example.com", "password": "your_password"}\n'
            "  Or set EGA_USERNAME and EGA_PASSWORD environment variables.",
            file=sys.stderr,
        )
        return downloaded

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as cf:
        json.dump(creds, cf)
        cred_path = cf.name

    try:
        # List files
        result = subprocess.run(
            ["pyega3", "-cf", cred_path, "files", accession],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"  pyega3 error: {result.stderr.strip()}", file=sys.stderr)
            return downloaded

        files = []
        for line in result.stdout.strip().split("\n"):
            line = line.strip()
            if line and re.match(r"EGAF\d+", line):
                parts = line.split()
                if len(parts) >= 1:
                    files.append(parts[0])

        if not files:
            print(f"  No files found in dataset {accession}")
            return downloaded

        print(f"  Found {len(files)} file(s)")

        # Download
        result = subprocess.run(
            ["pyega3", "-cf", cred_path, "fetch", accession,
             "--output-dir", str(output_dir)],
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            for f in output_dir.iterdir():
                if f.is_file() and f not in downloaded:
                    downloaded.append(f)
                    print(f"    Downloaded: {f.name}")
        else:
            print(f"  pyega3 error: {result.stderr.strip()}", file=sys.stderr)

    finally:
        Path(cred_path).unlink(missing_ok=True)

    return downloaded


# ---- Main ----

def main():
    parser = argparse.ArgumentParser(
        description="Download scRNA-seq data from public repositories (GEO, ArrayExpress, Zenodo, CELLxGENE, etc.)"
    )
    parser.add_argument("accession", help="Accession ID (e.g., GSE139129, E-MTAB-13502, \"Zenodo 7942968\", cellxgene:<uuid>)")
    parser.add_argument(
        "--output-dir", "-o", default=None,
        help="Output directory (default: raw_data/{accession}/)",
    )
    parser.add_argument(
        "--skip-download", action="store_true",
        help="Skip download, only extract existing archives",
    )
    args = parser.parse_args()

    accession = args.accession

    # Detect provider
    provider_name, download_fn = detect_provider(accession)
    print(f"=== {provider_name}: {accession} ===")

    output_dir = Path(args.output_dir) if args.output_dir else Path(
        f"raw_data/{accession}"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        downloaded = download_fn(accession, output_dir)
        if not downloaded:
            print(
                f"ERROR: Could not download any files for {accession}. "
                f"Check the accession ID and try again.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Post-download extraction
        _post_download_extract(output_dir)
    else:
        print("  Skipping download (--skip-download)")

    # Report
    files = sorted(
        f for f in output_dir.iterdir()
        if f.is_file() and not f.name.endswith(".tar")
    )
    print(f"\nDone. {len(files)} file(s) in {output_dir}:")
    for f in files:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
