#!/usr/bin/env python3
"""
paper-to-markdown — Convert PDF files to Markdown using the MinerU API.

Usage:
  paper_to_markdown.py validate              Check that MINERU_API_KEY is set and valid
  paper_to_markdown.py convert <path>        Convert PDF(s) at <path> (file or directory) to Markdown

Output is always placed under ./papers/<pmid_or_stem>/
"""

import os
import re
import sys
import time
import zipfile
import shutil
import requests
from pathlib import Path

API_BASE = "https://mineru.net/api/v4"
POLL_INTERVAL = 10  # seconds
MAX_WAIT = 1800  # 30 minutes
PAPERS_DIR = "papers"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _load_dotenv():
    """Load MINERU_API_KEY from .env file if not already set in environment."""
    if os.environ.get("MINERU_API_KEY"):
        return
    d = Path.cwd()
    while True:
        env_file = d / ".env"
        if env_file.exists():
            try:
                with open(env_file, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("#") or "=" not in line:
                            continue
                        k, _, v = line.partition("=")
                        if k.strip() == "MINERU_API_KEY":
                            val = v.strip().strip('"').strip("'")
                            if val:
                                os.environ["MINERU_API_KEY"] = val
                            return
            except OSError:
                pass
        if d.parent == d:
            break
        d = d.parent


def get_api_key():
    _load_dotenv()
    return os.environ.get("MINERU_API_KEY", "")


def _headers(token):
    return {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {token}",
    }


# ---------------------------------------------------------------------------
# PMID extraction
# ---------------------------------------------------------------------------

# PMID is an 8-digit number (rarely 7)
_PMID_RE = re.compile(r"\bPMID\s*[:=-]?\s*(\d{7,8})\b", re.IGNORECASE)
_PUBMED_ID_RE = re.compile(r"PubMed\s*(?:ID)?\s*[:=-]?\s*(\d{7,8})\b", re.IGNORECASE)
_PUBMED_URL_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d{7,8})")
_PUBMED_BRACKET_RE = re.compile(r"\[PubMed\s*[:—-]?\s*(\d{7,8})\]")
_DOI_RE = re.compile(r"\b(10\.\d{4,}/[^\s\]]+)")
_PMC_RE = re.compile(r"\bPMCID\s*[:=-]?\s*(PMC\d+)\b", re.IGNORECASE)

_PMID_PATTERNS = [_PMID_RE, _PUBMED_ID_RE, _PUBMED_URL_RE, _PUBMED_BRACKET_RE]

_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def _lookup_pmid_by_doi(doi):
    """Look up PMID from DOI via NCBI E-utilities."""
    try:
        r = requests.get(
            _ESEARCH_URL,
            params={"db": "pubmed", "term": f"{doi}[doi]", "retmode": "json"},
            timeout=10,
        )
        if r.status_code == 200:
            data = r.json()
            ids = data.get("esearchresult", {}).get("idlist", [])
            if ids:
                return ids[0]
    except Exception:
        pass
    return None


def extract_pmid(md_path):
    """Try to find a PMID in the markdown file. Returns the PMID string or None."""
    try:
        with open(md_path, "r", encoding="utf-8") as f:
            content = ""
            for i, line in enumerate(f):
                if i >= 500:
                    break
                content += line
    except Exception:
        return None

    # 1) direct PMID patterns in text
    for pattern in _PMID_PATTERNS:
        m = pattern.search(content)
        if m:
            return m.group(1)

    # 2) DOI → PMID via NCBI API
    m = _DOI_RE.search(content)
    if m:
        doi = m.group(1).rstrip(".")
        pmid = _lookup_pmid_by_doi(doi)
        if pmid:
            return pmid

    return None


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

def validate():
    token = get_api_key()
    if not token:
        print("ERROR: MINERU_API_KEY is not set.")
        print()
        print("To configure:")
        print("  1. Register at https://mineru.net")
        print("  2. Get your API token from the dashboard")
        print("  3. Export it in your shell:")
        print("       export MINERU_API_KEY=\"<your-jwt-token>\"")
        print("     Or add that line to ~/.zshrc / ~/.bashrc to persist.")
        sys.exit(1)

    try:
        r = requests.get(
            f"{API_BASE}/extract-results/batch/00000000-0000-0000-0000-000000000000",
            headers=_headers(token),
            timeout=15,
        )
        if r.status_code in (401, 403):
            print("ERROR: MINERU_API_KEY is invalid or expired.")
            print()
            print("To fix:")
            print("  1. Go to https://mineru.net and log in")
            print("  2. Regenerate your API token from the dashboard")
            print("  3. Update the environment variable:")
            print("       export MINERU_API_KEY=\"<new-token>\"")
            sys.exit(1)
        print("MINERU_API_KEY is valid.")
    except requests.exceptions.ConnectionError:
        print("ERROR: Cannot reach mineru.net. Check your network.")
        sys.exit(1)
    except requests.exceptions.Timeout:
        print("ERROR: Connection to mineru.net timed out.")
        sys.exit(1)


# ---------------------------------------------------------------------------
# collect PDFs
# ---------------------------------------------------------------------------

def collect_pdfs(path):
    p = Path(path).resolve()
    if not p.exists():
        print(f"ERROR: path does not exist: {path}")
        sys.exit(1)

    if p.is_file():
        if p.suffix.lower() != ".pdf":
            print(f"ERROR: file is not a PDF: {path}")
            sys.exit(1)
        return [p]

    pdfs = sorted([f for f in p.rglob("*.pdf") if f.is_file()])
    pdfs += sorted([f for f in p.rglob("*.PDF") if f.is_file()])
    seen = set()
    uniq = []
    for f in pdfs:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    uniq.sort(key=lambda x: str(x).lower())

    if not uniq:
        print(f"ERROR: no PDF files found under: {path}")
        sys.exit(1)

    return uniq


# ---------------------------------------------------------------------------
# upload
# ---------------------------------------------------------------------------

def request_upload_urls(pdf_paths, token):
    files_data = [{"name": f.name, "data_id": f.stem} for f in pdf_paths]

    r = requests.post(
        f"{API_BASE}/file-urls/batch",
        headers=_headers(token),
        json={"files": files_data, "model_version": "vlm"},
        timeout=60,
    )
    if r.status_code in (401, 403):
        print("ERROR: API key rejected during upload request.")
        print("Run 'paper_to_markdown.py validate' to check your key.")
        sys.exit(1)
    if r.status_code != 200:
        print(f"ERROR: upload request failed (HTTP {r.status_code})")
        print(r.text)
        sys.exit(1)

    body = r.json()
    if body.get("code") != 0:
        print(f"ERROR: API returned error: {body.get('msg')}")
        sys.exit(1)

    batch_id = body["data"]["batch_id"]
    upload_urls = body["data"]["file_urls"]
    return batch_id, upload_urls


def upload_files(pdf_paths, upload_urls):
    success = 0
    failed = 0
    for i, (file_path, url) in enumerate(zip(pdf_paths, upload_urls), 1):
        print(f"  [{i}/{len(pdf_paths)}] Uploading {file_path.name} ...", end=" ")
        try:
            with open(file_path, "rb") as f:
                r = requests.put(url, data=f, timeout=300)
            if r.status_code == 200:
                print("done")
                success += 1
            else:
                print(f"FAILED (HTTP {r.status_code})")
                failed += 1
        except Exception as e:
            print(f"FAILED ({e})")
            failed += 1
    return success, failed


# ---------------------------------------------------------------------------
# poll
# ---------------------------------------------------------------------------

def poll_batch(batch_id, token):
    url = f"{API_BASE}/extract-results/batch/{batch_id}"
    deadline = time.time() + MAX_WAIT

    print(f"  Polling batch {batch_id} ...")
    while time.time() < deadline:
        r = requests.get(url, headers=_headers(token), timeout=30)
        if r.status_code != 200:
            print(f"  WARNING: poll returned HTTP {r.status_code}, retrying...")
            time.sleep(POLL_INTERVAL)
            continue

        body = r.json()
        if body.get("code") != 0:
            print(f"  WARNING: API error: {body.get('msg')}, retrying...")
            time.sleep(POLL_INTERVAL)
            continue

        results = body["data"].get("extract_result", [])
        states = {item.get("state") for item in results}

        if all(s == "done" for s in states):
            print("  All files done.")
            return results
        if "failed" in states:
            failed = [item["file_name"] for item in results if item.get("state") == "failed"]
            print(f"  WARNING: {len(failed)} file(s) failed: {failed}")

        pending = sum(1 for s in states if s != "done")
        print(f"  {pending} file(s) still processing, waiting {POLL_INTERVAL}s ...")
        time.sleep(POLL_INTERVAL)

    print("ERROR: timed out waiting for batch completion.")
    sys.exit(1)


# ---------------------------------------------------------------------------
# download & organise
# ---------------------------------------------------------------------------

def download_and_organise(extract_results, token):
    """Download zip, extract, find PMID, place in papers/<pmid_or_stem>/."""
    papers_root = Path(PAPERS_DIR)
    papers_root.mkdir(parents=True, exist_ok=True)

    success = 0
    failed = 0
    pending = 0

    for item in extract_results:
        data_id = item.get("data_id")
        file_name = item.get("file_name", data_id)
        state = item.get("state")
        zip_url = item.get("full_zip_url", "")
        err_msg = item.get("err_msg", "")

        print(f"  [{file_name}] state={state}", end="")

        if state != "done":
            print(f" — skipping ({err_msg})" if err_msg else " — skipping")
            pending += 1
            continue

        if not zip_url:
            print(" — no zip URL")
            failed += 1
            continue

        # use a staging directory under papers/
        stage_dir = papers_root / f".{data_id}_stage"
        stage_dir.mkdir(parents=True, exist_ok=True)
        zip_path = stage_dir / f"{data_id}.zip"

        try:
            # download zip
            print(" — downloading ...", end=" ")
            r = requests.get(zip_url, stream=True, timeout=300)
            if r.status_code != 200:
                print(f"FAILED (HTTP {r.status_code})")
                failed += 1
                shutil.rmtree(stage_dir, ignore_errors=True)
                continue

            with open(zip_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            # extract into stage dir
            extract_dir = stage_dir / "extracted"
            extract_dir.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(extract_dir)

            # find the main .md file
            md_files = list(extract_dir.rglob("*.md"))
            if not md_files:
                print("WARNING: no .md found in zip")
                failed += 1
                shutil.rmtree(stage_dir, ignore_errors=True)
                continue

            main_md = None
            for md in md_files:
                if md.parent == extract_dir or "auto" in md.parent.name.lower():
                    main_md = md
                    break
            if main_md is None:
                main_md = md_files[0]

            # extract PMID from the markdown content
            pmid = extract_pmid(main_md)
            name = pmid if pmid else data_id
            print(f"— PMID={pmid if pmid else 'N/A'}", end=" ")

            # final output directory: papers/<pmid_or_stem>/
            final_dir = papers_root / name
            # handle duplicate: if dir already exists, append suffix
            if final_dir.exists():
                for n in range(2, 100):
                    alt = papers_root / f"{name}_{n}"
                    if not alt.exists():
                        name = f"{name}_{n}"
                        final_dir = alt
                        break
            final_dir.mkdir(parents=True, exist_ok=True)

            # move .md → final_dir/<name>.md
            target_md = final_dir / f"{name}.md"
            shutil.move(str(main_md), str(target_md))

            # move images/
            img_src = main_md.parent / "images"
            img_dst = final_dir / "images"
            if img_src.exists() and img_src.is_dir():
                if img_dst.exists():
                    shutil.rmtree(img_dst)
                shutil.move(str(img_src), str(img_dst))

            print(f"→ papers/{name}/")
            success += 1

        except Exception as e:
            print(f"FAILED ({e})")
            failed += 1
        finally:
            shutil.rmtree(stage_dir, ignore_errors=True)

    return success, failed, pending


# ---------------------------------------------------------------------------
# convert
# ---------------------------------------------------------------------------

def convert(path):
    token = get_api_key()
    if not token:
        print("ERROR: MINERU_API_KEY is not set.")
        print("Run 'paper_to_markdown.py validate' for setup instructions.")
        sys.exit(1)

    pdf_paths = collect_pdfs(path)
    print(f"Found {len(pdf_paths)} PDF(s) to convert.")

    # Step 1 — request upload URLs
    print("Requesting upload URLs ...")
    batch_id, upload_urls = request_upload_urls(pdf_paths, token)
    print(f"  batch_id: {batch_id}")

    # Step 2 — upload
    print("Uploading files ...")
    up_ok, up_fail = upload_files(pdf_paths, upload_urls)
    if up_fail:
        print(f"WARNING: {up_fail} upload(s) failed; these will be skipped at download.")
    if up_ok == 0:
        print("ERROR: no files uploaded successfully.")
        sys.exit(1)

    # Step 3 — poll
    results = poll_batch(batch_id, token)

    # Step 4 — download, extract PMID, place in papers/
    print("Downloading results ...")
    dl_ok, dl_fail, dl_pending = download_and_organise(results, token)

    # Summary
    print()
    print("=" * 60)
    print("Summary:")
    print(f"  Total PDFs:      {len(pdf_paths)}")
    print(f"  Converted:       {dl_ok}")
    print(f"  Failed:          {dl_fail}")
    print(f"  Still pending:   {dl_pending}")
    print(f"  Output:          {PAPERS_DIR}/")
    print("=" * 60)


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    command = sys.argv[1]

    if command == "validate":
        validate()
    elif command == "convert":
        if len(sys.argv) < 3:
            print("Usage: paper_to_markdown.py convert <file-or-directory>")
            sys.exit(1)
        convert(sys.argv[2])
    else:
        print(f"Unknown command: {command}")
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
