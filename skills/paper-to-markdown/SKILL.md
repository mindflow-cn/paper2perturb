---
name: paper-to-markdown
description: "Convert one or more scientific-paper PDFs to Markdown and extracted figures with the MinerU API. Use when a Paper2Perturb workflow starts from PDF files or needs searchable paper text and figure assets."
---

# Paper to Markdown

Convert PDF files to Markdown using the MinerU API.

## Usage

```
/paper-to-markdown <path>
```

`<path>`: a PDF file or a directory containing PDF files.

- **File mode**: converts that single PDF into `papers/<pmid_or_stem>/`.
- **Directory mode**: converts every `.pdf` file under that directory into its own directory under `papers/`.

## Prerequisites

This skill requires a valid MinerU API key. Two ways to configure:

**Option 1 — `.env` file (recommended):** Add your key to the project's `.env` file:

```
MINERU_API_KEY=your-jwt-token-here
```

**Option 2 — environment variable:**

```bash
export MINERU_API_KEY="your-jwt-token-here"
```

To obtain a key, register at https://mineru.net and get your API token from the dashboard.

## How to get an API key

1. Go to https://mineru.net and register an account.
2. Log in and navigate to the API / Token management page.
3. Copy your JWT token.
4. Either add it to `.env` or export it as shown above.

## Workflow

### Step 0: Validate API key

Run the validation check:

```bash
python3 "$SKILL_DIR/scripts/paper_to_markdown.py" validate
```

If the key is missing or invalid, the script prints instructions on how to configure it. **Stop here** — do not proceed until the key is configured.

### Step 1: Convert PDF(s)

```bash
python3 "$SKILL_DIR/scripts/paper_to_markdown.py" convert <path>
```

The script:
1. Scans the path for PDF files.
2. Requests batch upload URLs from MinerU.
3. Uploads each PDF.
4. Polls for extraction completion.
5. Downloads and extracts results as `<stem>.md` files.

Output Markdown and image files are placed under `papers/<pmid_or_stem>/`.

### Step 2: Report

Print a summary of:
- How many PDFs were processed
- Where the .md files were written
- Any failures and their reasons

## Notes

- Large PDFs may take several minutes to process. The script polls every 10 seconds and will time out after 30 minutes.
- Only `.pdf` files (case-insensitive) are processed.
- Images extracted from PDFs are placed in an `images/` subdirectory next to each Markdown file.
