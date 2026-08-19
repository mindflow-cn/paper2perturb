# Paper2Perturb

**A data curation agent for [SimuCella](https://mindflow-cn.github.io/simucella/).**

Paper2Perturb turns scientific papers and public single-cell datasets into
structured, evidence-linked perturbation metadata and standardized h5ad files.
It packages the workflow as reusable agent skills for Codex and Claude Code.

[中文文档](README_CN.md)

## What it does

```text
paper PDF
  -> searchable Markdown and figures
  -> perturbation metadata and paper evidence
  -> public scRNA-seq data download
  -> cell-type annotation
  -> standardized h5ad files
  -> metadata and h5ad validation
```

Paper2Perturb currently focuses on human single-cell small-molecule
perturbation studies. The extraction workflow records eligible experiments in
`result.xlsx` and routes unsupported experimental designs to
`result_excluded.xlsx`.

## Skills

| Skill | Responsibility |
|---|---|
| `extract-perturbation` | Extract study, drug, cell, condition, target-gene, direction, and evidence metadata. |
| `build-h5ad` | Download public scRNA-seq data and build standardized h5ad and JSON outputs. |
| `validate-metadata` | Run field/schema checks and verify extracted claims against the paper. |
| `validate-h5ad` | Check consistency across h5ad, `test_case.json`, and spreadsheet metadata. |

`validate-metadata` intentionally combines deterministic validation and
paper-evidence review. These are two stages of one metadata audit rather than
separate skills.

## Quick start

Requirements: Python 3.10 or newer, Codex or Claude Code, a MinerU API key for
PDF conversion, and a Qwen API key for figure understanding.

```bash
git clone https://github.com/mindflow-cn/paper2perturb.git
cd paper2perturb
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
```

Configure MinerU and Qwen in `.env`:

```bash
MINERU_API_KEY=your-mineru-jwt-token
Qwen_API_KEY=your-dashscope-api-key
```

`MINERU_API_KEY` is used for automatic PDF conversion. `Qwen_API_KEY` is used
for figure understanding during perturbation extraction. You may alternatively
export `DASHSCOPE_API_KEY`; it takes precedence over `Qwen_API_KEY`.

Install the skills into a working project:

```bash
./scripts/install-skills.sh /path/to/working-project
```

The installer creates links in both `.agents/skills/` and `.claude/skills/`.
You can then invoke a skill explicitly, for example:

```text
Use $extract-perturbation on papers/example.pdf.
Use $build-h5ad with result.xlsx and GSE139129.
Use $validate-metadata to audit PMID 34591417 against the local paper.
```

Agent clients may also trigger skills automatically from natural-language
requests.

Read each skill's `SKILL.md` for its complete workflow, inputs, and outputs.

## Repository layout

```text
paper2perturb/
├── skills/                  # Agent skills and bundled implementations
├── scripts/                 # Project maintenance and installation tools
├── README.md
├── README_CN.md
├── requirements.txt
└── .env.example
```

Runtime outputs such as papers, downloaded matrices, h5ad files, spreadsheets,
logs, and API credentials are intentionally excluded from version control.

## Contributing

Keep each skill self-contained. Put operational instructions in `SKILL.md`,
reusable code in `scripts/`, and detailed domain rules in `references/`. Run the
validation commands in `AGENTS.md` before submitting a change. The basic
project check is:

```bash
python3 scripts/validate_project.py
```

## License

Apache License 2.0. See [LICENSE](LICENSE).
