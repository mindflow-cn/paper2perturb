---
name: extract-perturbation
description: "Extract structured single-cell small-molecule perturbation metadata and gene-direction evidence from a scientific paper into result.xlsx. Use with a PDF or a Markdown paper directory when curating Paper2Perturb datasets; combine deterministic API helpers with semantic paper and figure review."
---

# Extract Perturbation

Extract benchmark and case data from a paper (PDF or Markdown) and write to `result.xlsx`.

## Usage

```
/extract-perturbation <path>
```

`<path>`: a PDF file **or** a case directory.

- **PDF mode**: If `<path>` is a `.pdf` file, first convert it to Markdown using `paper-to-markdown`, then proceed with `./papers/<pmid_or_stem>/`.
- **Directory mode**: If `<path>` is a directory, proceeds directly with extraction. The directory must contain a `.md` file and `images/` subdirectory.

Example: `./paper_md/test_case/case1` (contains `1.md` + `images/`)

## Test case structure

```
./paper_md/test_case/
├── case1/
│   ├── 1.md              # Paper markdown with image references
│   └── images/           # Figure images (*.jpg, *.png)
│       ├── ffd81ee9...jpg
│       └── ...
├── case2/
│   └── ...
```

## Workflow

### PDF input preprocessing

If the input `<path>` is a `.pdf` file (not a directory), first convert it to Markdown using paper-to-markdown:

1. Run paper-to-markdown:
   ```bash
   python3 skills/paper-to-markdown/scripts/paper_to_markdown.py convert <path>
   ```
2. The script outputs to `./papers/<pmid>/` (or `./papers/<stem>/` if PMID is not found). This directory contains the Markdown file and `images/`.
3. Use this output directory as `<case_dir>` for all subsequent steps.
4. **If paper-to-markdown fails** (API key missing, network error, etc.), report the error and **stop**.

If `<path>` is already a directory, use it directly as `<case_dir>` and proceed below.

### Step 0: Eligibility check (MANDATORY — do NOT skip)

**Before any extraction work, first determine whether this paper is eligible.** Read the full Markdown file and assess against the two mandatory criteria below. Report your decision to the user before proceeding.

#### Criterion 1: scRNA-seq drug perturbation data with public accession

The paper must:
- Contain **scRNA-seq data** (not just bulk RNA-seq, ATAC-seq, or proteomics alone)
- Have a **drug perturbation** experimental design (the scRNA-seq must include drug-treated conditions, not just untreated samples)
- Provide a **publicly accessible data accession** for the scRNA-seq data. The `dataset_accession` value MUST use one of the following formats so that `/build-h5ad` download.py can auto-detect the provider:

  | Repository | Format | Example |
  |-----------|--------|---------|
  | **GEO** | `GSE\d+` | `GSE139129` |
  | **ArrayExpress** | `E-MTAB-\d+`, `E-GEOD-\d+`, `E-MEXP-\d+` | `E-MTAB-13502` |
  | **Zenodo** | `Zenodo \d+` or `10.5281/zenodo.\d+` (DOI) | `Zenodo 7942968` |
  | **CZ CELLxGENE Discover** | `cellxgene:<collection_uuid>` | `cellxgene:c2879de0-affc` |
  | **GitHub** | `GitHub: owner/repo` | `GitHub: federicogiorgi/panobinostat` |
  | **ENA** (PRJEB), **DDBJ**, **Biosino** (OEP), **figshare**, etc. | Use repository prefix: `figshare: \d+` | _not yet supported by download.py_ |

  Repositories marked as unsupported should still be recorded in `dataset_accession` with a clear repository prefix, but they cannot be processed by `/build-h5ad`.

**If NO:** Report: "This paper lacks scRNA-seq drug perturbation data with a publicly accessible accession. Not eligible for benchmark extraction." **Stop here.**

#### Criterion 2: Gene targets with explicit direction in paper text

The paper must contain at least one gene/protein target where the **up/down regulation direction is explicitly stated** in the paper text (Results section or figure captions), and the data supporting the claim was produced by **wet-lab experiments performed in this paper**.

Wet-lab data includes (but is not limited to):
- **scRNA-seq**, **bulk RNA-seq**, **PRO-seq** — sequencing-based assays performed by the authors
- **Flow cytometry** for specific protein markers
- **Western blot**, **smFISH**, **Immunofluorescence (IF)**, **RT-qPCR**
- **Knockdown / overexpression** of individual genes followed by functional assays
- **ChIP-seq**, **ATAC-seq** — epigenomic assays performed by the authors

**What does NOT count:**
- Purely computational predictions, meta-analyses, or public database mining without the paper's own wet-lab data
- Pathway enrichment results (GSEA, KEGG, GO, Hallmark gene sets) — these describe programs, not specific gene targets
- Genes mentioned only in passing (e.g. in the introduction/discussion) without data in this paper
- Pooled genome-wide CRISPR screens that only report hit lists without naming individual gene-level direction

**How to check:** Scan the Results section and figure captions for explicit direction statements about specific genes. Example: "Gene X was upregulated after treatment (Fig. Y)" or "we observed decreased expression of Gene Z". **Count the gene targets with explicit direction.**

**If NO targets with explicit direction found:** Report: "This paper lacks clearly stated gene-level targets with direction. Not eligible for benchmark extraction." **Stop here.**

**If YES:** Report the number and identity of targets found, then proceed to check Criterion 3.

#### Criterion 3: Human species only — non-human data is discarded

The benchmark framework only accepts **human** scRNA-seq data. Data from other species (mouse, rat, zebrafish, etc.) is not eligible and should be **discarded entirely** (not saved to any output file).

**How to check:**
- Look at the paper's Methods or Results section for species information
- Check the `taxon` field from the `geo` API helper — it should indicate "Homo sapiens" or "human"
- Check the GEO/ArrayExpress record's organism field
- The paper may state the species explicitly (e.g. "mouse model", "murine", "rat")

**If non-human (mouse, rat, zebrafish, etc.):** Report: "This paper's data is from [species], not human. Not eligible for benchmark extraction." **Stop here.** Do not proceed to extraction or data writing.

**If human:** Proceed to Step 1.

#### Eligibility report format

Before proceeding, output a brief report:

```
## Eligibility Check: [Paper Title]

**Criterion 1 — scRNA-seq drug perturbation data:** [PASS/FAIL]
- Data accessions found: [list accession numbers and repositories, or "none"]
- Drug perturbation: [description or "N/A"]

**Criterion 2 — Gene targets with explicit direction:** [PASS/FAIL]
- Targets with explicit direction: [list specific genes and source, or "none"]

**Criterion 3 — Human species:** [PASS/FAIL]
- Species detected: [Homo sapiens / Mus musculus / etc.]

**Decision:** [PROCEED → result.xlsx / PROCEED → result_excluded.xlsx / STOP]
[Explain routing decision]
```

**Important:** If ANY of the three criteria fails, you MUST stop. Do not proceed to API calls, figure extraction, or any data writing. Report the decision to the user and end the skill execution.

- **Criterion 1 or 2 fails:** Paper lacks required scRNA-seq data or clearly stated gene targets — STOP.
- **Criterion 3 fails (non-human species):** Data is from non-human species — STOP. Animal data is discarded entirely and NOT saved to either result.xlsx or result_excluded.xlsx.

### Data routing: result.xlsx vs result_excluded.xlsx

Some extracted data does not belong in the main benchmark and should be written to `result_excluded.xlsx` instead of `result.xlsx`. Both files share the same schema (same column headers), so the write helper works identically for either target — just pass `result_excluded.xlsx` as the `result_path`.

The following types of data should be routed to `result_excluded.xlsx`:

#### Filter 1: Data not measured by scRNA-seq (bulk, RT-qPCR, etc.)

If the paper's supporting data for a drug condition is **not scRNA-seq** — including bulk RNA-seq, microarray, bulk proteomics, or targeted assays like RT-qPCR — route it to the excluded Paper2Perturb output.

- **How to detect**:
  - The paper describes "bulk RNA-seq", "RNA-seq of tissue", "microarray analysis", or the GEO platform is a bulk array (e.g. Affymetrix, Illumina BeadChip).
  - A drug treatment condition is **only validated by RT-qPCR** (or western blot, ELISA, etc.) and was never included in the scRNA-seq experiment. Check: does the paper show UMAP/clustering/DEG analysis for this condition? If the condition only appears in bar charts of qPCR data, it likely lacks scRNA-seq.
  - **Per-drug check (CRITICAL)**: Before creating a benchmark for a drug, verify that the paper actually performed scRNA-seq on cells treated with that specific drug. Do not assume all drugs mentioned in the paper were sequenced. Check the Methods section ("Sample preparation for scRNA-seq") and the Results section for scRNA-seq analysis of each drug condition.
- **Action**: Write to `result_excluded.xlsx`. Note the measurement method in `notes` (e.g. "bulk RT-qPCR only, no scRNA-seq for this condition").

#### Filter 2: Non-small-molecule drugs

If the perturbation agent is **not a small molecule** (e.g. monoclonal antibodies, peptides, proteins, siRNA, cytokines, growth factors, cell therapy), it falls outside the small-molecule drug perturbation benchmark scope.

- **How to detect**:
  - The SMILES lookup via the `smiles` helper returns `""` for biologics
  - Drug names ending in "-mab" (monoclonal antibodies), "-cept" (receptor-Fc fusions), or other biologic naming patterns
  - The paper explicitly describes the agent as "antibody", "peptide", "recombinant protein", "siRNA", "shRNA", etc.
  - Examples: nivolumab, pembrolizumab, trastuzumab, etanercept, bevacizumab
- **Action**: Write to `result_excluded.xlsx`. Note the agent class in `notes`.

#### Filter 3: Multiple drugs simultaneously (combination therapy)

If the paper studies **multiple small-molecule drugs applied simultaneously** to the same cell line (combination therapy rather than single-drug perturbation), the benchmark framework cannot isolate individual drug effects.

- **How to detect**:
  - The paper describes "combination treatment", "co-treatment", "dual therapy", or "combined with"
  - Multiple drug names appear in the same treatment group label (e.g. "DrugA+DrugB")
  - The perturbation design compares DrugA+DrugB vs control, without single-drug arms
- **Action**: Write to `result_excluded.xlsx`. Set `perturbation_scope: "combination"` and `perturbation_name` to list all drugs.
- **Exception**: If the paper includes **both** single-drug arms and combination arms, extract the single-drug arms to `result.xlsx` and the combination arms to `result_excluded.xlsx`.

#### Filter 4: Drug resistance experiments (perturb_var = resistance_status)

If the primary comparison is between **resistant vs sensitive** cell populations (not treated vs untreated), the experimental design measures resistance status rather than drug perturbation response.

- **How to detect**:
  - The paper compares "resistant" vs "parental" or "sensitive" cell lines
  - The experimental groups are labeled "R" (resistant), "S" (sensitive), "parental", "WT"
  - No drug treatment is applied in the scRNA-seq experiment — cells are profiled at baseline to compare intrinsic resistance states
  - `perturb_var` would be `resistance_status`
- **Action**: Write to `result_excluded.xlsx`. The resistance comparison is a fundamentally different experimental design from drug perturbation.
- **Exception**: If the paper profiles **both** drug-treated cells and untreated resistant/sensitive comparisons, extract the drug-treated arms to `result.xlsx` and the resistance-status comparisons to `result_excluded.xlsx`.

#### Filter 5: Non-human species (from Criterion 3)

Non-human data (mouse, rat, zebrafish, etc.) is **discarded entirely** — it is NOT saved to result_excluded.xlsx or any other output file. See Step 0 Criterion 3.

#### Routing summary

| Condition | Destination |
|-----------|-------------|
| scRNA-seq + small molecule + single drug + human + treated vs control | `result.xlsx` |
| Bulk data, RT-qPCR-only conditions (no scRNA-seq for that drug) | `result_excluded.xlsx` |
| Non-small-molecule agent (mAb, peptide, siRNA, etc.) | `result_excluded.xlsx` |
| Combination therapy (multiple drugs simultaneously) | `result_excluded.xlsx` |
| Resistance experiment (perturb_var=resistance_status) | `result_excluded.xlsx` |
| Non-human species | **Discard** (do not save) |

**When a single paper contains both eligible and ineligible data**, extract both sets separately: eligible benchmarks/cases to `result.xlsx`, ineligible ones to `result_excluded.xlsx`.

### Step 1: Read the paper and extract identifiers

Use the Read tool to read the **full** Markdown file. From the text, extract:
- **PMID** (look for citations like `[PubMed: 38987605]`, DOI, or paper title for API lookup)
- **All data accessions** — scan the Data Availability section for public repository IDs. **Format each accession for `/build-h5ad` compatibility:**
  - GEO: `GSE\d+` (e.g. `GSE139129`)
  - ArrayExpress: `E-MTAB-\d+`, `E-GEOD-\d+`, `E-MEXP-\d+`
  - Zenodo: `Zenodo \d+` or `10.5281/zenodo.\d+`
  - CELLxGENE: `cellxgene:<uuid>`
  - GitHub: `GitHub: owner/repo`
  - Other repositories (figshare, Biosino OEP, ENA PRJEB, etc.): use `RepositoryName: accession` format
- **Drug names** mentioned in the study
- **Data availability section** — note which accessions map to which experiments
- **Cell counts from figure captions** — e.g. Extended Data Fig. 6e-i gives per-condition cell N values

### Step 1a: Extract all figure/table captions to file (MANDATORY)

**Before any case-level gene extraction, extract every figure and table caption from the paper markdown into a standalone file.** This file becomes the single source of truth for all downstream caption lookups — no ad-hoc grep needed later.

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" extract-captions ./papers/<case_dir>/<file>.md ./papers/<case_dir>/captions.md
```

This produces `captions.md` containing only the figure and table captions, each prefixed with its marker line for easy reference (e.g. `### Figure 1. ...`, `### Table S3. ...`). The extraction identifies captions by their marker patterns (`Figure \d+`, `Table S?\d+`, `Supplemental Figure \d+`, `Extended Data Fig\. \d+`) and captures the full caption text until the next marker, image reference, or section boundary.

**If the helper's extraction misses captions or the paper uses non-standard naming**, fall back to manually copying the captions into `captions.md` by reading the markdown and collecting all figure/table caption blocks.

During Step 3 case extraction, **read `captions.md` directly** whenever you need to cross-reference a figure or table caption — no more ad-hoc grep.

### Step 1b: Read figure images for gene-level extraction

**This step is critical for case-level target_genes accuracy.** The paper text typically only names a few representative genes. The actual gene lists for each cell line's modules are visible only in the figures.

#### Optimized workflow: overview first, then targeted extraction

Paper images directories often contain 100+ sub-panel images, most of which are irrelevant (confluency curves, cell cycle plots, experimental diagrams, etc.). Blindly reading sub-panels one by one wastes significant time. Instead:

**Phase 1 — Run `figure_overview` on all candidate images (run in parallel):**

Figures in the markdown typically have sub-panel images placed **above** the caption, arranged in sequential panel order (a, b, c, ...). The combined multi-panel figure (if present) is usually the last or one of the images. However, not all figures have a combined image — some only have individual sub-panel images.

**Critical: images after the caption may belong to the next figure.** The previous figure's last sub-panel image is the last one before its caption. The image after the caption often starts the next figure.

**Step 1a — Collect candidate images (scan ABOVE the caption):**

For each figure referenced in `source_location`, find ALL images between this figure's caption and the **previous** figure's caption (or the previous figure's last content). These images, placed above the caption in panel order, are the figure's sub-panels. Also check the first image after the caption — it may be a combined figure for this figure, or it may belong to the next figure (verify with overview).

**Step 1b — Run `figure_overview` on ALL candidates in parallel:**

```bash
SKILL_DIR=skills/extract-perturbation

# Run overview on every candidate image — combined or individual — IN PARALLEL
python3 "$SKILL_DIR/scripts/figure_extract.py" <candidate_image_1> --mode figure_overview &
python3 "$SKILL_DIR/scripts/figure_extract.py" <candidate_image_2> --mode figure_overview &
python3 "$SKILL_DIR/scripts/figure_extract.py" <candidate_image_3> --mode figure_overview &
# ... for ALL candidate images
```

The overview output tells you, for each image:
- Whether it is a combined figure or a single panel
- What the panel shows (UMAP, dotplot, heatmap, violin plot, confluency curve, etc.)
- Which cell line(s) and condition(s) it corresponds to

**Step 1c — Classify each image:**

From the overview results, classify each image:
- **DEG dotplot/heatmap with gene labels** → candidate for `--mode genes`
- **Violin/bar plot with per-condition cell counts** → candidate for `--mode cell_counts`
- **UMAP, confluency curve, experimental diagram, etc.** → skip
- **Belongs to a different figure** → ignore (reassign to the correct figure)

**Step 1d — Cross-check caption against discovered panels:**

After classifying all candidate images, verify that every cell line mentioned in the caption has a corresponding DEG panel among the discovered images. If a cell line's dot plot was not found among the candidates, broaden the search to nearby images outside the initial candidate range before concluding it's missing.

**Phase 2 — Targeted gene extraction (only on confirmed DEG panels):**

Only run `--mode genes` on images that Phase 1 classified as DEG dotplots/heatmaps. Skip images that are:
- Confluency / growth curves
- Cell cycle phase plots
- Experimental design / schematic diagrams
- UMAP/tSNE without gene annotations
- Western blots / microscopy images

```bash
# Only run on panels the overview confirmed as having gene info
python3 "$SKILL_DIR/scripts/figure_extract.py" <dotplot_subpanel_image> --mode genes --json
```

The `--mode genes` output includes **direction and conditions**, not just gene names.

**Critical: LLM must still cross-reference with paper text.** Qwen-VL may misread some gene symbols (e.g. PAX3→PAX8, WIT1→WT1) and condition labels. The LLM corrects these against the paper's figure caption and known gene symbols before writing to xlsx. Gene symbols that cannot be confidently corrected should be omitted rather than guessed.

**When figures are NOT dotplots/heatmaps with readable gene labels (violin plots, pathway bar charts, enrichment plots):** `figure_extract.py --mode genes` cannot extract gene names from these figure types. In this case, the **figure caption is the primary supplementary gene source**. The body text typically names only 2-3 representative genes (often with "e.g."), while the figure caption names additional genes shown in individual panels. **Always read the relevant caption from `captions.md` for gene symbols and merge with the body text list.** See "e.g. is non-exhaustive" rule in Step 3.

**Phase 3 — Cell counts from violin plots:**

For cell count estimates, only read the violin/bar plot images that the overview identified:

```bash
python3 "$SKILL_DIR/scripts/figure_extract.py" <violin_subpanel_image> --mode cell_counts --json
```

#### Key rules

- **Overview ALL candidate images, not just combined figures.** If a combined figure image is missing, the individual sub-panel images near the caption may still be present. Run `figure_overview` on every image near the relevant figure caption — combined or individual — to discover which DEG panels are available.
- **Never assume a combined figure exists.** Some markdowns only contain scattered individual sub-panel images. The workflow must handle both formats.
- **Run overviews in parallel.** If the paper has 3-4 figures with DEG data, launch all overview `figure_extract.py` calls simultaneously via background bash jobs.
- **Never read the same image twice.** Track which image files you've already processed.
- Target genes MUST come from the specific figure for the specific cell line. Do not copy gene lists across cell lines even if they share a module name. Each figure shows a different subset of genes for that cell line.
- For cell counts, use `--mode cell_counts` on violin plot images identified by the overview. Sum per-condition N values when the paper provides them.
- **Do not assume all sub-panels mentioned in the caption exist as images.** Markdowns may only include a subset of figure panels. When a needed panel image is genuinely missing, acknowledge the limitation rather than silently substituting text-inferred gene lists.

### Step 2: Call Python helpers for deterministic data

Run these via Bash (all independent, call in parallel where possible):

The helper script is at `skills/extract-perturbation/scripts/extract_helpers.py`. Use the `SKILL_DIR` variable for portability:

```bash
SKILL_DIR=skills/extract-perturbation
python3 "$SKILL_DIR/scripts/extract_helpers.py" pubmed "{PMID}"
```
→ Returns JSON: `{"pmid", "doi", "title", "year", "journal"}`

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" geo "{GSE}"
```
→ Run once per GEO accession. Returns JSON: `{"accession", "title", "platform", "n_samples", "bioproject", "taxon", "subseries"}`. **Only works for GEO (GSE) accessions.** For other repositories (ArrayExpress, ENA, Biosino, figshare, Zenodo, etc.), extract metadata manually from the paper's Data Availability section and the repository page.

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" smiles "{DRUG_NAME}"
```
→ Run once per drug. Returns SMILES string or empty string for non-small-molecule agents.

### Step 3: LLM semantic extraction

Using the full paper text and the API results from Step 2, extract structured data as JSON.

#### Benchmark granularity

One benchmark = one unique combination of **(cell type, drug, dataset accession)**.

- **Split by cell type within the same sample.** When a paper's scRNA-seq data identifies multiple distinct cell types or subtypes (e.g. different neuronal subtypes, different immune cell populations, different epithelial subclusters), and the paper reports cell-type-specific drug response gene expression changes, create a **separate benchmark for each cell type + drug combination**. Do NOT merge all cell types from one sample into a single benchmark — different cell types have different baseline expression profiles and different drug responses, and they should be evaluated independently. Each benchmark's `cell_context` and `description` should specify the cell type, not just the broader cell line or tissue.
- **Do NOT over-split by sequencing platform.** If the same cell type + drug combination spans two dataset accessions that used different platforms (e.g. inDrop vs 10x), merge them into one benchmark using the primary accession as `dataset_accession`. Only split when the accessions genuinely represent independent experiments with different biological conditions.
- **PDX / in vivo models**: include only when the paper provides dedicated scRNA-seq data with clear drug-perturbation design and its own dataset accession.
- **secondary_accession**: prefer the BioProject ID (PRJNAxxxxxx) when available (from GEO API or ENA). For non-GEO repositories (Biosino, figshare, Zenodo, etc.), use the project-level or study-level accession if one exists, or leave empty if not applicable. Do not use the superseries (GSExxxxxx) as secondary_accession.
- **Response-stratified splitting**: When the paper reports opposite drug effects for different patient subgroups (e.g. "gene X decreased in responders, increased in non-responders"), do NOT pool all patients into one benchmark. Create separate benchmarks per subgroup (e.g. `_R` for responders, `_NR` for non-responders), each with its own cases and correct expected direction. Verify the deposited data contains subgroup labels before splitting. See Claim Testability Check 4.

#### Benchmark field guidance

**dose_design**: `single` means each experimental condition uses one dose level, even in a dose-escalation design where different adapted populations were generated at different doses. `multi` means the same cells were simultaneously exposed to multiple doses in a single experiment (e.g. a dose-response plate).

**cell_context**: Use a concise human-readable description that includes the **specific cell type**, without specific mutation details. Examples: `human <cell_type_name> from <disease_context> cell line` (not just `human <disease_context> cell line`), `human <mutation> <cancer_type> cell line`.

**disease**: Use the generic disease name. Examples: `<generic_disease_name>` (not `<mutation>-mutant <generic_disease_name>`), `<mutation> <cancer_type>`.

**paper_title**: Remove trailing period from PubMed title.

**description**: One-sentence summary that names the **specific cell type**. Example: `<cell_type> from <cell_line> treated with <drug>.</`

**cell_type_original**: The original cell type name used in the paper. If `source_type` is `cell_line`, leave this column empty.

**cell_type_standard**: Standardized cell type name using the CellTypist nomenclature system. If `source_type` is `cell_line`, leave this column empty.

**cell_type_markers**: Marker genes selected based on tissue, dataset context, and cell_type information, used to identify or validate the cell type annotation. May include positive and negative markers, e.g. `CD3E+, CD8A+, CD4-`. If `source_type` is `cell_line`, or the benchmark does not involve cell type annotation, leave empty.

#### Case extraction: wet-lab gene targets with explicit direction

**Extract cases for gene targets whose up/down direction is explicitly stated in the paper text, supported by the paper's own wet-lab data.** The benchmark framework requires that the paper's own experimental data (scRNA-seq, RNA-seq, flow cytometry, Western blot, etc.) demonstrates the claimed direction for the named genes.

**What qualifies as a case target:**
- Genes whose **up/down direction is explicitly stated** in the paper's Results or figure captions (e.g. "Gene X was upregulated", "we observed decreased expression of Gene Y")
- The data supporting the direction claim comes from **wet-lab experiments performed in this paper** — scRNA-seq, bulk RNA-seq, PRO-seq, ChIP-seq, flow cytometry, Western blot, smFISH, RT-qPCR, etc.
- Genes shown in figure panels (heatmaps, dot plots, volcano plots, feature plots) where the paper explicitly names them and states their direction

**What does NOT qualify (do NOT extract as cases):**
- Pathway enrichment results (GSEA, hallmark gene sets, KEGG, GO terms) — these describe broad biological programs, not specific gene targets
- Genes only listed in supplementary tables without explicit direction in the text
- Purely computational predictions or public database analyses without the paper's own wet-lab data
- Genes the paper mentions in passing (e.g. in introduction/discussion) without showing own data

**How to identify targets:** Scan the Results section and figure captions for genes named with explicit direction. The best targets appear in figures (Fig. 5C heatmaps, volcano plots, dot plots) where the paper visually shows expression changes and names the genes in the caption or surrounding text. Use `figure_extract.py --mode genes` to extract gene lists from figures, then cross-reference with the paper text to confirm direction.

**How to ensure complete cell type coverage:** After identifying all drug-response gene targets, scan the Results text for phrases indicating that an effect spans multiple cell types: "in both X and Y," "common to X," "observed across X and Y," "in X and Y neurons," etc. For every such phrase, create parallel cases under EVERY cell type's benchmark — not just the first one named. The same `original_statement` and gene list apply to all mentioned cell types that share the same drug condition and direction.

#### MANDATORY — Caption cross-reference before finalizing target_genes (CRITICAL)

**For every case, before writing to xlsx, you MUST execute the three-source gene audit below.** This is not optional — it is a hard gate. You may not call `write-xlsx` or `write_to_xlsx` until every case row has passed this audit.

The body text, figure captions, and figure_extract output are three complementary gene sources, each with different coverage:
- **Body text**: names 2–3 representative genes per finding (often with "e.g.") — **incomplete by design**
- **Figure captions**: the authoritative written record of every gene shown in each panel — **always read**
- **Figure_extract output** (overview / --mode genes): visual gene detection from the image itself — **may catch genes the caption describes without naming**

All three sources must be cross-referenced before finalizing `target_genes`. Relying on only one or two of the three is the primary cause of missing genes.

**The captions file already exists at `./papers/<case_dir>/captions.md`** because Step 1a generated it. Use it directly.

##### Three-source gene audit (execute per case row — output required)

**Step A — Collect genes from all three sources:**

For each case row, produce an explicit accounting:

```
Case: [benchmark] / [cell type] / [relation]

Source 1 — Body text genes:
  (list every gene symbol found in the body text passages that support this case)
  e.g. GENE_A, GENE_B, GENE_C

Source 2 — Caption genes:
  (read captions from captions.md for EVERY figure/table referenced in source_location)
  e.g. Fig. X caption names: GENE_A, GENE_B, GENE_D
  (gene symbols follow standard nomenclature: uppercase letters + digits, optional hyphens;
   also check for protein names/aliases that map to gene symbols)

Source 3 — Figure_extract genes (if figure_extract.py was run on relevant panels):
  (list every gene detected by figure_overview or --mode genes for the panels relevant to this case)
  e.g. Fig. X panel overview detected: GENE_A, GENE_B, GENE_E
```

**Step B — Compute the union and resolve differences:**

```
Union of all three sources: GENE_A, GENE_B, GENE_C, GENE_D, GENE_E

Genes already in target_genes: GENE_A, GENE_B, GENE_C

NEW genes to add (with rationale):
  - GENE_D: named in Fig. X caption, same direction/conditions → ADD
  - GENE_E: detected by figure_overview on Fig. X panel, same direction/conditions,
    symbol verified against known gene nomenclature → ADD

Genes NOT added (with reason):
  - (any gene from sources 2/3 that is ambiguous or has unclear direction)
```

**Step C — Each gene found in source 2 or 3 that is NOT in source 1 MUST be accounted for.**

- **Source 2 genes (caption-named) are NON-NEGOTIABLE.** Every gene explicitly named in the figure/table caption MUST appear in `target_genes`. You may NOT exclude a caption-named gene with "gene set size," "kept focused," or any other justification. The caption is the authors' own written record of what the figure shows — if they named it, it belongs in the case. The only valid reason to exclude a caption-named gene is if it was measured under different experimental conditions (different drug/dose/cell type) and does not share the case's merge key.
- **Source 3 genes (figure_extract) must be accounted for.** Either add to `target_genes` (after verifying the symbol against known gene nomenclature and confirming same direction/conditions) or document the specific reason for exclusion in `notes`. Valid exclusion reasons: symbol could not be confidently verified, belongs to a different condition group, or figure_extract clearly misread the label.
- No gene from any source may be silently dropped.

**Step D — Append audit summary to case `notes`:**

```
"Caption cross-check: [N] genes added from caption ([list]), [M] genes confirmed from figure_extract. [K] genes excluded ([list]) with reason."
```

If no new genes were found: `"Caption cross-check: all caption/figure_extract genes already in target_genes — no additions."`

##### Key rules

- **"e.g." in body text is a red flag.** When the body text uses "e.g.", "such as", or "including" before a gene list, the list is guaranteed incomplete. The figure caption and figure_extract output MUST be consulted to recover the full gene set. A caption cross-check that finds zero new genes for an "e.g." list is almost certainly incomplete — re-examine the caption more carefully. **Every gene named in the caption that matches the case's conditions MUST be added to `target_genes` regardless of how many body-text genes are already present.**
- **Every figure/table in `source_location` must have its caption read.** If `source_location` references "Figures 3B, 3D, 4C", you must read the captions for Figure 3 AND Figure 4 from `captions.md`.
- **Figure_extract output is not redundant with captions.** Qwen-VL may detect a gene symbol in a violin plot title that the caption describes only as "cholesterol genes" without naming individually. Conversely, captions may name genes that Qwen-VL misreads. The two sources complement each other — use both.
- **Protein names ≠ gene symbols.** Captions and figure_extract may report protein names (e.g. "CHOP") that differ from the official gene symbol (e.g. "DDIT3"). Always map to the official gene symbol before adding to `target_genes`.
- **This audit applies regardless of whether `figure_extract.py` was used.** For cases where no images were processed, sources 2 and 3 collapse into one (caption only), but the audit must still be performed.

**Non-coding RNA targets:** Always check for lncRNAs and other non-coding RNAs. These are often discussed in separate sentences or paragraphs and are easy to overlook.

**Gene sets for multi-gene modules:** When the paper groups multiple genes together (e.g. stress response module, lineage markers), and they share the same direction under the same conditions, merge them into a single case row with `target_genes` as a JSON array. Do NOT create separate rows per gene when the merge key (benchmark + relation + conditions) is identical.

**Pathway-level enrichment is NOT a case:** EMT, cholesterol metabolism, drug metabolism, cell cycle — these are biological programs identified by enrichment analysis. They are not specific, testable predictions. Do NOT create `gene_set` cases for broad pathway enrichment results.

**Gene set size:** Keep gene sets focused (typically 4-10 genes). A 20+ gene set from a broad pathway definition is not a specific target — it's a pathway-level observation. Narrow to the core genes that are actually shown in the specific figure. **This rule never justifies removing a gene that is explicitly named in the figure/table caption. Caption-named genes always stay in `target_genes` regardless of total count.**

**Beware biological priors:** Do not apply general biological heuristics (e.g. "drug tolerance → dedifferentiation → lineage markers down") without verifying the direction in the paper's actual data. Different cell lines and perturbagens produce different cell-type-specific responses. Always cross-check the direction against the specific figure for the specific cell line, not against what is typical in the field.

#### Claim testability checks (MANDATORY — before creating each case)

Before creating a case, verify the claim is actually testable in the pseudobulk per-cell mean expression framework. The following claim types are NOT testable as simple per-cell expression direction cases. Check each case against all three.

##### Check 1: Proportion/frequency claims — NOT testable

If the paper describes changes in terms of **"proportion of X+ cells"**, **"frequency of Y+ cells"**, or **"percentage of Z+ cells"**, the claim is about cell composition, not per-cell mean expression. Pseudobulk per-cell expression cannot test proportion changes — these are orthogonal dimensions.

- **How to detect**: Look for "proportion", "frequency", "percentage", "%", "fraction of cells expressing", "ratio of X+ to X- cells" in the Results text and figure captions. Bar charts showing "% positive cells" or pie charts are proportion claims.
- **Example**: Paper claims "proportion of <GENE_A>+ and <GENE_B>+ cells increased after treatment" (Fig. X). This is a proportion change, not a per-cell expression change.
- **Action**: Do NOT create cases for proportion/frequency claims. Document skipped claims in `notes` with reason: "proportion claim, not testable in pseudobulk per-cell expression framework".

##### Check 2: Evidence modality mismatch

If the claim's supporting evidence comes from a **different measurement modality** than the scRNA-seq data (e.g. bulk microarray, bulk RT-qPCR, Western blot of total tissue lysate), verify that the scRNA-seq data captures the relevant cell population.

- **How to detect**: For each case, trace the evidence in `source_location` back to the paper. Ask: (a) What measurement method produced this evidence? (b) What cell population does that method sample? (c) Does the scRNA-seq data capture the same cell population?
- **Common mismatch patterns**:
  - Bulk tissue microarray/RT-qPCR where the gene is expressed in cell types NOT captured by scRNA-seq (e.g. keratinocytes when scRNA-seq only captured immune cells)
  - Whole-tissue Western blot where the protein is primarily in stromal cells but scRNA-seq is tumor-enriched
  - Serum/plasma biomarkers where the measurement is circulating protein, not cellular mRNA
- **Example**: A gene's evidence comes from bulk tissue microarray/RT-PCR where the gene is mainly expressed in a cell type NOT captured by the scRNA-seq data (e.g. stromal cells when scRNA-seq only captured immune cells). The evidence modality and data modality measure different populations.
- **Action**: If the evidence modality cannot be reproduced in the scRNA-seq data due to cell population mismatch, do NOT create the case. If there is partial overlap but uncertainty, flag in `notes`: "evidence from [modality], may not be captured in scRNA-seq [cell population]".

##### Check 3: Response-stratified claims — requires benchmark splitting

If the paper reports **different directions for different patient subgroups** (responders vs non-responders, sensitive vs resistant, etc.), the claim direction depends on the subgroup. Pooling all patients into one benchmark mixes opposite signals.

- **How to detect**: Look for phrases like "decreased in X and increased in Y", "opposite trends in X vs Y", "stratified by response", "responders showed... while non-responders showed...", "in X patients... but in Y patients...". Check if figures have separate panels or faceted views for different patient subgroups.
- **Example**: Paper states "<GENE_A>, <GENE_B>, and <GENE_C> decreased in responders and increased in non-responders following treatment" (Fig. X). Direction depends on response status — pooled analysis shows whichever subgroup dominates.
- **Action**: Do NOT create a single benchmark pooling all patients. Instead, create **separate benchmarks for each subgroup** (e.g. `<PMID>_01` for responders, `<PMID>_02` for non-responders). Each gets its own cases with the correct expected direction per subgroup. Verify that the data contains the subgroup labels needed to split — check barcode metadata or sample annotations.
- **When subgroup labels are missing from data**: If the paper defines subgroups (R/NR) but the deposited data lacks the subgroup annotations needed to split cells, flag in `notes`: "data lacks [subgroup] annotations needed to split benchmark — requires manual curation of clinical metadata".

##### Testability report in notes

For every case, the `notes` field should briefly confirm testability or flag the specific issue. Examples:
- `"expression claim, testable in pseudobulk"` (normal case)
- `"proportion claim — skipped; paper reports %<GENE>+ cells not per-cell expression"`
- `"evidence from bulk RT-qPCR (Fig. 3E); scRNA-seq captures same cell population — testable"`
- `"response-stratified: R and NR have opposite directions — benchmark split into _R and _NR"`

#### Case granularity and coverage

- **Every benchmark must have at least one case.** If a benchmark has no extractable cases, reconsider whether it should be a benchmark.
- **One case row = one unique (benchmark + conditions + relation + original_text_context) combination.** The case granularity is per cell type + drug + experimental conditions, NOT per individual gene. Multiple genes that share the same `relation`, `dose_groups`, and `time_groups` within the same benchmark are merged into a single case row with `target_genes` as a JSON array of all genes. For example, if GENE_A, GENE_B, and GENE_C are all DOWN under the same dose conditions, they go in one row: `target_genes: ["GENE_A","GENE_B","GENE_C"]`.
- **Singleton genes are merged when conditions match**: genes validated individually that share the same direction and conditions go into the same case row. Do NOT create separate rows per gene when they share identical `perturb_var`, `control`, `relation`, `dose_groups`, and `time_groups`.
- **Different relations or conditions = separate rows**: if two genes have different directions (UP vs DOWN) or are compared under different conditions (different `dose_groups` or `time_groups`), they remain in separate case rows. The merge key is `(benchmark_id, perturb_var, control, relation, dose_groups, time_groups)`.
- **Original text proximity determines merge vs split (CRITICAL)**: Even when genes share the same merge key, check how they are discussed in the original paper. If the paper mentions their up/down regulation **together in the same sentence or same sentence cluster** (consecutive sentences about the same finding), they represent a single biological conclusion → **merge into one case**. If the two genes are discussed in **separate, distant parts of the paper** (different paragraphs, different sections, different figure discussions) with no textual connection, they represent **different biological conclusions** → **split into separate cases**, one per biological context. The `original_statement` field will naturally differ between merged vs split cases: merged cases quote the single passage covering all genes together; split cases each quote their respective distant passages.
- **Gene sets are for multi-gene modules**: use `target_type: "gene_set"` when the paper groups multiple genes into a named module/program that are validated together or shown as a coherent set in a specific figure. When independently validated singleton genes are merged due to shared conditions, use `target_type: "gene"` with `gene_set_source: "singleton_target"` and `target_id` concatenating the gene names (e.g. `GENE_A_GENE_B_GENE_C`).
- **Stress response modules must always be separate cases** — even when the paper discusses two stress modules in the same sentence or same figure (e.g. a UPR/ER-stress module and an antioxidant/redox module). Do NOT merge them into one combined case. Each module has a distinct `target_id` and `gene_set_source`.
- **Only create stress module cases when explicitly shown in the figure**: if a supporting cell line's figure only shows lineage/surface markers without stress module genes, do NOT create stress module cases for that cell line. Do not assume all cell lines show all module types.
- For cell lines appearing in multiple benchmarks, create parallel cases for each benchmark — do not only create cases for the "main" cell line.
- **Cross-model validation coverage**: for supporting benchmarks (non-primary cell lines), create the full set of parallel cases. No minimum case limit — some benchmarks naturally have only 1 case.
- **"Common to multiple cell types" requires parallel cases in ALL mentioned cell types (CRITICAL)**: When the paper describes a drug response gene set as occurring "in both cell type A and cell type B," "common to cell type X," or any similar phrasing that indicates the same direction under the same drug in multiple cell types, you MUST create parallel cases under EVERY mentioned cell type's benchmark. Do NOT only create cases for the first-mentioned or "primary" cell type. The original_statement quoting "in both DAn1 and DAn2" belongs to both benchmarks — include it in both. After drafting all cases, run a completeness check: for each drug, list which cell types have cases, and verify this matches the set of cell types the paper discusses for that drug.
- **EMT handling**: the EMT/mesenchymal phenotype is captured as a SINGLE case per benchmark — not two. For the primary benchmark, use `target_type: "cell_state"` with `target_id: "EMT_score"` and include the actual gene list in `target_genes`. For supporting benchmarks, use `target_type: "gene_set"` with `target_id: "paper_defined_EMT_stemness_markers"`. Never create a separate `cell_state` EMT_score case with a generic placeholder like `["Hallmark_EMT_gene_set"]`.

#### Case field guidance

**relation**: Default to `UP` / `DOWN`. Only use `SERIES_UP` / `SERIES_DOWN` when the paper **explicitly and quantitatively** demonstrates strict monotonic progression across dose or time (e.g. "expression progressively increases with each dose step"). In practice, `SERIES_*` is rarely appropriate — when in doubt, use `UP` / `DOWN`.

**comparison**: Describe the actual comparison design. Example: `dose-escalation trend relative to C`. Do not use generic labels like `treated_vs_control`.

**gene_set_source**: Cite the specific figure or data source. Examples: `singleton_target`, `Extended Data Fig. 4f paper-defined marker set`, `Fig. 2a,b paper-defined module A`.

**quantitative_support_type**: Use specific types as a JSON list. Valid types: `expression_comparison`, `cluster_DE_support`, `module_score`, `DE_analysis`, `fold_change`, `p_value`, `smFISH_validation`, `protein_validation`, `bulk_validation`, `knockdown_validation`, `reporter_validation`. Include both computational evidence types AND experimental validation types when the paper provides them.

**target_genes**: Always a JSON array of gene symbols. Contains ALL genes that share the same direction and experimental conditions within this case row. `["GENE1"]` for a single gene, `["GENE1","GENE2","GENE3"]` for multiple genes with the same relation and conditions. **Critical**: extract the gene list from the **specific figure referenced in `source_location`**, not from a global module definition elsewhere in the paper. The same module may have different visible genes in different cell lines/figures. **MANDATORY**: before finalizing, run the caption cross-reference procedure (see Step 3) — read the figure/table captions from `captions.md` for each source referenced in `source_location` to find additional gene symbols not present in the body text. The body text + figure captions together define the complete gene set for the case.

**target_id**: For single-gene cases, use the gene symbol (e.g. `GENE1`). For multi-gene cases where genes share the same direction and conditions, concatenate gene names with underscores (e.g. `GENE_A_GENE_B_GENE_C`). For gene sets/modules, use a descriptive name prefixed with `paper_defined_` (e.g. `paper_defined_stress_module`, `paper_defined_lineage_markers`). For cell states, use a score name (e.g. `EMT_score`).

**cell_type**: Include the cell line name and disease context. Example: `<CELL_LINE> <disease_context> cell line`, not just `<CELL_LINE>`.

**dose_groups**: **Always write actual dose values regardless of `perturb_var`.** Even when `perturb_var` is `time` or `resistance_status`, if the experiment involves drug treatment at specific doses, record those doses here. List only the dose groups actually compared in this case, not all doses in the experiment. **Format as dose value with unit** (e.g. `["1.0 uM", "2.0 uM", "4.0 uM"]`), not paper labels (e.g. `["T1", "T2", "T4"]`). **Use a narrow dose range** — select the 2-3 highest doses where the expression change is clearest, not the full dose ladder. If the experiment genuinely has no dose component (e.g. purely untreated comparison), use `[]`.

**original_statement**: 原封不动地将论文原话复制过来，不要同义转述或者做任何删改。Copy the paper's exact wording from the figure caption or results text character by character — no paraphrasing, no rewriting, no deletion or modification of any kind. This is NOT a field for synthesized summaries or templated sentences with gene names swapped in. Use this field to anchor the case to specific paper evidence. When multiple genes are merged into one case row, include the quote(s) that cover all genes, separated by newlines if from different parts of the paper.

**time_groups**: For non-time-response cases, use `"[]"` (empty JSON array string). For time-course experiments, each case should focus on the specific time point(s) where the expression change is most pronounced. **All time_groups values must be expressed in day units.** Convert hours to days (e.g., 24h → `"1 day"`, 72h → `"3 day"`), weeks to days (e.g., 2 weeks → `"14 day"`). A gene that peaks at day 2 should have `time_groups: ["2 day"]`, not `["D1","D2","D4","D9","D11"]`. Use `time_mode: "single"` with a narrow time window. Only use `time_mode: "multi"` when the paper explicitly demonstrates sustained change across all listed time points with no temporal specificity.

**paper_time_index**: For non-time-response cases, use `"[]"`.

#### Output JSON structure

```json
{
  "benchmarks": [
    {
      "dataset_accession": "GSEXXXXX",
      "secondary_accession": "PRJNAXXXXXX",
      "species": "human",
      "cell_type_original": "",
      "cell_type_standard": "",
      "cell_type_markers": "",
      "tissue": "<tissue>",
      "source_type": "cell_line",
      "cell_context": "human <disease_context> cell line",
      "disease": "<disease_name>",
      "platform": "<platform>",
      "perturbation_type": "drug",
      "perturbation_scope": "single",
      "perturbation_name": "<drug_name>",
      "description": "<one-sentence summary>",
      "control_type": "untreated",
      "dose_design": "single",
      "time_design": "multi"
    }
  ],
  "cases": [
    {
      "benchmark_index": 0,
      "source_location": "<Figure reference>",
      "original_statement": "<verbatim quote from paper>",
      "drug": "<drug_name>",
      "cell_type": "<cell_line> <disease> cell line",
      "target_genes": ["<GENE1>", "<GENE2>", "<GENE3>"],
      "target_type": "gene",
      "target_id": "<GENE1>_<GENE2>_<GENE3>",
      "perturb_var": "dose",
      "control": "C",
      "dose_groups": ["<dose1>", "<dose2>", "<dose3>"],
      "time_groups": [],
      "relation": "DOWN",
      "comparison": "dose-escalation trend relative to C",
      "gene_set_source": "singleton_target",
      "has_quantitative_support": true,
      "quantitative_support_type": ["expression_comparison", "cluster_DE_support"],
      "quantitative_support_detail": "<figure> shows <description of change>",
      "response_timescale": "<description of time/dose scale>",
      "experiment_design": "<summary of experimental design>",
      "is_dose_response": true,
      "paper_dose_index": ["<label1>", "<label2>", "<label3>"],
      "is_time_response": false,
      "time_mode": "single",
      "paper_time_index": [],
      "notes": ""
    }
  ]
}
```

`benchmark_index` in cases refers to the 0-based index into the `benchmarks` array — used in Step 4 to link benchmark_id.

### Step 4: Generate IDs, merge, and write to xlsx

**`write-xlsx` auto-creates the file if it doesn't exist and appends rows.** No manual init step is needed.

For each benchmark, generate a benchmark_id:

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" gen-benchmark-id "{PMID}" result.xlsx
```

Then for each case, generate a test_id:

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" gen-test-id "{BENCHMARK_ID}" result.xlsx
```

Merge all data into final rows:
- Benchmark rows = pubmed metadata + smiles + geo metadata + LLM semantic fields + generated benchmark_id
- Case rows = LLM semantic fields + benchmark_id (from linked benchmark) + test_id + pmid + dataset_accession

Write to xlsx (appends to existing data if file exists, auto-creates if not). **Choose the target file based on the data routing rules above:**

```bash
# Eligible data → result.xlsx
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result.xlsx '{JSON_WITH_BENCHMARKS_AND_CASES}'

# Ineligible data → result_excluded.xlsx
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result_excluded.xlsx '{JSON_WITH_BENCHMARKS_AND_CASES}'
```

If a single paper contains both eligible and ineligible data, call `write-xlsx` **twice** — once for each target file with the appropriate subset of benchmarks and cases.

The JSON argument must be: `{"benchmarks": [...], "cases": [...]}` where each dict's keys match the column headers. Both files share the same schema.

**Important**: When calling write-xlsx, pipe the JSON via stdin or use a temp file if the JSON is too long for a command-line argument:

```bash
python3 -c "
import json, sys
sys.path.insert(0, 'skills/extract-perturbation/scripts')
from extract_helpers import write_to_xlsx
data = json.load(open('/tmp/extract_data.json'))
write_to_xlsx('result.xlsx', data['benchmarks'], data['cases'])
"
```

For `result_excluded.xlsx`, change the path accordingly. The `write_to_xlsx` function works identically for any path.

### Step 5: Confirm

Print a summary:
- Paper title and PMID
- Number of benchmarks added to `result.xlsx`
- Number of cases added to `result.xlsx`
- Number of benchmarks written to `result_excluded.xlsx` (if any)
- Number of cases written to `result_excluded.xlsx` (if any)
- List of benchmark_ids created (with destination file)
- For data routed to `result_excluded.xlsx`, note which filter(s) triggered the routing

## Field Schema Reference

### Benchmark-level fields

| Field | Type | Source | Description |
|-------|------|--------|-------------|
| benchmark_id | string | **generated** | `<PMID>_<Index:02d>` via `gen-benchmark-id` (Index starts at 01 per PMID) |
| pmid | string | **paper** | 来源文献 PMID |
| paper_doi | string | **pubmed API** | 来源文献 DOI |
| paper_title | string | **pubmed API** | 论文标题 |
| year | int | **pubmed API** | 发表年份 |
| journal | string | **pubmed API** | 期刊名 |
| dataset_accession | string | **paper** | 主数据编号，必须使用 `/build-h5ad` 可识别的格式：GEO→`GSE\d+`；ArrayExpress→`E-MTAB-\d+` 等；Zenodo→`Zenodo \d+`；CELLxGENE→`cellxgene:<uuid>`；GitHub→`GitHub: owner/repo`；其他→`仓库名: 编号`。优先用 subseries 或子数据集 |
| secondary_accession | string | **paper / API** | BioProject ID（PRJNAxxxxxx）或其他项目级编号；非 GEO 仓库可从论文或仓库页面获取；不用 superseries |
| species | enum | **LLM** | human / mouse / other |
| cell_type_original | string | **LLM** | 论文中的细胞类型原始命名（若source_type是cell_line，则无需填写此列） |
| cell_type_standard | string | **LLM** | 细胞类型标准化命名，使用CellTypist命名体系（若source_type是cell_line，则无需填写此列） |
| cell_type_markers | string | **LLM** | 结合 tissue、dataset context 和 cell_type 信息选择的 marker genes，用于识别或验证该细胞类型注释。可包含阳性和阴性 marker，例如 CD3E+, CD8A+, CD4-。若 source_type 为 cell_line，或该 benchmark 不涉及细胞类型注释，则无需填写。 |
| tissue | string | **LLM** | 组织来源，如 Skin、Lung、Ovary |
| source_type | string | **LLM** | cell_line / organoid / primary_culture / patient_sample / PDX / co_culture_model |
| cell_context | string | **LLM** | 简洁的样本背景，如 "human <mutation> <cancer_type> cell line" |
| disease | string | **LLM** | 疾病或模型背景 |
| platform | string | **LLM** / paper | 技术平台，如 10x Chromium、inDrop |
| perturbation_type | string | **LLM** | 扰动类型，如 drug |
| perturbation_scope | enum | **LLM** | single / combination / both |
| perturbation_name | string | **LLM** | 具体扰动名称 |
| smiles | string | **smiles API** | 药物 SMILES（非小分子填空字符串） |
| description | string | **LLM** | 一句话 benchmark 描述 |
| control_type | string | **LLM** | 对照类型，如 untreated、DMSO、vehicle |
| dose_design | enum | **LLM** | single（每个条件一个剂量）/ multi（同时多剂量） |
| time_design | enum | **LLM** | single / multi |

### Case-level fields

| Field | Type | Source | Description |
|-------|------|--------|-------------|
| test_id | string | **generated** | `<benchmark_id>_<Test_Index:02d>` via `gen-test-id` or `assign-ids` (Test_Index starts at 01 per benchmark). |
| benchmark_id | string | **generated** | 关联 benchmark 主键 |
| pmid | string | **paper** | 来源文献 PMID |
| dataset_accession | string | **from benchmark** | 数据集编号 |
| source_location | string | **LLM** | 原文图号或段落位置，如 "Extended Data Fig. 4f" |
| original_statement | string | **LLM** | 原文原话，一字不改地复制。Must be a verbatim quote, not a paraphrase or synthesized sentence. |
| drug | string | **LLM** | 药物或处理名称 |
| cell_type | string | **LLM** | 细胞类型含疾病背景，如 "<cell_line> <mutation> <cancer_type> cell line" |
| target_genes | list[string] | **LLM** | 基因列表，JSON 数组格式。Contains ALL genes sharing the same direction and conditions. Single: `["GENE1"]`; multi: `["GENE1","GENE2","GENE3"]`. |
| target_type | enum | **LLM** | gene / gene_set / pathway / cell_state |
| target_id | string | **LLM** | 单基因用基因名；多基因用基因名拼接（如 `GENE_A_GENE_B_GENE_C`）；gene_set 用 paper_defined_ 前缀 |
| perturb_var | enum | **LLM** | 主要比较维度。**仅允许三个值：`dose` / `time` / `resistance_status`**。禁止填 `treatment`、`drug`、`condition` 等其他值。多剂量实验用 `dose`；治疗前后比较用 `time`；耐药/敏感状态比较用 `resistance_status` |
| control | string | **LLM** | 对照组定义，如 "C"、"untreated" |
| dose_groups | list[string] | **LLM** | 剂量组标签，JSON 格式，如 ["1.0 uM","2.0 uM","4.0 uM"] |
| time_groups | list[string] | **LLM** | 时间组，必须以天(day)为单位，如 ["1 day"]、["14 day"] |
| relation | string | **LLM** | UP / DOWN / SERIES_UP / SERIES_DOWN |
| comparison | string | **LLM** | 描述实际比较设计，如 "dose-escalation trend relative to C" |
| gene_set_source | string | **LLM** | singleton_target / 具体图表引用（如 "Extended Data Fig. 4f paper-defined marker set"） |
| has_quantitative_support | bool | **LLM** | 是否有定量支持 |
| quantitative_support_type | list[string] | **LLM** | JSON 格式，如 ["expression_comparison","cluster_DE_support"] |
| quantitative_support_detail | string | **LLM** | 定量支持的具体说明 |
| response_timescale | string | **LLM** | 原文中的时间尺度描述 |
| experiment_design | string | **LLM** | 该 case 对应的实验设计摘要 |
| is_dose_response | bool | **LLM** | 是否为多剂量 case |
| paper_dose_index | list[string] | **LLM** | 论文原始剂量标签，JSON 格式 |
| is_time_response | bool | **LLM** | 是否为时间趋势 case |
| time_mode | enum | **LLM** | single / multi |
| paper_time_index | list[string] | **LLM** | 原始时间标签，JSON 格式 |
| notes | string | **LLM** | 备注。Must include a brief testability confirmation or flag specific issues (proportion claim, modality mismatch, response stratification). See Claim Testability Checks. |

## Common Pitfalls

### Not checking data routing criteria (CRITICAL)
After extraction, always check the 5 filter conditions before writing to xlsx. It's easy to overlook that a paper uses bulk RNA-seq, or that the "drug" is actually a monoclonal antibody. The most common routing mistakes:
- Extracting bulk RNA-seq data without noticing it's not scRNA-seq (check GEO platform field)
- Treating monoclonal antibodies (names ending in -mab) as small molecules — SMILES lookup returns empty string
- Writing combination therapy data to result.xlsx when the paper only studied DrugA+DrugB, not single drugs
- Including resistance-status comparisons (resistant vs parental lines profiled without drug treatment) as drug perturbation cases
- Writing non-human species data to result.xlsx or result_excluded.xlsx (animal data should be discarded)
- **Including drug conditions that lack scRNA-seq data** — a drug may be mentioned in the paper but only measured by RT-qPCR or western blot, never by scRNA-seq. Verify each drug condition has scRNA-seq analysis (UMAP/clustering/DEG) before creating a benchmark.

### Not splitting benchmarks by distinct cell subtypes (CRITICAL)
When a paper's scRNA-seq data contains multiple distinct cell types (e.g. different neuronal subtypes, immune cell populations, epithelial subclusters), and the paper reports **cell-type-specific drug responses**, create a separate benchmark for each cell type + drug combination. Different cell types have different baseline expression profiles and respond differently to the same drug — pooling them into one benchmark dilutes the signal. The paper will typically report DEGs separately per cell type, and cases reference different gene lists per cell type. If the paper says "only 17% of DEGs were common to both cell types," separate benchmarks are mandatory.

### Creating cases for only one cell type when the paper describes effects in multiple cell types (CRITICAL)
When the paper uses phrases like "in both cell type A and B," "common to cell type X," or "observed across cell types Y and Z," you MUST create parallel cases under EVERY mentioned cell type's benchmark. Typical failure mode: the paper says "upregulation of stress response genes in both DAn1 and DAn2 neurons," and the extractor creates cases only for DAn1. The DAn2 benchmark ends up with zero cases despite the paper explicitly stating the same effect occurs there. **After drafting all cases, run a per-drug completeness check**: for each drug, enumerate which cell types have cases and which cell types the paper discusses for that drug. Any cell type discussed but missing cases is an error.

### Wrong case-to-benchmark mapping when multiple benchmarks share the same drug (CRITICAL)
When benchmarks are split by cell type, multiple benchmarks may share the same `perturbation_name` (e.g. two rotenone benchmarks for different cell subtypes). **Always use the two-phase write workflow** (see "ID generation" above): write benchmarks first, assign IDs, read back the assigned `benchmark_id` for each benchmark, then write cases with the correct `benchmark_id` directly. This eliminates the need to match cases by drug name — the LLM already knows which case belongs to which benchmark from the extraction step.

**Automated safeguards:**
- Use the two-phase write workflow — NOT `propagate-ids` — to assign benchmark_ids to cases.
- Always run `validate-mapping` after ID assignment as a final safety check.
- `propagate-ids` is a recovery tool for legacy data; do not use it in normal extraction.

### Skipping the eligibility check (CRITICAL)
The most costly error — spending time on full API lookups, figure extraction, and semantic analysis only to discover at the end that the paper has no scRNA-seq data or no clearly stated gene targets. **Always run Step 0 first.** If either criterion fails, stop immediately and report to the user. Do not proceed to API calls or figure extraction "just to see what's there."

### Extracting genes without explicit paper direction (CRITICAL)
Only extract genes whose up/down direction is explicitly stated in the paper text. Do NOT infer direction from general biological knowledge or assume that all genes in a figure follow the same trend. The paper must name the gene and state its direction (e.g. "X was upregulated") in the Results or figure caption.

### Applying biological priors over paper data
Do NOT assume expression direction based on general biological heuristics. What is "typical" for a given drug-cell line combination in the literature may not hold in the specific paper's data. Always verify direction against the paper's actual figures and text for the specific cell line and condition, not against your knowledge of the field. Cell-type-specific responses can diverge significantly from general patterns.

### Missing non-coding RNA targets
Papers often validate lncRNAs and other non-coding RNAs using the same assays as coding genes (smFISH). These are discussed in separate sentences or paragraphs and are easy to overlook. Always scan validation experiment descriptions for non-coding RNA targets.

### Extracting pathway-level enrichment as cases
GSEA results, hallmark gene sets, KEGG pathways, and GO term enrichment describe broad biological programs observed in the data — they are NOT specific gene targets. Do NOT create `gene_set` cases for "EMT pathway", "cholesterol metabolism", "drug metabolism", "cell cycle" etc. These are descriptive observations, not testable predictions.

### Gene set size inflation
Gene sets with 20+ genes are typically pathway-level observations, not focused gene targets. If a gene set is that large, verify that each gene is individually shown in the referenced figure — otherwise narrow to the core 4-10 genes. **This rule never justifies removing a gene that is explicitly named in the figure/table caption. Caption-named genes are always included regardless of total gene count.**

### relation: over-using SERIES_UP / SERIES_DOWN
Most dose-escalation papers describe trends qualitatively ("decreases along the resistance continuum"). This is NOT sufficient justification for `SERIES_DOWN`. Reserve `SERIES_*` only when every adjacent dose pair shows the same direction with statistical significance. **Default to `UP` / `DOWN`.**

### Merging stress modules into one case
These are distinct biological programs. Even when mentioned in the same figure caption sentence, **always create separate cases** with separate `target_genes` lists specific to that cell line.

### Splitting genes that share the same conditions into separate rows (CRITICAL)
The case granularity is per (cell type + drug + conditions + relation + original_text_context), NOT per gene. When multiple genes share the same `relation`, `dose_groups`, and `time_groups` within the same benchmark, merge them into a single case row with all genes in `target_genes`. Do NOT create separate rows for each gene when they have identical merge keys. This bloats the Case-level sheet and obscures the fact that the genes form a coherent set validated under the same experimental conditions. For example, if GENE_A, GENE_B, and GENE_C are all DOWN under the same dose-escalation conditions, create ONE row with `target_genes: ["GENE_A","GENE_B","GENE_C"]`, not three separate rows.

**Exception — different biological conclusions**: If genes share the same merge key but their regulation is discussed in separate, distant parts of the paper (different paragraphs/sections/figures) with no textual connection, they represent different biological conclusions and should be split into separate cases. The `original_statement` field will naturally reflect different source passages. Merging two unrelated biological findings just because they happen to share conditions obscures the distinct conclusions the paper is making.

### Paraphrasing or synthesizing original_statement (CRITICAL)
`original_statement` must be a verbatim quote copied character-by-character from the paper. Do NOT paraphrase, summarize, rewrite, or synthesize the text. Do NOT create templated sentences where only the gene name changes between rows (e.g. "XX gene expression decreased along the dose-escalation continuum" — this is synthesis, not a quote). Do NOT delete or modify any part of the original text. If the paper says "GENE_A, GENE_B, and GENE_C regulons were uniformly downregulated", copy that exact sentence — do not shorten it, restructure it, or adapt it. When the paper's original wording doesn't name every gene in the case's `target_genes`, either: (a) find and quote the passages that do name them, or (b) note in `notes` which genes are implied by the context. Never fill gaps by inventing text.

### Skipping caption cross-reference for "e.g." gene lists (CRITICAL — see mandatory three-source audit in Step 3)

**Failure mode**: Body text says "pathway X biosynthesis (e.g., <GENE_A>, <GENE_B>, <GENE_C>)" for a drug-treated condition. The figure caption for the same panel explicitly names an additional gene <GENE_D> shown in a violin plot with a log2 fold change value. The figure_overview output on that violin plot image also reports <GENE_D>. But the LLM compiles `target_genes` using only the body text "e.g." list — missing <GENE_D>.

**Why this happens**: The LLM reads the caption and runs figure_overview, but when building the final JSON it defaults to the body text gene list without executing the gene-by-gene comparison step.

**Prevention**: The three-source gene audit in Step 3 forces an explicit diff between body text, caption, and figure_extract gene lists before writing. If any gene from source 2 or 3 is not in the final `target_genes`, the audit output makes the omission visible. A caption cross-check that finds zero new genes for a body text "e.g." list is a red flag — re-examine.

When the paper body text uses "e.g.", "such as", or "including" before a gene list, those genes are **representative examples, not the full set**. The figure caption for the same panel often names additional genes that were omitted from the body text. **Every case with an "e.g." in its original_statement MUST produce at least one new gene from the caption or figure_extract, or the audit must explicitly justify why none were found.**

### Finding caption genes but excluding them (CRITICAL)

**Failure mode**: The three-source audit correctly identifies genes in the figure caption that are missing from `target_genes`. But instead of adding them, the LLM writes a note like "SQLE found in Fig 4C caption, not added to keep gene set focused" and moves on. The audit becomes a documentation exercise rather than a corrective gate.

**Why this happens**: The old Step C allowed "document why it was excluded" as an alternative to adding the gene. The LLM treats "gene set size" as a valid exclusion reason even for author-named genes.

**Prevention**: Caption-named genes are NON-NEGOTIABLE (see Step C). If a gene appears in a figure/table caption with matching conditions/direction, it goes into `target_genes`. No exception. "Gene set size" or "keeping it focused" are never valid reasons to drop a caption-named gene. The only valid exclusion for a caption-named gene is a genuine condition mismatch (different drug, dose, or cell type than the case row).

### Using global module gene lists instead of figure-specific genes (CRITICAL)
This is the most common and most damaging error. Do NOT copy-paste the same `target_genes` from a global module definition. Instead, extract the genes **visible in the specific figure referenced by `source_location`**. The same biological module (e.g. TF_A targets, TF_B targets) shows COMPLETELY DIFFERENT genes in different cell lines/figures — each dotplot/heatmap displays a distinct subset of genes for that specific cell line. **Always run `figure_extract.py` on the exact image file to get the gene list.** Never assume that genes listed in the paper text are the complete set — the text only names 2-3 representative examples per module, while the actual figure shows many more.

### Assuming combined figure images always exist (CRITICAL)
Not all figures in the markdown have a combined multi-panel image. Some figures may only have a few individual sub-panel images scattered near the caption. When you run `figure_overview` on an image near a figure caption expecting a combined figure but it turns out to be a different figure, do NOT skip the entire figure's image extraction. Instead, run `figure_overview` on ALL individual images near that caption to discover which sub-panels are actually available. A supporting cell line's DEG dot plot may exist as a standalone image even when the combined figure is missing. Missing a standalone image silently causes gene lists to be inferred from paper text instead of extracted from the figure — leading to inaccurate, cross-cell-line-identical lists.

### Missing cross-model validation cases
Supporting cell line benchmarks should have the same case structure as the primary benchmark: lineage/dedifferentiation DOWN, EMT UP, stress response modules UP. Don't reduce coverage just because the supporting data appears in supplementary figures.

### Missing cell-state cases
When the paper computes signature scores (EMT, stemness, etc.), extract a `target_type: "cell_state"` case even if the score is not the main focus.

### secondary_accession: use project-level BioProject when available
For GEO datasets, all benchmarks from the same paper share the **superseries BioProject** (from the parent GSE record's `bioproject` field). Do NOT use per-subseries BioProjects. For non-GEO repositories (Biosino, ArrayExpress, ENA, etc.), use the project-level or study-level accession if one exists. The superseries BioProject PRJNA is the correct identifier for GEO data.

### paper_title: remove trailing period
PubMed titles often end with a period. Strip it before writing to xlsx.

### perturb_var: do NOT use "treatment" (CRITICAL)
`perturb_var` only accepts three values: **`dose`**, **`time`**, or **`resistance_status`**. Never use `"treatment"`, `"drug"`, `"condition"`, or any other value. For pre-vs-post treatment comparisons where the main dimension is the time elapsed since treatment initiation, use `"time"` (with `time_groups` and `is_time_response: true`). For multi-dose experiments, use `"dose"` (with `is_dose_response: true`). For comparing resistant vs sensitive populations where the primary axis is resistance state, use `"resistance_status"`. **Note:** `dose_groups` is always written with actual dose values regardless of `perturb_var` choice.

### dose_groups: format as dose value + unit
Use dose values with units (e.g. `["1.0 uM", "2.0 uM", "4.0 uM"]`), NOT paper labels (e.g. `["T1", "T2", "T4"]`).

### disease/cell_context: generic names, no mutation details
Use generic disease names: `<generic_disease_name>` (not `<mutation>-mutant <generic_disease_name>`). For cell_context: `human <disease_context> cell line` (not `human <mutation>-mutant <disease_context> cell line`).

### time_groups: use peak time point in day units
For time-course experiments, each case should focus on the specific time point(s) where the expression change is most pronounced, not the entire time range. **All time_groups values must be expressed in day units.** Convert hours to days (e.g., 24h → `"1 day"`, 72h → `"3 day"`), weeks to days (e.g., 2 weeks → `"14 day"`). Use paper labels like "Post"/"Pre" only when the paper does not report exact time duration. A gene that peaks at day 2 should have `time_groups: ["2 day"]`, not `["D1","D2","D4","D9","D11"]`. This makes each case a specific, testable prediction. Use `time_mode: "single"` with a narrow time window. Only use `time_mode: "multi"` when the paper explicitly demonstrates sustained change across all listed time points with no temporal specificity.

### source_location: cite all relevant figures
When evidence for a case spans multiple figures (e.g. main figure for primary cell line + extended data figure for supporting cell lines), cite all relevant figures. Format: `Extended Data Fig. 3e,j and 4f,l; Supplementary Tables 2-5`. Do not cite only a single figure when the paper shows the same finding in multiple panels.

### Missing non-obvious case types
- **Immune-modulatory markers**: For certain cancer types (e.g. melanoma), the immune evasion / MHC pathway may be modulated during drug adaptation. Check for class II MHC markers and related immune-modulatory genes that may be downregulated along the adaptation continuum.
- **Surface marker switching**: For certain cancer types (e.g. lung cancer), check for surface marker switching as part of the dedifferentiation phenotype (e.g. one surface marker down, another surface marker up). These should be separate singleton gene cases (`target_type: "gene"`).

### Extracting proportion/frequency claims as expression cases (CRITICAL)
When the paper reports changes in the **proportion** or **frequency** of cells expressing a marker (e.g. "% <GENE>+ cells increased"), do NOT create a per-cell expression case. Proportion and per-cell mean expression are orthogonal dimensions — a population can have more <GENE>+ cells while each cell expresses less <GENE>. Pseudobulk per-cell expression cannot test proportion changes. Check for keywords: "proportion", "frequency", "percentage", "%", "fraction of cells". See Claim Testability Check 1.

### Extracting claims supported only by a different measurement modality (CRITICAL)
When a claim's evidence comes from bulk microarray, bulk RT-qPCR, or whole-tissue Western blot, verify that the scRNA-seq data captures the same cell population. A gene expressed in keratinocytes cannot be validated with scRNA-seq that only captured immune cells, even if both come from "skin." For each case, trace `source_location` back to the paper and check: what measurement method? What cell population? Does scRNA-seq capture it? See Claim Testability Check 2.

### Pooling response-stratified claims into a single benchmark (CRITICAL)
When the paper reports **opposite directions for different patient subgroups** (e.g. "decreased in responders, increased in non-responders"), do NOT create one benchmark pooling all patients. The pooled pseudobulk direction depends arbitrarily on which subgroup has more cells. Instead, split into separate benchmarks per subgroup (e.g. `_R` and `_NR` suffixes). Verify the data contains subgroup labels. See Claim Testability Check 3.

### Boolean types: Python True/False, not strings
Use Python `True`/`False` for boolean fields (`has_quantitative_support`, `is_dose_response`, `is_time_response`). Do NOT use strings `"TRUE"`/`"FALSE"`.

### Empty array/list fields
For non-time-response cases, `time_groups` and `paper_time_index` should be empty arrays (`[]`). Do not omit them — the xlsx column must be filled.

## ID generation (two-phase write)

When extracting multiple benchmarks from one paper, the LLM already knows which cases belong to which benchmark (via `benchmark_index` in the extraction JSON). Use a **two-phase write** that preserves this linkage directly — no guessing, no post-hoc propagation needed.

### Phase 1: Write benchmarks, assign IDs, read back

```bash
# Step 1: Write only benchmark rows (empty benchmark_id, test_id)
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result.xlsx '{"benchmarks": [...], "cases": []}'

# Step 2: Assign benchmark IDs (fills per-PMID sequential indices)
python3 "$SKILL_DIR/scripts/extract_helpers.py" assign-ids result.xlsx

# Step 3: Read back assigned IDs
python3 "$SKILL_DIR/scripts/extract_helpers.py" read-benchmark-ids result.xlsx
```

`read-benchmark-ids` outputs a JSON array of all benchmark rows with their `benchmark_id`, `perturbation_name`, and `cell_context`. The new rows are appended at the end — map them back to your extraction JSON by order: `benchmark_index` 0 → the 1st new row in the output, 1 → 2nd, etc. Each case's `benchmark_index` now resolves to a concrete `benchmark_id`.

### Phase 2: Write cases with benchmark_ids, generate test_ids

```bash
# Step 4: Write case rows with the now-known benchmark_id
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result.xlsx '{"benchmarks": [], "cases": [...]}'

# Step 5: Generate test_ids
python3 "$SKILL_DIR/scripts/extract_helpers.py" assign-ids result.xlsx
```

### Verify

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" validate-mapping result.xlsx
```

Fix any reported errors before proceeding. With two-phase writes, `validate-mapping` should typically pass on the first run.

### For `result_excluded.xlsx`

Same two-phase workflow, applied independently:

```bash
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result_excluded.xlsx '{"benchmarks": [...], "cases": []}'
python3 "$SKILL_DIR/scripts/extract_helpers.py" assign-ids result_excluded.xlsx
python3 "$SKILL_DIR/scripts/extract_helpers.py" read-benchmark-ids result_excluded.xlsx
# --- map IDs to cases ---
python3 "$SKILL_DIR/scripts/extract_helpers.py" write-xlsx result_excluded.xlsx '{"benchmarks": [], "cases": [...]}'
python3 "$SKILL_DIR/scripts/extract_helpers.py" assign-ids result_excluded.xlsx
python3 "$SKILL_DIR/scripts/extract_helpers.py" validate-mapping result_excluded.xlsx
```

### `propagate-ids` is a recovery tool only

Do NOT use `propagate-ids` in the normal extraction workflow. It exists only to repair legacy xlsx files where case-to-benchmark mapping was lost. The two-phase write above makes it unnecessary — cases are written with the correct `benchmark_id` from the start.

## Important notes

- Read the **entire** paper — data availability and methods are often at the end.
- One paper may yield multiple benchmarks (different cell types, different drugs, different dataset accessions).
- When a dataset has a superseries/project with subseries/sub-datasets, prefer using the most granular accession as `dataset_accession` for each benchmark, with the project-level ID as `secondary_accession`.
- For non-GEO repositories, ensure the data is publicly downloadable before proceeding. `dataset_accession` must use the format required by `/build-h5ad` download.py: Zenodo → `Zenodo \d+`; CELLxGENE → `cellxgene:<uuid>`; GitHub → `GitHub: owner/repo`. For repositories not yet supported by download.py (Biosino, ENA, figshare, etc.), use `RepositoryName: accession` format and note the limitation.
- List-type fields (target_genes, dose_groups, etc.) are serialized as JSON strings by the write helper.
- `write_to_xlsx` auto-creates the target file (whether `result.xlsx` or `result_excluded.xlsx`) if it doesn't exist and appends data. No manual init needed.
- Column headers are defined as `BENCHMARK_HEADERS` and `CASE_HEADERS` constants in `extract_helpers.py` — the single source of truth for the xlsx schema. Both `result.xlsx` and `result_excluded.xlsx` use identical schemas.
- The `smiles` helper queries PubChem by drug name and falls back to ChEMBL — no need for a separate search step. An empty return value is a signal that the agent may be a biologic (see Filter 2 in Data Routing).
