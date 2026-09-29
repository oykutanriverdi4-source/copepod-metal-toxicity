# Copepod Metal Toxicity

Reproducibility repository for the MSc thesis:

**Comparative Acute Metal Toxicity in Marine and Estuarine Copepods: A Quantitative Assessment of Relative LC50 Patterns Across Biological and Experimental Contexts**

This repository contains the cleaned analytical workflow used to construct the final dataset, evaluate evidence support, fit the retained statistical models, perform within-Reference and robustness analyses, assess matched-metal covariance, and regenerate the thesis-facing figures.

## Scope

The workflow focuses on mortality-based acute LC50 records for marine and estuarine copepods. The final quantitative comparison uses concentrations harmonized to micromoles of target metal per litre where the concentration basis could be resolved.

The public repository is limited to analyses retained in the final thesis. Developmental exploratory work that was not retained in the thesis, including PCA and simulation work, is intentionally not part of the public workflow.

## Core workflow

The reproducible analytical sequence is:

1. prepare the ECOTOX dataset;
2. integrate source-verified supplementary literature records;
3. audit support for the planned comparisons;
4. fit the pooled 96-h Cu-Cd-Zn model;
5. compare metals within References;
6. analyse repeated exposure durations;
7. assess salinity and temperature evidence;
8. assess pH evidence and confounding;
9. compare developmental stages within References;
10. audit support for the developmental-stage model;
11. fit the developmental-stage model;
12. audit support for the taxonomic-order model;
13. fit the taxonomic-order model;
14. analyse matched-metal covariance and correlation;
15. audit sex-specific evidence;
16. regenerate the thesis-facing figure suite.

The master runner is:

```r
source("scripts/00_run_core_pipeline.R")
```

Open `copepod-metal-toxicity.Rproj` first so that all relative paths resolve from the repository root.

## Frozen numerical benchmarks

A successful core run reproduces the following workflow checks:

- 353 Results in the combined ECOTOX + supplementary-evidence dataset;
- 304 Results with harmonized positive molar LC50 values;
- 131 Results in the pooled 96-h Cu-Cd-Zn model;
- 69 Results in the 96-h Cu-Cd Adult/Nauplii model;
- 101 Results in the 96-h Cu-Cd Calanoida/Harpacticoida model.

These checks are implemented in the individual scripts and again in `scripts/00_run_core_pipeline.R`.

## Software

The final thesis workflow was run in **R 4.4.2** on Windows 11 x64. Principal modelling and inference packages include `lme4`, `lmerTest`, `emmeans`, and `pbkrtest`; figures use `ggplot2`/`tidyverse`. Data import and workflow utilities include `readxl`, `readr`, `dplyr`, `tidyr`, `stringr`, `tibble`, `purrr`, `writexl`, and `openxlsx`.

Each major analysis writes a `sessionInfo.txt` file so that the package environment used for that run is recorded with the outputs.

## Repository structure

```text
copepod-metal-toxicity/
├── README.md
├── .gitignore
├── copepod-metal-toxicity.Rproj
├── data/
│   ├── raw/
│   │   └── ecotox/
│   ├── curated/
│   │   ├── ecotox/
│   │   ├── wos/
│   │   └── source_verification/
│   ├── processed/
│   └── external/
│       ├── wos/
│       └── envirotox/
├── scripts/
│   ├── 00_run_core_pipeline.R
│   ├── analysis/
│   ├── figures/
│   └── database_audits/
│       ├── wos/
│       └── envirotox/
├── outputs/
├── results/
│   └── figures/
└── documentation/
```

`outputs/` contains computational and QA outputs. `results/figures/` contains the reader-facing thesis figures regenerated from the frozen analysis outputs.

## Data provenance

The main analytical data lineage is:

```text
original ECOTOX export
        ↓
ECOTOX screening + harmonization
        ↓
source-verified supplementary literature records
        ↓
combined dataset
        ↓
analysis-specific subsets
        ↓
models, within-Reference comparisons, robustness and figures
```

See `documentation/data_provenance.md` for the role of each retained data file.

## External database audits

The core analysis pipeline does **not** require raw Web of Science or EnviroTox exports. The WoS export is kept local because it comes from a licensed bibliographic database. EnviroTox is publicly accessible, but this repository conservatively keeps the downloaded raw export local and instead provides the exact search settings, audit scripts, and compact audit summaries.

Optional audit scripts are provided under:

```text
scripts/database_audits/wos/
scripts/database_audits/envirotox/
```

To reproduce those audits, obtain your own export from the relevant database and place it under the documented `data/external/` path. The exact search settings and screening flows are documented in:

- `documentation/wos_search.md`
- `documentation/envirotox_search.md`

NORMAN and CAFE were supplementary coverage checks in the thesis and did not contribute separate quantitative records to the combined dataset.

## Figures

Reader-facing figures are generated with:

```r
source("scripts/figures/01_make_main_thesis_figures.R")
```

PNG and PDF files are written to:

```text
results/figures/
```

Plot data, figure objects, provenance and QA information are retained under:

```text
outputs/16_make_main_thesis_figures/
```

The figure script does not refit the scientific models; it reconstructs the final figures from frozen analytical outputs.

## Reproducibility boundary

The scripts preserve the final thesis workflow, including the retained eligibility and harmonization decisions, model formulas, planned contrasts, Reference and Species random effects, LORO/LOSO checks, environmental complete-case sensitivity, within-Reference matching rules, and Reference-level covariance bootstrap.

For Web of Science, users reproducing the audit should obtain their own export through their licensed access. For EnviroTox, users can reproduce the documented public search and export the corresponding records; the raw downloaded workbook is kept outside this repository as a conservative redistribution choice.

## Citation and archive

The archived thesis reproducibility release is available on Zenodo:

**Tanrıverdi, Ö. (2026). _Comparative Acute Metal Toxicity in Marine and Estuarine Copepods: Reproducibility Code and Data_ (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23031360**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23031360.svg)](https://doi.org/10.5281/zenodo.23031360)

- **Archived release:** `v1.0.0`
- **Version-specific DOI:** `10.5281/zenodo.23031360`
- **Zenodo record:** https://zenodo.org/records/23031360
- **GitHub repository:** https://github.com/oykutanriverdi4-source/copepod-metal-toxicity

The Zenodo record is the frozen archival copy of release `v1.0.0`. The GitHub `main` branch may contain later documentation or metadata updates.
