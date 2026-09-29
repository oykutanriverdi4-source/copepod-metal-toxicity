# Web of Science search and screening workflow

## Search

**Database:** Web of Science Core Collection  
**Search date:** 5 September 2026  
**Field:** Topic

Exact query:

```text
copepod* AND ("LC50" OR "LC 50" OR "median lethal concentration") AND
(copper OR cadmium OR zinc OR nickel OR lead OR silver OR chromium OR mercury)
```

The Topic field included title, abstract, author keywords, and Keywords Plus. The subsequent screening decision used titles, abstracts, and author keywords; Keywords Plus was not used as a screening criterion.

## Search flow

The final thesis workflow recorded the following sequence:

```text
117 WoS bibliographic records
        ↓ preliminary title/abstract/author-keyword screening
75 candidate publications
        ↓ comparison with the ECOTOX retrieval
24 automatic exact DOI/title overlaps
+ 4 high-confidence fuzzy matches verified manually
= 28 ECOTOX-overlapping publications
        ↓
47 publications not matched to the ECOTOX retrieval
        ↓ original-publication review and harmonization
17 publications contributed 57 quantitative LC50 Results
```

The four manually verified fuzzy overlaps were retained as explicit human adjudications rather than converted into a new automated matching rule.

## Repository scripts

```text
scripts/database_audits/wos/
├── 01_screen_wos_records.R
├── 02_check_wos_ecotox_overlap.R
└── 03_prepare_wos_unique_reference_review.R
```

### 01 — Preliminary screening

Reads the local raw WoS export and reproduces the 117-record screening step.

Expected local input:

```text
data/external/wos/wos_raw_export.xls
```

The script checks the frozen benchmarks:

- 117 raw records;
- 75 candidate publications (`Potentially eligible` + `Manual review required`).

### 02 — ECOTOX overlap audit

Compares the 75 candidate publications with the ECOTOX reference sampling frame.

Matching hierarchy:

1. DOI exact match, where available;
2. normalized-title exact match;
3. fuzzy-title suggestion for manual verification only.

Fuzzy matches are not treated as automatic matches.

### 03 — Unique-reference triage

Organizes the 47 publications not matched to ECOTOX for structured original-publication review. It does not itself make the final 17-publication / 57-Result inclusion decision.

## Why the raw WoS export is not in the public repository

The raw export is a provider-supplied bibliographic export obtained through database access. It is intentionally kept as a local external input rather than redistributed in this repository.

The public archive therefore provides:

- the exact query;
- search date and search field;
- the screening and overlap scripts;
- compact count summaries and workflow documentation;
- the curated source-verified records used in the final analytical dataset.

A user with appropriate Web of Science access can run the documented search, export the records, place the export at the expected local path, and rerun the audit scripts.

## Interpretation

Web of Science was used as an evidence-identification route. Quantitative LC50 values included in the combined dataset were extracted from and checked against the original publications rather than treated as values supplied by the WoS bibliographic export.
