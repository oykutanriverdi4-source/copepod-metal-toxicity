# EnviroTox search and overlap-audit workflow

## Search

**Database:** EnviroTox version 2.0.0  
**Search date:** 5 September 2026  
**Interface:** Advanced Search

The following conditions were combined with AND:

```text
Heavy Metals contains "1"
Trophic Level contains "INVERT"
Test type = "A"
Test statistic = "LC50"
```

## Preliminary screening flow

The final thesis workflow recorded:

```text
3,375 raw EnviroTox test records
        ↓ copepod orders: Calanoida, Harpacticoida, Cyclopoida
106 records
        ↓ eight prespecified metals
76 records
        ↓ exclude 22 records classified as Freshwater
54 Saltwater candidate records
```

The 54 Saltwater candidates represented 17 species and 23 source citations.

The EnviroTox `Medium` field was used as a database-level screening proxy for marine/estuarine relevance. Experimental conditions were considered during subsequent source review.

## ECOTOX overlap audit

The 54 candidates were compared with the ECOTOX copepod retrieval.

The Result-level comparison considered:

- source/reference identity;
- species;
- metal;
- exposure duration;
- LC50 value after conversion to mg/L where needed.

Numerical equality was evaluated using base R `all.equal()` with its default tolerance and the larger absolute value as the comparison scale.

Final automatic audit:

```text
54 EnviroTox candidates
├── 38 exact ECOTOX Result overlaps
└── 16 unresolved candidates
    ├── 12 same source + species + metal, but no matching duration
    ├──  2 same source/species/metal/duration, but different LC50
    ├──  1 source-only match
    └──  1 source not matched
```

The 16 unresolved records were retained for original-source review rather than automatically counted as additional evidence.

## Source-review outcome used in the thesis

Review of the unresolved candidates established that no separate EnviroTox-derived quantitative records should be added to the combined dataset.

Important source-specific outcomes included:

- the zinc record attributed to Nipper et al. (1993) was confirmed as already represented in ECOTOX;
- the copper record for *Calanus plumchrus* attributed to Reeve et al. (1976) was treated as a duplicate after resolving the exposure-duration discrepancy and considering the correction reported by Heuschele et al. (2022);
- the copper LC50 from Sosnowski et al. (1979) was traced to an existing ECOTOX record that had been excluded during initial screening and was reinstated after review of the original publication;
- the Luoma et al. (1983) record had an exposure duration of 20–24 h rather than the 96 h listed in EnviroTox, so eligibility under the thesis 24–96 h window could not be established;
- the remaining unresolved candidates could not be confirmed as additional eligible Results because of unresolved differences in duration, concentration, source attribution, or unavailable original-source evidence.

## Repository scripts

```text
scripts/database_audits/envirotox/
├── 01_screen_envirotox_records.R
├── 02_check_envirotox_ecotox_overlap.R
└── 03_prepare_envirotox_manual_review.R
```

The third script prepares an editable manual-review workbook. It intentionally reports `WORKBOOK_PREPARED` rather than claiming that source review has been completed automatically.

## Local input and public-repository boundary

Expected local raw input:

```text
data/external/envirotox/envirotox_raw_export.xlsx
```

EnviroTox is publicly accessible. In this repository, the downloaded raw EnviroTox workbook and detailed row-level audit workbooks are nevertheless kept local as a conservative redistribution choice. The exact search settings, scripts, compact screening counts, matching-method documentation, and final analytical provenance are retained so that the audit can be reproduced from a fresh EnviroTox export.

## Final analytical role

EnviroTox functioned as a coverage and overlap audit. It did not contribute a separate set of quantitative records to the final combined dataset, although it contributed to source checking and to the recovery of one ECOTOX record.
