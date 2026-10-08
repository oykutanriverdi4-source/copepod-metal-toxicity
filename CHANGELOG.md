# Change log

## Unreleased — source and reporting correction, 7 October 2026

- Correct O'Brien (1988) C6, Result E_R115201, to Adult in the harmonized working master; retain the reported Cu LC50 value of 12 micromoles per litre.
- Correct the 13 O'Brien harmonization descriptions and distinguish culture temperature from a separately verified numerical test temperature.
- Update the expected stage support from 69 to 70 records; retain the historical ECOTOX membership check separately from the current corrected analysis.
- Repair two Order interaction-extreme reports so that they select the actual minimum and maximum rather than the two smallest values.
- Keep singularity and other recorded numerical alerts separately identifiable; preserve raw optimizer messages.
- Add a Reference-FE numerical status report and explicit point-direction versus interval-inclusion counts for the existing pooled deletion checks.
- Clarify that direction-only stability summaries do not establish statistical support. No model formula or planned contrast is changed by these reporting repairs.
- Supply a checked-run launcher, reporting tests and a file/input manifest. The figure-generation script is unchanged in this correction: it reads the regenerated outputs.

**Release status:** the integrated correction has not yet been executed in the user's R environment or published. The archive identity in CITATION.cff still identifies the previous v1.0.0 release; do not describe the present working tree as that exact archive. This change log does not close the broader scientific source/quality audit or thesis-text revision.
