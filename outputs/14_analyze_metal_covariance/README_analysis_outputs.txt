Start with 02_tables/metal_covariance_main_scope_overview.csv.
This run writes to outputs/14_analyze_metal_covariance.
RNG settings are explicitly fixed and recorded in 01_audit/metal_covariance_run_settings.csv.
For identical inputs and computation environment, bootstrap sampling no longer inherits the session RNG kind.
Complementary, exploratory covariance analysis; no p-values calculated.
A point is one Reference_ID + Source_Context_ID; within-metal values use geometric means.
Each reference has equal total weight; contexts within a reference share that weight.
Covariance uses normalized weighted moments, not an unbiased sample-covariance denominator.
Bootstrap bounds are 2.5% and 97.5% percentiles for equal-reference correlation r.
They are not intervals for covariance components. Valid resample counts are reported.
Bootstrap and leave-one-reference-out apply to approved_ledger only.
The three-reference bootstrap floor is a reporting rule, not proof of sufficient support.
Very few references can yield unstable or degenerate bootstrap distributions.
A reference with one context contributes no within-reference variation.
Two contexts give a degenerate correlation; one reference supports no between-reference inference.
Between/total covariance may be negative or above one; it is not variance explained.
Strict verification is applied to the entire context, as in the original script.
The order-label sensitivity excludes missing Species/Order as well as equal labels, as before.
Source weighting does not establish independence. Associations do not demonstrate mixture effects.
Pairs with at least two contexts are plotted; all pairs remain in the Word/support tables.
Source names/titles are in 01_audit/metal_covariance_reference_index.csv; author-year labels are not inferred.
Source-verification decisions are frozen in the curated ledgers and are not changed here.
This plan followed earlier exploratory results; it is not a preregistered analysis.
This public-repository script preserves the final thesis covariance calculations and reproducibility settings.
