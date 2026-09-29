# ============================================================
# ESS THESIS REPRODUCIBILITY WORKFLOW
# 06_compare_exposure_duration.R
# ============================================================
# PURPOSE
# Reproduce the source-verified within-test exposure-duration
# comparisons reported in Sections 3.7.3 and 4.4.1 of the thesis.
#
# SCIENTIFIC DESIGN
# - Exposure duration is the factor allowed to vary.
# - Different observation times are compared only when source
#   verification establishes that they belong to the same acute
#   toxicity-test trajectory.
# - Repeated Results at the same observation time are reduced by
#   geometric mean before comparison.
# - The primary source-verified set contains 19 duration series.
# - Longest-versus-shortest comparisons are complemented by:
#     * six fixed duration-pair summaries among 24/48/72/96 h,
#     * complete-profile repeated-measures analyses for
#       48/72/96 h and 24/48/72/96 h,
#     * a stricter sensitivity excluding flagged extrapolated
#       endpoint LC50 values,
#     * an extended sensitivity including two secondary-supported
#       ECOTOX series with less certain same-test membership.
#
# IMPORTANT
# The source-adjudication decisions and statistical calculations
# are preserved from the final thesis analysis. Repository cleanup
# changes paths, file organization, labels, and documentation only.
# ============================================================


# ============================================================
# 0. PACKAGES + PATHS
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

options(width = 220)

input_file <- file.path(
  "data", "processed", "combined_ecotox_wos_harmonized.csv"
)

output_root <- file.path(
  "outputs",
  "06_compare_exposure_duration"
)

verification_dir <- file.path(output_root, "01_source_verification")
table_dir        <- file.path(output_root, "02_tables")
figure_dir       <- file.path(output_root, "03_figures")
object_dir       <- file.path(output_root, "04_objects")

for (d in c(output_root, verification_dir, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 1. LOAD COMBINED HARMONIZED DATA
# ============================================================

if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "Reference_Number",
  "Test_Number",
  "Result_Number",
  "Species",
  "Metal",
  "Lifestage",
  "Duration_days",
  "Temperature",
  "Salinity",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Harmonization_Status",
  "LC50_umol_L",
  "Source_Origin",
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "External_Stage_Sex",
  "External_Source_Note"
)

stopifnot(all(required_columns %in% names(dat)))

# Frozen benchmarks for the canonical combined dataset.
stopifnot(
  nrow(dat) == 353,
  n_distinct(dat$Reference_ID) == 84,
  sum(dat$Source_Origin == "ECOTOX") == 296,
  sum(dat$Source_Origin == "WoS_supplemental") == 57,
  sum(dat$Harmonization_Status == "HARMONIZED", na.rm = TRUE) == 304
)

dat <- dat %>%
  mutate(
    Lifestage = na_if(trimws(as.character(Lifestage)), ""),
    External_Stage_Sex = na_if(trimws(as.character(External_Stage_Sex)), ""),
    Stage_Label = coalesce(External_Stage_Sex, Lifestage)
  )

harm <- dat %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0,
    !is.na(Duration_days)
  )


# ============================================================
# 2. REPRODUCE ORIGINAL ECOTOX METADATA-DERIVED DISCOVERY SET
# ============================================================
# This reconstructs the metadata-derived ECOTOX candidate scan used
# before source adjudication and confirms that the frozen ECOTOX
# discovery baseline has not drifted.

ecotox <- harm %>%
  filter(Source_Origin == "ECOTOX")

duration_strict <- ecotox %>%
  filter(
    !is.na(Lifestage),
    !is.na(Temperature),
    !is.na(Salinity)
  )

duration_context_ecotox <- duration_strict %>%
  group_by(
    Reference_Number,
    Species,
    Metal,
    Lifestage,
    Duration_days,
    Temperature,
    Salinity,
    Exposure_Type,
    Conc_Type,
    Basis_Interpretation
  ) %>%
  summarise(
    n_results = n(),
    n_tests = n_distinct(Test_Number),
    LC50_duration_context = exp(mean(log(LC50_umol_L))),
    .groups = "drop"
  )

series_key_ecotox <- c(
  "Reference_Number",
  "Species",
  "Metal",
  "Lifestage",
  "Temperature",
  "Salinity",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation"
)

duration_series_ecotox <- duration_context_ecotox %>%
  group_by(across(all_of(series_key_ecotox))) %>%
  filter(n_distinct(Duration_days) >= 2) %>%
  ungroup()

summarize_series <- function(x, grouping_vars) {
  x %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      n_durations = n_distinct(Duration_days),
      Durations = paste(sort(unique(Duration_days)), collapse = " | "),
      shortest_duration = min(Duration_days),
      longest_duration = max(Duration_days),
      LC50_shortest = LC50_duration_context[which.min(Duration_days)],
      LC50_longest = LC50_duration_context[which.max(Duration_days)],
      .groups = "drop"
    ) %>%
    mutate(
      ratio_long_short = LC50_longest / LC50_shortest,
      ln_ratio_long_short = log(ratio_long_short),
      direction = case_when(
        ln_ratio_long_short < 0 ~ "Lower LC50 at longer duration",
        ln_ratio_long_short > 0 ~ "Higher LC50 at longer duration",
        TRUE ~ "No change"
      ),
      percent_change = (ratio_long_short - 1) * 100
    )
}

duration_series_summary_ecotox <- summarize_series(
  duration_series_ecotox,
  series_key_ecotox
) %>%
  arrange(Reference_Number, Species, Metal, Lifestage)

discovery_total_ecotox <- duration_series_summary_ecotox %>%
  summarise(
    total_series = n(),
    unique_references = n_distinct(Reference_Number),
    lower_at_longer_duration = sum(ln_ratio_long_short < 0),
    higher_at_longer_duration = sum(ln_ratio_long_short > 0)
  )

stopifnot(
  discovery_total_ecotox$total_series == 7,
  discovery_total_ecotox$unique_references == 6
)


# ============================================================
# 3. ORIGINAL ECOTOX SOURCE-ADJUDICATION LEDGER
# ============================================================

ecotox_source_ledger <- tribble(
  ~Reference_Number, ~Species, ~Metal, ~Lifestage,
  ~Source_Verification_Status, ~Primary_Verified_Use, ~Extended_Use,
  ~Verification_Note,
  
  3744, "Pseudodiaptomus coronatus", "Zn", "Adult",
  "SECONDARY_SUPPORTED_SAME_SERIES_UNRESOLVED", "NO", "YES",
  paste(
    "Original Lussier & Cardin (1985) Zn memorandum unavailable.",
    "EPA secondary documents support provenance/endpoints, but original",
    "same-series membership of 72 h and 96 h was not independently verified."
  ),
  
  8445, "Acartia tonsa", "Hg", "F2 generation",
  "SOURCE_VERIFIED_SAME_SERIES", "YES", "YES",
  paste(
    "Sosnowski & Gentile (1978), Table 3 reports Hg LC50 values",
    "at 24, 48, 72 and 96 h within the same static 96-h F2 assay."
  ),
  
  8445, "Acartia tonsa", "Hg", "F6 generation",
  "SOURCE_VERIFIED_SAME_SERIES", "YES", "YES",
  paste(
    "Sosnowski & Gentile (1978), Table 3 reports Hg LC50 values",
    "at 24, 48, 72 and 96 h within the same static 96-h F6 assay."
  ),
  
  14601, "Acartia clausi", "Cd", "Adult",
  "SECONDARY_SUPPORTED_SAME_SERIES_UNRESOLVED", "NO", "YES",
  paste(
    "Original Lussier & Cardin (1985) Cd memorandum unavailable.",
    "EPA secondary documents support provenance/endpoints, but original",
    "same-series membership of 72 h and 96 h was not independently verified."
  ),
  
  119554, "Tisbe battagliai", "Cu", "Gestation",
  "SOURCE_VERIFIED_SAME_SERIES", "YES", "YES",
  paste(
    "Diz et al. (2009) follows gravid adult females within the same Cu assay",
    "and reports adult LC50 values at 48 and 72 h."
  ),
  
  156333, "Nitocra spinipes", "Cu", "Adult",
  "EXCLUDE_SEPARATE_EXPERIMENTAL_EXERCISES", "NO", "NO",
  paste(
    "Ward et al. (2011): the 72-h and 96-h estimates came from separate",
    "experimental exercises and are not one duration trajectory."
  ),
  
  166160, "Tigriopus japonicus", "Cu", "Adult",
  "SOURCE_VERIFIED_SAME_SERIES", "YES", "YES",
  paste(
    "Bao et al. (2013) conducted one 96-h semi-static acute Cu assay;",
    "Table 2 reports the 72-h and 96-h LC50 values from that assay."
  )
)

stopifnot(
  nrow(ecotox_source_ledger) == 7,
  !anyDuplicated(
    ecotox_source_ledger %>%
      select(Reference_Number, Species, Metal, Lifestage)
  )
)

ecotox_series_adjudicated <- duration_series_summary_ecotox %>%
  left_join(
    ecotox_source_ledger,
    by = c("Reference_Number", "Species", "Metal", "Lifestage")
  )

stopifnot(
  nrow(ecotox_series_adjudicated) == 7,
  !any(is.na(ecotox_series_adjudicated$Source_Verification_Status))
)


# ============================================================
# 4. WoS SOURCE-VERIFIED SAME-SERIES LEDGER
# ============================================================
# Test_ID was assigned during source verification. For the rows below,
# repeated durations within a Test_ID are confirmed to belong to the
# same acute assay trajectory in the original paper.

wos_series_ledger <- tribble(
  ~Reference_ID, ~Test_ID, ~Species, ~Metal, ~Stage_Label,
  ~Endpoint_Extrapolation_Flag, ~Verification_Note,
  
  "W_002", "W_002_T01", "Pseudodiaptomus annandalei", "Cd", "Adult male",
  "NO", "Kadiene et al. (2017): male 72-h and 96-h Cd LC50s from the same acute assay.",
  
  "W_002", "W_002_T02", "Pseudodiaptomus annandalei", "Cd", "Adult female",
  "YES", "Kadiene et al. (2017): female 72-h and 96-h Cd LC50s from the same assay; 72-h LC50 exceeds the highest tested concentration.",
  
  "W_002", "W_002_T03", "Eurytemora affinis", "Cd", "Adult female",
  "NO", "Kadiene et al. (2017): female 48/72/96-h Cd LC50s from the same acute assay.",
  
  "W_002", "W_002_T04", "Eurytemora affinis", "Cd", "Adult male",
  "NO", "Kadiene et al. (2017): male 72-h and 96-h Cd LC50s from the same acute assay.",
  
  "W_003", "W_003_T01", "Calanus finmarchicus", "Hg", "Copepodite V (CV)",
  "YES", "Tollefsen et al. (2017): 24/48/72/96-h Hg LC50s from the same static assay; 24-h LC50 exceeds the tested range.",
  
  "W_004", "W_004_T01", "Pseudodiaptomus marinus", "Ni", "Adult",
  "YES", "Tlili et al. (2016): 24/48/72/96-h Ni LC50s from the same assay; 24-h LC50 slightly exceeds the highest tested concentration.",
  
  "W_008", "W_008_T01", "Calanus glacialis", "Zn", "Copepodite IV–V",
  "NO", "Cherkashin (2020): 24-h and 48-h Zn LC50s from the same acute exposure series.",
  
  "W_008", "W_008_T02", "Neocalanus plumchrus", "Zn", "Copepodite IV–V",
  "NO", "Cherkashin (2020): 24-h and 48-h Zn LC50s from the same acute exposure series.",
  
  "W_010", "W_010_T01", "Pseudodiaptomus marinus", "Hg", "Adult",
  "YES", "Tlili et al. (2019): 24/48/72/96-h Hg LC50s from the same assay; 24-h and 48-h estimates exceed the highest tested treatment.",
  
  "W_013", "W_013_T01", "Eurytemora affinis", "Cd", "Adult female",
  "NO", "Zidour et al. (2019): retained female Cd 48-h and 72-h endpoints belong to the same assay; duplicate 96-h value was removed before merge.",
  
  "W_013", "W_013_T02", "Eurytemora affinis", "Cu", "Adult female",
  "NO", "Zidour et al. (2019): female 72-h and 96-h Cu LC50s from the same assay.",
  
  "W_013", "W_013_T04", "Eurytemora affinis", "Cu", "Adult male",
  "NO", "Zidour et al. (2019): male 48/72/96-h Cu LC50s from the same assay.",
  
  "W_013", "W_013_T06", "Eurytemora affinis", "Ni", "Adult male",
  "NO", "Zidour et al. (2019): male 72-h and 96-h Ni LC50s from the same assay.",
  
  "W_015", "W_015_T01", "Calanus finmarchicus", "Hg", "Copepodite V (CV)",
  "NO", "Øverjordet et al. (2014): 48-h and 96-h Hg LC50s from the same acute assay.",
  
  "W_015", "W_015_T02", "Calanus glacialis", "Hg", "Copepodite V (CV)",
  "NO", "Øverjordet et al. (2014): 48-h and 96-h Hg LC50s from the same acute assay."
)

stopifnot(
  nrow(wos_series_ledger) == 15,
  !anyDuplicated(wos_series_ledger$Test_ID)
)


# ============================================================
# 5. RECONSTRUCT WoS SAME-SERIES ENDPOINTS
# ============================================================

wos <- harm %>%
  filter(Source_Origin == "WoS_supplemental")

wos_duration_context <- wos %>%
  inner_join(
    wos_series_ledger %>%
      select(
        Reference_ID,
        Test_ID,
        Ledger_Species = Species,
        Ledger_Metal = Metal,
        Ledger_Stage_Label = Stage_Label,
        Endpoint_Extrapolation_Flag,
        Verification_Note
      ),
    by = c("Reference_ID", "Test_ID")
  ) %>%
  filter(
    Species == Ledger_Species,
    Metal == Ledger_Metal
  ) %>%
  group_by(
    Reference_ID,
    Test_ID,
    Species,
    Metal,
    Stage_Label,
    Endpoint_Extrapolation_Flag,
    Verification_Note,
    Duration_days
  ) %>%
  summarise(
    n_results = n(),
    LC50_duration_context = exp(mean(log(LC50_umol_L))),
    Result_IDs = paste(Result_ID, collapse = " | "),
    .groups = "drop"
  )

wos_series_summary <- wos_duration_context %>%
  group_by(
    Reference_ID,
    Test_ID,
    Species,
    Metal,
    Stage_Label,
    Endpoint_Extrapolation_Flag,
    Verification_Note
  ) %>%
  filter(n_distinct(Duration_days) >= 2) %>%
  summarise(
    n_durations = n_distinct(Duration_days),
    Durations = paste(sort(unique(Duration_days)), collapse = " | "),
    shortest_duration = min(Duration_days),
    longest_duration = max(Duration_days),
    LC50_shortest = LC50_duration_context[which.min(Duration_days)],
    LC50_longest = LC50_duration_context[which.max(Duration_days)],
    .groups = "drop"
  ) %>%
  mutate(
    ratio_long_short = LC50_longest / LC50_shortest,
    ln_ratio_long_short = log(ratio_long_short),
    direction = case_when(
      ln_ratio_long_short < 0 ~ "Lower LC50 at longer duration",
      ln_ratio_long_short > 0 ~ "Higher LC50 at longer duration",
      TRUE ~ "No change"
    ),
    percent_change = (ratio_long_short - 1) * 100,
    Source_Origin = "WoS_supplemental",
    Source_Verification_Status = "SOURCE_VERIFIED_SAME_SERIES",
    Primary_Verified_Use = "YES",
    Extended_Use = "YES"
  )

stopifnot(
  nrow(wos_series_summary) == 15,
  n_distinct(wos_series_summary$Reference_ID) == 7
)


# ============================================================
# 6. RECONSTRUCT ECOTOX VERIFIED + EXTENDED SERIES
# ============================================================

ecotox_verified_summary <- ecotox_series_adjudicated %>%
  filter(Primary_Verified_Use == "YES") %>%
  transmute(
    Reference_ID = paste0("E_", Reference_Number),
    Test_ID = NA_character_,
    Species,
    Metal,
    Stage_Label = Lifestage,
    Endpoint_Extrapolation_Flag = "NO",
    Verification_Note,
    n_durations,
    Durations,
    shortest_duration,
    longest_duration,
    LC50_shortest,
    LC50_longest,
    ratio_long_short,
    ln_ratio_long_short,
    direction,
    percent_change,
    Source_Origin = "ECOTOX",
    Source_Verification_Status,
    Primary_Verified_Use,
    Extended_Use
  )

ecotox_extended_only <- ecotox_series_adjudicated %>%
  filter(
    Primary_Verified_Use == "NO",
    Extended_Use == "YES"
  ) %>%
  transmute(
    Reference_ID = paste0("E_", Reference_Number),
    Test_ID = NA_character_,
    Species,
    Metal,
    Stage_Label = Lifestage,
    Endpoint_Extrapolation_Flag = "NO",
    Verification_Note,
    n_durations,
    Durations,
    shortest_duration,
    longest_duration,
    LC50_shortest,
    LC50_longest,
    ratio_long_short,
    ln_ratio_long_short,
    direction,
    percent_change,
    Source_Origin = "ECOTOX",
    Source_Verification_Status,
    Primary_Verified_Use,
    Extended_Use
  )

stopifnot(
  nrow(ecotox_verified_summary) == 4,
  nrow(ecotox_extended_only) == 2
)


# ============================================================
# 7. COMBINED PRIMARY SOURCE-VERIFIED DURATION SET
# ============================================================

duration_verified_summary <- bind_rows(
  ecotox_verified_summary,
  wos_series_summary %>%
    select(names(ecotox_verified_summary))
) %>%
  arrange(Metal, Reference_ID, Species, Stage_Label)

duration_verified_total <- duration_verified_summary %>%
  summarise(
    total_series = n(),
    unique_references = n_distinct(Reference_ID),
    unique_species = n_distinct(Species),
    lower_at_longer_duration = sum(ln_ratio_long_short < 0),
    higher_at_longer_duration = sum(ln_ratio_long_short > 0),
    no_change = sum(ln_ratio_long_short == 0)
  )

stopifnot(
  duration_verified_total$total_series == 19,
  duration_verified_total$unique_references == 10,
  duration_verified_total$unique_species == 9,
  duration_verified_total$lower_at_longer_duration == 18,
  duration_verified_total$higher_at_longer_duration == 1
)

direction_by_metal <- duration_verified_summary %>%
  group_by(Metal) %>%
  summarise(
    Series = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    Lower_at_longer = sum(ln_ratio_long_short < 0),
    Higher_at_longer = sum(ln_ratio_long_short > 0),
    Median_ratio_long_short = median(ratio_long_short),
    Min_ratio = min(ratio_long_short),
    Max_ratio = max(ratio_long_short),
    .groups = "drop"
  ) %>%
  arrange(Metal)

source_contribution <- duration_verified_summary %>%
  group_by(Source_Origin) %>%
  summarise(
    Series = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    Lower_at_longer = sum(ln_ratio_long_short < 0),
    Higher_at_longer = sum(ln_ratio_long_short > 0),
    .groups = "drop"
  )


# ============================================================
# 8. EXTENDED SENSITIVITY SET
# ============================================================

duration_extended_summary <- bind_rows(
  duration_verified_summary,
  ecotox_extended_only
) %>%
  arrange(Metal, Reference_ID, Species, Stage_Label)

duration_extended_total <- duration_extended_summary %>%
  summarise(
    total_series = n(),
    unique_references = n_distinct(Reference_ID),
    unique_species = n_distinct(Species),
    lower_at_longer_duration = sum(ln_ratio_long_short < 0),
    higher_at_longer_duration = sum(ln_ratio_long_short > 0),
    no_change = sum(ln_ratio_long_short == 0)
  )

stopifnot(
  duration_extended_total$total_series == 21,
  duration_extended_total$unique_references == 12,
  duration_extended_total$unique_species == 11,
  duration_extended_total$lower_at_longer_duration == 20,
  duration_extended_total$higher_at_longer_duration == 1
)


# ============================================================
# 9. STRICTER SENSITIVITY: REMOVE SERIES WITH EXTRAPOLATED ENDPOINTS
# ============================================================
# This does not replace the main source-verified set. It asks whether
# the overall direction depends on series where an endpoint LC50 was
# extrapolated beyond the tested concentration range.

duration_non_extrapolated <- duration_verified_summary %>%
  filter(Endpoint_Extrapolation_Flag != "YES")

duration_non_extrapolated_total <- duration_non_extrapolated %>%
  summarise(
    total_series = n(),
    unique_references = n_distinct(Reference_ID),
    unique_species = n_distinct(Species),
    lower_at_longer_duration = sum(ln_ratio_long_short < 0),
    higher_at_longer_duration = sum(ln_ratio_long_short > 0),
    no_change = sum(ln_ratio_long_short == 0)
  )

stopifnot(
  duration_non_extrapolated_total$total_series == 15,
  duration_non_extrapolated_total$lower_at_longer_duration == 14,
  duration_non_extrapolated_total$higher_at_longer_duration == 1
)


# ============================================================
# 10. EXTRACT ENDPOINT-LEVEL PLOT DATA
# ============================================================

# ECOTOX verified endpoint rows
ecotox_verified_keys <- ecotox_source_ledger %>%
  filter(Primary_Verified_Use == "YES") %>%
  select(Reference_Number, Species, Metal, Lifestage)

ecotox_verified_rows <- duration_series_ecotox %>%
  semi_join(
    ecotox_verified_keys,
    by = c("Reference_Number", "Species", "Metal", "Lifestage")
  ) %>%
  transmute(
    Reference_ID = paste0("E_", Reference_Number),
    Test_ID = NA_character_,
    Species,
    Metal,
    Stage_Label = Lifestage,
    Duration_days,
    LC50_duration_context,
    Source_Origin = "ECOTOX",
    Endpoint_Extrapolation_Flag = "NO"
  )

# WoS verified endpoint rows
wos_verified_rows <- wos_duration_context %>%
  transmute(
    Reference_ID,
    Test_ID,
    Species,
    Metal,
    Stage_Label,
    Duration_days,
    LC50_duration_context,
    Source_Origin = "WoS_supplemental",
    Endpoint_Extrapolation_Flag
  )

duration_verified_rows <- bind_rows(
  ecotox_verified_rows,
  wos_verified_rows
) %>%
  mutate(
    series_id = if_else(
      Source_Origin == "ECOTOX",
      paste(Reference_ID, Species, Metal, Stage_Label, sep = " | "),
      paste(Reference_ID, Test_ID, Species, Metal, Stage_Label, sep = " | ")
    )
  ) %>%
  group_by(series_id) %>%
  arrange(Duration_days, .by_group = TRUE) %>%
  mutate(
    LC50_baseline = first(LC50_duration_context),
    ln_relative_LC50 = log(LC50_duration_context / LC50_baseline)
  ) %>%
  ungroup()

stopifnot(
  n_distinct(duration_verified_rows$series_id) == 19
)



# ============================================================
# 10A. SIX FIXED DURATION PAIRS WITHIN VERIFIED SERIES
# ============================================================
# Ratio = later LC50 / earlier LC50; below 1 means lower LC50 later.
# Do not pool different time contrasts. ALL_METALS is a descriptive
# available-evidence average, not a universal or metal-adjusted time effect.
# Average series log ratios per reference before averaging references.
# Conservative shared-origin sensitivity groups W_002/W_013 for inference
# only; it never joins endpoints from different assays into a new series.
# Extrapolation time keys below come from the existing source ledger notes.
# NO means no flagged extrapolation in that ledger, not proof of precision.
# Pointwise t CIs assume independent, approximately normal unit log ratios.
# Exact sign tests test direction among non-tied units, not mean magnitude.
# Holm family: six time pairs per metal/scope/unit scheme. CIs are unadjusted.
# R stats: t.test, binom.test, p.adjust; Holm (1979), SJS 6:65-70.

duration_pair_rows <- duration_verified_rows %>% mutate(
  Flagged_time = (Test_ID == "W_002_T02" & Duration_days == 3) |
    (Test_ID %in% c("W_003_T01", "W_004_T01") & Duration_days == 1) |
    (Test_ID == "W_010_T01" & Duration_days %in% c(1, 2)),
  Flagged_time = coalesce(Flagged_time, FALSE))
flag_check <- duration_pair_rows %>% group_by(series_id) %>% summarise(
  Has_flagged_time = any(Flagged_time),
  Ledger_flag = any(Endpoint_Extrapolation_Flag == "YES"), .groups = "drop")
stopifnot(all(flag_check$Has_flagged_time == flag_check$Ledger_flag))
time_specs <- combn(c(1, 2, 3, 4), 2, simplify = FALSE)
pair_list <- list()
for (sid in unique(duration_pair_rows$series_id)) {
  d <- filter(duration_pair_rows, series_id == sid)
  stopifnot(!anyDuplicated(d$Duration_days), n_distinct(d$Reference_ID) == 1,
    n_distinct(d$Metal) == 1, n_distinct(d$Species) == 1,
    all(is.finite(d$LC50_duration_context)), all(d$LC50_duration_context > 0))
  for (tt in time_specs) {
    ix <- match(tt, d$Duration_days)
    if (anyNA(ix)) next
    pair_list[[length(pair_list) + 1L]] <- tibble(
      Series_ID = sid, Reference_ID = d$Reference_ID[1], Metal = d$Metal[1],
      Species = d$Species[1], Earlier_h = tt[1] * 24, Later_h = tt[2] * 24,
      Comparison = paste(tt * 24, collapse = "-"),
      Earlier_LC50 = d$LC50_duration_context[ix[1]],
      Later_LC50 = d$LC50_duration_context[ix[2]],
      Log_ratio = log(Later_LC50 / Earlier_LC50),
      Ratio_later_earlier = Later_LC50 / Earlier_LC50,
      Pair_has_flagged_time = any(d$Flagged_time[ix]),
      Series_has_flagged_time = any(d$Flagged_time))
  }
}
duration_pair_contexts <- bind_rows(pair_list)
summary_list <- list(); unit_list <- list(); omission_list <- list()
for (scope in c("ALL_VERIFIED", "PAIR_WITHOUT_FLAGGED_TIMES", "WHOLE_SERIES_WITHOUT_FLAGS")) {
  z <- duration_pair_contexts
  if (scope == "PAIR_WITHOUT_FLAGGED_TIMES") z <- filter(z, !Pair_has_flagged_time)
  if (scope == "WHOLE_SERIES_WITHOUT_FLAGS") z <- filter(z, !Series_has_flagged_time)
  for (scheme in c("REFERENCE", "SHARED_ORIGIN_SENSITIVITY")) {
    z <- z %>% mutate(Unit_ID = if_else(scheme == "SHARED_ORIGIN_SENSITIVITY" &
      Reference_ID %in% c("W_002", "W_013"), "W002_W013", Reference_ID))
    for (metal in c("ALL_METALS", sort(unique(duration_pair_contexts$Metal)))) {
      for (tt in time_specs) {
        label <- paste(tt * 24, collapse = "-")
        d <- filter(z, Comparison == label)
        if (metal != "ALL_METALS") d <- filter(d, Metal == metal)
        u <- d %>% group_by(Unit_ID) %>% summarise(Series = n(),
          Mean_log_ratio = mean(Log_ratio), Ratio_later_earlier = exp(Mean_log_ratio),
          .groups = "drop") %>% mutate(Scope = scope, Unit_scheme = scheme,
            Metal = metal, Comparison = label, .before = 1)
        unit_list[[length(unit_list) + 1L]] <- u
        x <- u$Mean_log_ratio; k <- length(x)
        estimate <- if (k) exp(mean(x)) else NA_real_
        lo <- hi <- NA_real_
        status <- if (k < 2) "FEWER_THAN_TWO_UNITS" else "ZERO_VARIANCE"
        if (k >= 2 && is.finite(sd(x)) && sd(x) > 0) {
          margin <- qt(0.975, k-1) * sd(x) / sqrt(k)
          lo <- exp(mean(x)-margin); hi <- exp(mean(x)+margin)
          status <- "POINTWISE_T_INTERVAL_ASSUMPTION_DEPENDENT"
        }
        neg <- sum(x < -1e-12); pos <- sum(x > 1e-12)
        pv <- if (k >= 2 && neg + pos > 0)
          binom.test(neg, neg + pos, p = 0.5)$p.value else NA_real_
        omin <- omax <- NA_real_
        if (k >= 2) {
          o <- vapply(seq_len(k), function(j) exp(mean(x[-j])), numeric(1))
          omin <- min(o); omax <- max(o)
          omission_list[[length(omission_list) + 1L]] <- tibble(Scope = scope,
            Unit_scheme = scheme, Metal = metal, Comparison = label,
            Omitted_unit = u$Unit_ID, Remaining_units = k-1,
            Ratio_later_earlier = o)
        }
        summary_list[[length(summary_list) + 1L]] <- tibble(Scope = scope,
          Unit_scheme = scheme, Metal = metal, Comparison = label,
          Earlier_h = tt[1]*24, Later_h = tt[2]*24, Matched_series = nrow(d),
          References = n_distinct(d$Reference_ID), Units = k,
          Ratio_later_earlier = estimate, Percent_change = 100*(estimate-1),
          CI95_low = lo, CI95_high = hi, Interval_status = status,
          Lower_later_units = neg, Higher_later_units = pos, Tied_units = k-neg-pos,
          Sign_p = pv, Leave_one_unit_out_min = omin, Leave_one_unit_out_max = omax,
          Support_note = if (k == 0) "NO_MATCHED_PAIR" else
            if (k == 1) "ONE_UNIT_DESCRIPTIVE_ONLY" else "REPORT_SMALL_UNIT_COUNTS")
      }
    }
  }
}
duration_pair_effect_summary <- bind_rows(summary_list) %>%
  group_by(Scope, Unit_scheme, Metal) %>%
  mutate(Sign_p_Holm = p.adjust(Sign_p, "holm", n = 6L)) %>% ungroup()
duration_pair_unit_effects <- bind_rows(unit_list)
duration_pair_omissions <- bind_rows(omission_list)
write_csv(duration_pair_contexts, file.path(table_dir, "duration_pair_contexts.csv"))
write_csv(duration_pair_effect_summary, file.path(table_dir, "duration_pair_effect_summary.csv"))
write_csv(duration_pair_unit_effects, file.path(table_dir, "duration_pair_unit_effects.csv"))
write_csv(duration_pair_omissions, file.path(table_dir, "duration_pair_leave_one_unit_out.csv"))
cat("\n--- FIXED TIME PAIRS: DESCRIPTIVE ALL-METAL SUMMARY ---\n")
print(duration_pair_effect_summary %>% filter(Metal == "ALL_METALS",
  Unit_scheme == "REFERENCE"), n = Inf)

# ============================================================
# 10B. COMPLETE-PROFILE REPEATED MEASURES (LOG LC50)
# ============================================================
# Question: do time means differ within complete assay profiles?
# LC50 estimates are repeated, not individual animals.
# Select complete SERIES first; never stitch times across series.
# Average complete-series logs within reference; references have equal weight.
# Shared-origin sensitivity averages reference profiles within origin.
# Pooled metals describe this supported mixture, not a universal time effect.
# ANOVA: balanced one-factor repeated measures, Greenhouse-Geisser correction.
# Friedman: stats::friedman.test, asymptotic chi-square p (small n caution).
# Sources: R stats aov / friedman.test documentation;
# car R/Anova.R Greenhouse-Geisser implementation (trace covariance formula).
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/aov.html
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/friedman.test.html
# https://rdrr.io/cran/car/src/R/Anova.R
# No omnibus test here establishes monotonicity or identifies a time pair.

rm_test_profile <- function(y) {
  n <- nrow(y); k <- ncol(y)
  out <- tibble(ANOVA_F = NA_real_, DF_time = k-1,
    DF_error = (n-1)*(k-1), GG_epsilon = NA_real_, ANOVA_p_GG = NA_real_,
    Friedman_chisq = NA_real_, Friedman_p_asymptotic = NA_real_,
    ANOVA_status = "FEWER_THAN_TWO_UNITS",
    Friedman_status = "FEWER_THAN_TWO_UNITS")
  if (n < 2) return(out)
  stopifnot(all(is.finite(y)))
  residual <- sweep(sweep(y, 1, rowMeans(y)), 2, colMeans(y)) + mean(y)
  ss_error <- sum(residual^2)
  ss_time <- n * sum((colMeans(y)-mean(y))^2)
  out$ANOVA_status <- "DEGENERATE_RESIDUAL_VARIANCE"
  if (ss_error > 1e-14 * max(1, sum((y-mean(y))^2))) {
    out$ANOVA_F <- (ss_time/(k-1))/(ss_error/((n-1)*(k-1)))
    h <- diag(k) - matrix(1/k, k, k)
    sc <- h %*% cov(y) %*% h
    eps <- sum(diag(sc))^2 / ((k-1)*sum(sc^2))
    out$GG_epsilon <- max(1/(k-1), min(1, eps))
    out$ANOVA_p_GG <- pf(out$ANOVA_F, (k-1)*out$GG_epsilon,
      (n-1)*(k-1)*out$GG_epsilon, lower.tail = FALSE)
    out$ANOVA_status <- "COMPUTED_SMALL_N_ASSUMPTION_DEPENDENT"
  }
  fr <- suppressWarnings(tryCatch(stats::friedman.test(y), error = function(e) NULL))
  out$Friedman_status <- "DEGENERATE_RANKS_OR_TEST_ERROR"
  if (!is.null(fr) && is.finite(fr$p.value)) {
    out$Friedman_chisq <- unname(fr$statistic)
    out$Friedman_p_asymptotic <- fr$p.value
    out$Friedman_status <- "ASYMPTOTIC_P_SMALL_N_CAUTION"
  }
  out
}

rm_summaries <- list(); rm_profiles <- list(); rm_eligibility <- list()
rm_specs <- list(c(2,3,4), c(1,2,3,4))
rm_scopes <- c("ALL_VERIFIED", "ANALYZED_TIMES_WITHOUT_FLAGS", "WHOLE_SERIES_WITHOUT_FLAGS")
rm_metals <- c("ALL_METALS", sort(unique(duration_pair_rows$Metal)))
for (tt in rm_specs) {
  label <- paste(tt*24, collapse = "-")
  for (scope in rm_scopes) {
    eligible_ids <- character()
    for (sid in unique(duration_pair_rows$series_id)) {
      d <- filter(duration_pair_rows, series_id == sid)
      complete <- all(tt %in% d$Duration_days)
      time_flag <- any(d$Flagged_time[d$Duration_days %in% tt])
      series_flag <- any(d$Flagged_time)
      reason <- if (!complete) "MISSING_REQUIRED_TIME" else
        if (scope == "ANALYZED_TIMES_WITHOUT_FLAGS" && time_flag) "FLAGGED_ANALYZED_TIME" else
        if (scope == "WHOLE_SERIES_WITHOUT_FLAGS" && series_flag) "FLAGGED_SERIES" else "INCLUDED"
      if (reason == "INCLUDED") eligible_ids <- c(eligible_ids, sid)
      rm_eligibility[[length(rm_eligibility)+1L]] <- tibble(Time_set = label,
        Scope = scope, Series_ID = sid, Reference_ID = d$Reference_ID[1],
        Metal = d$Metal[1], Eligibility = reason)
    }
    z <- filter(duration_pair_rows, series_id %in% eligible_ids, Duration_days %in% tt)
    for (metal in rm_metals) {
      d <- z
      if (metal != "ALL_METALS") d <- filter(d, Metal == metal)
      refs <- d %>% group_by(Reference_ID, Duration_days) %>% summarise(
        Mean_log_LC50 = mean(log(LC50_duration_context)), .groups = "drop")
      for (scheme in c("REFERENCE", "SHARED_ORIGIN_SENSITIVITY")) {
        u <- refs %>% mutate(Unit_ID = if_else(scheme == "SHARED_ORIGIN_SENSITIVITY" &
          Reference_ID %in% c("W_002", "W_013"), "W002_W013", Reference_ID)) %>%
          group_by(Unit_ID, Duration_days) %>% summarise(
            Mean_log_LC50 = mean(Mean_log_LC50), .groups = "drop")
        ids <- sort(unique(u$Unit_ID))
        y <- matrix(NA_real_, length(ids), length(tt), dimnames = list(ids, as.character(tt*24)))
        if (nrow(u)) {
          y[cbind(match(u$Unit_ID, ids), match(u$Duration_days, tt))] <- u$Mean_log_LC50
          stopifnot(nrow(u) == length(ids)*length(tt), !anyNA(y))
        }
        rm_profiles[[length(rm_profiles)+1L]] <- u %>% mutate(Time_set = label,
          Scope = scope, Unit_scheme = scheme, Metal = metal,
          Time_h = Duration_days*24, .before = 1)
        rm_summaries[[length(rm_summaries)+1L]] <- bind_cols(tibble(Time_set = label,
          Scope = scope, Unit_scheme = scheme, Metal = metal,
          Complete_series = n_distinct(d$series_id), References = n_distinct(d$Reference_ID),
          Units = length(ids), Metals_present = paste(sort(unique(d$Metal)), collapse = ";")),
          rm_test_profile(y))
      }
    }
  }
}
duration_repeated_measures_summary <- bind_rows(rm_summaries) %>%
  group_by(Scope, Unit_scheme) %>% mutate(
    ANOVA_p_GG_Holm = p.adjust(ANOVA_p_GG, "holm", n = length(rm_specs)*length(rm_metals)),
    Friedman_p_Holm = p.adjust(Friedman_p_asymptotic, "holm", n = length(rm_specs)*length(rm_metals))) %>%
  ungroup()
# Holm family: both time sets and all pooled/metal-specific tests within each
# sensitivity scope and unit scheme; ANOVA and Friedman reported separately.
duration_repeated_measures_profiles <- bind_rows(rm_profiles)
duration_repeated_measures_eligibility <- bind_rows(rm_eligibility)
write_csv(duration_repeated_measures_summary, file.path(table_dir, "duration_repeated_measures_summary.csv"))
write_csv(duration_repeated_measures_profiles, file.path(table_dir, "duration_repeated_measures_profiles.csv"))
write_csv(duration_repeated_measures_eligibility, file.path(table_dir, "duration_repeated_measures_eligibility.csv"))
cat("\n--- COMPLETE-PROFILE REPEATED MEASURES: ALL METALS ---\n")
print(filter(duration_repeated_measures_summary, Metal == "ALL_METALS", Unit_scheme == "REFERENCE"), n = Inf)

# ============================================================
# 11. FIGURES
# ============================================================

duration_verified_plot <- ggplot(
  duration_verified_rows,
  aes(
    x = Duration_days,
    y = ln_relative_LC50,
    group = series_id,
    linetype = Source_Origin
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.4) +
  facet_wrap(~ Metal, ncol = 3) +
  scale_x_continuous(breaks = 1:4) +
  labs(
    x = "Exposure duration (days)",
    y = "ln(relative LC50)",
    linetype = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.justification = "left"
  )

# A second figure restricted to non-extrapolated series.
non_extrapolated_ids <- duration_non_extrapolated %>%
  mutate(
    series_id = if_else(
      Source_Origin == "ECOTOX",
      paste(Reference_ID, Species, Metal, Stage_Label, sep = " | "),
      paste(Reference_ID, Test_ID, Species, Metal, Stage_Label, sep = " | ")
    )
  ) %>%
  pull(series_id)

duration_non_extrapolated_rows <- duration_verified_rows %>%
  filter(series_id %in% non_extrapolated_ids)

duration_non_extrapolated_plot <- ggplot(
  duration_non_extrapolated_rows,
  aes(
    x = Duration_days,
    y = ln_relative_LC50,
    group = series_id,
    linetype = Source_Origin
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.4) +
  facet_wrap(~ Metal, ncol = 3) +
  scale_x_continuous(breaks = 1:4) +
  labs(
    x = "Exposure duration (days)",
    y = "ln(relative LC50)",
    linetype = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    legend.justification = "left"
  )


# ============================================================
# 12. WRITE TABLES + QA OUTPUTS
# ============================================================

write_csv(
  ecotox_source_ledger,
  file.path(verification_dir, "ECOTOX_duration_source_verification_ledger.csv")
)

write_csv(
  wos_series_ledger,
  file.path(verification_dir, "WoS_duration_source_verification_ledger.csv")
)

write_csv(
  ecotox_series_adjudicated,
  file.path(table_dir, "ECOTOX_duration_candidate_series_adjudicated.csv")
)

write_csv(
  wos_series_summary,
  file.path(table_dir, "WoS_source_verified_duration_series.csv")
)

write_csv(
  duration_verified_summary,
  file.path(table_dir, "duration_verified_series_combined.csv")
)

write_csv(
  duration_verified_total,
  file.path(table_dir, "duration_verified_total_combined.csv")
)

write_csv(
  direction_by_metal,
  file.path(table_dir, "duration_direction_by_metal.csv")
)

write_csv(
  source_contribution,
  file.path(table_dir, "duration_source_contribution.csv")
)

write_csv(
  duration_extended_summary,
  file.path(table_dir, "duration_extended_sensitivity_series_combined.csv")
)

write_csv(
  duration_extended_total,
  file.path(table_dir, "duration_extended_sensitivity_total_combined.csv")
)

write_csv(
  duration_non_extrapolated,
  file.path(table_dir, "duration_non_extrapolated_sensitivity_series.csv")
)

write_csv(
  duration_non_extrapolated_total,
  file.path(table_dir, "duration_non_extrapolated_sensitivity_total.csv")
)

write_csv(
  duration_verified_rows,
  file.path(table_dir, "duration_verified_plot_data_combined.csv")
)

# Explicitly save the one counter-direction series for audit visibility.
counter_direction_series <- duration_verified_summary %>%
  filter(ln_ratio_long_short > 0)

write_csv(
  counter_direction_series,
  file.path(table_dir, "duration_counter_direction_series.csv")
)

# Evidence-set synthesis table used to document the effect of source verification.
duration_ecotox_vs_combined <- tibble(
  Evidence_Set = c(
    "Original ECOTOX source-verified",
    "Combined ECOTOX + WoS source-verified",
    "Combined strict non-extrapolated",
    "Combined extended sensitivity"
  ),
  Series = c(
    4L,
    duration_verified_total$total_series,
    duration_non_extrapolated_total$total_series,
    duration_extended_total$total_series
  ),
  References = c(
    3L,
    duration_verified_total$unique_references,
    duration_non_extrapolated_total$unique_references,
    duration_extended_total$unique_references
  ),
  Lower_at_longer = c(
    4L,
    duration_verified_total$lower_at_longer_duration,
    duration_non_extrapolated_total$lower_at_longer_duration,
    duration_extended_total$lower_at_longer_duration
  ),
  Higher_at_longer = c(
    0L,
    duration_verified_total$higher_at_longer_duration,
    duration_non_extrapolated_total$higher_at_longer_duration,
    duration_extended_total$higher_at_longer_duration
  )
)

write_csv(
  duration_ecotox_vs_combined,
  file.path(table_dir, "duration_ECOTOX_vs_combined_summary.csv")
)


# ============================================================
# 13. SAVE FIGURES
# ============================================================

ggsave(
  file.path(figure_dir, "duration_within_verified_combined.png"),
  duration_verified_plot,
  width = 8.2,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figure_dir, "duration_within_verified_combined.pdf"),
  duration_verified_plot,
  width = 8.2,
  height = 5.2,
  units = "in"
)

ggsave(
  file.path(figure_dir, "duration_within_non_extrapolated_sensitivity.png"),
  duration_non_extrapolated_plot,
  width = 8.2,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figure_dir, "duration_within_non_extrapolated_sensitivity.pdf"),
  duration_non_extrapolated_plot,
  width = 8.2,
  height = 5.2,
  units = "in"
)


# ============================================================
# 14. SAVE OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    ecotox_discovery = duration_series_summary_ecotox,
    ecotox_source_ledger = ecotox_source_ledger,
    wos_series_ledger = wos_series_ledger,
    duration_verified_summary = duration_verified_summary,
    duration_extended_summary = duration_extended_summary,
    duration_non_extrapolated = duration_non_extrapolated,
    duration_verified_rows = duration_verified_rows,
    duration_pair_contexts = duration_pair_contexts,
    duration_pair_effect_summary = duration_pair_effect_summary,
    duration_pair_unit_effects = duration_pair_unit_effects,
    duration_pair_omissions = duration_pair_omissions,
    duration_repeated_measures_summary = duration_repeated_measures_summary,
    duration_repeated_measures_profiles = duration_repeated_measures_profiles,
    duration_repeated_measures_eligibility = duration_repeated_measures_eligibility,
    duration_verified_plot = duration_verified_plot,
    duration_non_extrapolated_plot = duration_non_extrapolated_plot
  ),
  file.path(object_dir, "duration_within_verified_combined_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)


# ============================================================
# 15. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("06 WITHIN-TEST EXPOSURE-DURATION COMPARISONS\n")
cat("============================================================\n")

cat("\n--- ORIGINAL ECOTOX METADATA DISCOVERY ---\n")
print(discovery_total_ecotox)

cat("\n--- COMBINED SOURCE-VERIFIED SET ---\n")
print(duration_verified_total)

cat("\n--- DIRECTION BY METAL ---\n")
print(direction_by_metal)

cat("\n--- SOURCE CONTRIBUTION ---\n")
print(source_contribution)

cat("\n--- COUNTER-DIRECTION SERIES ---\n")
print(
  counter_direction_series %>%
    select(
      Reference_ID,
      Species,
      Metal,
      Stage_Label,
      shortest_duration,
      longest_duration,
      LC50_shortest,
      LC50_longest,
      ratio_long_short,
      percent_change
    )
)

cat("\n--- STRICT NON-EXTRAPOLATED SENSITIVITY ---\n")
print(duration_non_extrapolated_total)

cat("\n--- EXTENDED SENSITIVITY SET ---\n")
print(duration_extended_total)

cat(
  "\nInterpretive summary:\n",
  "The source-verified ECOTOX subset shows 4/4 series with lower LC50 at longer duration.\n",
  "The combined source-verified set shows 18/19 series in that direction.\n",
  "The single counter-direction series is retained and reported.\n",
  "A stricter sensitivity excluding series with flagged extrapolated endpoint LC50s is also exported.\n",
  sep = ""
)

manifest <- tibble(
  Role = c(
    "Input: canonical combined dataset",
    "Output: primary duration-series summary",
    "Output: fixed duration-pair summary",
    "Output: complete-profile repeated-measures summary",
    "Output: primary duration figure",
    "Output: R analysis objects"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "duration_verified_series_combined.csv"),
    file.path(table_dir, "duration_pair_effect_summary.csv"),
    file.path(table_dir, "duration_repeated_measures_summary.csv"),
    file.path(figure_dir, "duration_within_verified_combined.png"),
    file.path(object_dir, "duration_within_verified_combined_objects.rds")
  )
) %>%
  mutate(
    Exists = file.exists(Path),
    MD5 = if_else(
      Exists,
      unname(tools::md5sum(Path)),
      NA_character_
    )
  )

write_csv(
  manifest,
  file.path(output_root, "input_output_manifest.csv")
)

writeLines(
  c(
    "STATUS: 06 COMPARE EXPOSURE DURATION = PASS",
    paste0(
      "Primary source-verified duration series: ",
      duration_verified_total$total_series
    ),
    paste0(
      "Lower LC50 at longer duration: ",
      duration_verified_total$lower_at_longer_duration,
      "/",
      duration_verified_total$total_series
    )
  ),
  con = file.path(output_root, "RUN_COMPLETE.txt")
)

required_outputs <- c(
  file.path(verification_dir, "ECOTOX_duration_source_verification_ledger.csv"),
  file.path(verification_dir, "WoS_duration_source_verification_ledger.csv"),
  file.path(table_dir, "duration_verified_series_combined.csv"),
  file.path(table_dir, "duration_direction_by_metal.csv"),
  file.path(table_dir, "duration_counter_direction_series.csv"),
  file.path(table_dir, "duration_ECOTOX_vs_combined_summary.csv"),
  file.path(figure_dir, "duration_within_verified_combined.png"),
  file.path(object_dir, "duration_within_verified_combined_objects.rds"),
  file.path(output_root, "input_output_manifest.csv"),
  file.path(output_root, "sessionInfo.txt"),
  file.path(output_root, "RUN_COMPLETE.txt")
)

stopifnot(all(file.exists(required_outputs)))

cat("\nSTATUS: 06 COMPARE EXPOSURE DURATION = PASS\n")

# ============================================================
# END
# ============================================================
