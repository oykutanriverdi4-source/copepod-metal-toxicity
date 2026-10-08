# ============================================================
# COPEPOD METAL TOXICITY THESIS
# 09_compare_stages_within_references.R
# ============================================================
# STATUS: REPOSITORY VERSION / SOURCE-AUDITED / COMBINED ECOTOX + WoS
#
# PURPOSE
# Preserve the original exact within-context developmental-stage
# analysis and add only newly source-verified comparable WoS evidence.
#
# IMPORTANT
# - This script does NOT fit the broader Metal x Stage model.
# - It preserves the original ECOTOX source decisions and exact-context
#   matching logic from the original ECOTOX within-stage workflow.
# - WoS rows are NOT admitted by an automatic "same Reference" rule.
#   They enter only through an explicit source-verified ledger.
# - Charry et al. (2019) is NOT a direct stage series because the
#   nauplii and adult LC50 endpoints use different durations.
# - Kadiene et al. (2019) IS a direct stage series: nauplii and
#   copepodids were tested in the same 96-h Cd experiment under the
#   same temperature/salinity conditions.
#
# EXPECTED FINAL ARCHITECTURE
# Original ECOTOX: 9 series / 6 References / 5 Species
# Combined:        10 series / 7 References / 6 Species
# ============================================================


# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

required_packages <- c("tidyverse", "gt")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(gt)

options(width = 220)

output_root <- file.path(
  "outputs",
  "09_compare_stages_within_references"
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

input_file <- file.path(
  "data",
  "processed",
  "combined_ecotox_wos_harmonized.csv"
)

if (!file.exists(input_file)) {
  stop(
    "Missing required input file: ", input_file,
    "\nRun 02_integrate_wos_verified_records.R first."
  )
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
  "Duration_days",
  "Lifestage",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Harmonization_Status",
  "LC50_umol_L",
  "Temperature",
  "Salinity",
  "Temperature_Raw_Mean",
  "Temperature_Raw_Min",
  "Temperature_Raw_Max",
  "Temperature_Raw_Units",
  "Salinity_Raw_Mean",
  "Salinity_Raw_Min",
  "Salinity_Raw_Max",
  "Salinity_Raw_Units",
  "Source_Origin",
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "External_Stage_Sex"
)

stopifnot(all(required_columns %in% names(dat)))

stopifnot(
  nrow(dat) == 353,
  n_distinct(dat$Reference_ID) == 84,
  sum(dat$Source_Origin == "ECOTOX") == 296,
  sum(dat$Source_Origin == "WoS_supplemental") == 57
)

harm <- dat %>%
  mutate(
    Lifestage = na_if(trimws(as.character(Lifestage)), ""),
    External_Stage_Sex = na_if(trimws(as.character(External_Stage_Sex)), "")
  ) %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0
  )

ecotox <- harm %>%
  filter(Source_Origin == "ECOTOX")

wos <- harm %>%
  filter(Source_Origin == "WoS_supplemental")


# ============================================================
# 2. PRESERVE ORIGINAL ECOTOX SOURCE-LEVEL DECISIONS
# ============================================================

stage_source_decisions <- tribble(
  ~Reference_Number, ~Metal, ~Final_stage_use,
  ~Stage_verification_status, ~Stage_verification_note,
  
  2977, "Cu", "YES", "SOURCE_COMMON_PROTOCOL_CULTURE_TEMP_ONLY",
  "O'Brien et al. 1988 pp.60-62: common experimental protocol and SOW medium; 16 +/- 1 C is reported in the culture paragraph, not separately as numeric test temperature. C6 is Adult. This comparison is protocol-supported, not a verified numeric-temperature match.",
  
  19281, "Cd", "YES", "SOURCE_VERIFIED_ENV",
  "Forget et al. (1998); Temperature 20 C and salinity 35",
  
  14137, "Cu", "YES", "SOURCE_VERIFIED_RAW_RANGE",
  "Adult-Nauplii contrast; RAW/source environmental range recovered",
  
  14137, "Cd", "YES", "SOURCE_VERIFIED_RAW_RANGE",
  "Adult-Nauplii contrast; RAW/source environmental range recovered",
  
  11097, "Cu", "YES", "SOURCE_STAGE_SPLIT",
  "Five source-defined developmental contexts",
  
  11097, "Cd", "YES", "SOURCE_STAGE_SPLIT",
  "Five source-defined developmental contexts",
  
  14474, "Cu", "NO", "LIKELY_REREPORT",
  "Do not count as independent developmental-stage evidence",
  
  6045, "Hg", "NO", "SOURCE_VALUE_HOLD",
  "Hg source LC50 values remain unresolved"
)

stage_source_labels <- tribble(
  ~Reference_Number, ~Test_Number, ~Result_Number,
  ~Stage_label_source, ~Stage_rank_source, ~Stage_Source_Context_ID,
  
  11097, 1135989, 90609, "Nauplius_1d",           2.1, "11097_NAUPLIUS_1D",
  11097, 1135990, 90610, "Nauplius_5d",           2.2, "11097_NAUPLIUS_5D",
  11097, 1135991, 90611, "Copepodid_10d",         3.1, "11097_COPEPODID_10D",
  11097, 1135992, 90612, "Adult_ovigerous_band",  4.1, "11097_OVIG_BAND",
  11097, 1135993, 90613, "Adult_ovigerous_sac",   4.2, "11097_OVIG_SAC",
  
  11097, 1135994, 90926, "Nauplius_1d",           2.1, "11097_NAUPLIUS_1D",
  11097, 1135995, 90927, "Nauplius_5d",           2.2, "11097_NAUPLIUS_5D",
  11097, 1135996, 90928, "Copepodid_10d",         3.1, "11097_COPEPODID_10D",
  11097, 1135997, 90929, "Adult_ovigerous_band",  4.1, "11097_OVIG_BAND",
  11097, 1135998, 90930, "Adult_ovigerous_sac",   4.2, "11097_OVIG_SAC"
)

# Source review 2026-10-07: the O'Brien key below records a common
# protocol, NOT a measured or explicitly restated numeric test temperature.
# Culture temperature stays in the source-review ledger; no 16 C value is
# copied into the main dataset or the temperature display for that source.
stage_env_source <- tribble(
  ~Reference_Number, ~Metal,
  ~Temperature_source_key, ~Salinity_source_key,
  ~Temperature_source_display, ~Salinity_source_display,
  
  2977, "Cu",
  "SOURCE|REF2977|COMMON_PROTOCOL|TEST_TEMP_NOT_RESTATED",
  "SOURCE|REF2977|SAL=35",
  NA_real_, 35,
  
  19281, "Cd",
  "SOURCE|REF19281|TEMP=20C",
  "SOURCE|REF19281|SAL=35",
  20, 35
)

write_csv(
  stage_source_decisions,
  file.path(verification_dir, "ECOTOX_stage_source_decisions.csv")
)

write_csv(
  stage_source_labels,
  file.path(verification_dir, "ECOTOX_stage_source_recovered_labels.csv")
)

write_csv(
  stage_env_source,
  file.path(verification_dir, "ECOTOX_stage_source_environmental_recovery.csv")
)


# ============================================================
# 3. HELPER FUNCTIONS FOR EXACT ENVIRONMENTAL PROVENANCE
# ============================================================

format_env_number <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    format(
      signif(as.numeric(x), digits = 12),
      trim = TRUE,
      scientific = FALSE
    )
  )
}

format_env_unit <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "NA"
  toupper(x)
}


# ============================================================
# 4. RECONSTRUCT ORIGINAL ECOTOX EXACT STAGE PIPELINE
# ============================================================

source_label_missing <- stage_source_labels %>%
  anti_join(
    ecotox,
    by = c("Reference_Number", "Test_Number", "Result_Number")
  )

stopifnot(nrow(source_label_missing) == 0)

stage_prepared_ecotox <- ecotox %>%
  left_join(
    stage_source_decisions,
    by = c("Reference_Number", "Metal")
  ) %>%
  left_join(
    stage_source_labels,
    by = c("Reference_Number", "Test_Number", "Result_Number")
  ) %>%
  left_join(
    stage_env_source,
    by = c("Reference_Number", "Metal")
  ) %>%
  mutate(
    Final_stage_use = coalesce(Final_stage_use, "YES"),
    Stage_verification_status = coalesce(
      Stage_verification_status,
      "METADATA_DERIVED"
    ),
    
    Stage_label_final = coalesce(
      Stage_label_source,
      Lifestage
    ),
    
    Stage_rank_final = coalesce(
      Stage_rank_source,
      case_when(
        Lifestage == "Egg" ~ 1,
        Lifestage == "Nauplii" ~ 2,
        Lifestage %in% c("Copepodite", "Copepodid") ~ 3,
        Lifestage == "Adult" ~ 4,
        TRUE ~ NA_real_
      )
    ),
    
    Developmental_stage_use = case_when(
      !is.na(Stage_label_source) ~ "YES",
      Lifestage %in% c(
        "Egg",
        "Nauplii",
        "Copepodite",
        "Copepodid",
        "Adult"
      ) ~ "YES",
      TRUE ~ "NO"
    ),
    
    Temperature_raw_key = case_when(
      !is.na(Temperature_Raw_Mean) |
        !is.na(Temperature_Raw_Min) |
        !is.na(Temperature_Raw_Max) ~
        paste0(
          "RAW",
          "|MEAN=", format_env_number(Temperature_Raw_Mean),
          "|MIN=", format_env_number(Temperature_Raw_Min),
          "|MAX=", format_env_number(Temperature_Raw_Max),
          "|UNIT=", format_env_unit(Temperature_Raw_Units)
        ),
      TRUE ~ NA_character_
    ),
    
    Salinity_raw_key = case_when(
      !is.na(Salinity_Raw_Mean) |
        !is.na(Salinity_Raw_Min) |
        !is.na(Salinity_Raw_Max) ~
        paste0(
          "RAW",
          "|MEAN=", format_env_number(Salinity_Raw_Mean),
          "|MIN=", format_env_number(Salinity_Raw_Min),
          "|MAX=", format_env_number(Salinity_Raw_Max),
          "|UNIT=", format_env_unit(Salinity_Raw_Units)
        ),
      TRUE ~ NA_character_
    ),
    
    Temperature_match_key = coalesce(
      Temperature_source_key,
      Temperature_raw_key
    ),
    
    Salinity_match_key = coalesce(
      Salinity_source_key,
      Salinity_raw_key
    ),
    
    Temperature_display = coalesce(
      Temperature_source_display,
      Temperature
    ),
    
    Salinity_display = coalesce(
      Salinity_source_display,
      Salinity
    )
  )

stage_exact_rows_ecotox <- stage_prepared_ecotox %>%
  filter(
    Final_stage_use == "YES",
    Developmental_stage_use == "YES",
    !is.na(Stage_label_final),
    !is.na(Stage_rank_final),
    !is.na(Temperature_match_key),
    !is.na(Salinity_match_key)
  )

stage_context_ecotox <- stage_exact_rows_ecotox %>%
  group_by(
    Reference_Number,
    Species,
    Metal,
    Duration_days,
    Temperature_match_key,
    Salinity_match_key,
    Exposure_Type,
    Conc_Type,
    Basis_Interpretation,
    Stage_label_final,
    Stage_rank_final
  ) %>%
  summarise(
    n_results = n(),
    n_tests = n_distinct(Test_Number),
    LC50_stage_context_umol_L =
      exp(mean(log(LC50_umol_L))),
    Temperature_display = first(Temperature_display),
    Salinity_display = first(Salinity_display),
    .groups = "drop"
  )

ecotox_series_key <- c(
  "Reference_Number",
  "Species",
  "Metal",
  "Duration_days",
  "Temperature_match_key",
  "Salinity_match_key",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation"
)

stage_series_ecotox <- stage_context_ecotox %>%
  group_by(across(all_of(ecotox_series_key))) %>%
  summarise(
    n_stages = n_distinct(Stage_label_final),
    Stages = paste(
      Stage_label_final[order(Stage_rank_final)],
      collapse = " | "
    ),
    Temperature = first(Temperature_display),
    Salinity = first(Salinity_display),
    .groups = "drop"
  ) %>%
  filter(n_stages >= 2) %>%
  arrange(Metal, Reference_Number)

ecotox_total <- stage_series_ecotox %>%
  summarise(
    Stage_series = n(),
    References = n_distinct(Reference_Number),
    Species = n_distinct(Species)
  )

stopifnot(
  ecotox_total$Stage_series == 9,
  ecotox_total$References == 6,
  ecotox_total$Species == 5
)

stage_output_ecotox <- stage_context_ecotox %>%
  semi_join(
    stage_series_ecotox %>%
      select(all_of(ecotox_series_key)),
    by = ecotox_series_key
  ) %>%
  mutate(
    Reference_ID = paste0("E_", Reference_Number),
    Source_Origin = "ECOTOX",
    Source_Verification_Status = "SOURCE_AUDITED_OR_EXACT_METADATA"
  )


# ============================================================
# 5. WoS SOURCE-VERIFIED STAGE LEDGER
# ============================================================
# Kadiene et al. (2019):
# - same species
# - same metal
# - same 96-h acute experiment
# - same 26 C
# - same salinity 15
# - nauplii vs copepodids
#
# Source:
# Methods: 25 organisms per beaker, same Cd concentration series,
# triplicate, unfed, 26 C, 96 h.
# Results: 96-h LC50 = 40.3 ug/L Cd for nauplii and 120.4 ug/L Cd
# for copepodids.
#
# Charry et al. (2019) is explicitly not included here because stage
# and duration differ (48-h nauplii vs 96-h adults).

wos_stage_source_ledger <- tribble(
  ~Reference_ID, ~Species, ~Metal, ~Duration_days,
  ~Expected_Temperature, ~Expected_Salinity,
  ~Expected_Stages, ~Final_stage_use,
  ~Stage_verification_status, ~Stage_verification_note,
  
  "W_009",
  "Pseudodiaptomus annandalei",
  "Cd",
  4,
  26,
  15,
  "Nauplii | Copepodid",
  "YES",
  "SOURCE_VERIFIED_SAME_EXPERIMENT",
  paste(
    "Kadiene et al. (2019): nauplii and copepodids were tested",
    "in the same 96-h Cd acute-toxicity experiment at 26 C and salinity 15;",
    "reported LC50s were 40.3 and 120.4 ug Cd/L, respectively."
  )
)

write_csv(
  wos_stage_source_ledger,
  file.path(verification_dir, "WoS_stage_source_verification_ledger.csv")
)


# ============================================================
# 6. EXTRACT + VERIFY W_009
# ============================================================

w009_raw <- wos %>%
  filter(
    Reference_ID == "W_009",
    Species == "Pseudodiaptomus annandalei",
    Metal == "Cd",
    Duration_days == 4
  )

# We expect exactly the two acute stage-specific Results retained
# during source verification.
stopifnot(nrow(w009_raw) == 2)

# Preserve any combined metadata but source-fill temperature/salinity
# if those fields are missing. Contradictory non-missing values fail QA.
if (any(!is.na(w009_raw$Temperature))) {
  stopifnot(all(abs(w009_raw$Temperature[!is.na(w009_raw$Temperature)] - 26) < 0.01))
}
if (any(!is.na(w009_raw$Salinity))) {
  stopifnot(all(abs(w009_raw$Salinity[!is.na(w009_raw$Salinity)] - 15) < 0.01))
}

stage_output_wos <- w009_raw %>%
  mutate(
    Stage_label_final = case_when(
      str_detect(
        str_to_lower(coalesce(External_Stage_Sex, Lifestage, "")),
        "naupl"
      ) ~ "Nauplii",
      
      str_detect(
        str_to_lower(coalesce(External_Stage_Sex, Lifestage, "")),
        "copepod"
      ) ~ "Copepodid",
      
      TRUE ~ NA_character_
    ),
    
    Stage_rank_final = case_when(
      Stage_label_final == "Nauplii" ~ 2,
      Stage_label_final == "Copepodid" ~ 3,
      TRUE ~ NA_real_
    ),
    
    Temperature_display = coalesce(Temperature, 26),
    Salinity_display = coalesce(Salinity, 15),
    
    Temperature_match_key =
      "SOURCE|W_009|KADIENE2019|TEMP=26C",
    
    Salinity_match_key =
      "SOURCE|W_009|KADIENE2019|SAL=15",
    
    Source_Verification_Status =
      "SOURCE_VERIFIED_SAME_EXPERIMENT"
  ) %>%
  filter(!is.na(Stage_label_final)) %>%
  group_by(
    Reference_ID,
    Species,
    Metal,
    Duration_days,
    Temperature_match_key,
    Salinity_match_key,
    Stage_label_final,
    Stage_rank_final,
    Source_Origin,
    Source_Verification_Status
  ) %>%
  summarise(
    n_results = n(),
    n_tests = n_distinct(Test_ID),
    LC50_stage_context_umol_L =
      exp(mean(log(LC50_umol_L))),
    Temperature_display = first(Temperature_display),
    Salinity_display = first(Salinity_display),
    .groups = "drop"
  )

stopifnot(
  nrow(stage_output_wos) == 2,
  setequal(stage_output_wos$Stage_label_final, c("Nauplii", "Copepodid"))
)

# Source-value regression check after molar harmonization.
w009_value_check <- stage_output_wos %>%
  select(Stage_label_final, LC50_stage_context_umol_L) %>%
  pivot_wider(
    names_from = Stage_label_final,
    values_from = LC50_stage_context_umol_L
  )

stopifnot(
  abs(w009_value_check$Nauplii - (40.3 / 112.414)) < 0.002,
  abs(w009_value_check$Copepodid - (120.4 / 112.414)) < 0.002
)


# ============================================================
# 7. STANDARDIZE ECOTOX + WoS OUTPUT STRUCTURE
# ============================================================

stage_output_ecotox_std <- stage_output_ecotox %>%
  transmute(
    Reference_ID,
    Reference_Number = as.character(Reference_Number),
    Species,
    Metal,
    Duration_days,
    Temperature_match_key,
    Salinity_match_key,
    Stage_label_final,
    Stage_rank_final,
    n_results,
    n_tests,
    LC50_stage_context_umol_L,
    Temperature_display,
    Salinity_display,
    Source_Origin,
    Source_Verification_Status
  )

stage_output_wos_std <- stage_output_wos %>%
  transmute(
    Reference_ID,
    Reference_Number = NA_character_,
    Species,
    Metal,
    Duration_days,
    Temperature_match_key,
    Salinity_match_key,
    Stage_label_final,
    Stage_rank_final,
    n_results,
    n_tests,
    LC50_stage_context_umol_L,
    Temperature_display,
    Salinity_display,
    Source_Origin,
    Source_Verification_Status
  )

stage_output_data <- bind_rows(
  stage_output_ecotox_std,
  stage_output_wos_std
)


# ============================================================
# 8. IDENTIFY FINAL COMBINED STAGE SERIES
# ============================================================

combined_series_key <- c(
  "Reference_ID",
  "Species",
  "Metal",
  "Duration_days",
  "Temperature_match_key",
  "Salinity_match_key"
)

stage_series_final <- stage_output_data %>%
  group_by(across(all_of(combined_series_key))) %>%
  summarise(
    n_stages = n_distinct(Stage_label_final),
    Stages = paste(
      Stage_label_final[order(Stage_rank_final)],
      collapse = " | "
    ),
    Temperature = first(Temperature_display),
    Salinity = first(Salinity_display),
    Source_Origin = first(Source_Origin),
    Source_Verification_Status = first(Source_Verification_Status),
    .groups = "drop"
  ) %>%
  filter(n_stages >= 2) %>%
  arrange(Metal, Reference_ID)

stage_output_data <- stage_output_data %>%
  semi_join(
    stage_series_final %>%
      select(all_of(combined_series_key)),
    by = combined_series_key
  ) %>%
  group_by(across(all_of(combined_series_key))) %>%
  mutate(Series_ID = cur_group_id()) %>%
  ungroup()

combined_total <- stage_series_final %>%
  summarise(
    Stage_series = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species)
  )

stopifnot(
  combined_total$Stage_series == 10,
  combined_total$References == 7,
  combined_total$Species == 6,
  n_distinct(stage_output_data$Series_ID) == 10
)


# ============================================================
# 9. SUPPORT + SOURCE CONTRIBUTION
# ============================================================

stage_support_combined <- stage_series_final %>%
  group_by(Metal) %>%
  summarise(
    Stage_series = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(desc(Stage_series), Metal)

stage_source_contribution <- stage_series_final %>%
  group_by(Source_Origin) %>%
  summarise(
    Stage_series = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

stage_ECOTOX_vs_combined <- tibble(
  Evidence_Set = c(
    "Original ECOTOX exact/source-audited",
    "Combined ECOTOX + WoS source-verified"
  ),
  Stage_series = c(
    ecotox_total$Stage_series,
    combined_total$Stage_series
  ),
  References = c(
    ecotox_total$References,
    combined_total$References
  ),
  Species = c(
    ecotox_total$Species,
    combined_total$Species
  )
)


# ============================================================
# 10. CLEAN DISPLAY LABELS
# ============================================================

stage_output_data <- stage_output_data %>%
  mutate(
    Stage_display = case_when(
      Stage_label_final == "Egg" ~ "Egg",
      Stage_label_final == "Nauplii" ~ "Nauplii",
      Stage_label_final == "Nauplius_1d" ~ "Nauplius 1 d",
      Stage_label_final == "Nauplius_5d" ~ "Nauplius 5 d",
      Stage_label_final == "Copepodite" ~ "Copepodite",
      Stage_label_final == "Copepodid" ~ "Copepodid",
      Stage_label_final == "Copepodid_10d" ~ "Copepodid 10 d",
      Stage_label_final == "Adult" ~ "Adult",
      Stage_label_final == "Adult_ovigerous_band" ~ "Adult ovigerous band",
      Stage_label_final == "Adult_ovigerous_sac" ~ "Adult ovigerous sac",
      TRUE ~ str_replace_all(Stage_label_final, "_", " ")
    )
  )


# ============================================================
# 11. FIRST-TO-LAST DIRECTION SUMMARY
# ============================================================
# Descriptive triangulation only; not a pooled causal stage model.

stage_direction_summary <- stage_output_data %>%
  arrange(Series_ID, Stage_rank_final) %>%
  group_by(
    Series_ID,
    Reference_ID,
    Species,
    Metal,
    Duration_days,
    Source_Origin
  ) %>%
  summarise(
    earliest_stage = first(Stage_display),
    latest_stage = last(Stage_display),
    LC50_earliest = first(LC50_stage_context_umol_L),
    LC50_latest = last(LC50_stage_context_umol_L),
    ratio_latest_earliest = LC50_latest / LC50_earliest,
    ln_ratio_latest_earliest = log(ratio_latest_earliest),
    direction = case_when(
      ratio_latest_earliest > 1 ~
        "Later stage higher LC50 / earlier stage more sensitive",
      ratio_latest_earliest < 1 ~
        "Later stage lower LC50 / later stage more sensitive",
      TRUE ~ "No change"
    ),
    .groups = "drop"
  )

wos_new_stage_direction <- stage_direction_summary %>%
  filter(Source_Origin == "WoS_supplemental")

stopifnot(
  nrow(wos_new_stage_direction) == 1,
  wos_new_stage_direction$ratio_latest_earliest > 2.9,
  wos_new_stage_direction$ratio_latest_earliest < 3.1
)



# ============================================================
# 11A. ALL STAGE PAIRS WITHIN THE EXISTING ACCEPTED SERIES
# ============================================================
# Preserve source-defined ages and reproductive groups. Only the synonyms
# Copepodite and Copepodid share the comparison label Copepodid.
# A and B follow the original display rank; this does NOT establish that
# ovigerous band/sac are successive developmental stages.
# Ratio B/A > 1 means A has the lower LC50, not a mortality-percent ratio.
# Source verification is inherited, not upgraded by this calculation.
# Exact-metadata matching is not automatically full-text verification.
# No comparison crosses Series_ID, species, metal or exposure duration.
# Within a metal/pair/scope, average log ratios per reference first.
# References then have equal weight. No cross-metal pooling is performed.
# t intervals assume independent, approximately normal reference log ratios.
# Exact sign tests address direction among non-ties, not mean magnitude.
# Small reference counts remain limited evidence even if a CI excludes 1.
# Sources: R stats t.test / binom.test / p.adjust; Holm (1979), SJS 6:65-70.

stage_pair_list <- list()
for (sid in unique(stage_output_data$Series_ID)) {
  d <- stage_output_data %>% filter(Series_ID == sid) %>% arrange(Stage_rank_final)
  stopifnot(n_distinct(d$Reference_ID) == 1L, n_distinct(d$Species) == 1L,
            n_distinct(d$Metal) == 1L, n_distinct(d$Duration_days) == 1L,
            !anyDuplicated(d$Stage_label_final),
            all(is.finite(d$LC50_stage_context_umol_L)),
            all(d$LC50_stage_context_umol_L > 0))
  if (nrow(d) < 2L) next
  for (ij in combn(seq_len(nrow(d)), 2, simplify = FALSE)) {
    i <- ij[1]; j <- ij[2]
    a <- d$Stage_label_final[i]; z <- d$Stage_label_final[j]
    if (a == "Copepodite") a <- "Copepodid"
    if (z == "Copepodite") z <- "Copepodid"
    stopifnot(a != z)
    stage_pair_list[[length(stage_pair_list) + 1L]] <- tibble(
      Series_ID = sid, Reference_ID = d$Reference_ID[1], Species = d$Species[1],
      Metal = d$Metal[1], Duration_days = d$Duration_days[1],
      Source_Origin = d$Source_Origin[1],
      Source_Verification_Status = d$Source_Verification_Status[1],
      Stage_A_original = d$Stage_label_final[i], Stage_B_original = d$Stage_label_final[j],
      Stage_A = a, Stage_B = z, Comparison = paste(z, "vs", a),
      LC50_A = d$LC50_stage_context_umol_L[i], LC50_B = d$LC50_stage_context_umol_L[j],
      Log_ratio_B_over_A = log(LC50_B / LC50_A), Ratio_B_over_A = LC50_B / LC50_A)
  }
}
stage_pair_contexts <- bind_rows(stage_pair_list)
expected_stage_pairs <- stage_output_data %>% count(Series_ID) %>%
  summarise(total = sum(n * (n - 1) / 2)) %>% pull(total)
stopifnot(nrow(stage_pair_contexts) == expected_stage_pairs)
stage_pair_specs <- stage_pair_contexts %>% distinct(Stage_A, Stage_B, Comparison)
stage_pair_references <- list(); stage_pair_summaries <- list(); stage_pair_loro_runs <- list()
for (scope in c("ALL_DURATIONS", "96H")) {
  scope_data <- stage_pair_contexts
  if (scope == "96H") scope_data <- filter(scope_data, Duration_days == 4)
  for (metal in sort(unique(stage_output_data$Metal))) {
    for (row in seq_len(nrow(stage_pair_specs))) {
      spec <- stage_pair_specs[row, ]
      d <- scope_data %>% filter(Metal == metal, Comparison == spec$Comparison)
      r <- d %>% group_by(Reference_ID) %>% summarise(
        Series = n(), Mean_log_ratio = mean(Log_ratio_B_over_A),
        Ratio_B_over_A = exp(Mean_log_ratio), .groups = "drop") %>%
        mutate(Scope = scope, Metal = metal, Comparison = spec$Comparison, .before = 1)
      stage_pair_references[[length(stage_pair_references) + 1L]] <- r
      x <- r$Mean_log_ratio; k <- length(x)
      ratio <- if (k) exp(mean(x)) else NA_real_
      low <- high <- NA_real_
      status <- if (k < 2L) "FEWER_THAN_TWO_REFERENCES" else "ZERO_VARIANCE"
      if (k >= 2L && is.finite(sd(x)) && sd(x) > 0) {
        margin <- qt(0.975, k - 1L) * sd(x) / sqrt(k)
        low <- exp(mean(x) - margin); high <- exp(mean(x) + margin)
        status <- "POINTWISE_T_INTERVAL_ASSUMPTION_DEPENDENT"
      }
      a_lower <- sum(x > 1e-12); b_lower <- sum(x < -1e-12)
      ties <- k - a_lower - b_lower
      sign_p <- if (k >= 2L && a_lower + b_lower > 0L)
        binom.test(a_lower, a_lower + b_lower, p = 0.5)$p.value else NA_real_
      omit_min <- omit_max <- NA_real_
      if (k >= 2L) {
        o <- vapply(seq_len(k), function(j) exp(mean(x[-j])), numeric(1))
        omit_min <- min(o); omit_max <- max(o)
        stage_pair_loro_runs[[length(stage_pair_loro_runs) + 1L]] <- tibble(
          Scope = scope, Metal = metal, Comparison = spec$Comparison,
          Omitted_reference = r$Reference_ID, Remaining_references = k - 1L,
          Ratio_B_over_A = o)
      }
      stage_pair_summaries[[length(stage_pair_summaries) + 1L]] <- tibble(
        Scope = scope, Metal = metal, Stage_A = spec$Stage_A, Stage_B = spec$Stage_B,
        Comparison = spec$Comparison, Matched_series = nrow(d), References = k,
        Ratio_B_over_A = ratio, CI95_low = low, CI95_high = high, Interval_status = status,
        A_lower_references = a_lower, B_lower_references = b_lower, Tied_references = ties,
        Sign_p = sign_p, LORO_min = omit_min,
        LORO_max = omit_max,
        Support_note = if (k == 0L) "NO_MATCHED_PAIR" else
          if (k == 1L) "ONE_REFERENCE_DESCRIPTIVE_ONLY" else "LIMITED_REFERENCE_LEVEL_EVIDENCE")
    }
  }
}
stage_pair_reference_effects <- bind_rows(stage_pair_references)
# One test family per scope: every metal x stage-pair combination in this
# inventory. Untested combinations remain NA; no selection by p value.
stage_pair_effect_summary <- bind_rows(stage_pair_summaries) %>%
  group_by(Scope) %>% mutate(Sign_p_Holm = p.adjust(Sign_p, "holm", n = n())) %>% ungroup()
stage_pair_loro <- bind_rows(stage_pair_loro_runs)
write_csv(stage_pair_contexts, file.path(table_dir, "stage_pair_contexts.csv"))
write_csv(stage_pair_reference_effects, file.path(table_dir, "stage_pair_reference_effects.csv"))
write_csv(stage_pair_effect_summary, file.path(table_dir, "stage_pair_effect_summary.csv"))
write_csv(stage_pair_loro, file.path(table_dir, "stage_pair_loro.csv"))
cat("\n--- STAGE PAIR EFFECTS: B/A LC50 RATIOS ---\n")
print(stage_pair_effect_summary %>% filter(References > 0), n = Inf)

# ============================================================
# 12. DETAILED ONE-ROW-PER-SERIES TABLE
# ============================================================

stage_table_data <- stage_output_data %>%
  arrange(Series_ID, Stage_rank_final) %>%
  group_by(
    Series_ID,
    Reference_ID,
    Species,
    Metal,
    Duration_days,
    Source_Origin
  ) %>%
  summarise(
    Temperature = first(Temperature_display),
    Salinity = first(Salinity_display),
    n_stages = n_distinct(Stage_label_final),
    Stage_profile = paste0(
      Stage_display,
      ": ",
      format(
        signif(LC50_stage_context_umol_L, 4),
        trim = TRUE,
        scientific = FALSE
      ),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Duration_h = Duration_days * 24
  ) %>%
  arrange(Metal, Reference_ID)


# ============================================================
# 13. VISUAL TABLE
# ============================================================

stage_visual_table <- stage_table_data %>%
  select(
    Reference_ID,
    Source_Origin,
    Metal,
    Species,
    Duration_h,
    Temperature,
    Salinity,
    Stage_profile
  ) %>%
  gt() %>%
  tab_header(
    title = "Exact within-context developmental-stage series",
    subtitle = "10 source-audited series from 7 References"
  ) %>%
  cols_label(
    Reference_ID = "Reference",
    Source_Origin = "Source",
    Metal = "Metal",
    Species = "Species",
    Duration_h = "Duration (h)",
    Temperature = "Temperature (C)",
    Salinity = "Salinity",
    Stage_profile = "Developmental stage: LC50 (umol/L)"
  ) %>%
  fmt_number(
    columns = c(Temperature, Salinity),
    decimals = 2
  ) %>%
  cols_align(
    align = "center",
    columns = c(
      Reference_ID,
      Source_Origin,
      Metal,
      Duration_h,
      Temperature,
      Salinity
    )
  ) %>%
  cols_align(
    align = "left",
    columns = c(Species, Stage_profile)
  ) %>%
  opt_row_striping()

gtsave(
  stage_visual_table,
  file.path(table_dir, "stage_within_combined_table.html")
)


# ============================================================
# 14. QA / INTERMEDIATE DEVELOPMENTAL-STAGE LC50 FIGURE
# ============================================================
# Thesis-facing figures are produced later by the figure-suite script.
# Preserve the original design principle:
# - one panel = one exact comparable stage context
# - no pooled regression
# - no connecting line (stages are comparable source-defined groups,
#   not repeated measurements of the same individual organisms)

stage_plot_data <- stage_output_data %>%
  mutate(
    Series_panel = paste0(
      Reference_ID,
      " | ",
      Metal,
      "\n",
      Species
    )
  )

stage_within_plot <- ggplot(
  stage_plot_data,
  aes(
    x = Stage_rank_final,
    y = LC50_stage_context_umol_L
  )
) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = Stage_display),
    hjust = -0.05,
    vjust = -0.7,
    size = 2.8,
    check_overlap = TRUE
  ) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c(
      "Egg",
      "Nauplii",
      "Copepodid",
      "Adult"
    )
  ) +
  facet_wrap(
    ~ Series_panel,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    x = "Developmental stage",
    y = expression(LC[50] ~ (mu * mol ~ L^{-1}))
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 8.7
    )
  )

ggsave(
  file.path(figure_dir, "stage_within_combined_plot.png"),
  stage_within_plot,
  width = 11,
  height = 9,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(figure_dir, "stage_within_combined_plot.pdf"),
  stage_within_plot,
  width = 11,
  height = 9
)


# ============================================================
# 15. WRITE NUMERICAL OUTPUTS
# ============================================================

write_csv(
  stage_support_combined,
  file.path(table_dir, "stage_within_combined_support.csv")
)

write_csv(
  stage_source_contribution,
  file.path(table_dir, "stage_within_source_contribution.csv")
)

write_csv(
  stage_ECOTOX_vs_combined,
  file.path(table_dir, "stage_within_ECOTOX_vs_combined_summary.csv")
)

write_csv(
  stage_series_final,
  file.path(table_dir, "stage_within_combined_series.csv")
)

write_csv(
  stage_table_data,
  file.path(table_dir, "stage_within_combined_table_data.csv")
)

write_csv(
  stage_direction_summary,
  file.path(table_dir, "stage_within_first_to_last_direction_combined.csv")
)

write_csv(
  wos_new_stage_direction,
  file.path(table_dir, "new_WoS_stage_within_direction.csv")
)

write_csv(
  stage_output_data,
  file.path(table_dir, "stage_within_combined_plot_data.csv")
)


# ============================================================
# 16. SAVE ANALYSIS OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    ECOTOX_source_decisions = stage_source_decisions,
    ECOTOX_source_labels = stage_source_labels,
    ECOTOX_environmental_recovery = stage_env_source,
    WoS_stage_source_ledger = wos_stage_source_ledger,
    stage_series_final = stage_series_final,
    stage_output_data = stage_output_data,
    stage_table_data = stage_table_data,
    stage_direction_summary = stage_direction_summary,
    stage_pair_contexts = stage_pair_contexts,
    stage_pair_effect_summary = stage_pair_effect_summary,
    stage_pair_reference_effects = stage_pair_reference_effects,
    stage_pair_loro = stage_pair_loro,
    stage_within_plot = stage_within_plot
  ),
  file.path(object_dir, "stage_within_reference_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)


# ============================================================
# 17. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("09 WITHIN-REFERENCE DEVELOPMENTAL-STAGE COMPARISONS\n")
cat("============================================================\n")

cat("\n--- ORIGINAL ECOTOX ---\n")
print(ecotox_total)

cat("\n--- COMBINED ---\n")
print(combined_total)

cat("\n--- SUPPORT BY METAL ---\n")
print(stage_support_combined, n = Inf)

cat("\n--- SOURCE CONTRIBUTION ---\n")
print(stage_source_contribution)

cat("\n--- NEW WoS DIRECT STAGE SERIES ---\n")
print(wos_new_stage_direction, width = Inf)

cat(
  "\nInterpretive note:\n",
  "The original 9 exact ECOTOX developmental-stage series are preserved.\n",
  "One new source-verified WoS series is added: Kadiene et al. (2019),\n",
  "Pseudodiaptomus annandalei, Cd, 96 h, Nauplii versus Copepodid.\n",
  "The later-stage LC50 is approximately three times the Nauplii LC50,\n",
  "supporting greater acute sensitivity of the earlier stage in this context.\n",
  "This remains descriptive/direct evidence and does not replace the broader\n",
  "96-h Cu-Cd Adult/Nauplii interaction model evaluated in the later stage-support and stage-model scripts.\n",
  sep = ""
)

required_outputs <- c(
  file.path(
    verification_dir,
    "WoS_stage_source_verification_ledger.csv"
  ),
  file.path(
    table_dir,
    "stage_within_combined_support.csv"
  ),
  file.path(
    table_dir,
    "stage_within_combined_series.csv"
  ),
  file.path(
    table_dir,
    "new_WoS_stage_within_direction.csv"
  ),
  file.path(
    figure_dir,
    "stage_within_combined_plot.png"
  ),
  file.path(
    object_dir,
    "stage_within_reference_objects.rds"
  )
)

stopifnot(all(file.exists(required_outputs)))

# Reproducibility manifest and completion marker
manifest <- tibble(
  Role = c(
    "Input: canonical combined dataset",
    "Output: combined stage series",
    "Output: stage-pair effect summary",
    "Output: stage-pair LORO",
    "Output: first-to-last direction summary",
    "Output: QA stage figure",
    "Output: R analysis objects"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "stage_within_combined_series.csv"),
    file.path(table_dir, "stage_pair_effect_summary.csv"),
    file.path(table_dir, "stage_pair_loro.csv"),
    file.path(table_dir, "stage_within_first_to_last_direction_combined.csv"),
    file.path(figure_dir, "stage_within_combined_plot.png"),
    file.path(object_dir, "stage_within_reference_objects.rds")
  )
) %>%
  mutate(
    Exists = file.exists(Path),
    MD5 = ifelse(Exists, unname(tools::md5sum(Path)), NA_character_)
  )

write_csv(manifest, file.path(output_root, "input_output_manifest.csv"))

writeLines(
  c(
    "STATUS: 09 COMPARE STAGES WITHIN REFERENCES = PASS",
    paste0("Stage series: ", combined_total$Stage_series),
    paste0("References: ", combined_total$References),
    paste0("Species: ", combined_total$Species)
  ),
  file.path(output_root, "RUN_COMPLETE.txt")
)

cat("\nSTATUS: 09 COMPARE STAGES WITHIN REFERENCES = PASS\n")
# ============================================================
# END
# ============================================================
