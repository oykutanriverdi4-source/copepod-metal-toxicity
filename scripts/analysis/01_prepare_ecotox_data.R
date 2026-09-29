# ============================================================
# ECOTOX THESIS REPRODUCIBILITY ARCHIVE
# 01_prepare_ecotox_data.R
# ============================================================
# STATUS: CORE DATA-PREPARATION SCRIPT
#
# PURPOSE
# Reconstruct and verify the frozen 295-row ECOTOX baseline, then apply
# the documented Sosnowski (1979) source-review correction to produce
# the 296-row ECOTOX dataset used before Web of Science integration.
#
# INPUTS
#   data/raw/ecotox/ecotox_raw_export.xlsx
#   data/curated/ecotox/ecotox_working_master_harmonized.xlsx
#
# PRIMARY OUTPUT
#   data/processed/ecotox_harmonized.csv
#
# QA / PROVENANCE OUTPUTS
#   outputs/01_prepare_ecotox_data/
#
# IMPORTANT REPRODUCIBILITY PRINCIPLE
# R does not rediscover decisions that required reading original sources.
# Source-dependent decisions (eligibility, recovered point estimates,
# concentration-basis interpretations, unresolved Cr basis, etc.) are
# taken from the curated Working Master ledger and are then applied and
# checked deterministically by this script.
#
# PIPELINE
# ECOTOX raw export
#   -> programmatic copepod taxonomy screen
#   -> verify against documented 410-row screening ledger
#   -> apply documented eligibility/source decisions
#   -> reconstruct frozen 295-row Primary baseline
#   -> recompute molar harmonization from documented basis rules
#   -> recover selected environmental means from raw min/max values
#   -> apply documented Sosnowski (1979) correction
#   -> write 296-row ECOTOX dataset for downstream WoS integration
#
# THIS SCRIPT DOES NOT
# - perform the Web of Science search or screening
# - select the 96-h modelling domain
# - fit toxicity models
# ============================================================


# ============================================================
# 0. PACKAGES, INPUTS, AND OUTPUT LOCATIONS
# ============================================================

library(readxl)
library(readr)
library(dplyr)
library(tibble)

options(width = 200)

# This script uses paths relative to the repository root.
# Do not place a user-specific setwd() call in this file.
raw_file <- file.path(
  "data", "raw", "ecotox",
  "ecotox_raw_export.xlsx"
)

working_master_file <- file.path(
  "data", "curated", "ecotox",
  "ecotox_working_master_harmonized.xlsx"
)

processed_data_dir <- file.path("data", "processed")
canonical_output_file <- file.path(
  processed_data_dir,
  "ecotox_harmonized.csv"
)

output_dir <- file.path(
  "outputs",
  "01_prepare_ecotox_data"
)

qa_dir <- file.path(output_dir, "01_QA")
ledger_dir <- file.path(output_dir, "02_documented_decisions")
lookup_dir <- file.path(output_dir, "03_metadata_lookups")

for (d in c(processed_data_dir, output_dir, qa_dir, ledger_dir, lookup_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Fail early with an informative message if the repository is not being
# run from its root directory or the archived inputs are missing.
required_inputs <- c(raw_file, working_master_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    paste0(
      "Required input file(s) not found:\n  - ",
      paste(missing_inputs, collapse = "\n  - "),
      "\n\nRun this script from the repository root and confirm the data/ folder structure."
    )
  )
}

# Preserve a pre-existing processed output as an optional rerun regression
# reference before writing a new copy.
existing_target <- NULL

if (file.exists(canonical_output_file)) {
  existing_target <- read_csv(
    canonical_output_file,
    show_col_types = FALSE
  )
}

# The documented Sosnowski correction is applied after reconstruction of
# the frozen 295-row baseline, so compare that baseline separately.
if (!is.null(existing_target)) {
  existing_target <- existing_target %>%
    filter(as.character(Result_Number) != "42366")
}

cat("\n============================================\n")
cat("ECOTOX DATA PREPARATION + HARMONIZATION\n")
cat("============================================\n")
cat("RAW file:", raw_file, "\n")
cat("Working Master:", working_master_file, "\n")


# ============================================================
# 1. READ IMMUTABLE RAW EXPORT
# ============================================================

raw <- suppressWarnings(
  read_excel(
    raw_file,
    sheet = "Aquatic-Export",
    col_types = "text",
    .name_repair = "minimal"
  )
) %>%
  mutate(
    # Excel row number, including the header row at row 1.
    Raw_Row = row_number() + 1L,
    .before = 1
  )

# Frozen raw benchmark.
stopifnot(nrow(raw) == 2299)

required_raw_columns <- c(
  "Test Number",
  "Result Number",
  "Reference Number",
  "Chemical Analysis",
  "Species Order",
  "Title",
  "Temperature Mean",
  "Temperature Min",
  "Temperature Max",
  "Temperature Units",
  "Salinity Mean",
  "Salinity Min",
  "Salinity Max",
  "Salinity Units",
  "pH Mean",
  "pH Min",
  "pH Max",
  "Organic Carbon Mean",
  "Organic Carbon Min",
  "Organic Carbon Max",
  "Organic Carbon Units"
)

stopifnot(
  all(required_raw_columns %in% names(raw))
)

raw_lookup <- raw %>%
  transmute(
    Raw_Row,
    Test_Number = `Test Number`,
    Result_Number = `Result Number`,
    Reference_Number = `Reference Number`,
    Order_Raw = `Species Order`,
    Chemical_Analysis = `Chemical Analysis`,
    Title_Raw = Title,
    
    Temperature_Raw_Mean = `Temperature Mean`,
    Temperature_Raw_Min = `Temperature Min`,
    Temperature_Raw_Max = `Temperature Max`,
    Temperature_Raw_Units = `Temperature Units`,
    
    Salinity_Raw_Mean = `Salinity Mean`,
    Salinity_Raw_Min = `Salinity Min`,
    Salinity_Raw_Max = `Salinity Max`,
    Salinity_Raw_Units = `Salinity Units`,
    
    pH_Raw_Mean = `pH Mean`,
    pH_Raw_Min = `pH Min`,
    pH_Raw_Max = `pH Max`,
    # The archived ECOTOX export has no separate pH-units field.
    pH_Raw_Units = NA_character_,
    
    Organic_Carbon_Raw_Mean = `Organic Carbon Mean`,
    Organic_Carbon_Raw_Min = `Organic Carbon Min`,
    Organic_Carbon_Raw_Max = `Organic Carbon Max`,
    Organic_Carbon_Raw_Units = `Organic Carbon Units`
  )

make_record_key <- function(test, result, reference) {
  paste(
    as.character(test),
    as.character(result),
    as.character(reference),
    sep = "|"
  )
}

# RAW was deliberately imported as text so mixed-type ECOTOX columns
# cannot be silently coerced. Convert only the numeric metadata fields
# that are actually used downstream.
parse_raw_numeric <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  suppressWarnings(as.numeric(x))
}

raw_lookup <- raw_lookup %>%
  mutate(
    across(
      c(
        Temperature_Raw_Mean, Temperature_Raw_Min, Temperature_Raw_Max,
        Salinity_Raw_Mean, Salinity_Raw_Min, Salinity_Raw_Max,
        pH_Raw_Mean, pH_Raw_Min, pH_Raw_Max,
        Organic_Carbon_Raw_Mean, Organic_Carbon_Raw_Min, Organic_Carbon_Raw_Max
      ),
      parse_raw_numeric
    ),
    record_key = make_record_key(
      Test_Number,
      Result_Number,
      Reference_Number
    )
  )

stopifnot(
  n_distinct(raw_lookup$record_key) == nrow(raw_lookup)
)


# ============================================================
# 2. REPRODUCE TAXONOMIC COPEPOD SCREEN FROM RAW
# ============================================================

copepod_orders <- c(
  "Calanoida",
  "Harpacticoida",
  "Cyclopoida"
)

raw_copepods <- raw_lookup %>%
  filter(Order_Raw %in% copepod_orders)

stopifnot(
  nrow(raw_copepods) == 410,
  n_distinct(raw_copepods$Test_Number) == 401,
  n_distinct(raw_copepods$Reference_Number) == 72,
  n_distinct(raw_copepods$Result_Number) == 410
)

raw_screen_summary <- tibble(
  Metric = c(
    "Raw ECOTOX Results",
    "Definite copepod Results",
    "Definite copepod Tests",
    "Definite copepod References"
  ),
  Value = c(
    nrow(raw),
    nrow(raw_copepods),
    n_distinct(raw_copepods$Test_Number),
    n_distinct(raw_copepods$Reference_Number)
  )
)

write_csv(
  raw_screen_summary,
  file.path(qa_dir, "00_raw_and_copepod_screen_summary.csv")
)


# ============================================================
# 3. LOAD DOCUMENTED 410-ROW SCREENING / SOURCE-DECISION LEDGER
# ============================================================

screening_ledger <- read_excel(
  working_master_file,
  sheet = "Screening_410"
) %>%
  mutate(
    record_key = make_record_key(
      Test_Number,
      Result_Number,
      Reference_Number
    )
  )

stopifnot(
  nrow(screening_ledger) == 410,
  n_distinct(screening_ledger$record_key) == 410
)

# The documented ledger must refer to exactly the same 410 raw
# copepod Result records selected programmatically above.
raw_only_keys <- setdiff(
  raw_copepods$record_key,
  screening_ledger$record_key
)

ledger_only_keys <- setdiff(
  screening_ledger$record_key,
  raw_copepods$record_key
)

key_set_audit <- tibble(
  Audit = c(
    "Raw copepod keys absent from ledger",
    "Ledger keys absent from raw copepod screen"
  ),
  Count = c(
    length(raw_only_keys),
    length(ledger_only_keys)
  )
)

write_csv(
  key_set_audit,
  file.path(qa_dir, "01_screening_key_set_audit.csv")
)

if (length(raw_only_keys) > 0 || length(ledger_only_keys) > 0) {
  write_csv(
    tibble(record_key = raw_only_keys),
    file.path(qa_dir, "01a_raw_only_keys.csv")
  )
  write_csv(
    tibble(record_key = ledger_only_keys),
    file.path(qa_dir, "01b_ledger_only_keys.csv")
  )
  stop("The documented Screening_410 ledger does not match the raw 410-row copepod screen.")
}

# Check that row identity and taxonomy also agree.
ledger_raw_identity_audit <- screening_ledger %>%
  select(
    record_key,
    Ledger_Raw_Row = Raw_Row,
    Ledger_Order = Order
  ) %>%
  left_join(
    raw_copepods %>%
      select(
        record_key,
        Raw_Raw_Row = Raw_Row,
        Raw_Order = Order_Raw
      ),
    by = "record_key"
  ) %>%
  mutate(
    Raw_Row_match = Ledger_Raw_Row == Raw_Raw_Row,
    Order_match = Ledger_Order == Raw_Order
  )

stopifnot(
  all(ledger_raw_identity_audit$Raw_Row_match),
  all(ledger_raw_identity_audit$Order_match)
)


# ============================================================
# 4. APPLY DOCUMENTED ELIGIBILITY / SOURCE DECISIONS
# ============================================================

screening_disposition <- screening_ledger %>%
  count(
    Eligibility_Status,
    name = "Results"
  ) %>%
  arrange(Eligibility_Status)

expected_disposition <- tibble(
  Eligibility_Status = c(
    "PRIMARY_KEEP",
    "PRIMARY_RECOVERED",
    "SECONDARY_MODIFIER",
    "EXCLUDE_PRIMARY",
    "PENDING_PRIMARY_RECOVERY"
  ),
  Expected_Results = c(
    284L,
    11L,
    56L,
    41L,
    18L
  )
)

screening_disposition_check <- expected_disposition %>%
  left_join(
    screening_disposition,
    by = "Eligibility_Status"
  ) %>%
  mutate(
    Results = if_else(is.na(Results), 0L, as.integer(Results)),
    Match = Results == Expected_Results
  )

write_csv(
  screening_disposition_check,
  file.path(qa_dir, "02_screening_disposition_summary.csv")
)

stopifnot(
  all(screening_disposition_check$Match)
)

primary_295 <- screening_ledger %>%
  filter(
    Eligibility_Status %in% c(
      "PRIMARY_KEEP",
      "PRIMARY_RECOVERED"
    )
  ) %>%
  arrange(Raw_Row)

stopifnot(
  nrow(primary_295) == 295,
  n_distinct(primary_295$Test_Number) == 288,
  n_distinct(primary_295$Reference_Number) == 66,
  sum(primary_295$Eligibility_Status == "PRIMARY_RECOVERED") == 11
)

primary_summary <- tibble(
  Metric = c(
    "Primary Results",
    "Primary Tests",
    "Primary References",
    "Direct Primary Results",
    "Source-recovered Primary Results"
  ),
  Value = c(
    nrow(primary_295),
    n_distinct(primary_295$Test_Number),
    n_distinct(primary_295$Reference_Number),
    sum(primary_295$Eligibility_Status == "PRIMARY_KEEP"),
    sum(primary_295$Eligibility_Status == "PRIMARY_RECOVERED")
  )
)

write_csv(
  primary_summary,
  file.path(qa_dir, "03_primary_295_summary.csv")
)

# Export the source/eligibility ledger that R actually applied.
write_csv(
  screening_ledger %>%
    select(
      Raw_Row,
      Test_Number,
      Result_Number,
      Reference_Number,
      Metal,
      Species,
      Order,
      Lifestage,
      Analysis_Point,
      Analysis_Units,
      Eligibility_Status,
      Eligibility_Reason
    ),
  file.path(
    ledger_dir,
    "screening_and_source_decision_ledger_410.csv"
  )
)


# ============================================================
# 5. LOAD DOCUMENTED HARMONIZATION-BASIS DECISIONS
# ============================================================
#
# Numeric molar values are recomputed below.
# This sheet is used only for documented source/basis decisions,
# metal roles and an independent regression target for the numeric
# harmonization calculation.
# ============================================================

harmonization_ledger_full <- read_excel(
  working_master_file,
  sheet = "Harmonized_295"
) %>%
  mutate(
    record_key = make_record_key(
      Test_Number,
      Result_Number,
      Reference_Number
    )
  )

stopifnot(
  nrow(harmonization_ledger_full) == 295,
  n_distinct(harmonization_ledger_full$record_key) == 295,
  setequal(
    primary_295$record_key,
    harmonization_ledger_full$record_key
  )
)

harmonization_decisions <- harmonization_ledger_full %>%
  select(
    record_key,
    Atomic_Weight_g_mol,
    Basis_Interpretation,
    Harmonization_Status,
    Harmonization_Rule,
    Metal_Role,
    Core_Pairwise_Ready,
    Ledger_LC50_target_mg_L = LC50_target_mg_L,
    Ledger_LC50_umol_L = LC50_umol_L,
    Ledger_ln_LC50_umol_L = ln_LC50_umol_L
  )

write_csv(
  harmonization_decisions %>%
    select(
      -starts_with("Ledger_")
    ),
  file.path(
    ledger_dir,
    "harmonization_basis_decision_ledger_295.csv"
  )
)


# ============================================================
# 6. RECOMPUTE TARGET-METAL MOLAR LC50
# ============================================================

harmonized_calc <- primary_295 %>%
  left_join(
    harmonization_decisions,
    by = "record_key"
  ) %>%
  mutate(
    Analysis_Point_numeric = as.numeric(Analysis_Point),
    Analysis_Units_clean = trimws(as.character(Analysis_Units)),
    
    # Target-metal mass concentration in mg/L.
    LC50_target_mg_L_calc = case_when(
      Basis_Interpretation == "PENDING_SOURCE_BASIS" ~ NA_real_,
      
      # O'Brien Ref. 2977 source table is on a 10^-6 M scale:
      # the numeric value is directly umol Cu/L.
      Basis_Interpretation == "SOURCE_VERIFIED_MICROMOLAR" ~
        Analysis_Point_numeric * Atomic_Weight_g_mol / 1000,
      
      Analysis_Units_clean %in% c("ug/L", "ppb") ~
        Analysis_Point_numeric / 1000,
      
      Analysis_Units_clean %in% c("mg/L", "ppm") ~
        Analysis_Point_numeric,
      
      Analysis_Units_clean == "nM" ~
        (Analysis_Point_numeric / 1000) * Atomic_Weight_g_mol / 1000,
      
      TRUE ~ NA_real_
    ),
    
    # Target-metal concentration in umol/L.
    LC50_umol_L_calc = case_when(
      Basis_Interpretation == "PENDING_SOURCE_BASIS" ~ NA_real_,
      
      Basis_Interpretation == "SOURCE_VERIFIED_MICROMOLAR" ~
        Analysis_Point_numeric,
      
      Analysis_Units_clean %in% c("ug/L", "ppb") ~
        (Analysis_Point_numeric / 1000) * 1000 / Atomic_Weight_g_mol,
      
      Analysis_Units_clean %in% c("mg/L", "ppm") ~
        Analysis_Point_numeric * 1000 / Atomic_Weight_g_mol,
      
      Analysis_Units_clean == "nM" ~
        Analysis_Point_numeric / 1000,
      
      TRUE ~ NA_real_
    ),
    
    ln_LC50_umol_L_calc = case_when(
      !is.na(LC50_umol_L_calc) & LC50_umol_L_calc > 0 ~
        log(LC50_umol_L_calc),
      TRUE ~ NA_real_
    )
  )

# Every harmonized row must have a supported unit/basis rule.
unsupported_harmonized_units <- harmonized_calc %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    (
      is.na(LC50_target_mg_L_calc) |
        is.na(LC50_umol_L_calc) |
        is.na(ln_LC50_umol_L_calc)
    )
  )

if (nrow(unsupported_harmonized_units) > 0) {
  write_csv(
    unsupported_harmonized_units,
    file.path(qa_dir, "04a_unsupported_harmonized_units.csv")
  )
  stop("At least one HARMONIZED record could not be reconstructed from its documented unit/basis rule.")
}

near_equal <- function(a, b, tol = 1e-10) {
  both_na <- is.na(a) & is.na(b)
  both_present <- !is.na(a) & !is.na(b)
  
  out <- both_na
  out[both_present] <-
    abs(a[both_present] - b[both_present]) <=
    tol * pmax(
      1,
      abs(a[both_present]),
      abs(b[both_present])
    )
  
  out
}

harmonization_numeric_audit <- harmonized_calc %>%
  transmute(
    record_key,
    Test_Number,
    Result_Number,
    Reference_Number,
    Metal,
    Analysis_Point,
    Analysis_Units,
    Basis_Interpretation,
    
    LC50_target_mg_L_calc,
    Ledger_LC50_target_mg_L,
    target_mg_match = near_equal(
      LC50_target_mg_L_calc,
      Ledger_LC50_target_mg_L
    ),
    
    LC50_umol_L_calc,
    Ledger_LC50_umol_L,
    umol_match = near_equal(
      LC50_umol_L_calc,
      Ledger_LC50_umol_L
    ),
    
    ln_LC50_umol_L_calc,
    Ledger_ln_LC50_umol_L,
    ln_match = near_equal(
      ln_LC50_umol_L_calc,
      Ledger_ln_LC50_umol_L
    )
  )

harmonization_mismatches <- harmonization_numeric_audit %>%
  filter(
    !target_mg_match |
      !umol_match |
      !ln_match
  )

write_csv(
  harmonization_numeric_audit,
  file.path(qa_dir, "04_harmonization_numeric_audit_all_rows.csv")
)

if (nrow(harmonization_mismatches) > 0) {
  write_csv(
    harmonization_mismatches,
    file.path(qa_dir, "04b_harmonization_numeric_mismatches.csv")
  )
  stop("Recomputed molar harmonization does not match the documented Harmonized_295 ledger.")
}

harmonized_295 <- harmonized_calc %>%
  mutate(
    LC50_target_mg_L = LC50_target_mg_L_calc,
    LC50_umol_L = LC50_umol_L_calc,
    ln_LC50_umol_L = ln_LC50_umol_L_calc
  )

stopifnot(
  sum(harmonized_295$Harmonization_Status == "HARMONIZED") == 246,
  sum(harmonized_295$Harmonization_Status == "PENDING_CR_SOURCE_BASIS") == 49,
  sum(
    harmonized_295$Metal %in% c("Cu", "Cd", "Zn") &
      harmonized_295$Harmonization_Status == "HARMONIZED"
  ) == 186
)

harmonization_summary <- harmonized_295 %>%
  count(
    Metal,
    Harmonization_Status,
    name = "Results"
  ) %>%
  arrange(Metal, Harmonization_Status)

write_csv(
  harmonization_summary,
  file.path(qa_dir, "05_harmonization_status_by_metal.csv")
)


# ============================================================
# 7. JOIN RAW ENVIRONMENTAL METADATA + FULL RAW TITLE
# ============================================================

raw_primary_metadata <- raw_lookup %>%
  filter(record_key %in% primary_295$record_key) %>%
  select(
    record_key,
    Title_Raw,
    Chemical_Analysis,
    Temperature_Raw_Mean,
    Temperature_Raw_Min,
    Temperature_Raw_Max,
    Temperature_Raw_Units,
    Salinity_Raw_Mean,
    Salinity_Raw_Min,
    Salinity_Raw_Max,
    Salinity_Raw_Units,
    pH_Raw_Mean,
    pH_Raw_Min,
    pH_Raw_Max,
    pH_Raw_Units,
    Organic_Carbon_Raw_Mean,
    Organic_Carbon_Raw_Min,
    Organic_Carbon_Raw_Max,
    Organic_Carbon_Raw_Units
  )

stopifnot(
  nrow(raw_primary_metadata) == 295,
  n_distinct(raw_primary_metadata$record_key) == 295
)

# This lookup allows later scripts to use raw protocol metadata
# without reopening the 244-column RAW workbook.
write_csv(
  harmonized_295 %>%
    select(
      record_key,
      Test_Number,
      Result_Number,
      Reference_Number,
      Metal
    ) %>%
    left_join(
      raw_primary_metadata %>%
        select(
          record_key,
          Chemical_Analysis,
          starts_with("Temperature_Raw_"),
          starts_with("Salinity_Raw_"),
          starts_with("pH_Raw_"),
          starts_with("Organic_Carbon_Raw_")
        ),
      by = "record_key"
    ) %>%
    select(-record_key),
  file.path(
    lookup_dir,
    "primary_raw_protocol_environment_lookup.csv"
  )
)


# ============================================================
# 8. ENVIRONMENTAL METADATA REPAIR
# ============================================================
#
# Historical data-preparation rule used for the final R-ready CSV:
# - preserve an existing cleaned mean when present;
# - otherwise, if RAW min and max are both available, use their
#   midpoint;
# - otherwise retain missingness.
#
# No environmental value is imputed from another test/reference.
# ============================================================

prepared <- harmonized_295 %>%
  left_join(
    raw_primary_metadata,
    by = "record_key"
  )

# In the current frozen data, there are no cases where the cleaned
# mean is missing but a RAW mean is present. If this changes later,
# stop rather than silently invent a new repair rule.
stopifnot(
  sum(is.na(prepared$Temperature) & !is.na(prepared$Temperature_Raw_Mean)) == 0,
  sum(is.na(prepared$Salinity) & !is.na(prepared$Salinity_Raw_Mean)) == 0,
  sum(is.na(prepared$pH) & !is.na(prepared$pH_Raw_Mean)) == 0,
  sum(is.na(prepared$Organic_Carbon) & !is.na(prepared$Organic_Carbon_Raw_Mean)) == 0
)

prepared <- prepared %>%
  mutate(
    Temperature_original = Temperature,
    Salinity_original = Salinity,
    pH_original = pH,
    Organic_Carbon_original = Organic_Carbon,
    
    Temperature = case_when(
      !is.na(Temperature_original) ~ Temperature_original,
      !is.na(Temperature_Raw_Min) & !is.na(Temperature_Raw_Max) ~
        (Temperature_Raw_Min + Temperature_Raw_Max) / 2,
      TRUE ~ NA_real_
    ),
    Temperature_Repair_Status = case_when(
      !is.na(Temperature_original) ~ "EXISTING_MEAN_PRESERVED",
      !is.na(Temperature_Raw_Min) & !is.na(Temperature_Raw_Max) ~
        "RAW_RANGE_MIDPOINT_RECOVERED",
      TRUE ~ "MISSING"
    ),
    
    Salinity = case_when(
      !is.na(Salinity_original) ~ Salinity_original,
      !is.na(Salinity_Raw_Min) & !is.na(Salinity_Raw_Max) ~
        (Salinity_Raw_Min + Salinity_Raw_Max) / 2,
      TRUE ~ NA_real_
    ),
    Salinity_Repair_Status = case_when(
      !is.na(Salinity_original) ~ "EXISTING_MEAN_PRESERVED",
      !is.na(Salinity_Raw_Min) & !is.na(Salinity_Raw_Max) ~
        "RAW_RANGE_MIDPOINT_RECOVERED",
      TRUE ~ "MISSING"
    ),
    
    pH = case_when(
      !is.na(pH_original) ~ pH_original,
      !is.na(pH_Raw_Min) & !is.na(pH_Raw_Max) ~
        (pH_Raw_Min + pH_Raw_Max) / 2,
      TRUE ~ NA_real_
    ),
    pH_Repair_Status = case_when(
      !is.na(pH_original) ~ "EXISTING_MEAN_PRESERVED",
      !is.na(pH_Raw_Min) & !is.na(pH_Raw_Max) ~
        "RAW_RANGE_MIDPOINT_RECOVERED",
      TRUE ~ "MISSING"
    ),
    
    Organic_Carbon = case_when(
      !is.na(Organic_Carbon_original) ~ Organic_Carbon_original,
      !is.na(Organic_Carbon_Raw_Min) & !is.na(Organic_Carbon_Raw_Max) ~
        (Organic_Carbon_Raw_Min + Organic_Carbon_Raw_Max) / 2,
      TRUE ~ NA_real_
    ),
    Organic_Carbon_Repair_Status = case_when(
      !is.na(Organic_Carbon_original) ~ "EXISTING_MEAN_PRESERVED",
      !is.na(Organic_Carbon_Raw_Min) & !is.na(Organic_Carbon_Raw_Max) ~
        "RAW_RANGE_MIDPOINT_RECOVERED",
      TRUE ~ "MISSING"
    ),
    
    # Use the untruncated RAW title in the canonical CSV.
    Title = Title_Raw
  )

# Frozen repair benchmarks from the final analysis-ready dataset.
stopifnot(
  sum(prepared$Temperature_Repair_Status == "EXISTING_MEAN_PRESERVED") == 233,
  sum(prepared$Temperature_Repair_Status == "RAW_RANGE_MIDPOINT_RECOVERED") == 30,
  sum(prepared$Temperature_Repair_Status == "MISSING") == 32,
  
  sum(prepared$Salinity_Repair_Status == "EXISTING_MEAN_PRESERVED") == 221,
  sum(prepared$Salinity_Repair_Status == "RAW_RANGE_MIDPOINT_RECOVERED") == 42,
  sum(prepared$Salinity_Repair_Status == "MISSING") == 32,
  
  sum(prepared$pH_Repair_Status == "EXISTING_MEAN_PRESERVED") == 35,
  sum(prepared$pH_Repair_Status == "RAW_RANGE_MIDPOINT_RECOVERED") == 33,
  sum(prepared$pH_Repair_Status == "MISSING") == 227,
  
  sum(prepared$Organic_Carbon_Repair_Status == "EXISTING_MEAN_PRESERVED") == 4,
  sum(prepared$Organic_Carbon_Repair_Status == "RAW_RANGE_MIDPOINT_RECOVERED") == 3,
  sum(prepared$Organic_Carbon_Repair_Status == "MISSING") == 288
)

environment_repair_summary <- bind_rows(
  prepared %>%
    count(Temperature_Repair_Status, name = "Results") %>%
    transmute(
      Variable = "Temperature",
      Repair_Status = Temperature_Repair_Status,
      Results
    ),
  prepared %>%
    count(Salinity_Repair_Status, name = "Results") %>%
    transmute(
      Variable = "Salinity",
      Repair_Status = Salinity_Repair_Status,
      Results
    ),
  prepared %>%
    count(pH_Repair_Status, name = "Results") %>%
    transmute(
      Variable = "pH",
      Repair_Status = pH_Repair_Status,
      Results
    ),
  prepared %>%
    count(Organic_Carbon_Repair_Status, name = "Results") %>%
    transmute(
      Variable = "Organic_Carbon",
      Repair_Status = Organic_Carbon_Repair_Status,
      Results
    )
)

write_csv(
  environment_repair_summary,
  file.path(qa_dir, "06_environment_repair_summary.csv")
)


# ============================================================
# 9. BUILD FINAL 55-COLUMN ANALYSIS-READY DATASET
# ============================================================

final_data <- prepared %>%
  arrange(Raw_Row) %>%
  select(
    # Original cleaned/source-adjudicated analytical fields
    Raw_Row,
    Test_Number,
    Result_Number,
    Reference_Number,
    Metal,
    Chemical_Name,
    Species,
    Order,
    Lifestage,
    Duration_days,
    Exposure_Type,
    Conc_Type,
    Mean_Op,
    Mean,
    Min,
    Max,
    Units,
    Analysis_Point,
    Analysis_Units,
    Temperature,
    Salinity,
    pH,
    Organic_Carbon,
    Title,
    Eligibility_Status,
    Eligibility_Reason,
    
    # Harmonization fields
    Atomic_Weight_g_mol,
    Basis_Interpretation,
    LC50_target_mg_L,
    LC50_umol_L,
    ln_LC50_umol_L,
    Harmonization_Status,
    Harmonization_Rule,
    Metal_Role,
    Core_Pairwise_Ready,
    
    # RAW environmental audit / repair provenance
    Temperature_Raw_Mean,
    Temperature_Raw_Min,
    Temperature_Raw_Max,
    Temperature_Raw_Units,
    Temperature_Repair_Status,
    
    Salinity_Raw_Mean,
    Salinity_Raw_Min,
    Salinity_Raw_Max,
    Salinity_Raw_Units,
    Salinity_Repair_Status,
    
    pH_Raw_Mean,
    pH_Raw_Min,
    pH_Raw_Max,
    pH_Raw_Units,
    pH_Repair_Status,
    
    Organic_Carbon_Raw_Mean,
    Organic_Carbon_Raw_Min,
    Organic_Carbon_Raw_Max,
    Organic_Carbon_Raw_Units,
    Organic_Carbon_Repair_Status
  )

stopifnot(
  nrow(final_data) == 295,
  ncol(final_data) == 55,
  n_distinct(
    make_record_key(
      final_data$Test_Number,
      final_data$Result_Number,
      final_data$Reference_Number
    )
  ) == 295
)


# ============================================================
# 10. OPTIONAL REGRESSION TEST AGAINST EXISTING CANONICAL CSV
# ============================================================
#
# The existing CSV is NOT an input to the reconstruction.
# If present, it is used only as a regression reference to prove
# that Script 01 reconstructs the dataset used by
# the completed thesis analyses.
# ============================================================

reference_comparison_summary <- tibble(
  Column = character(),
  Mismatch_Count = integer(),
  Max_Abs_Difference = double()
)

reference_mismatch_details <- tibble(
  Row = integer(),
  Column = character(),
  Rebuilt = character(),
  Existing = character()
)

if (!is.null(existing_target)) {
  
  if (!identical(names(final_data), names(existing_target))) {
    stop(
      paste0(
        "Existing canonical CSV has a different column schema.\n",
        "Rebuilt columns: ", paste(names(final_data), collapse = " | "), "\n",
        "Existing columns: ", paste(names(existing_target), collapse = " | ")
      )
    )
  }
  
  if (nrow(final_data) != nrow(existing_target)) {
    stop("Existing canonical CSV has a different row count.")
  }
  
  comparison_rows <- vector("list", length(names(final_data)))
  detail_rows <- list()
  detail_counter <- 0L
  
  for (j in seq_along(names(final_data))) {
    
    col_name <- names(final_data)[[j]]
    a <- final_data[[j]]
    b <- existing_target[[j]]
    
    if (is.numeric(a) && is.numeric(b)) {
      match_vec <- near_equal(a, b, tol = 1e-9)
      abs_diff <- abs(a - b)
      abs_diff[is.na(abs_diff)] <- 0
      max_diff <- if (length(abs_diff) == 0) 0 else max(abs_diff)
    } else {
      a_chr <- ifelse(is.na(a), "<NA>", as.character(a))
      b_chr <- ifelse(is.na(b), "<NA>", as.character(b))
      match_vec <- a_chr == b_chr
      max_diff <- NA_real_
    }
    
    mismatch_idx <- which(!match_vec)
    
    comparison_rows[[j]] <- tibble(
      Column = col_name,
      Mismatch_Count = length(mismatch_idx),
      Max_Abs_Difference = max_diff
    )
    
    if (length(mismatch_idx) > 0) {
      for (ii in mismatch_idx) {
        detail_counter <- detail_counter + 1L
        detail_rows[[detail_counter]] <- tibble(
          Row = ii,
          Column = col_name,
          Rebuilt = ifelse(is.na(a[[ii]]), "<NA>", as.character(a[[ii]])),
          Existing = ifelse(is.na(b[[ii]]), "<NA>", as.character(b[[ii]]))
        )
      }
    }
  }
  
  reference_comparison_summary <- bind_rows(comparison_rows)
  
  if (length(detail_rows) > 0) {
    reference_mismatch_details <- bind_rows(detail_rows)
  }
  
  write_csv(
    reference_comparison_summary,
    file.path(qa_dir, "07_existing_canonical_csv_comparison.csv")
  )
  
  write_csv(
    reference_mismatch_details,
    file.path(qa_dir, "07a_existing_canonical_csv_mismatch_details.csv")
  )
  
  total_mismatches <- sum(reference_comparison_summary$Mismatch_Count)
  
  cat("\nExisting canonical CSV regression mismatches:", total_mismatches, "\n")
  
  if (total_mismatches > 0) {
    stop(
      paste0(
        "Script 01 did not exactly reconstruct the existing processed ECOTOX CSV. ",
        "Inspect outputs/01_prepare_ecotox_data/01_QA/07*.csv before overwriting."
      )
    )
  }
}


# ============================================================
# 10b. DOCUMENTED SOURCE-REVIEW CORRECTION: SOSNOWSKI (1979)
# The preceding QA reproduces the historical 295-row baseline unchanged.
# Source review, pp. 449-452: natural populations sampled on different dates;
# the former deliberate-modifier exclusion is unsupported for this result.
# 72 h Cu LC50 = 9 ug/L; general 24-96 h pool only, not the 96 h main model.
# The corrected row is appended only after the frozen 295-row baseline has passed QA.

sos <- screening_ledger %>% filter(as.character(Result_Number) == "42366")
stopifnot(nrow(sos) == 1L,
          as.character(sos$Reference_Number) == "8446",
          as.character(sos$Test_Number) == "1111947",
          sos$Species == "Acartia tonsa", sos$Metal == "Cu",
          as.numeric(sos$Analysis_Point) == 9,
          as.character(sos$Analysis_Units) == "ug/L",
          sos$Eligibility_Status == "SECONDARY_MODIFIER",
          !any(as.character(final_data$Result_Number) == "42366"))
sos_raw <- raw_lookup %>% filter(record_key == sos$record_key)
stopifnot(nrow(sos_raw) == 1L)

# A typed empty row preserves the established 55-column schema.
sos_row <- final_data[NA_integer_, ]
for (nm in intersect(names(sos), names(sos_row))) {
  sos_row[[nm]] <- sos[[nm]]
}
for (nm in intersect(names(sos_raw), names(sos_row))) {
  if (grepl("_Raw_", nm)) sos_row[[nm]] <- sos_raw[[nm]]
}
sos_row$Title <- sos_raw$Title_Raw
sos_row$Duration_days <- 3
sos_row$Exposure_Type <- "Flow-through"
sos_row$Conc_Type <- "Total"
sos_row$Temperature <- 20
sos_row$Salinity <- 30
sos_row$Eligibility_Status <- "PRIMARY_KEEP"
sos_row$Eligibility_Reason <- paste(
  "Source review: Sosnowski (1979), pp. 449-452; natural populations,",
  "not a deliberate modifier experiment; waterborne Cu LC50, 72 h."
)
sos_row$Atomic_Weight_g_mol <- 63.546
sos_row$Basis_Interpretation <- "ECOTOX_STANDARDIZED_TARGET_METAL"
sos_row$LC50_target_mg_L <- 9 / 1000
sos_row$LC50_umol_L <- sos_row$LC50_target_mg_L * 1000 / 63.546
sos_row$ln_LC50_umol_L <- log(sos_row$LC50_umol_L)
sos_row$Harmonization_Status <- "HARMONIZED"
sos_row$Harmonization_Rule <- "Source-confirmed 9 ug Cu/L converted to mg Cu/L and umol Cu/L."
sos_row$Metal_Role <- "CORE"
sos_row$Core_Pairwise_Ready <- "YES"
sos_row$Temperature_Repair_Status <- "SOURCE_VERIFIED"
sos_row$Salinity_Repair_Status <- "SOURCE_VERIFIED"
# Do not infer unverified pH or organic-carbon values.
sos_row$pH <- NA_real_
sos_row$Organic_Carbon <- NA_real_
sos_row$pH_Repair_Status <- "MISSING"
sos_row$Organic_Carbon_Repair_Status <- "MISSING"

baseline_data <- final_data
final_data <- bind_rows(final_data, sos_row) %>% arrange(Raw_Row)
stopifnot(nrow(final_data) == 296L,
          n_distinct(final_data$Result_Number) == 296L,
          sum(final_data$Harmonization_Status == "HARMONIZED") == 247L,
          identical(final_data %>% filter(as.character(Result_Number) != "42366"),
                    baseline_data),
          identical(final_data %>% filter(Duration_days == 4),
                    baseline_data %>% filter(Duration_days == 4)))

# Preserve both the original decision and the applied correction.
write_csv(tibble(
  Result_Number = 42366, Reference_Number = 8446,
  Previous_Status = as.character(sos$Eligibility_Status),
  Previous_Reason = as.character(sos$Eligibility_Reason),
  Applied_Status = sos_row$Eligibility_Status,
  Applied_Reason = sos_row$Eligibility_Reason,
  Duration_days = 3, LC50_target_mg_L = 0.009
), file.path(ledger_dir, "Sosnowski_1979_source_review_correction.csv"))
write_csv(sos_row, file.path(ledger_dir, "Sosnowski_1979_applied_row.csv"), na = "")
write_csv(screening_ledger, file.path(ledger_dir, "screening_410_baseline_before_correction.csv"))
screening_ledger$Eligibility_Status[as.character(screening_ledger$Result_Number) == "42366"] <- "PRIMARY_KEEP"
screening_ledger$Eligibility_Reason[as.character(screening_ledger$Result_Number) == "42366"] <- sos_row$Eligibility_Reason
write_csv(screening_ledger, file.path(ledger_dir, "screening_and_source_decision_ledger_410.csv"))
write_csv(final_data %>% count(Eligibility_Status, Harmonization_Status, name = "Results"),
          file.path(qa_dir, "08_final_after_source_review_summary.csv"))
cat("\nSource-review correction applied: Result 42366; final ECOTOX rows = 296.\n")
cat("Earlier QA summaries describe the unchanged 295-row baseline.\n")

# ============================================================
# 11. WRITE REBUILT + CANONICAL DATASETS
# ============================================================

rebuilt_output_file <- file.path(
  output_dir,
  "ecotox_harmonized_rebuilt.csv"
)

write_csv(
  final_data,
  rebuilt_output_file,
  na = ""
)

# Only reached after all frozen QA and optional regression checks
# have passed.
write_csv(
  final_data,
  canonical_output_file,
  na = ""
)


# ============================================================
# 12. SAVE REPRODUCIBILITY INFORMATION
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

input_manifest <- tibble(
  Role = c(
    "Immutable raw ECOTOX export",
    "Documented human/source-decision ledger",
    "Processed ECOTOX dataset for downstream WoS integration"
  ),
  File = c(
    raw_file,
    working_master_file,
    canonical_output_file
  )
)

write_csv(
  input_manifest,
  file.path(output_dir, "input_output_manifest.csv")
)


# ============================================================
# 13. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n\n")
cat("============================================\n")
cat("FINAL DATA-PREPARATION QA SUMMARY\n")
cat("============================================\n")

cat("\n--- RAW / TAXONOMY ---\n")
print(raw_screen_summary)

cat("\n--- SCREENING DISPOSITION ---\n")
print(screening_disposition_check)

cat("\n--- PRIMARY DATASET ---\n")
print(primary_summary)

cat("\n--- HARMONIZATION STATUS ---\n")
print(harmonization_summary)

cat("\n--- ENVIRONMENTAL REPAIR ---\n")
print(environment_repair_summary)

cat("\n--- FINAL DATASET ---\n")
cat("Rows:", nrow(final_data), "\n")
cat("Columns:", ncol(final_data), "\n")
cat("Harmonized quantitative Results:", sum(final_data$Harmonization_Status == "HARMONIZED"), "\n")
cat("Pending Cr source-basis Results:", sum(final_data$Harmonization_Status == "PENDING_CR_SOURCE_BASIS"), "\n")

if (!is.null(existing_target)) {
  cat("Existing canonical baseline (excluding corrected Result 42366) reproduced exactly: TRUE\n")
}

cat("\nProcessed ECOTOX dataset written to:\n")
cat(canonical_output_file, "\n")

cat("\nReproducibility package written to:\n")
cat(output_dir, "\n")

cat("\nSTATUS: 01 PREPARE ECOTOX DATA = PASS\n")

# ============================================================
# END
# ==================================