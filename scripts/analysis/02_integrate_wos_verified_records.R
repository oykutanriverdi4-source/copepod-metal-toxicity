# ============================================================
# 02_integrate_wos_verified_records.R
# ============================================================
# PURPOSE
# Integrate the source-verified Web of Science (WoS) LC50 records
# with the reconstructed ECOTOX dataset used in this thesis.
#
# SCIENTIFIC SCOPE
# - Script 01 remains the authoritative ECOTOX reconstruction.
# - WoS records are read only from the curated source-verification
#   workbook produced after publication-level review.
# - Only rows marked Include_in_combined == "YES" are integrated.
# - Reported LC50 values are not silently corrected in this script.
# - ECOTOX native Reference/Test/Result identifiers are preserved.
# - Cross-source analyses use project-wide identifiers:
#       ECOTOX Reference_ID: E_<ECOTOX Reference Number>
#       WoS    Reference_ID: W_001, W_002, ...
# - WoS Result_ID values are rebuilt deterministically from Test_ID.
# - Metal/domain selection is not performed here. That is handled
#   in the subsequent support/comparability audit.
#
# INPUTS
#   data/processed/ecotox_harmonized.csv
#   data/curated/wos/wos_source_verified_records.xlsx
#
# CORE OUTPUT
#   data/processed/combined_ecotox_wos_harmonized.csv
#
# QA / AUDIT OUTPUTS
#   outputs/02_integrate_wos_verified_records/
# ============================================================


# ============================================================
# 0. PACKAGES + PATHS
# ============================================================

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

options(width = 220)

# ---- Inputs ----
ecotox_file <- file.path(
  "data", "processed",
  "ecotox_harmonized.csv"
)

wos_file <- file.path(
  "data", "curated", "wos",
  "wos_source_verified_records.xlsx"
)

# ---- Core processed output ----
processed_data_dir <- file.path("data", "processed")
combined_csv <- file.path(
  processed_data_dir,
  "combined_ecotox_wos_harmonized.csv"
)

# ---- QA / audit outputs ----
output_dir <- file.path(
  "outputs",
  "02_integrate_wos_verified_records"
)

qa_dir <- file.path(output_dir, "01_QA")
data_dir <- file.path(output_dir, "02_data")

for (d in c(processed_data_dir, output_dir, qa_dir, data_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required_inputs <- c(
  "ECOTOX processed dataset" = ecotox_file,
  "WoS source-verification workbook" = wos_file
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
  stop(
    paste0(
      "Required input file(s) not found:\n  - ",
      paste(unname(missing_inputs), collapse = "\n  - "),
      "\n\nRun Script 01 first and confirm the repository data/ folder structure."
    )
  )
}

# Confirm the two curated WoS sheets required by this script exist.
wos_sheets <- excel_sheets(wos_file)
required_wos_sheets <- c("Source_Audit", "Results")
missing_wos_sheets <- setdiff(required_wos_sheets, wos_sheets)

if (length(missing_wos_sheets) > 0L) {
  stop(
    paste0(
      "WoS source-verification workbook is missing required sheet(s): ",
      paste(missing_wos_sheets, collapse = " | ")
    )
  )
}


# ============================================================
# 1. LOAD THE SCRIPT-01 ECOTOX DATASET
# ============================================================

ecotox <- read_csv(
  ecotox_file,
  show_col_types = FALSE
)

# Final source-reviewed Script-01 benchmark:
# frozen 295-row baseline + reviewed Sosnowski Result 42366.
stopifnot(
  nrow(ecotox) == 296,
  n_distinct(ecotox$Test_Number) == 289,
  n_distinct(ecotox$Reference_Number) == 67
)

# Require the reviewed Sosnowski result, not merely an increased row count.
sosnowski_input <- ecotox %>%
  filter(as.character(Result_Number) == "42366")

if (nrow(sosnowski_input) != 1L) {
  stop(
    paste0(
      "Expected exactly one source-reviewed Sosnowski Result 42366. ",
      "Run the approved Script 01 first."
    )
  )
}

stopifnot(
  as.character(sosnowski_input$Reference_Number) == "8446",
  as.character(sosnowski_input$Test_Number) == "1111947",
  sosnowski_input$Species == "Acartia tonsa",
  sosnowski_input$Metal == "Cu",
  sosnowski_input$Duration_days == 3,
  sosnowski_input$Eligibility_Status == "PRIMARY_KEEP",
  sosnowski_input$Harmonization_Status == "HARMONIZED",
  abs(sosnowski_input$LC50_target_mg_L - 0.009) < 1e-12
)

# Full 55-column ECOTOX schema expected from the approved Script 01.
required_ecotox_columns <- c(
  "Raw_Row",
  "Test_Number",
  "Result_Number",
  "Reference_Number",
  "Metal",
  "Chemical_Name",
  "Species",
  "Order",
  "Lifestage",
  "Duration_days",
  "Exposure_Type",
  "Conc_Type",
  "Mean_Op",
  "Mean",
  "Min",
  "Max",
  "Units",
  "Analysis_Point",
  "Analysis_Units",
  "Temperature",
  "Salinity",
  "pH",
  "Organic_Carbon",
  "Title",
  "Eligibility_Status",
  "Eligibility_Reason",
  "Atomic_Weight_g_mol",
  "Basis_Interpretation",
  "LC50_target_mg_L",
  "LC50_umol_L",
  "ln_LC50_umol_L",
  "Harmonization_Status",
  "Harmonization_Rule",
  "Metal_Role",
  "Core_Pairwise_Ready",
  "Temperature_Raw_Mean",
  "Temperature_Raw_Min",
  "Temperature_Raw_Max",
  "Temperature_Raw_Units",
  "Temperature_Repair_Status",
  "Salinity_Raw_Mean",
  "Salinity_Raw_Min",
  "Salinity_Raw_Max",
  "Salinity_Raw_Units",
  "Salinity_Repair_Status",
  "pH_Raw_Mean",
  "pH_Raw_Min",
  "pH_Raw_Max",
  "pH_Raw_Units",
  "pH_Repair_Status",
  "Organic_Carbon_Raw_Mean",
  "Organic_Carbon_Raw_Min",
  "Organic_Carbon_Raw_Max",
  "Organic_Carbon_Raw_Units",
  "Organic_Carbon_Repair_Status"
)

missing_ecotox_columns <- setdiff(
  required_ecotox_columns,
  names(ecotox)
)

if (length(missing_ecotox_columns) > 0L) {
  stop(
    paste0(
      "The ECOTOX input is not the approved 55-column Script-01 dataset.\n",
      "Missing columns: ",
      paste(missing_ecotox_columns, collapse = " | ")
    )
  )
}


# ============================================================
# 2. ADD PROJECT-WIDE IDENTIFIERS TO ECOTOX
# ============================================================
# Native ECOTOX numeric identifiers remain untouched.
# ============================================================

ecotox_aug <- ecotox %>%
  mutate(
    Source_Origin = "ECOTOX",
    Reference_ID = paste0("E_", Reference_Number),
    Test_ID = paste0("E_T", Test_Number),
    Result_ID = paste0("E_R", Result_Number),
    External_Audit_ID = NA_character_,
    External_DOI = NA_character_,
    External_Reference_Label = NA_character_,
    External_Stage_Sex = NA_character_,
    External_Concentration_Type_Detail = NA_character_,
    External_Source_Location = NA_character_,
    External_Source_Note = NA_character_
  )

stopifnot(
  n_distinct(ecotox_aug$Reference_ID) == 67,
  n_distinct(ecotox_aug$Test_ID) == 289,
  n_distinct(ecotox_aug$Result_ID) == 296
)


# ============================================================
# 3. LOAD THE CURATED WoS SOURCE-VERIFICATION WORKBOOK
# ============================================================

wos_results_raw <- read_excel(
  wos_file,
  sheet = "Results"
)

wos_source_audit <- read_excel(
  wos_file,
  sheet = "Source_Audit"
)

required_wos_columns <- c(
  "Audit_ID",
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "Reference",
  "Species",
  "Order",
  "Stage_Sex",
  "Metal",
  "Duration_h",
  "LC50_reported",
  "Unit",
  "Chemical",
  "Concentration_basis",
  "Concentration_type",
  "Temp_C",
  "Salinity",
  "Include_in_combined",
  "Extrapolation_flag",
  "Duplicate_flag",
  "Source_location",
  "Note"
)

missing_wos_columns <- setdiff(
  required_wos_columns,
  names(wos_results_raw)
)

if (length(missing_wos_columns) > 0L) {
  stop(
    paste0(
      "WoS Results sheet is missing required columns: ",
      paste(missing_wos_columns, collapse = " | ")
    )
  )
}

# Source_Audit is used to attach DOI and final publication role.
required_audit_columns <- c(
  "Audit_ID",
  "Reference",
  "DOI",
  "Final_role"
)

missing_audit_columns <- setdiff(
  required_audit_columns,
  names(wos_source_audit)
)

if (length(missing_audit_columns) > 0L) {
  stop(
    paste0(
      "WoS Source_Audit sheet is missing required columns: ",
      paste(missing_audit_columns, collapse = " | ")
    )
  )
}


# ============================================================
# 4. FREEZE THE FINAL SOURCE-VERIFIED WoS KEEP SET
# ============================================================
# These QA benchmarks describe the source-verified dataset used
# in the final thesis. They should change only if the documented
# source audit is intentionally revised.
# ============================================================

wos_keep <- wos_results_raw %>%
  mutate(
    Include_in_combined = toupper(trimws(as.character(Include_in_combined))),
    Duplicate_flag = toupper(trimws(as.character(Duplicate_flag))),
    Extrapolation_flag = toupper(trimws(as.character(Extrapolation_flag)))
  ) %>%
  filter(Include_in_combined == "YES") %>%
  left_join(
    wos_source_audit %>%
      select(
        Audit_ID,
        DOI,
        Final_role
      ),
    by = "Audit_ID"
  )

stopifnot(
  nrow(wos_keep) == 57,
  n_distinct(wos_keep$Reference_ID) == 17,
  all(wos_keep$Duplicate_flag != "YES"),
  all(!is.na(wos_keep$LC50_reported)),
  all(wos_keep$LC50_reported > 0),
  all(wos_keep$Metal %in% c("Cu", "Cd", "Zn", "Ni", "Pb", "Ag", "Cr", "Hg"))
)

# Explicitly retained Yang et al. (2022) Zn record.
stopifnot(
  sum(wos_keep$Reference_ID == "W_017") == 1L,
  sum(
    wos_keep$DOI == "10.1016/j.chemosphere.2022.134099",
    na.rm = TRUE
  ) == 1L
)

# Final source-verification yield by metal.
expected_wos_by_metal <- tribble(
  ~Metal, ~Expected_Results,
  "Ag",  1L,
  "Cd", 15L,
  "Cu",  8L,
  "Hg", 12L,
  "Ni", 12L,
  "Zn",  9L
)

wos_by_metal_check <- wos_keep %>%
  count(Metal, name = "Results") %>%
  full_join(expected_wos_by_metal, by = "Metal") %>%
  mutate(
    Results = replace_na(Results, 0L),
    Expected_Results = replace_na(Expected_Results, 0L),
    Match = Results == Expected_Results
  ) %>%
  arrange(Metal)

write_csv(
  wos_by_metal_check,
  file.path(qa_dir, "01_wos_keep_results_by_metal.csv")
)

stopifnot(all(wos_by_metal_check$Match))


# ============================================================
# 5. STANDARDIZE WoS IDENTIFIERS + BASIC METADATA
# ============================================================
# Reference_ID and Test_ID come from the human-reviewed ledger.
# Result_ID is rebuilt deterministically as:
#   <Test_ID>_R01, <Test_ID>_R02, ...
# ============================================================

normalize_unit <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[\u00B5\u03BC]", "u")
  x <- str_to_lower(x)
  x <- str_replace_all(x, "\\s+", "")

  case_when(
    str_detect(x, "^ug") ~ "ug/L",
    str_detect(x, "^mg") ~ "mg/L",
    TRUE ~ NA_character_
  )
}

normalize_lifestage <- function(x) {
  x_low <- str_to_lower(as.character(x))

  case_when(
    str_detect(x_low, "naupli") ~ "Nauplii",
    str_detect(x_low, "copepodite|copepodid") ~ "Copepodite",
    str_detect(x_low, "adult|female|male") ~ "Adult",
    TRUE ~ NA_character_
  )
}

normalize_conc_type <- function(x) {
  x_low <- str_to_lower(as.character(x))

  case_when(
    str_detect(x_low, "nominal") ~ "Nominal",
    str_detect(x_low, "measured") ~ "Measured",
    TRUE ~ "Unspecified"
  )
}

wos_keep <- wos_keep %>%
  arrange(
    Reference_ID,
    Test_ID,
    Duration_h,
    Metal,
    Species,
    Stage_Sex
  ) %>%
  group_by(Test_ID) %>%
  mutate(
    Result_ID_final = paste0(
      Test_ID,
      "_R",
      sprintf("%02d", row_number())
    )
  ) %>%
  ungroup() %>%
  mutate(
    Unit_normalized = normalize_unit(Unit),
    Lifestage_final = normalize_lifestage(Stage_Sex),
    Conc_Type_final = normalize_conc_type(Concentration_type),
    Duration_days_final = as.numeric(Duration_h) / 24,
    Temp_C = as.numeric(Temp_C),
    Salinity = as.numeric(Salinity),
    LC50_reported = as.numeric(LC50_reported)
  )

unrecognized_units <- wos_keep %>%
  filter(is.na(Unit_normalized)) %>%
  select(Audit_ID, Reference_ID, Species, Metal, Unit)

if (nrow(unrecognized_units) > 0L) {
  print(unrecognized_units, n = Inf)
  write_csv(
    unrecognized_units,
    file.path(qa_dir, "wos_unrecognized_units.csv")
  )
  stop(
    paste0(
      "Unrecognized WoS concentration units. ",
      "Inspect wos_unrecognized_units.csv."
    )
  )
}

stopifnot(
  n_distinct(wos_keep$Result_ID_final) == nrow(wos_keep),
  all(!is.na(wos_keep$Test_ID)),
  all(!is.na(wos_keep$Reference_ID)),
  all(!is.na(wos_keep$Unit_normalized)),

  # The Yang et al. ZnCl2 arm does not explicitly specify life stage.
  # Keep that value missing; all other retained stage values are required.
  all(
    !is.na(wos_keep$Lifestage_final) |
      (
        wos_keep$Reference_ID == "W_017" &
          wos_keep$DOI == "10.1016/j.chemosphere.2022.134099"
      )
  ),
  all(wos_keep$Duration_days_final %in% c(1, 2, 3, 4))
)

write_csv(
  wos_keep %>%
    select(
      Audit_ID,
      Reference_ID,
      Test_ID,
      Original_Result_ID = Result_ID,
      Final_Result_ID = Result_ID_final,
      Reference,
      Species,
      Metal,
      Duration_h
    ),
  file.path(qa_dir, "02_wos_id_crosswalk.csv")
)


# ============================================================
# 6. HARMONIZE THE RETAINED WoS LC50 VALUES
# ============================================================

atomic_weights <- tribble(
  ~Metal, ~Atomic_Weight_g_mol,
  "Ag", 107.8682,
  "Cd", 112.414,
  "Cr", 51.9961,
  "Cu", 63.546,
  "Hg", 200.592,
  "Ni", 58.6934,
  "Pb", 207.2,
  "Zn", 65.38
)

# Molecular weight used only where the source-reported LC50 is
# explicitly represented as HgCl2 compound mass.
MW_HgCl2_g_mol <- 271.496

wos_harmonized_calc <- wos_keep %>%
  left_join(atomic_weights, by = "Metal") %>%
  mutate(
    Reported_mass_mg_L = case_when(
      Unit_normalized == "ug/L" ~ LC50_reported / 1000,
      Unit_normalized == "mg/L" ~ LC50_reported,
      TRUE ~ NA_real_
    ),

    Basis_Class = case_when(
      str_detect(
        str_to_lower(as.character(Concentration_basis)),
        "compound-mass"
      ) ~ "COMPOUND_MASS",

      str_detect(
        str_to_lower(as.character(Concentration_basis)),
        "elemental"
      ) ~ "ELEMENTAL_METAL_MASS",

      TRUE ~ "UNRESOLVED"
    ),

    # Target-metal mass concentration (mg metal/L).
    LC50_target_mg_L_calc = case_when(
      Basis_Class == "ELEMENTAL_METAL_MASS" ~
        Reported_mass_mg_L,

      Basis_Class == "COMPOUND_MASS" & Metal == "Hg" ~
        Reported_mass_mg_L * Atomic_Weight_g_mol / MW_HgCl2_g_mol,

      TRUE ~ NA_real_
    ),

    # Target-metal molar concentration (umol metal/L).
    LC50_umol_L_calc = case_when(
      Basis_Class == "ELEMENTAL_METAL_MASS" ~
        LC50_target_mg_L_calc * 1000 / Atomic_Weight_g_mol,

      # HgCl2 contains one Hg atom per molecule, so umol Hg/L
      # equals umol HgCl2/L for the reported compound-mass LC50.
      Basis_Class == "COMPOUND_MASS" & Metal == "Hg" ~
        Reported_mass_mg_L * 1000 / MW_HgCl2_g_mol,

      TRUE ~ NA_real_
    ),

    ln_LC50_umol_L_calc = case_when(
      !is.na(LC50_umol_L_calc) & LC50_umol_L_calc > 0 ~
        log(LC50_umol_L_calc),
      TRUE ~ NA_real_
    ),

    Basis_Interpretation_final = case_when(
      Basis_Class == "ELEMENTAL_METAL_MASS" ~
        "WOS_SOURCE_VERIFIED_ELEMENTAL_METAL_MASS",
      Basis_Class == "COMPOUND_MASS" & Metal == "Hg" ~
        "WOS_SOURCE_VERIFIED_HGCL2_COMPOUND_MASS",
      TRUE ~
        "WOS_UNRESOLVED_BASIS"
    ),

    Harmonization_Rule_final = case_when(
      Basis_Class == "ELEMENTAL_METAL_MASS" ~
        "Reported metal-mass LC50 converted to umol/L using target-metal atomic weight",
      Basis_Class == "COMPOUND_MASS" & Metal == "Hg" ~
        "Reported HgCl2-mass LC50 converted using HgCl2 molecular weight (1 Hg per molecule)",
      TRUE ~
        "UNRESOLVED"
    )
  )

# No retained WoS row enters the quantitative augmentation unless
# harmonization can be reproduced deterministically.
unresolved_keep_rows <- wos_harmonized_calc %>%
  filter(
    is.na(Atomic_Weight_g_mol) |
      Basis_Class == "UNRESOLVED" |
      is.na(LC50_target_mg_L_calc) |
      is.na(LC50_umol_L_calc) |
      is.na(ln_LC50_umol_L_calc) |
      LC50_umol_L_calc <= 0
  )

write_csv(
  unresolved_keep_rows,
  file.path(qa_dir, "03_unresolved_wos_harmonization_rows.csv")
)

if (nrow(unresolved_keep_rows) > 0L) {
  stop(
    paste0(
      "At least one WoS row marked Include_in_combined=YES could not be harmonized. ",
      "Inspect 03_unresolved_wos_harmonization_rows.csv."
    )
  )
}


# ============================================================
# 7. MAP WoS RECORDS TO THE SCRIPT-01 ECOTOX SCHEMA
# ============================================================
# Only metadata explicitly captured during source verification are
# populated. Unextracted fields remain NA rather than being inferred.
# ============================================================

wos_55 <- wos_harmonized_calc %>%
  transmute(
    Raw_Row = NA_real_,
    Test_Number = NA_real_,
    Result_Number = NA_real_,
    Reference_Number = NA_real_,
    Metal = as.character(Metal),
    Chemical_Name = as.character(Chemical),
    Species = as.character(Species),
    Order = as.character(Order),
    Lifestage = Lifestage_final,
    Duration_days = Duration_days_final,
    Exposure_Type = "Waterborne",
    Conc_Type = Conc_Type_final,
    Mean_Op = "=",
    Mean = LC50_reported,
    Min = NA_real_,
    Max = NA_real_,
    Units = Unit_normalized,
    Analysis_Point = LC50_reported,
    Analysis_Units = Unit_normalized,
    Temperature = Temp_C,
    Salinity = Salinity,
    pH = NA_real_,
    Organic_Carbon = NA_real_,
    Title = NA_character_,
    Eligibility_Status = "WOS_SOURCE_VERIFIED_KEEP",
    Eligibility_Reason = "Supplementary WoS record; source verified and retained for quantitative augmentation",
    Atomic_Weight_g_mol = Atomic_Weight_g_mol,
    Basis_Interpretation = Basis_Interpretation_final,
    LC50_target_mg_L = LC50_target_mg_L_calc,
    LC50_umol_L = LC50_umol_L_calc,
    ln_LC50_umol_L = ln_LC50_umol_L_calc,
    Harmonization_Status = "HARMONIZED",
    Harmonization_Rule = Harmonization_Rule_final,
    Metal_Role = "PENDING_COMBINED_SUPPORT_AUDIT",
    Core_Pairwise_Ready = "PENDING_COMBINED_SUPPORT_AUDIT",

    # These are direct source-reported WoS values, not repaired
    # ECOTOX environmental fields.
    Temperature_Raw_Mean = Temp_C,
    Temperature_Raw_Min = NA_real_,
    Temperature_Raw_Max = NA_real_,
    Temperature_Raw_Units = if_else(!is.na(Temp_C), "C", NA_character_),
    Temperature_Repair_Status = if_else(
      !is.na(Temp_C),
      "SOURCE_REPORTED_WOS",
      "NOT_REPORTED_WOS"
    ),

    Salinity_Raw_Mean = Salinity,
    Salinity_Raw_Min = NA_real_,
    Salinity_Raw_Max = NA_real_,
    Salinity_Raw_Units = if_else(
      !is.na(Salinity),
      "PSU/ppt as reported",
      NA_character_
    ),
    Salinity_Repair_Status = if_else(
      !is.na(Salinity),
      "SOURCE_REPORTED_WOS",
      "NOT_REPORTED_WOS"
    ),

    pH_Raw_Mean = NA_real_,
    pH_Raw_Min = NA_real_,
    pH_Raw_Max = NA_real_,
    pH_Raw_Units = NA_character_,
    pH_Repair_Status = "NOT_EXTRACTED_WOS",

    Organic_Carbon_Raw_Mean = NA_real_,
    Organic_Carbon_Raw_Min = NA_real_,
    Organic_Carbon_Raw_Max = NA_real_,
    Organic_Carbon_Raw_Units = NA_character_,
    Organic_Carbon_Repair_Status = "NOT_EXTRACTED_WOS"
  )

stopifnot(
  identical(names(wos_55), required_ecotox_columns),
  nrow(wos_55) == 57,
  all(wos_55$Harmonization_Status == "HARMONIZED")
)


# ============================================================
# 8. ADD WoS PROVENANCE FIELDS
# ============================================================

wos_aug <- bind_cols(
  wos_55,
  wos_harmonized_calc %>%
    transmute(
      Source_Origin = "WoS_supplemental",
      Reference_ID = as.character(Reference_ID),
      Test_ID = as.character(Test_ID),
      Result_ID = as.character(Result_ID_final),
      External_Audit_ID = as.character(Audit_ID),
      External_DOI = as.character(DOI),
      External_Reference_Label = as.character(Reference),
      External_Stage_Sex = as.character(Stage_Sex),
      External_Concentration_Type_Detail = as.character(Concentration_type),
      External_Source_Location = as.character(Source_location),
      External_Source_Note = as.character(Note)
    )
)

stopifnot(
  n_distinct(wos_aug$Reference_ID) == 17,
  n_distinct(wos_aug$Test_ID) >= 17,
  n_distinct(wos_aug$Result_ID) == 57
)


# ============================================================
# 9. ALIGN THE TWO SOURCES AND CREATE THE COMBINED DATASET
# ============================================================

provenance_columns <- c(
  "Source_Origin",
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "External_Audit_ID",
  "External_DOI",
  "External_Reference_Label",
  "External_Stage_Sex",
  "External_Concentration_Type_Detail",
  "External_Source_Location",
  "External_Source_Note"
)

combined_column_order <- c(
  required_ecotox_columns,
  provenance_columns
)

ecotox_aug <- ecotox_aug %>%
  select(all_of(combined_column_order))

wos_aug <- wos_aug %>%
  select(all_of(combined_column_order))

combined <- bind_rows(
  ecotox_aug,
  wos_aug
)

stopifnot(
  nrow(combined) == 353,
  n_distinct(combined$Result_ID) == 353,
  sum(combined$Source_Origin == "ECOTOX") == 296,
  sum(combined$Source_Origin == "WoS_supplemental") == 57,
  sum(combined$Harmonization_Status == "HARMONIZED") == 304
)

# Corrected ECOTOX Sosnowski result must survive the merge exactly once.
stopifnot(sum(combined$Result_ID == "E_R42366") == 1L)

write_csv(
  combined %>% filter(Result_ID == "E_R42366"),
  file.path(qa_dir, "sosnowski_1979_merged_record.csv"),
  na = ""
)

# Confirm identity, units and intentionally missing metadata for the
# explicitly added Yang et al. (2022) Zn result.
yang_check <- combined %>%
  filter(Reference_ID == "W_017")

stopifnot(
  nrow(yang_check) == 1L,
  yang_check$Result_ID == "W_017_T01_R01",
  yang_check$External_DOI == "10.1016/j.chemosphere.2022.134099",
  yang_check$Species == "Tigriopus japonicus",
  yang_check$Metal == "Zn",
  yang_check$Duration_days == 4,
  abs(yang_check$LC50_target_mg_L - 7.03) < 1e-10,
  abs(yang_check$LC50_umol_L - 7.03 * 1000 / 65.38) < 1e-8,
  is.na(yang_check$Lifestage),
  is.na(yang_check$Temperature),
  is.na(yang_check$Salinity),
  is.na(yang_check$pH),
  yang_check$Conc_Type == "Unspecified"
)

write_csv(
  yang_check,
  file.path(qa_dir, "yang_2022_added_record.csv"),
  na = ""
)


# ============================================================
# 10. CROSS-SOURCE RESULT-LEVEL DUPLICATE QA
# ============================================================
# Flagging only: no row is automatically excluded here.
# Same Species x Metal x Duration and <=1% relative LC50 difference
# are surfaced for provenance review.
# ============================================================

possible_cross_source_duplicates <- wos_aug %>%
  select(
    WOS_Reference_ID = Reference_ID,
    WOS_Result_ID = Result_ID,
    WOS_Reference = External_Reference_Label,
    Species,
    Metal,
    Duration_days,
    WOS_LC50_umol_L = LC50_umol_L
  ) %>%
  inner_join(
    ecotox_aug %>%
      select(
        ECOTOX_Reference_ID = Reference_ID,
        ECOTOX_Result_ID = Result_ID,
        ECOTOX_Reference_Number = Reference_Number,
        Species,
        Metal,
        Duration_days,
        ECOTOX_LC50_umol_L = LC50_umol_L
      ) %>%
      filter(!is.na(ECOTOX_LC50_umol_L)),
    by = c("Species", "Metal", "Duration_days"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    Relative_Difference = abs(
      WOS_LC50_umol_L - ECOTOX_LC50_umol_L
    ) / pmax(
      abs(WOS_LC50_umol_L),
      abs(ECOTOX_LC50_umol_L)
    )
  ) %>%
  filter(Relative_Difference <= 0.01) %>%
  arrange(
    Relative_Difference,
    Species,
    Metal,
    Duration_days
  )

write_csv(
  possible_cross_source_duplicates,
  file.path(qa_dir, "04_possible_cross_source_result_duplicates.csv")
)


# ============================================================
# 11. DESCRIPTIVE MERGE-SUPPORT SNAPSHOT
# ============================================================
# This is merge QA only. It does not replace the formal support
# and comparability audit performed by Script 03.
# ============================================================

merge_support_snapshot <- combined %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0
  ) %>%
  group_by(
    Source_Origin,
    Metal,
    Duration_days
  ) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(
    Metal,
    Duration_days,
    Source_Origin
  )

write_csv(
  merge_support_snapshot,
  file.path(qa_dir, "05_merge_support_snapshot.csv")
)

merge_summary <- tibble(
  Metric = c(
    "ECOTOX rows",
    "WoS retained rows",
    "Combined rows",
    "ECOTOX References",
    "WoS References",
    "Combined References",
    "Combined harmonized quantitative rows"
  ),
  Value = c(
    nrow(ecotox_aug),
    nrow(wos_aug),
    nrow(combined),
    n_distinct(ecotox_aug$Reference_ID),
    n_distinct(wos_aug$Reference_ID),
    n_distinct(combined$Reference_ID),
    sum(combined$Harmonization_Status == "HARMONIZED")
  )
)

write_csv(
  merge_summary,
  file.path(qa_dir, "06_merge_summary.csv")
)


# ============================================================
# 12. WRITE DATA OUTPUTS
# ============================================================

wos_csv <- file.path(
  data_dir,
  "wos_harmonized_verified.csv"
)

combined_xlsx <- file.path(
  data_dir,
  "combined_ecotox_wos_harmonized.xlsx"
)

# Curated WoS augmentation retained as an audit-friendly derived file.
write_csv(
  wos_aug,
  wos_csv,
  na = ""
)

# This is the canonical machine-readable dataset used by all
# subsequent core analysis scripts.
write_csv(
  combined,
  combined_csv,
  na = ""
)

# Excel is a convenience/audit export. The CSV above remains the
# canonical machine-readable downstream input.
excel_sheets <- list(
  Combined = combined,
  WoS_Harmonized = wos_aug,
  WoS_Source_Audit = wos_source_audit,
  WoS_All_Audit_Rows = wos_results_raw,
  Merge_QA = merge_summary,
  Possible_Duplicates = possible_cross_source_duplicates
)

xlsx_written <- FALSE

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    excel_sheets,
    path = combined_xlsx
  )
  xlsx_written <- TRUE
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(
    excel_sheets,
    file = combined_xlsx,
    overwrite = TRUE
  )
  xlsx_written <- TRUE
} else {
  warning(
    paste0(
      "CSV outputs were written successfully, but no XLSX writer is installed. ",
      "Install 'writexl' or 'openxlsx' and rerun this script only if the ",
      "convenience audit workbook is required."
    )
  )
}


# ============================================================
# 13. SAVE SESSION + INPUT/OUTPUT MANIFEST
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo.txt")
)

manifest <- tibble(
  Role = c(
    "Script-01 ECOTOX processed dataset",
    "Curated WoS source-verification workbook",
    "Harmonized WoS quantitative augmentation",
    "Combined analysis-ready dataset (canonical CSV)",
    "Combined audit workbook (optional XLSX)"
  ),
  File = c(
    ecotox_file,
    wos_file,
    wos_csv,
    combined_csv,
    combined_xlsx
  ),
  Created = c(
    FALSE,
    FALSE,
    file.exists(wos_csv),
    file.exists(combined_csv),
    xlsx_written && file.exists(combined_xlsx)
  )
)

write_csv(
  manifest,
  file.path(output_dir, "input_output_manifest.csv")
)


# ============================================================
# 14. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n============================================\n")
cat("02 - WoS VERIFIED RECORD INTEGRATION COMPLETE\n")
cat("============================================\n\n")

print(merge_summary)

cat("\nWoS retained Results by metal:\n")
print(
  wos_aug %>%
    count(Metal, sort = TRUE)
)

cat("\nWoS retained Results by metal and duration:\n")
print(
  wos_aug %>%
    count(Metal, Duration_days, sort = TRUE)
)

cat(
  "\nPossible near-identical ECOTOX/WoS result matches requiring provenance review: ",
  nrow(possible_cross_source_duplicates),
  "\n",
  sep = ""
)

cat("\nCanonical combined CSV:\n  ", combined_csv, "\n", sep = "")

if (xlsx_written) {
  cat("Optional audit workbook:\n  ", combined_xlsx, "\n", sep = "")
} else {
  cat("Optional audit workbook: not written (no XLSX writer installed).\n")
}

cat("\nFinal combined checks:\n")
cat("  Total rows: ", nrow(combined), "\n", sep = "")
cat(
  "  Harmonized quantitative Results: ",
  sum(combined$Harmonization_Status == "HARMONIZED"),
  "\n",
  sep = ""
)
cat("  ECOTOX rows: ", sum(combined$Source_Origin == "ECOTOX"), "\n", sep = "")
cat(
  "  WoS supplemental rows: ",
  sum(combined$Source_Origin == "WoS_supplemental"),
  "\n",
  sep = ""
)

cat("\nSTATUS: 02 INTEGRATE WOS VERIFIED RECORDS = PASS\n")
