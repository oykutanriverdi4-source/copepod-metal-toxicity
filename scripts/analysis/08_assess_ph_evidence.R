# ============================================================
# ESS THESIS REPRODUCIBILITY WORKFLOW
# 08_assess_ph_evidence.R
# ============================================================
# Purpose
#   Assess quantitative pH metadata and repeated-pH contexts in the
#   final combined harmonized dataset, including the source-verified
#   second-pass pH extraction for the 17 retained WoS References.
#
# This script:
#   1) retains repaired ECOTOX pH metadata from the combined dataset;
#   2) maps source-verified pH metadata to the 57 retained WoS Results;
#   3) builds one unified pH metadata layer across 304 harmonized Results;
#   4) identifies within-Reference repeated-pH candidate contexts;
#   5) checks whether pH variation can be separated from salinity and
#      temperature variation; and
#   6) keeps Wei et al. (2021) as direct CO2-driven acidification evidence,
#      not as a numerical pH coefficient.
#
# No pooled, developmental-stage, or taxonomic-order model is refitted.
# No LC50 value or harmonization decision is changed.
# ============================================================

# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)
library(readxl)

options(width = 220)

output_root <- file.path(
  "outputs",
  "08_assess_ph_evidence"
)

audit_dir  <- file.path(output_root, "01_audit")
table_dir  <- file.path(output_root, "02_tables")
figure_dir <- file.path(output_root, "03_figures")
object_dir <- file.path(output_root, "04_objects")

for (d in c(output_root, audit_dir, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 1. INPUT FILES
# ============================================================

combined_file <- file.path(
  "data",
  "processed",
  "combined_ecotox_wos_harmonized.csv"
)

wos_ledger_file <- file.path(
  "data",
  "curated",
  "wos",
  "wos_source_verified_records.xlsx"
)

stopifnot(
  file.exists(combined_file),
  file.exists(wos_ledger_file)
)


# ============================================================
# 2. LOAD COMBINED HARMONIZED DATA
# ============================================================

dat <- read_csv(
  combined_file,
  show_col_types = FALSE
)

required_columns <- c(
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "Reference_Number",
  "Source_Origin",
  "Species",
  "Metal",
  "Lifestage",
  "External_Stage_Sex",
  "Duration_days",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Harmonization_Status",
  "LC50_umol_L",
  "Temperature",
  "Salinity",
  "pH"
)

stopifnot(all(required_columns %in% names(dat)))

harm <- dat %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0,
    !is.na(Duration_days)
  ) %>%
  mutate(
    Lifestage =
      na_if(trimws(as.character(Lifestage)), ""),
    
    External_Stage_Sex =
      na_if(trimws(as.character(External_Stage_Sex)), ""),
    
    Stage_context = coalesce(
      External_Stage_Sex,
      Lifestage,
      "NOT_REPORTED"
    ),
    
    Exposure_context = coalesce(
      na_if(trimws(as.character(Exposure_Type)), ""),
      "NOT_REPORTED"
    ),
    
    Conc_context = coalesce(
      na_if(trimws(as.character(Conc_Type)), ""),
      "NOT_REPORTED"
    ),
    
    Basis_context = coalesce(
      na_if(trimws(as.character(Basis_Interpretation)), ""),
      "NOT_REPORTED"
    )
  )

stopifnot(
  nrow(harm) == 304,
  sum(harm$Source_Origin == "ECOTOX") == 247,
  sum(harm$Source_Origin == "WoS_supplemental") == 57
)


# ============================================================
# 3. LOAD FROZEN WoS SOURCE + RESULT LEDGERS
# ============================================================

source_audit <- read_excel(
  wos_ledger_file,
  sheet = "Source_Audit"
)

results_all <- read_excel(
  wos_ledger_file,
  sheet = "Results"
)

required_source_cols <- c(
  "Audit_ID",
  "DOI"
)

required_result_cols <- c(
  "Audit_ID",
  "Reference_ID",
  "Result_ID",
  "Include_in_combined"
)

stopifnot(
  all(required_source_cols %in% names(source_audit)),
  all(required_result_cols %in% names(results_all))
)

results_keep <- results_all %>%
  filter(
    str_to_upper(
      str_trim(
        as.character(Include_in_combined)
      )
    ) == "YES"
  )

stopifnot(
  nrow(results_keep) == 57,
  n_distinct(results_keep$Reference_ID) == 17
)

normalize_doi <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_trim() %>%
    str_remove("^https?://(dx\\.)?doi\\.org/") %>%
    str_remove("^doi:\\s*") %>%
    na_if("")
}


# ============================================================
# 4. SOURCE-VERIFIED SECOND-PASS pH AUDIT FOR 17 WoS REFERENCES
# ============================================================

ph_audit <- tribble(
  ~Reference, ~DOI, ~Retained_Results, ~Metal,
  ~pH_Status, ~pH_Reported_Text, ~pH_Low, ~pH_High, ~pH_Central,
  ~pH_Scale, ~Deliberate_pH_Manipulation, ~Source_Location,
  ~Extraction_Note,
  
  "Yang et al. 2022",
  "10.1016/j.chemosphere.2022.134099",
  1, "Zn",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NOT_REPORTED",
  "Section 3.3: separate ZnCl2 LC50; Methods and supplementary information",
  "Exposure pH was not separately established for the ZnCl2 arm. Culture/TWP/leachate pH was not transferred to this record.",
  
  "Panneerselvam et al. 2018",
  "10.1007/s00128-018-2279-7",
  2, "Ni",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Methods: copepod culture / acute bioassay conditions",
  "Temperature and salinity are reported, but no numeric exposure-water pH was located in the main article.",
  
  "Kadiene et al. 2017",
  "10.1007/s10646-017-1848-6",
  9, "Cd",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Materials and methods: copepod cultures / cadmium exposure",
  "Salinity and species-specific temperatures are reported; no numeric test-water pH was located in the main article.",
  
  "Tollefsen et al. 2017",
  "10.1080/15287394.2017.1352198",
  4, "Hg",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Materials and methods: Experiments",
  "Acute test temperature and filtered seawater are described, but no numeric exposure-water pH was located in the main article.",
  
  "Tlili et al. 2016",
  "10.1016/j.chemosphere.2015.10.057",
  4, "Ni",
  "NUMERIC_PH_REPORTED",
  "8.4 ?? 0.2", 8.2, 8.6, 8.4,
  "Not reported", "NO",
  "Section 2.4 Lethal concentration determination experiments",
  "Acute lethality tests were conducted at pH 8.4??0.2. Stock-solution acidification is not treated as exposure pH.",
  
  "Lavorante et al. 2013",
  "10.1016/j.ecoenv.2013.05.010",
  1, "Zn",
  "MEASURED_BUT_NUMERIC_NOT_REPORTED",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Section 2.3 Description of test",
  "The paper states that pH was determined before the test, but no numeric test pH is given. Analytical-solution pH is not exposure pH.",
  
  "Farkas et al. 2020",
  "10.1016/j.aquatox.2020.105582",
  1, "Ag",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Section 2.2 Acute Ag toxicity test",
  "No numeric exposure-water pH was located in the main article.",
  
  "Vimercati et al. 2020",
  "10.3389/fpubh.2020.00192",
  3, "Zn",
  "NUMERIC_PH_REPORTED",
  "8 ?? 0.3", 7.7, 8.3, 8.0,
  "Not reported", "NO",
  "Table 2: Experimental conditions during acute toxicity tests",
  "For Tigriopus fulvus acute tests, pH was 8??0.3.",
  
  "Cherkashin 2020",
  "10.1134/S1063074020030037",
  4, "Zn",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Materials and methods",
  "No numeric pH for the copepod exposure was located.",
  
  "Kadiene et al. 2019",
  "10.1016/j.chemosphere.2019.05.220",
  2, "Cd",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Section 2.2.1 Acute exposure of P. annandalei to cadmium",
  "Salinity and temperature are reported, but no numeric test-water pH was located.",
  
  "Tlili et al. 2019",
  "10.1007/s41742-019-00202-y",
  4, "Hg",
  "NUMERIC_PH_REPORTED",
  "8.4 ?? 0.2", 8.2, 8.6, 8.4,
  "Not reported", "NO",
  "Materials and Methods: Lethal Concentration Experiments",
  "Acute lethality tests were run at pH 8.4??0.2.",
  
  "Charry et al. 2019",
  "10.1016/j.ecoenv.2019.03.022",
  2, "Cu",
  "NUMERIC_PH_REPORTED",
  "8 ?? 0.3", 7.7, 8.3, 8.0,
  "Not reported", "NO",
  "Results 3.3: validation parameters",
  "Validation parameters remained within pH 8??0.3.",
  
  "Yi et al. 2019",
  "10.1016/j.cbpc.2019.04.014",
  1, "Cd",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Sections 2.1???2.2 Chemicals and toxicity test",
  "No numeric exposure-water pH was located in the main article.",
  
  "Zidour et al. 2019",
  "10.1016/j.chemosphere.2018.12.148",
  11, "Cu/Cd/Ni",
  "NUMERIC_PH_REPORTED",
  "8.4 ?? 0.2", 8.2, 8.6, 8.4,
  "Not reported", "NO",
  "Section 2.3 Acute toxicity tests",
  "Culture water during the acute tests was pH 8.4??0.2.",
  
  "Biandolino et al. 2018",
  "10.1016/j.ecoenv.2017.08.041",
  1, "Cu",
  "NUMERIC_PH_REPORTED",
  "8 ?? 0.3", 7.7, 8.3, 8.0,
  "Not reported", "NO",
  "Table 2: Summary of experimental conditions",
  "Acute T. fulvus test conditions included pH 8??0.3.",
  
  "??verjordet et al. 2014",
  "10.1016/j.aquatox.2014.06.019",
  4, "Hg",
  "NO_NUMERIC_PH_FOUND",
  NA_character_, NA_real_, NA_real_, NA_real_,
  "Not reported", "NO",
  "Section 2.2 Acute toxicity testing",
  "No numeric exposure-water pH was located in the main article.",
  
  "Mohammed et al. 2010",
  "10.1007/s10646-010-0471-6",
  3, "Ni",
  "NUMERIC_PH_REPORTED",
  "7.90???8.25", 7.90, 8.25, 8.075,
  "Not reported", "NO",
  "Materials and methods: culture conditions + Acute toxicity test",
  "Culture water pH was 7.90???8.25; acute-test conditions were stated to be the same except for temperature and salinity."
) %>%
  mutate(
    DOI_key = normalize_doi(DOI)
  )

stopifnot(
  nrow(ph_audit) == 17,
  n_distinct(ph_audit$DOI_key) == 17
)


# ============================================================
# 5. MAP SECOND-PASS pH TO THE 57 WoS RESULTS
# ============================================================

source_key <- source_audit %>%
  transmute(
    Audit_ID,
    DOI_key = normalize_doi(DOI)
  ) %>%
  filter(!is.na(DOI_key))

wos_ph_result_level <- results_keep %>%
  left_join(
    source_key,
    by = "Audit_ID",
    relationship = "many-to-one"
  ) %>%
  left_join(
    ph_audit %>%
      select(
        DOI_key,
        pH_Status,
        pH_Reported_Text,
        pH_Low,
        pH_High,
        pH_Central,
        pH_Scale,
        Deliberate_pH_Manipulation,
        Source_Location,
        Extraction_Note
      ),
    by = "DOI_key",
    relationship = "many-to-one"
  ) %>%
  mutate(
    WoS_pH_numeric_known =
      pH_Status == "NUMERIC_PH_REPORTED" &
      !is.na(pH_Central)
  )

stopifnot(
  nrow(wos_ph_result_level) == 57,
  all(!is.na(wos_ph_result_level$pH_Status)),
  sum(wos_ph_result_level$WoS_pH_numeric_known) == 28,
  n_distinct(
    wos_ph_result_level$Reference_ID[
      wos_ph_result_level$WoS_pH_numeric_known
    ]
  ) == 7,
  all(if_else(wos_ph_result_level$Reference_ID == "W_017",
              wos_ph_result_level$Deliberate_pH_Manipulation == "NOT_REPORTED",
              wos_ph_result_level$Deliberate_pH_Manipulation == "NO"))
)

wos_ph_reference_level <- wos_ph_result_level %>%
  group_by(
    Reference_ID,
    DOI_key
  ) %>%
  summarise(
    Results = n(),
    pH_Status = first(pH_Status),
    pH_Reported_Text = first(pH_Reported_Text),
    pH_Low = first(pH_Low),
    pH_High = first(pH_High),
    pH_Central = first(pH_Central),
    pH_Scale = first(pH_Scale),
    Deliberate_pH_Manipulation =
      first(Deliberate_pH_Manipulation),
    Source_Location = first(Source_Location),
    Extraction_Note = first(Extraction_Note),
    .groups = "drop"
  )


# ============================================================
# 6. BUILD ONE UNIFIED pH METADATA LAYER FOR ALL 304 RESULTS
# ============================================================
# IMPORTANT ID NOTE
#
# Script 00b rebuilt WoS Result_ID values deterministically during the
# merge. Therefore the Result_ID values in the human audit workbook are
# not guaranteed to be identical to the final combined Result_ID values.
#
# The second-pass pH metadata extracted here are REFERENCE-LEVEL test-
# water metadata (constant/range within each retained quantitative
# Reference), not result-specific pH treatments.
#
# Therefore pH is mapped to the final combined dataset by universal
# Reference_ID, which is stable across the audit workbook and Script 00b.
# The original 16 references do not provide a deliberate pH contrast; Yang has no arm-specific pH.
# This mapping remains appropriate because none of the original 16
# quantitative WoS ADD References deliberately manipulated pH.

wos_ph_map <- wos_ph_reference_level %>%
  select(
    Reference_ID,
    WoS_pH_Status = pH_Status,
    WoS_pH_Reported_Text = pH_Reported_Text,
    WoS_pH_Low = pH_Low,
    WoS_pH_High = pH_High,
    WoS_pH_Central = pH_Central,
    WoS_pH_Scale = pH_Scale,
    WoS_pH_Source_Location = Source_Location,
    WoS_pH_Extraction_Note = Extraction_Note
  )

stopifnot(
  nrow(wos_ph_map) == 17,
  n_distinct(wos_ph_map$Reference_ID) == 17
)

harm_ph <- harm %>%
  left_join(
    wos_ph_map,
    by = "Reference_ID",
    relationship = "many-to-one"
  ) %>%
  mutate(
    pH_final = case_when(
      Source_Origin == "ECOTOX" ~ pH,
      Source_Origin == "WoS_supplemental" ~ WoS_pH_Central,
      TRUE ~ NA_real_
    ),
    
    pH_metadata_origin = case_when(
      Source_Origin == "ECOTOX" &
        !is.na(pH) ~
        "ECOTOX_REPORTED_OR_REPAIRED",
      
      Source_Origin == "WoS_supplemental" &
        !is.na(WoS_pH_Central) ~
        "WOS_SECOND_PASS_SOURCE_VERIFIED",
      
      Source_Origin == "WoS_supplemental" &
        is.na(WoS_pH_Central) ~
        "WOS_SECOND_PASS_NO_NUMERIC_PH",
      
      TRUE ~
        "NO_NUMERIC_PH"
    ),
    
    pH_scale_final = case_when(
      Source_Origin == "WoS_supplemental" ~
        WoS_pH_Scale,
      TRUE ~ NA_character_
    )
  )

stopifnot(
  nrow(harm_ph) == 304,
  sum(
    harm_ph$Source_Origin == "WoS_supplemental" &
      !is.na(harm_ph$pH_final)
  ) == 28,
  n_distinct(
    harm_ph$Reference_ID[
      harm_ph$Source_Origin == "WoS_supplemental" &
        !is.na(harm_ph$pH_final)
    ]
  ) == 7
)


# Yang must remain present, with no transferred numeric pH.
stopifnot(
  sum(harm_ph$Reference_ID == "W_017") == 1,
  all(is.na(harm_ph$pH_final[harm_ph$Reference_ID == "W_017"]))
)

# ============================================================
# 7. pH METADATA COVERAGE
# ============================================================

pH_metadata_coverage_by_source_and_metal <- harm_ph %>%
  group_by(Source_Origin, Metal) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    
    pH_known = sum(!is.na(pH_final)),
    pH_percent = 100 * mean(!is.na(pH_final)),
    
    .groups = "drop"
  ) %>%
  arrange(Source_Origin, Metal)

pH_metadata_coverage_overall <- harm_ph %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    
    pH_known = sum(!is.na(pH_final)),
    pH_percent = 100 * mean(!is.na(pH_final)),
    
    ECOTOX_pH_known =
      sum(
        Source_Origin == "ECOTOX" &
          !is.na(pH_final)
      ),
    
    WoS_pH_known =
      sum(
        Source_Origin == "WoS_supplemental" &
          !is.na(pH_final)
      )
  )

# WoS second-pass benchmark is frozen.
stopifnot(
  pH_metadata_coverage_overall$WoS_pH_known == 28
)


# ============================================================
# 8. FINAL WITHIN-REFERENCE pH CANDIDATE SCREEN
# ============================================================

base_context_vars <- c(
  "Reference_ID",
  "Source_Origin",
  "Species",
  "Metal",
  "Stage_context",
  "Duration_days",
  "Exposure_context",
  "Conc_context",
  "Basis_context"
)

pH_context_screen <- harm_ph %>%
  group_by(across(all_of(base_context_vars))) %>%
  summarise(
    Results = n(),
    
    pH_reported_n =
      sum(!is.na(pH_final)),
    
    pH_levels =
      n_distinct(pH_final, na.rm = TRUE),
    
    pH_values = paste(
      sort(unique(pH_final[!is.na(pH_final)])),
      collapse = " | "
    ),
    
    Temperature_levels =
      n_distinct(Temperature, na.rm = TRUE),
    
    Temperature_values = paste(
      sort(unique(Temperature[!is.na(Temperature)])),
      collapse = " | "
    ),
    
    Salinity_levels =
      n_distinct(Salinity, na.rm = TRUE),
    
    Salinity_values = paste(
      sort(unique(Salinity[!is.na(Salinity)])),
      collapse = " | "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    pH_candidate =
      pH_levels >= 2,
    
    Structural_status = case_when(
      pH_levels < 2 ~
        "NO_WITHIN_CONTEXT_PH_VARIATION",
      
      Temperature_levels <= 1 &
        Salinity_levels <= 1 ~
        "POTENTIALLY_ISOLATABLE_PH",
      
      Temperature_levels <= 1 &
        Salinity_levels >= 2 ~
        "PH_COVARIES_WITH_SALINITY",
      
      Temperature_levels >= 2 &
        Salinity_levels <= 1 ~
        "PH_COVARIES_WITH_TEMPERATURE",
      
      Temperature_levels >= 2 &
        Salinity_levels >= 2 ~
        "PH_COVARIES_WITH_MULTIPLE_ENVIRONMENTAL_FACTORS",
      
      TRUE ~
        "REVIEW"
    )
  )

pH_candidate_contexts <- pH_context_screen %>%
  filter(pH_candidate) %>%
  arrange(
    Structural_status,
    Reference_ID,
    Species,
    Metal
  )

pH_candidate_rows <- harm_ph %>%
  semi_join(
    pH_candidate_contexts %>%
      select(all_of(base_context_vars)),
    by = base_context_vars
  ) %>%
  select(
    Reference_ID,
    Source_Origin,
    Test_ID,
    Result_ID,
    Species,
    Metal,
    Stage_context,
    Duration_days,
    Exposure_context,
    Conc_context,
    Basis_context,
    Temperature,
    Salinity,
    pH_final,
    pH_metadata_origin,
    LC50_umol_L
  ) %>%
  arrange(
    Reference_ID,
    Species,
    Metal,
    pH_final
  )

# The pH metadata layer should still contain the same three
# ECOTOX repeated-pH contexts and no new WoS within-Reference pH
# treatment series.
stopifnot(
  nrow(pH_candidate_contexts) == 3,
  setequal(
    pH_candidate_contexts$Reference_ID,
    c("E_166133", "E_107036", "E_17219")
  ),
  all(pH_candidate_contexts$Salinity_levels >= 2),
  !any(
    pH_candidate_contexts$Structural_status ==
      "POTENTIALLY_ISOLATABLE_PH"
  ),
  !any(
    pH_candidate_contexts$Source_Origin ==
      "WoS_supplemental"
  )
)


# ============================================================
# 9. pH-SALINITY CONFOUNDING SUMMARY
# ============================================================

pH_salinity_confounded_contexts <- pH_candidate_rows %>%
  group_by(
    Reference_ID,
    Species,
    Metal,
    Stage_context,
    Duration_days,
    Exposure_context,
    Conc_context,
    Basis_context
  ) %>%
  summarise(
    n_points = n(),
    
    pH_min = min(pH_final),
    pH_max = max(pH_final),
    
    Salinity_min = min(Salinity),
    Salinity_max = max(Salinity),
    
    Temperature_levels =
      n_distinct(Temperature, na.rm = TRUE),
    
    Pearson_r_pH_salinity = cor(
      pH_final,
      Salinity,
      use = "complete.obs"
    ),
    
    Decision =
      "NOT_ISOLATABLE_AS_PH_EFFECT",
    
    Reason = paste(
      "pH changes together with salinity within this context;",
      "some contexts also contain minor temperature variation.",
      "The rows therefore cannot identify an independent pH effect."
    ),
    
    .groups = "drop"
  ) %>%
  arrange(Reference_ID)

stopifnot(
  nrow(pH_salinity_confounded_contexts) == 3,
  all(
    pH_salinity_confounded_contexts$Decision ==
      "NOT_ISOLATABLE_AS_PH_EFFECT"
  )
)


# ============================================================
# 10. WEI ET AL. (2021): DIRECT ACIDIFICATION EVIDENCE
# ============================================================
# This remains OUTSIDE the quantitative ADD master.
#
# Published acute comparison:
# 400 vs 1000 uatm pCO2
# 48-h Cd LC50: 12.03 vs 9.08 mg/L
#
# Numerical pH values + pH scale are NOT source-verified here.
# Therefore no pH value is inferred from pCO2.

wei2021_acidification <- tribble(
  ~Acidification_Comparison_ID,
  ~Reference_ID,
  ~Reference,
  ~DOI,
  ~Species,
  ~Metal,
  ~Duration_h,
  ~pCO2_uatm,
  ~pH_value,
  ~pH_scale,
  ~LC50_reported_mg_L,
  ~Evidence_role,
  ~Source_note,
  
  "PH_EXT_01",
  "ENV_WEI2021",
  "Wei et al. (2021)",
  "10.1016/j.marpolbul.2021.113145",
  "Tigriopus japonicus",
  "Cd",
  48,
  400,
  NA_real_,
  NA_character_,
  12.03,
  "DIRECT_ACIDIFICATION_NOT_NUMERIC_PH",
  paste(
    "Published source reports 400 uatm pCO2 and a 48-h Cd LC50 of 12.03 mg/L.",
    "A verified numerical pH value/scale is not used in this audit."
  ),
  
  "PH_EXT_01",
  "ENV_WEI2021",
  "Wei et al. (2021)",
  "10.1016/j.marpolbul.2021.113145",
  "Tigriopus japonicus",
  "Cd",
  48,
  1000,
  NA_real_,
  NA_character_,
  9.08,
  "DIRECT_ACIDIFICATION_NOT_NUMERIC_PH",
  paste(
    "Published source reports 1000 uatm pCO2 and a 48-h Cd LC50 of 9.08 mg/L.",
    "A verified numerical pH value/scale is not used in this audit."
  )
)

wei2021_acidification_summary <- wei2021_acidification %>%
  arrange(pCO2_uatm) %>%
  summarise(
    Comparison_ID =
      first(Acidification_Comparison_ID),
    
    Reference_ID =
      first(Reference_ID),
    
    Reference =
      first(Reference),
    
    Species =
      first(Species),
    
    Metal =
      first(Metal),
    
    Duration_h =
      first(Duration_h),
    
    low_pCO2_uatm =
      min(pCO2_uatm),
    
    high_pCO2_uatm =
      max(pCO2_uatm),
    
    LC50_low_pCO2_mg_L =
      LC50_reported_mg_L[which.min(pCO2_uatm)],
    
    LC50_high_pCO2_mg_L =
      LC50_reported_mg_L[which.max(pCO2_uatm)],
    
    ratio_highCO2_lowCO2 =
      LC50_high_pCO2_mg_L /
      LC50_low_pCO2_mg_L,
    
    ln_ratio_highCO2_lowCO2 =
      log(ratio_highCO2_lowCO2),
    
    Direction = case_when(
      ratio_highCO2_lowCO2 < 1 ~
        "Higher pCO2 / acidification -> lower LC50",
      
      ratio_highCO2_lowCO2 > 1 ~
        "Higher pCO2 / acidification -> higher LC50",
      
      TRUE ~
        "No change"
    ),
    
    pH_numeric_status =
      "NOT SOURCE-VERIFIED",
    
    Final_use =
      "DIRECT ACIDIFICATION CONTEXT; NOT A QUANTITATIVE pH COEFFICIENT"
  )


# ============================================================
# 11. pH DECISION
# ============================================================

isolatable_pH_contexts <- pH_candidate_contexts %>%
  filter(
    Structural_status ==
      "POTENTIALLY_ISOLATABLE_PH"
  )

pH_final_decision <- tibble(
  Evidence_layer = c(
    "Unified ECOTOX + WoS quantitative pH metadata",
    "Within-Reference repeated-pH contexts",
    "Quantitative WoS ADD references",
    "Supplementary deliberate acidification literature",
    "Formal pH model / coefficient"
  ),
  
  Support = c(
    paste0(
      pH_metadata_coverage_overall$pH_known,
      " / ",
      pH_metadata_coverage_overall$Results,
      " Results with numeric pH metadata; ",
      pH_metadata_coverage_overall$ECOTOX_pH_known,
      " ECOTOX + ",
      pH_metadata_coverage_overall$WoS_pH_known,
      " WoS"
    ),
    
    paste0(
      nrow(pH_candidate_contexts),
      " candidate contexts; ",
      nrow(isolatable_pH_contexts),
      " potentially isolatable"
    ),
    
    "7 / 17 References with numeric exposure pH; 28 / 57 retained Results; no verified deliberate pH manipulation; Yang ZnCl2 arm status not separately reported",
    
    "1 direct Reference (Wei et al. 2021)",
    
    "Unsupported"
  ),
  
  Decision = c(
    "METADATA_COVERAGE_IMPROVED_BUT_INCOMPLETE",
    "NO_ISOLATABLE_DIRECT_PH_CONTEXT",
    "NO_NEW_WITHIN_REFERENCE_PH_CONTRAST",
    "DIRECT_ACIDIFICATION_CONTEXT_ONLY",
    "DO_NOT_FIT"
  ),
  
  Interpretation = c(
    paste(
      "The WoS second pass improves pH metadata coverage substantially,",
      "but pH remains incompletely and non-uniformly reported across sources."
    ),
    
    paste(
      "Three repeated-pH ECOTOX contexts are present (E_166133, E_107036, E_17219),",
      "but salinity changes together with pH in all three;",
      "some also contain minor temperature variation."
    ),
    
    paste(
      "Numeric test-water pH was source-verifiable for seven quantitative WoS References,",
      "but it describes constant/range test conditions rather than a manipulated pH series."
    ),
    
    paste(
      "Higher pCO2 was associated with a lower 48-h Cd LC50,",
      "but the comparison is retained as CO2-driven acidification evidence;",
      "it is not converted into a pH-response slope."
    ),
    
    paste(
      "The evidence does not support a general pH coefficient.",
      "pH remains a relevant water-chemistry / carbonate-system context variable",
      "rather than an independently estimated predictor in this thesis."
    )
  )
)


# ============================================================
# 12. QA FIGURE: E_166133 CONFOUNDING EXAMPLE
# ============================================================

e166133_pH_context <- pH_candidate_rows %>%
  filter(Reference_ID == "E_166133") %>%
  arrange(Salinity)

e166133_confounded_plot <- ggplot(
  e166133_pH_context,
  aes(
    x = Salinity,
    y = pH_final
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(
    aes(
      label = paste0(
        "LC50=",
        round(LC50_umol_L, 2)
      )
    ),
    nudge_y = 0.035,
    size = 3.1
  ) +
  labs(
    title =
      "pH-salinity alignment in ECOTOX Reference E_166133",
    
    subtitle =
      "Confounding audit only: pH and salinity change together",
    
    x =
      "Salinity",
    
    y =
      "Reported pH",
    
    caption =
      "These observations are not treated as an independent pH effect."
  ) +
  theme_classic(base_size = 10)

ggsave(
  file.path(
    figure_dir,
    "E166133_pH_salinity_confounding_audit.png"
  ),
  e166133_confounded_plot,
  width = 7.8,
  height = 5.0,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 13. WRITE FINAL OUTPUTS
# ============================================================

write_csv(
  ph_audit,
  file.path(
    audit_dir,
    "wos_second_pass_ph_reference_audit.csv"
  ),
  na = ""
)

write_csv(
  wos_ph_result_level,
  file.path(
    audit_dir,
    "wos_second_pass_ph_result_level.csv"
  ),
  na = ""
)

write_csv(
  wos_ph_reference_level,
  file.path(
    audit_dir,
    "wos_second_pass_ph_reference_level.csv"
  ),
  na = ""
)

write_csv(
  harm_ph,
  file.path(
    audit_dir,
    "combined_304_with_ph_metadata_layer.csv"
  ),
  na = ""
)

write_csv(
  pH_metadata_coverage_by_source_and_metal,
  file.path(
    table_dir,
    "ph_metadata_coverage_by_source_and_metal.csv"
  )
)

write_csv(
  pH_metadata_coverage_overall,
  file.path(
    table_dir,
    "ph_metadata_coverage_overall.csv"
  )
)

write_csv(
  pH_context_screen,
  file.path(
    table_dir,
    "ph_all_context_screen.csv"
  )
)

write_csv(
  pH_candidate_contexts,
  file.path(
    table_dir,
    "ph_within_reference_candidate_contexts.csv"
  )
)

write_csv(
  pH_candidate_rows,
  file.path(
    table_dir,
    "ph_within_reference_candidate_rows.csv"
  )
)

write_csv(
  pH_salinity_confounded_contexts,
  file.path(
    table_dir,
    "ph_salinity_confounded_contexts.csv"
  )
)

write_csv(
  wei2021_acidification,
  file.path(
    table_dir,
    "wei2021_direct_acidification_values.csv"
  )
)

write_csv(
  wei2021_acidification_summary,
  file.path(
    table_dir,
    "wei2021_direct_acidification_summary.csv"
  )
)

write_csv(
  pH_final_decision,
  file.path(
    table_dir,
    "ph_final_decision_summary.csv"
  )
)


# ============================================================
# 14. SAVE OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    ph_audit =
      ph_audit,
    
    wos_ph_result_level =
      wos_ph_result_level,
    
    wos_ph_reference_level =
      wos_ph_reference_level,
    
    harm_ph =
      harm_ph,
    
    pH_metadata_coverage_overall =
      pH_metadata_coverage_overall,
    
    pH_candidate_contexts =
      pH_candidate_contexts,
    
    pH_candidate_rows =
      pH_candidate_rows,
    
    pH_salinity_confounded_contexts =
      pH_salinity_confounded_contexts,
    
    wei2021_acidification =
      wei2021_acidification,
    
    wei2021_acidification_summary =
      wei2021_acidification_summary,
    
    pH_final_decision =
      pH_final_decision,
    
    e166133_confounded_plot =
      e166133_confounded_plot
  ),
  file.path(
    object_dir,
    "ph_evidence_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_root,
    "sessionInfo.txt"
  )
)


# ============================================================
# 15. QA
# ============================================================

required_outputs <- c(
  file.path(
    table_dir,
    "ph_metadata_coverage_overall.csv"
  ),
  file.path(
    table_dir,
    "ph_within_reference_candidate_contexts.csv"
  ),
  file.path(
    table_dir,
    "ph_salinity_confounded_contexts.csv"
  ),
  file.path(
    table_dir,
    "wei2021_direct_acidification_summary.csv"
  ),
  file.path(
    table_dir,
    "ph_final_decision_summary.csv"
  ),
  file.path(
    object_dir,
    "ph_evidence_objects.rds"
  )
)

stopifnot(all(file.exists(required_outputs)))

# Reproducibility manifest
manifest <- tibble(
  role = c("input", "input", "output_object", "output_session"),
  path = c(
    combined_file,
    wos_ledger_file,
    file.path(object_dir, "ph_evidence_objects.rds"),
    file.path(output_root, "sessionInfo.txt")
  )
)
write_csv(manifest, file.path(output_root, "input_output_manifest.csv"))

writeLines(
  c(
    "STATUS: 08 ASSESS pH EVIDENCE = PASS",
    "No independent pH coefficient/model was fitted.",
    "Wei et al. (2021) retained as direct acidification context."
  ),
  con = file.path(output_root, "RUN_COMPLETE.txt")
)


# ============================================================
# 16. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("08 ASSESS pH EVIDENCE\n")
cat("============================================================\n")

cat("\n--- pH METADATA COVERAGE ---\n")
print(
  pH_metadata_coverage_overall,
  width = Inf
)

cat("\n--- WoS SECOND-PASS REFERENCE COVERAGE ---\n")
print(
  wos_ph_reference_level %>%
    count(pH_Status, name = "References"),
  n = Inf,
  width = Inf
)

cat("\n--- WITHIN-REFERENCE pH CANDIDATES ---\n")
print(
  pH_candidate_contexts,
  n = Inf,
  width = Inf
)

cat("\n--- pH-SALINITY CONFOUNDED CONTEXTS ---\n")
print(
  pH_salinity_confounded_contexts,
  n = Inf,
  width = Inf
)

cat("\n--- WEI 2021 DIRECT ACIDIFICATION ---\n")
print(
  wei2021_acidification_summary,
  width = Inf
)

cat("\n--- pH DECISION ---\n")
print(
  pH_final_decision,
  n = Inf,
  width = Inf
)

cat(
  "\nInterpretive rule:\n",
  "- WoS second-pass pH metadata are incorporated into this audit layer.\n",
  "- No quantitative WoS ADD Reference creates a manipulated within-Reference pH-LC50 contrast.\n",
  "- The three repeated-pH ECOTOX contexts remain confounded with salinity.\n",
  "- Wei et al. 2021 remains direct CO2-driven acidification evidence, not a pH slope.\n",
  "- No formal pH coefficient/model is supported.\n",
  sep = ""
)

cat(
  "\nSTATUS: 08 ASSESS pH EVIDENCE = PASS\n"
)

# ============================================================
# END
# ============================================================
