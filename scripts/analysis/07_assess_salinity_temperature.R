# ============================================================
# Copepod metal-toxicity thesis reproducibility repository
# 07_assess_salinity_temperature.R
# ============================================================
# PURPOSE
# Reproduce the final descriptive within-Reference salinity and
# temperature evidence used in the thesis.
#
# ANALYTICAL ROLE
# - No general salinity or temperature coefficient is fitted here.
# - Salinity is retained as direct within-Reference evidence with a
#   main evidence set plus two source-adjudicated extended comparisons.
# - Temperature is retained as three source-verified direct comparisons.
# - pH is assessed separately in 08_assess_ph_evidence.R.
#
# IMPORTANT
# The manually encoded salinity and publication-derived temperature
# values below are the frozen source-adjudicated values used for the
# final thesis analyses. They are not re-selected automatically from
# metadata at run time.
#
# Figures produced here are QA/intermediate figures. The polished
# thesis-facing figures are generated later by the final figure suite.
# ============================================================


# ============================================================
# 0. PACKAGES AND HELPERS
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

options(width = 220)

md5_or_na <- function(path) {
  if (file.exists(path)) {
    unname(tools::md5sum(path))
  } else {
    NA_character_
  }
}


# ============================================================
# 1. INPUTS AND OUTPUT DIRECTORIES
# ============================================================

input_file <- file.path(
  "data",
  "processed",
  "combined_ecotox_wos_harmonized.csv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Missing required input file: ", input_file,
      "\nRun 02_integrate_wos_verified_records.R first."
    )
  )
}

output_root <- file.path(
  "outputs",
  "07_assess_salinity_temperature"
)

audit_dir  <- file.path(output_root, "01_audit")
table_dir  <- file.path(output_root, "02_tables")
figure_dir <- file.path(output_root, "03_figures")
object_dir <- file.path(output_root, "04_objects")

for (d in c(output_root, audit_dir, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 2. LOAD FINAL COMBINED HARMONIZED DATA
# ============================================================

dat <- read_csv(
  input_file,
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
  "Chemical_Name"
)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Combined dataset is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

harm <- dat %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0
  )

stopifnot(nrow(harm) == 304)


# ============================================================
# 3. MAIN SALINITY DIRECT SET
# ============================================================
# Five source-adjudicated comparisons retained as the main salinity
# evidence set in the final thesis. These values are intentionally
# frozen rather than reconstructed through a broad metadata screen.

salinity_main <- tribble(
  ~Comparison_ID, ~Reference_ID, ~Reference_Number,
  ~Species, ~Metal, ~Salinity, ~LC50_umol_L,

  "SAL_MAIN_01", "E_166133", 166133,
  "Acartia tonsa", "Cu", 5, 0.657791,

  "SAL_MAIN_01", "E_166133", 166133,
  "Acartia tonsa", "Cu", 15, 1.060649,

  "SAL_MAIN_01", "E_166133", 166133,
  "Acartia tonsa", "Cu", 30, 1.710572,

  "SAL_MAIN_02", "E_107036", 107036,
  "Eurytemora affinis", "Cu", 5, 1.636610,

  "SAL_MAIN_02", "E_107036", 107036,
  "Eurytemora affinis", "Cu", 15, 1.063796,

  "SAL_MAIN_02", "E_107036", 107036,
  "Eurytemora affinis", "Cu", 25, 0.914298,

  "SAL_MAIN_03", "E_17219", 17219,
  "Eurytemora affinis", "Cd", 5, 0.459018,

  "SAL_MAIN_03", "E_17219", 17219,
  "Eurytemora affinis", "Cd", 15, 1.896561,

  "SAL_MAIN_03", "E_17219", 17219,
  "Eurytemora affinis", "Cd", 25, 0.737453,

  "SAL_MAIN_04", "E_6825", 6825,
  "Acartia lilljeborgi", "Zn", 25, 10.706638,

  "SAL_MAIN_04", "E_6825", 6825,
  "Acartia lilljeborgi", "Zn", 28, 8.412359,

  "SAL_MAIN_04", "E_6825", 6825,
  "Acartia lilljeborgi", "Zn", 32, 13.612726,

  "SAL_MAIN_05", "E_6825", 6825,
  "Temora stylifera", "Zn", 25, 0.458856,

  "SAL_MAIN_05", "E_6825", 6825,
  "Temora stylifera", "Zn", 28, 0.351790,

  "SAL_MAIN_05", "E_6825", 6825,
  "Temora stylifera", "Zn", 32, 0.474151
) %>%
  mutate(Evidence_tier = "MAIN")

stopifnot(
  n_distinct(salinity_main$Comparison_ID) == 5,
  n_distinct(salinity_main$Reference_ID) == 4,
  n_distinct(salinity_main$Species) == 4,
  all(salinity_main$LC50_umol_L > 0)
)


# ============================================================
# 4. EXTENDED SALINITY EVIDENCE: ECOTOX REFERENCE 2332
# ============================================================
# Ref. E_2332 is retained separately because temperature, pH and
# life-stage information are unavailable. For Zn, chemical form also
# differs across salinity levels: the salinity-7 ZnCl2 result is not
# combined with the salinity-15 and salinity-25 ZnSO4 results.

ref2332 <- harm %>%
  filter(
    Reference_ID == "E_2332",
    Species == "Nitocra spinipes",
    Duration_days == 4,
    Metal %in% c("Cd", "Zn")
  ) %>%
  select(
    Reference_ID,
    Reference_Number,
    Test_ID,
    Result_ID,
    Species,
    Metal,
    Chemical_Name,
    Duration_days,
    Exposure_Type,
    Conc_Type,
    Temperature,
    Salinity,
    LC50_umol_L
  ) %>%
  arrange(Metal, Salinity)

stopifnot(
  any(ref2332$Metal == "Cd" & ref2332$Salinity == 7),
  any(ref2332$Metal == "Cd" & ref2332$Salinity == 15),
  any(ref2332$Metal == "Zn" & ref2332$Salinity == 15),
  any(ref2332$Metal == "Zn" & ref2332$Salinity == 25)
)

salinity_extended <- bind_rows(
  ref2332 %>%
    filter(
      Metal == "Cd",
      Salinity %in% c(7, 15)
    ) %>%
    mutate(
      Comparison_ID = "SAL_EXT_01",
      Evidence_tier = "EXTENDED",
      Adjudication_note = paste(
        "CdCl2 at both salinity levels; retained as an extended",
        "two-level salinity comparison."
      )
    ),

  ref2332 %>%
    filter(
      Metal == "Zn",
      Salinity %in% c(15, 25),
      str_detect(
        str_to_lower(Chemical_Name),
        "sulfuric|sulf"
      )
    ) %>%
    mutate(
      Comparison_ID = "SAL_EXT_02",
      Evidence_tier = "EXTENDED",
      Adjudication_note = paste(
        "ZnSO4 at both retained salinity levels; the salinity-7",
        "ZnCl2 result is intentionally excluded from this comparison."
      )
    )
) %>%
  select(
    Comparison_ID,
    Evidence_tier,
    Reference_ID,
    Reference_Number,
    Test_ID,
    Result_ID,
    Species,
    Metal,
    Chemical_Name,
    Duration_days,
    Exposure_Type,
    Conc_Type,
    Temperature,
    Salinity,
    LC50_umol_L,
    Adjudication_note
  )

stopifnot(
  n_distinct(salinity_extended$Comparison_ID) == 2,
  n_distinct(salinity_extended$Reference_ID) == 1,
  n_distinct(salinity_extended$Species) == 1,
  all(table(salinity_extended$Comparison_ID) == 2)
)


# ============================================================
# 5. SALINITY DIRECTION SUMMARIES
# ============================================================

summarise_salinity_direction <- function(x, evidence_tier) {
  x %>%
    arrange(Comparison_ID, Salinity) %>%
    group_by(
      Comparison_ID,
      Reference_ID,
      Species,
      Metal
    ) %>%
    summarise(
      Evidence_tier = evidence_tier,
      n_levels = n_distinct(Salinity),
      min_salinity = min(Salinity),
      max_salinity = max(Salinity),
      LC50_at_min_salinity = LC50_umol_L[which.min(Salinity)],
      LC50_at_max_salinity = LC50_umol_L[which.max(Salinity)],
      ratio_highS_lowS = LC50_at_max_salinity / LC50_at_min_salinity,
      ln_ratio_highS_lowS = log(ratio_highS_lowS),
      direction = case_when(
        ratio_highS_lowS > 1 ~ "Higher salinity -> higher LC50",
        ratio_highS_lowS < 1 ~ "Higher salinity -> lower LC50",
        TRUE ~ "No change"
      ),
      .groups = "drop"
    )
}

salinity_main_direction <- summarise_salinity_direction(
  salinity_main,
  "MAIN"
)

salinity_extended_direction <- summarise_salinity_direction(
  salinity_extended,
  "EXTENDED"
)

salinity_all_direction <- bind_rows(
  salinity_main_direction,
  salinity_extended_direction
)

salinity_support_summary <- tibble(
  Evidence_set = c(
    "Main direct set",
    "Main + extended set"
  ),
  Comparisons = c(
    n_distinct(salinity_main$Comparison_ID),
    n_distinct(salinity_main$Comparison_ID) +
      n_distinct(salinity_extended$Comparison_ID)
  ),
  References = c(
    n_distinct(salinity_main$Reference_ID),
    n_distinct(c(
      salinity_main$Reference_ID,
      salinity_extended$Reference_ID
    ))
  ),
  Species = c(
    n_distinct(salinity_main$Species),
    n_distinct(c(
      salinity_main$Species,
      salinity_extended$Species
    ))
  )
)

stopifnot(
  salinity_support_summary$Comparisons[1] == 5,
  salinity_support_summary$References[1] == 4,
  salinity_support_summary$Comparisons[2] == 7,
  salinity_support_summary$References[2] == 5,
  salinity_support_summary$Species[2] == 5
)


# ============================================================
# 6. TEMPERATURE: ECOTOX E_9029 CLEAN PAIR
# ============================================================
# Ref. E_9029 contains three Cd temperature/salinity combinations.
# Only 13 C and 21 C are retained for the direct temperature contrast
# because both are at salinity 20. The 18 C result is at salinity 15.

temp_e9029 <- harm %>%
  filter(
    Reference_ID == "E_9029",
    Species == "Acartia tonsa",
    Metal == "Cd",
    Duration_days == 4,
    Salinity == 20,
    Temperature %in% c(13, 21)
  ) %>%
  transmute(
    Temperature_Comparison_ID = "TEMP_01",
    Reference_ID,
    Source_Origin = "ECOTOX",
    Reference_label = "E_9029",
    Species,
    Metal,
    Stage_context = coalesce(
      External_Stage_Sex,
      Lifestage,
      "Adult"
    ),
    Duration_h = Duration_days * 24,
    Temperature_C = Temperature,
    Salinity,
    LC50_umol_L,
    Include_clean_direct = "YES",
    Quality_note = paste(
      "Clean 13-vs-21 C pair at salinity 20; the 18 C / salinity 15",
      "result is excluded from this direct temperature contrast."
    )
  )

stopifnot(
  nrow(temp_e9029) == 2,
  setequal(temp_e9029$Temperature_C, c(13, 21))
)


# ============================================================
# 7. TEMPERATURE: LI ET AL. (2014)
# ============================================================
# Source-verified 96-h Cu LC50 values for adult Tigriopus japonicus:
# 4 C  ~20,000 ug Cu/L
# 15 C   9,300 ug Cu/L
# 25 C   2,500 ug Cu/L
# 32 C     243 ug Cu/L
# 38 C       1.1 ug Cu/L
#
# The 38 C endpoint is retained in the audit table but excluded from
# the clean direct set because control survival was only 23%.
# Values are reported on a Cu basis.

AW_Cu <- 63.546

temp_li2014_all <- tribble(
  ~Temperature_Comparison_ID,
  ~Reference_ID,
  ~Source_Origin,
  ~Reference_label,
  ~Species,
  ~Metal,
  ~Stage_context,
  ~Duration_h,
  ~Temperature_C,
  ~Salinity,
  ~LC50_reported_ug_L,
  ~Control_survival_percent,
  ~Include_clean_direct,
  ~Quality_note,

  "TEMP_02", "ENV_LI2014", "SOURCE_REVIEWED_PUBLICATION",
  "Li et al. (2014)",
  "Tigriopus japonicus", "Cu", "Adult",
  96, 4, NA_real_, 20000, 100, "YES",
  paste(
    "The publication describes the 4 C LC50 as around 20,000 ug Cu/L;",
    "surviving copepods were dormant."
  ),

  "TEMP_02", "ENV_LI2014", "SOURCE_REVIEWED_PUBLICATION",
  "Li et al. (2014)",
  "Tigriopus japonicus", "Cu", "Adult",
  96, 15, NA_real_, 9300, 99, "YES",
  NA_character_,

  "TEMP_02", "ENV_LI2014", "SOURCE_REVIEWED_PUBLICATION",
  "Li et al. (2014)",
  "Tigriopus japonicus", "Cu", "Adult",
  96, 25, NA_real_, 2500, 95, "YES",
  NA_character_,

  "TEMP_02", "ENV_LI2014", "SOURCE_REVIEWED_PUBLICATION",
  "Li et al. (2014)",
  "Tigriopus japonicus", "Cu", "Adult",
  96, 32, NA_real_, 243, 87, "YES",
  NA_character_,

  "TEMP_02", "ENV_LI2014", "SOURCE_REVIEWED_PUBLICATION",
  "Li et al. (2014)",
  "Tigriopus japonicus", "Cu", "Adult",
  96, 38, NA_real_, 1.1, 23, "NO",
  "Excluded from the clean direct set because control survival was 23%."
) %>%
  mutate(
    LC50_umol_L = LC50_reported_ug_L / AW_Cu
  )

temp_li2014_clean <- temp_li2014_all %>%
  filter(Include_clean_direct == "YES")

stopifnot(
  nrow(temp_li2014_all) == 5,
  nrow(temp_li2014_clean) == 4,
  max(temp_li2014_clean$Temperature_C) == 32
)


# ============================================================
# 8. TEMPERATURE: CHEN & DONG (2022)
# ============================================================
# The publication contains an internal numerical discrepancy for the
# 22 C Hg LC50: one passage states 1.098 mg/L, whereas the Results
# section states 1.020 (0.980-1.065) mg/L. The final direct ledger uses
# the Results-section value and records the discrepancy explicitly.
# The 25 C value is 0.521 (0.461-0.576) mg/L.

AW_Hg <- 200.59

temp_chen2022 <- tribble(
  ~Temperature_Comparison_ID,
  ~Reference_ID,
  ~Source_Origin,
  ~Reference_label,
  ~Species,
  ~Metal,
  ~Stage_context,
  ~Duration_h,
  ~Temperature_C,
  ~Salinity,
  ~LC50_reported_mg_L,
  ~Include_clean_direct,
  ~Quality_note,

  "TEMP_03", "ENV_CHEN2022", "SOURCE_REVIEWED_PUBLICATION",
  "Chen & Dong (2022)",
  "Tigriopus japonicus", "Hg", "Adult",
  24, 22, 29.5, 1.020, "YES",
  paste(
    "Results-section LC50 used. Another passage states 1.098 mg/L;",
    "the discrepancy is retained in this audit note."
  ),

  "TEMP_03", "ENV_CHEN2022", "SOURCE_REVIEWED_PUBLICATION",
  "Chen & Dong (2022)",
  "Tigriopus japonicus", "Hg", "Adult",
  24, 25, 29.5, 0.521, "YES",
  paste(
    "Two copepods died across the three 25 C control replicates;",
    "no control mortality was reported at 22 C."
  )
) %>%
  mutate(
    LC50_umol_L = LC50_reported_mg_L / AW_Hg * 1000
  )

stopifnot(nrow(temp_chen2022) == 2)


# ============================================================
# 9. FINAL TEMPERATURE DIRECT SET AND DIRECTION
# ============================================================

temperature_clean <- bind_rows(
  temp_e9029 %>%
    select(
      Temperature_Comparison_ID,
      Reference_ID,
      Source_Origin,
      Reference_label,
      Species,
      Metal,
      Stage_context,
      Duration_h,
      Temperature_C,
      Salinity,
      LC50_umol_L,
      Include_clean_direct,
      Quality_note
    ),

  temp_li2014_clean %>%
    select(
      Temperature_Comparison_ID,
      Reference_ID,
      Source_Origin,
      Reference_label,
      Species,
      Metal,
      Stage_context,
      Duration_h,
      Temperature_C,
      Salinity,
      LC50_umol_L,
      Include_clean_direct,
      Quality_note
    ),

  temp_chen2022 %>%
    select(
      Temperature_Comparison_ID,
      Reference_ID,
      Source_Origin,
      Reference_label,
      Species,
      Metal,
      Stage_context,
      Duration_h,
      Temperature_C,
      Salinity,
      LC50_umol_L,
      Include_clean_direct,
      Quality_note
    )
)

stopifnot(
  n_distinct(temperature_clean$Temperature_Comparison_ID) == 3,
  n_distinct(temperature_clean$Reference_ID) == 3,
  n_distinct(temperature_clean$Species) == 2,
  n_distinct(temperature_clean$Metal) == 3,
  all(temperature_clean$LC50_umol_L > 0)
)

temperature_direction <- temperature_clean %>%
  arrange(Temperature_Comparison_ID, Temperature_C) %>%
  group_by(
    Temperature_Comparison_ID,
    Reference_ID,
    Reference_label,
    Species,
    Metal,
    Duration_h
  ) %>%
  summarise(
    n_temperature_levels = n_distinct(Temperature_C),
    min_temperature = min(Temperature_C),
    max_temperature = max(Temperature_C),
    LC50_at_min_temperature = LC50_umol_L[which.min(Temperature_C)],
    LC50_at_max_temperature = LC50_umol_L[which.max(Temperature_C)],
    ratio_highT_lowT = LC50_at_max_temperature / LC50_at_min_temperature,
    ln_ratio_highT_lowT = log(ratio_highT_lowT),
    direction = case_when(
      ratio_highT_lowT < 1 ~ "Higher temperature -> lower LC50",
      ratio_highT_lowT > 1 ~ "Higher temperature -> higher LC50",
      TRUE ~ "No change"
    ),
    .groups = "drop"
  )

temperature_support_summary <- tibble(
  Comparisons = n_distinct(temperature_clean$Temperature_Comparison_ID),
  References = n_distinct(temperature_clean$Reference_ID),
  Species = n_distinct(temperature_clean$Species),
  Metals = n_distinct(temperature_clean$Metal),
  Lower_LC50_at_higher_temperature = sum(
    temperature_direction$ratio_highT_lowT < 1
  )
)

stopifnot(
  temperature_support_summary$Comparisons == 3,
  temperature_support_summary$References == 3,
  temperature_support_summary$Species == 2,
  temperature_support_summary$Metals == 3,
  temperature_support_summary$Lower_LC50_at_higher_temperature == 3
)


# ============================================================
# 10. SOURCE-ADJUDICATION / EXCLUSION LEDGER
# ============================================================

environmental_exclusion_ledger <- tribble(
  ~Reference_label,
  ~Factor,
  ~Final_role,
  ~Reason,

  "E_9029",
  "Temperature",
  "PARTIAL_EXCLUSION",
  paste(
    "The 18 C Cd result is excluded from the direct temperature pair",
    "because salinity is 15 rather than 20."
  ),

  "Li et al. (2014)",
  "Temperature",
  "PARTIAL_EXCLUSION",
  paste(
    "The 38 C Cu endpoint is retained in the audit but excluded from",
    "the clean direct set because control survival was 23%."
  ),

  "Holan et al. (2019)",
  "Temperature + salinity",
  "CONTEXT_ONLY",
  paste(
    "The copepod LC50 experiment lasted 7 d and is outside the thesis",
    "24-96 h direct environmental-comparison domain."
  )
)


# ============================================================
# 11. FINAL SALINITY/TEMPERATURE DECISION SUMMARY
# ============================================================

environmental_final_decision <- tribble(
  ~Component,
  ~Final_evidence_role,
  ~Support,
  ~Interpretation,

  "Salinity",
  "MAIN_DIRECT_PLUS_EXTENDED",
  "7 comparisons / 5 References / 5 Species",
  paste(
    "Directions are mixed or non-monotonic; no common salinity-toxicity",
    "direction or general coefficient is inferred."
  ),

  "Temperature",
  "DIRECT_DESCRIPTIVE",
  "3 comparisons / 3 References / 2 Species / 3 metals",
  paste(
    "All three retained comparisons show lower LC50 at the higher",
    "temperature, but support is too limited and heterogeneous for a",
    "general temperature coefficient."
  )
)


# ============================================================
# 12. QA FIGURES
# ============================================================
# These figures are retained for transparent QA and inspection. They
# are not intended to replace the polished thesis-facing figure suite.

salinity_plot_data <- bind_rows(
  salinity_main %>%
    transmute(
      Comparison_ID,
      Reference_ID,
      Species,
      Metal,
      Salinity,
      LC50_umol_L,
      Evidence_tier = "Main"
    ),
  salinity_extended %>%
    transmute(
      Comparison_ID,
      Reference_ID,
      Species,
      Metal,
      Salinity,
      LC50_umol_L,
      Evidence_tier = "Extended"
    )
) %>%
  mutate(
    Panel = paste0(Reference_ID, " | ", Metal, "\n", Species),
    Panel = factor(Panel, levels = unique(Panel))
  )

salinity_qa_plot <- ggplot(
  salinity_plot_data,
  aes(
    x = Salinity,
    y = LC50_umol_L,
    shape = Evidence_tier,
    group = Comparison_ID
  )
) +
  geom_line(alpha = 0.55) +
  geom_point(size = 2.8) +
  scale_y_log10() +
  facet_wrap(~ Panel, ncol = 3, scales = "free_x") +
  labs(
    title = "Within-Reference salinity comparisons",
    subtitle = "QA figure: main and extended source-adjudicated evidence",
    x = "Salinity",
    y = "LC50 (umol metal/L; log scale)",
    shape = "Evidence tier",
    caption = paste(
      "Extended evidence adds two E_2332 comparisons.",
      "The Zn comparison uses only comparable ZnSO4 rows at salinity 15 and 25."
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      lineheight = 1.08
    )
  )

ggsave(
  file.path(figure_dir, "salinity_within_reference_QA.png"),
  salinity_qa_plot,
  width = 10.8,
  height = 7.8,
  dpi = 300,
  bg = "white"
)

temperature_plot_data <- temperature_clean %>%
  mutate(
    Panel = paste0(
      Reference_label,
      " | ",
      Metal,
      "\n",
      Species
    )
  )

temperature_qa_plot <- ggplot(
  temperature_plot_data,
  aes(
    x = Temperature_C,
    y = LC50_umol_L,
    group = Temperature_Comparison_ID,
    shape = Source_Origin
  )
) +
  geom_line(alpha = 0.55) +
  geom_point(size = 2.8) +
  scale_y_log10() +
  facet_wrap(~ Panel, ncol = 3, scales = "free_x") +
  labs(
    title = "Within-study temperature comparisons",
    subtitle = "QA figure: three retained direct comparisons",
    x = "Temperature (C)",
    y = "LC50 (umol metal/L; log scale)",
    shape = "Source",
    caption = paste(
      "Direct descriptive comparisons only; no common temperature effect is fitted.",
      "Li et al. (2014) excludes the 38 C endpoint from the clean set because control survival was 23%."
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      lineheight = 1.08
    )
  )

ggsave(
  file.path(figure_dir, "temperature_within_reference_QA.png"),
  temperature_qa_plot,
  width = 10.5,
  height = 4.9,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 13. WRITE AUDIT TABLES AND ANALYTICAL TABLES
# ============================================================

write_csv(
  ref2332,
  file.path(
    audit_dir,
    "E2332_rows_reviewed_for_salinity_adjudication.csv"
  )
)

write_csv(
  temp_li2014_all,
  file.path(
    audit_dir,
    "temperature_Li2014_all_source_verified_points.csv"
  )
)

write_csv(
  environmental_exclusion_ledger,
  file.path(
    audit_dir,
    "salinity_temperature_exclusion_ledger.csv"
  )
)

write_csv(
  salinity_main,
  file.path(
    table_dir,
    "salinity_main_direct_values.csv"
  )
)

write_csv(
  salinity_extended,
  file.path(
    table_dir,
    "salinity_extended_E2332_values.csv"
  )
)

write_csv(
  salinity_all_direction,
  file.path(
    table_dir,
    "salinity_direction_summary_main_and_extended.csv"
  )
)

write_csv(
  salinity_support_summary,
  file.path(
    table_dir,
    "salinity_support_summary.csv"
  )
)

write_csv(
  temp_e9029,
  file.path(
    table_dir,
    "temperature_E9029_clean_pair.csv"
  )
)

write_csv(
  temp_li2014_clean,
  file.path(
    table_dir,
    "temperature_Li2014_clean_points.csv"
  )
)

write_csv(
  temp_chen2022,
  file.path(
    table_dir,
    "temperature_ChenDong2022_points.csv"
  )
)

write_csv(
  temperature_clean,
  file.path(
    table_dir,
    "temperature_final_direct_values.csv"
  )
)

write_csv(
  temperature_direction,
  file.path(
    table_dir,
    "temperature_direction_summary.csv"
  )
)

write_csv(
  temperature_support_summary,
  file.path(
    table_dir,
    "temperature_support_summary.csv"
  )
)

write_csv(
  environmental_final_decision,
  file.path(
    table_dir,
    "salinity_temperature_final_decision_summary.csv"
  )
)


# ============================================================
# 14. SAVE REUSABLE OBJECTS
# ============================================================

object_file <- file.path(
  object_dir,
  "salinity_temperature_analysis_objects.rds"
)

saveRDS(
  list(
    salinity_main = salinity_main,
    salinity_extended = salinity_extended,
    salinity_direction = salinity_all_direction,
    salinity_support_summary = salinity_support_summary,
    temperature_clean = temperature_clean,
    temperature_direction = temperature_direction,
    temperature_support_summary = temperature_support_summary,
    environmental_exclusion_ledger = environmental_exclusion_ledger,
    environmental_final_decision = environmental_final_decision,
    salinity_qa_plot = salinity_qa_plot,
    temperature_qa_plot = temperature_qa_plot
  ),
  object_file
)


# ============================================================
# 15. MANIFEST, SESSION INFO AND RUN RECORD
# ============================================================

primary_outputs <- c(
  file.path(table_dir, "salinity_support_summary.csv"),
  file.path(table_dir, "salinity_direction_summary_main_and_extended.csv"),
  file.path(table_dir, "temperature_support_summary.csv"),
  file.path(table_dir, "temperature_direction_summary.csv"),
  file.path(table_dir, "salinity_temperature_final_decision_summary.csv"),
  file.path(figure_dir, "salinity_within_reference_QA.png"),
  file.path(figure_dir, "temperature_within_reference_QA.png"),
  object_file
)

stopifnot(all(file.exists(primary_outputs)))

manifest <- tibble(
  Role = c(
    "Input: canonical combined dataset",
    "Output: salinity support summary",
    "Output: salinity direction summary",
    "Output: temperature support summary",
    "Output: temperature direction summary",
    "Output: final salinity/temperature decision summary",
    "Output: salinity QA figure",
    "Output: temperature QA figure",
    "Output: R analysis objects"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "salinity_support_summary.csv"),
    file.path(table_dir, "salinity_direction_summary_main_and_extended.csv"),
    file.path(table_dir, "temperature_support_summary.csv"),
    file.path(table_dir, "temperature_direction_summary.csv"),
    file.path(table_dir, "salinity_temperature_final_decision_summary.csv"),
    file.path(figure_dir, "salinity_within_reference_QA.png"),
    file.path(figure_dir, "temperature_within_reference_QA.png"),
    object_file
  )
) %>%
  mutate(
    Exists = file.exists(Path),
    MD5 = vapply(Path, md5_or_na, character(1))
  )

write_csv(
  manifest,
  file.path(output_root, "input_output_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)


# ============================================================
# 16. FINAL QA AND CONSOLE SUMMARY
# ============================================================

# Thesis-facing benchmark checks.
stopifnot(
  n_distinct(salinity_all_direction$Comparison_ID) == 7,
  n_distinct(salinity_all_direction$Reference_ID) == 5,
  n_distinct(c(salinity_main$Species, salinity_extended$Species)) == 5,
  nrow(temperature_direction) == 3,
  all(temperature_direction$ratio_highT_lowT < 1)
)

run_complete_file <- file.path(output_root, "RUN_COMPLETE.txt")

writeLines(
  c(
    "STATUS: 07 ASSESS SALINITY AND TEMPERATURE = PASS",
    paste0(
      "Salinity comparisons: ",
      n_distinct(salinity_all_direction$Comparison_ID),
      " across ",
      n_distinct(salinity_all_direction$Reference_ID),
      " References"
    ),
    paste0(
      "Temperature comparisons: ",
      nrow(temperature_direction),
      "; lower LC50 at higher temperature: ",
      sum(temperature_direction$ratio_highT_lowT < 1),
      "/",
      nrow(temperature_direction)
    )
  ),
  con = run_complete_file
)

cat("\n")
cat("============================================================\n")
cat("07 ASSESS SALINITY AND TEMPERATURE\n")
cat("============================================================\n")

cat("\n--- SALINITY SUPPORT ---\n")
print(salinity_support_summary, n = Inf, width = Inf)

cat("\n--- SALINITY DIRECTIONS ---\n")
print(salinity_all_direction, n = Inf, width = Inf)

cat("\n--- TEMPERATURE SUPPORT ---\n")
print(temperature_support_summary, width = Inf)

cat("\n--- TEMPERATURE DIRECTIONS ---\n")
print(temperature_direction, n = Inf, width = Inf)

cat("\n--- FINAL ANALYTICAL TREATMENT ---\n")
print(environmental_final_decision, n = Inf, width = Inf)

cat("\nSTATUS: 07 ASSESS SALINITY AND TEMPERATURE = PASS\n")
# ============================================================
# END
# ============================================================
