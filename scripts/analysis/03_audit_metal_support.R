# ============================================================
# 03_audit_metal_support.R
# ============================================================
# Purpose
# -------
# Audit evidence support and comparability before fitting the pooled
# metal model to the combined ECOTOX + source-verified WoS dataset.
#
# This script does not select metals from model coefficients or
# significance tests. It describes the evidence structure used to
# justify the thesis modelling domain, including duration support,
# cross-metal connectivity, publication/species concentration,
# biological coverage, environmental metadata and source origin.
#
# Input
# -----
# data/processed/combined_ecotox_wos_harmonized.csv
#
# Outputs
# -------
# outputs/03_audit_metal_support/
#   01_tables/   CSV audit tables
#   02_figures/  support-audit figures
#   03_objects/  reusable R objects for downstream analyses
#
# Reproducibility principle
# -------------------------
# The script checks frozen dataset benchmarks but does not alter
# eligibility, harmonisation or source-verification decisions.
# ============================================================

# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)

options(width = 220)

output_root <- file.path(
  "outputs",
  "03_audit_metal_support"
)

table_dir <- file.path(output_root, "01_tables")
figure_dir <- file.path(output_root, "02_figures")
object_dir <- file.path(output_root, "03_objects")

for (d in c(output_root, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 1. LOAD COMBINED DATA + MERGE QA BENCHMARKS
# ============================================================

input_file <- file.path(
  "data", "processed",
  "combined_ecotox_wos_harmonized.csv"
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Required combined ECOTOX + WoS dataset was not found:\n  ",
      input_file,
      "\n\nRun scripts/analysis/02_integrate_wos_verified_records.R first."
    )
  )
}

dat <- read_csv(input_file, show_col_types = FALSE)

required_columns <- c(
  "Metal",
  "Species",
  "Order",
  "Lifestage",
  "Duration_days",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Temperature",
  "Salinity",
  "pH",
  "Eligibility_Status",
  "Harmonization_Status",
  "LC50_umol_L",
  "ln_LC50_umol_L",
  "Source_Origin",
  "Reference_ID",
  "Test_ID",
  "Result_ID"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Combined dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Frozen merge benchmarks from the verified ECOTOX + WoS integration step.
stopifnot(
  nrow(dat) == 353,
  n_distinct(dat$Result_ID) == 353,
  n_distinct(dat$Reference_ID) == 84,
  sum(dat$Source_Origin == "ECOTOX") == 296,
  sum(dat$Source_Origin == "WoS_supplemental") == 57,
  sum(dat$Harmonization_Status == "HARMONIZED") == 304
)

preferred_metal_order <- c(
  "Cu", "Cd", "Zn", "Ag", "Cr", "Hg", "Ni", "Pb"
)

duration_levels <- 1:4


# ============================================================
# 2. DEFINE COMBINED CANDIDATE POOL + COMMON MOLAR POOL
# ============================================================

candidate_pool <- dat %>%
  filter(
    Eligibility_Status %in% c(
      "PRIMARY_KEEP",
      "PRIMARY_RECOVERED",
      "WOS_SOURCE_VERIFIED_KEEP"
    )
  )

common_pool <- candidate_pool %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0
  )

stopifnot(
  nrow(candidate_pool) == 353,
  nrow(common_pool) == 304,
  sum(candidate_pool$Harmonization_Status == "PENDING_CR_SOURCE_BASIS") == 49
)

candidate_metals <- preferred_metal_order[
  preferred_metal_order %in% unique(candidate_pool$Metal)
]

harmonized_metals <- preferred_metal_order[
  preferred_metal_order %in% unique(common_pool$Metal)
]

# Final expected harmonized support after source verification + merge.
expected_common_by_metal <- tribble(
  ~Metal, ~Expected_Results,
  "Cu", 113L,
  "Cd",  57L,
  "Zn",  49L,
  "Ag",  19L,
  "Hg",  37L,
  "Ni",  25L,
  "Pb",   4L
)

common_by_metal_qa <- common_pool %>%
  count(Metal, name = "Results") %>%
  full_join(expected_common_by_metal, by = "Metal") %>%
  mutate(Match = Results == Expected_Results) %>%
  arrange(factor(Metal, levels = harmonized_metals))

write_csv(
  common_by_metal_qa,
  file.path(table_dir, "00_common_pool_by_metal_QA.csv")
)

stopifnot(all(common_by_metal_qa$Match))

# Log integrity.
log_integrity <- common_pool %>%
  summarise(
    missing_log = sum(is.na(ln_LC50_umol_L)),
    max_abs_difference = max(
      abs(ln_LC50_umol_L - log(LC50_umol_L)),
      na.rm = TRUE
    )
  )

stopifnot(
  log_integrity$missing_log == 0,
  log_integrity$max_abs_difference < 1e-10
)


# ============================================================
# 3. OVERALL COMBINED SUPPORT BY METAL
# ============================================================

metal_harmonization_status <- candidate_pool %>%
  group_by(Metal) %>%
  summarise(
    Eligible_Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    Harmonized_Results = sum(Harmonization_Status == "HARMONIZED"),
    Pending_Basis_Results = sum(Harmonization_Status == "PENDING_CR_SOURCE_BASIS"),
    Common_Molar_Ready = all(Harmonization_Status == "HARMONIZED"),
    Harmonization_Statuses = paste(sort(unique(Harmonization_Status)), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    Metal = factor(Metal, levels = candidate_metals),
    Chemistry_Interpretation = case_when(
      Common_Molar_Ready ~ "Common target-metal molar LC50 available",
      Metal == "Cr" ~ "Common molar LC50 not assigned: source concentration basis unresolved",
      TRUE ~ "Review harmonization status"
    )
  ) %>%
  arrange(Metal) %>%
  mutate(Metal = as.character(Metal))

metal_support_common_pool <- common_pool %>%
  group_by(Metal) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  mutate(
    WoS_Share_pct = round(100 * WoS_Results / Results, 1),
    Metal = factor(Metal, levels = harmonized_metals)
  ) %>%
  arrange(Metal) %>%
  mutate(Metal = as.character(Metal))

write_csv(
  metal_harmonization_status,
  file.path(table_dir, "01_metal_harmonization_status.csv")
)

write_csv(
  metal_support_common_pool,
  file.path(table_dir, "02_common_pool_support_by_metal.csv")
)


# ============================================================
# 4. DURATION SUPPORT AUDIT
# ============================================================

# 96 h is not assumed here; all four acute durations are re-counted.
duration_overall <- tibble(Duration_days = duration_levels) %>%
  left_join(
    common_pool %>%
      group_by(Duration_days) %>%
      summarise(
        Results = n(),
        Tests = n_distinct(Test_ID),
        References = n_distinct(Reference_ID),
        Species = n_distinct(Species),
        Metals = n_distinct(Metal),
        ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
        WoS_Results = sum(Source_Origin == "WoS_supplemental"),
        .groups = "drop"
      ),
    by = "Duration_days"
  ) %>%
  mutate(
    across(
      c(Results, Tests, References, Species, Metals, ECOTOX_Results, WoS_Results),
      ~ replace_na(.x, 0L)
    ),
    Duration_h = 24 * Duration_days
  ) %>%
  select(Duration_days, Duration_h, everything())

# Frozen combined duration totals: original ECOTOX harmonized + retained WoS.
stopifnot(
  duration_overall %>% filter(Duration_days == 1) %>% pull(Results) == 27,
  duration_overall %>% filter(Duration_days == 2) %>% pull(Results) == 70,
  duration_overall %>% filter(Duration_days == 3) %>% pull(Results) == 30,
  duration_overall %>% filter(Duration_days == 4) %>% pull(Results) == 177
)

metal_duration_support <- expand_grid(
  Metal = harmonized_metals,
  Duration_days = duration_levels
) %>%
  left_join(
    common_pool %>%
      group_by(Metal, Duration_days) %>%
      summarise(
        Results = n(),
        Tests = n_distinct(Test_ID),
        References = n_distinct(Reference_ID),
        Species = n_distinct(Species),
        ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
        WoS_Results = sum(Source_Origin == "WoS_supplemental"),
        .groups = "drop"
      ),
    by = c("Metal", "Duration_days")
  ) %>%
  mutate(
    across(
      c(Results, Tests, References, Species, ECOTOX_Results, WoS_Results),
      ~ replace_na(.x, 0L)
    ),
    Duration_h = 24 * Duration_days,
    WoS_Share_pct = if_else(
      Results > 0,
      round(100 * WoS_Results / Results, 1),
      NA_real_
    )
  ) %>%
  arrange(
    Duration_days,
    factor(Metal, levels = harmonized_metals)
  )

candidate_metal_duration_support <- expand_grid(
  Metal = candidate_metals,
  Duration_days = duration_levels
) %>%
  left_join(
    candidate_pool %>%
      group_by(Metal, Duration_days) %>%
      summarise(
        Results = n(),
        Tests = n_distinct(Test_ID),
        References = n_distinct(Reference_ID),
        Species = n_distinct(Species),
        .groups = "drop"
      ),
    by = c("Metal", "Duration_days")
  ) %>%
  mutate(
    across(c(Results, Tests, References, Species), ~ replace_na(.x, 0L)),
    Duration_h = 24 * Duration_days
  ) %>%
  left_join(
    metal_harmonization_status %>%
      select(Metal, Common_Molar_Ready, Chemistry_Interpretation),
    by = "Metal"
  ) %>%
  arrange(
    Duration_days,
    factor(Metal, levels = candidate_metals)
  )

write_csv(
  duration_overall,
  file.path(table_dir, "03_duration_overall_support.csv")
)

write_csv(
  metal_duration_support,
  file.path(table_dir, "04_metal_x_duration_support_harmonized.csv")
)

write_csv(
  candidate_metal_duration_support,
  file.path(table_dir, "05_metal_x_duration_support_all_candidates.csv")
)


# ============================================================
# 5. GENERIC PAIRWISE CONNECTIVITY FUNCTION
# ============================================================

make_pair_connectivity <- function(data, metal_levels, duration_value = NULL) {
  
  if (!is.null(duration_value)) {
    data <- data %>% filter(Duration_days == duration_value)
  }
  
  metal_pairs <- combn(metal_levels, 2, simplify = FALSE)
  
  map_dfr(
    metal_pairs,
    function(pair) {
      A <- pair[1]
      B <- pair[2]
      
      dat_A <- data %>% filter(Metal == A)
      dat_B <- data %>% filter(Metal == B)
      
      refs_A <- unique(dat_A$Reference_ID[!is.na(dat_A$Reference_ID)])
      refs_B <- unique(dat_B$Reference_ID[!is.na(dat_B$Reference_ID)])
      
      spp_A <- unique(dat_A$Species[!is.na(dat_A$Species) & dat_A$Species != ""])
      spp_B <- unique(dat_B$Species[!is.na(dat_B$Species) & dat_B$Species != ""])
      
      cells_A <- dat_A %>%
        filter(!is.na(Reference_ID), !is.na(Species), Species != "") %>%
        distinct(Reference_ID, Species)
      
      cells_B <- dat_B %>%
        filter(!is.na(Reference_ID), !is.na(Species), Species != "") %>%
        distinct(Reference_ID, Species)
      
      shared_cells <- inner_join(
        cells_A,
        cells_B,
        by = c("Reference_ID", "Species")
      )
      
      tibble(
        Metal_A = A,
        Metal_B = B,
        A_Results = nrow(dat_A),
        A_References = n_distinct(dat_A$Reference_ID),
        A_Species = n_distinct(dat_A$Species),
        B_Results = nrow(dat_B),
        B_References = n_distinct(dat_B$Reference_ID),
        B_Species = n_distinct(dat_B$Species),
        Shared_References = length(intersect(refs_A, refs_B)),
        Shared_Species = length(intersect(spp_A, spp_B)),
        Shared_RefSpecies = nrow(shared_cells),
        Shared_Refs_within_RefSpecies = n_distinct(shared_cells$Reference_ID),
        Shared_Species_within_RefSpecies = n_distinct(shared_cells$Species)
      )
    }
  )
}

pair_connectivity_all <- make_pair_connectivity(
  common_pool,
  harmonized_metals,
  duration_value = NULL
)

pair_connectivity_96 <- make_pair_connectivity(
  common_pool,
  harmonized_metals,
  duration_value = 4
)

write_csv(
  pair_connectivity_all,
  file.path(table_dir, "06_pair_connectivity_all_durations.csv")
)

write_csv(
  pair_connectivity_96,
  file.path(table_dir, "07_pair_connectivity_96h.csv")
)


# ============================================================
# 6. 96-h DOMAIN SUPPORT
# ============================================================

pool_96 <- common_pool %>% filter(Duration_days == 4)
candidate_pool_96 <- candidate_pool %>% filter(Duration_days == 4)

stopifnot(
  nrow(pool_96) == 177,
  nrow(candidate_pool_96) == 200
)

# Frozen combined 96-h harmonized Result counts.
expected_96_by_metal <- tribble(
  ~Metal, ~Expected_Results_96,
  "Cu", 67L,
  "Cd", 37L,
  "Zn", 27L,
  "Ag", 18L,
  "Hg", 12L,
  "Ni", 12L,
  "Pb",  4L
)

metal_support_96 <- pool_96 %>%
  group_by(Metal) %>%
  summarise(
    Results_96 = n(),
    Tests_96 = n_distinct(Test_ID),
    References_96 = n_distinct(Reference_ID),
    Species_96 = n_distinct(Species),
    ECOTOX_Results_96 = sum(Source_Origin == "ECOTOX"),
    WoS_Results_96 = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  left_join(expected_96_by_metal, by = "Metal") %>%
  mutate(
    Count_QA = Results_96 == Expected_Results_96,
    WoS_Share_96_pct = round(100 * WoS_Results_96 / Results_96, 1)
  ) %>%
  arrange(factor(Metal, levels = harmonized_metals))

stopifnot(all(metal_support_96$Count_QA))

candidate_support_96 <- candidate_pool_96 %>%
  group_by(Metal) %>%
  summarise(
    Results_96 = n(),
    Tests_96 = n_distinct(Test_ID),
    References_96 = n_distinct(Reference_ID),
    Species_96 = n_distinct(Species),
    .groups = "drop"
  ) %>%
  left_join(
    metal_harmonization_status %>%
      select(Metal, Common_Molar_Ready, Chemistry_Interpretation),
    by = "Metal"
  ) %>%
  arrange(factor(Metal, levels = candidate_metals))

write_csv(
  metal_support_96,
  file.path(table_dir, "08_metal_support_96h_harmonized.csv")
)

write_csv(
  candidate_support_96,
  file.path(table_dir, "09_metal_support_96h_all_candidates.csv")
)


# ============================================================
# 7. SOURCE / SPECIES DOMINANCE AT 96 h
# ============================================================

ref_contribution_96 <- pool_96 %>%
  count(Metal, Reference_ID, name = "Results_from_reference") %>%
  group_by(Metal) %>%
  mutate(
    Metal_results = sum(Results_from_reference),
    Reference_share = Results_from_reference / Metal_results
  ) %>%
  ungroup()

ref_dominance_96 <- ref_contribution_96 %>%
  group_by(Metal) %>%
  summarise(
    Results = first(Metal_results),
    References = n_distinct(Reference_ID),
    Top_reference_results = max(Results_from_reference),
    Top_reference_share_pct = round(100 * max(Reference_share), 1),
    Top_2_reference_share_pct = round(
      100 * sum(sort(Reference_share, decreasing = TRUE)[seq_len(min(2, n()))]),
      1
    ),
    Single_result_references = sum(Results_from_reference == 1),
    .groups = "drop"
  )

species_contribution_96 <- pool_96 %>%
  count(Metal, Species, name = "Results_from_species") %>%
  group_by(Metal) %>%
  mutate(
    Metal_results = sum(Results_from_species),
    Species_share = Results_from_species / Metal_results
  ) %>%
  ungroup()

species_dominance_96 <- species_contribution_96 %>%
  group_by(Metal) %>%
  summarise(
    Results = first(Metal_results),
    Species = n_distinct(Species),
    Top_species_results = max(Results_from_species),
    Top_species_share_pct = round(100 * max(Species_share), 1),
    Top_2_species_share_pct = round(
      100 * sum(sort(Species_share, decreasing = TRUE)[seq_len(min(2, n()))]),
      1
    ),
    Single_result_species = sum(Results_from_species == 1),
    .groups = "drop"
  )

write_csv(
  ref_dominance_96,
  file.path(table_dir, "10_reference_dominance_96h.csv")
)

write_csv(
  species_dominance_96,
  file.path(table_dir, "11_species_dominance_96h.csv")
)


# ============================================================
# 8. STAGE SUPPORT SNAPSHOT AT 96 h
# ============================================================

stage_support_96 <- pool_96 %>%
  filter(!is.na(Lifestage), Lifestage != "") %>%
  group_by(Metal, Lifestage) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(
    factor(Metal, levels = harmonized_metals),
    Lifestage
  )

stage_missing_96 <- pool_96 %>%
  group_by(Metal) %>%
  summarise(
    Results = n(),
    Stage_known = sum(!is.na(Lifestage) & Lifestage != ""),
    Stage_missing = sum(is.na(Lifestage) | Lifestage == ""),
    Stage_known_pct = round(100 * Stage_known / Results, 1),
    .groups = "drop"
  )

write_csv(
  stage_support_96,
  file.path(table_dir, "12_stage_support_96h.csv")
)

write_csv(
  stage_missing_96,
  file.path(table_dir, "13_stage_metadata_completeness_96h.csv")
)


# ============================================================
# 9. ORDER SUPPORT SNAPSHOT AT 96 h
# ============================================================

order_support_96 <- pool_96 %>%
  filter(!is.na(Order), Order != "") %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(
    factor(Metal, levels = harmonized_metals),
    Order
  )

write_csv(
  order_support_96,
  file.path(table_dir, "14_order_support_96h.csv")
)


# ============================================================
# 10. PROTOCOL + ENVIRONMENTAL SUPPORT AT 96 h
# ============================================================

exposure_support_96 <- pool_96 %>%
  group_by(Metal, Exposure_Type) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

conc_type_support_96 <- pool_96 %>%
  group_by(Metal, Conc_Type) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

basis_support_96 <- pool_96 %>%
  group_by(Metal, Basis_Interpretation) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

safe_min <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

env_support_96 <- pool_96 %>%
  group_by(Metal) %>%
  summarise(
    Results = n(),
    Temp_known = sum(!is.na(Temperature)),
    Temp_known_pct = round(100 * Temp_known / Results, 1),
    Temp_min = safe_min(Temperature),
    Temp_max = safe_max(Temperature),
    Sal_known = sum(!is.na(Salinity)),
    Sal_known_pct = round(100 * Sal_known / Results, 1),
    Sal_min = safe_min(Salinity),
    Sal_max = safe_max(Salinity),
    Both_known = sum(!is.na(Temperature) & !is.na(Salinity)),
    Both_known_pct = round(100 * Both_known / Results, 1),
    pH_known = sum(!is.na(pH)),
    pH_known_pct = round(100 * pH_known / Results, 1),
    .groups = "drop"
  )

write_csv(
  exposure_support_96,
  file.path(table_dir, "15_exposure_type_support_96h.csv")
)

write_csv(
  conc_type_support_96,
  file.path(table_dir, "16_concentration_type_support_96h.csv")
)

write_csv(
  basis_support_96,
  file.path(table_dir, "17_basis_support_96h.csv")
)

write_csv(
  env_support_96,
  file.path(table_dir, "18_environmental_metadata_support_96h.csv")
)


# ============================================================
# 11. 96-h METAL-SELECTION EVIDENCE MATRIX
# ============================================================

ref_connectivity_summary_96 <- pair_connectivity_96 %>%
  select(Metal_A, Metal_B, Shared_References) %>%
  pivot_longer(c(Metal_A, Metal_B), names_to = "side", values_to = "Metal") %>%
  group_by(Metal) %>%
  summarise(
    Connected_metals_by_reference = sum(Shared_References > 0),
    Total_shared_references_across_pairs = sum(Shared_References),
    .groups = "drop"
  )

species_connectivity_summary_96 <- pair_connectivity_96 %>%
  select(Metal_A, Metal_B, Shared_Species) %>%
  pivot_longer(c(Metal_A, Metal_B), names_to = "side", values_to = "Metal") %>%
  group_by(Metal) %>%
  summarise(
    Connected_metals_by_species = sum(Shared_Species > 0),
    Total_shared_species_across_pairs = sum(Shared_Species),
    .groups = "drop"
  )

joint_connectivity_summary_96 <- pair_connectivity_96 %>%
  select(Metal_A, Metal_B, Shared_RefSpecies) %>%
  pivot_longer(c(Metal_A, Metal_B), names_to = "side", values_to = "Metal") %>%
  group_by(Metal) %>%
  summarise(
    Connected_metals_by_refspecies = sum(Shared_RefSpecies > 0),
    Total_shared_refspecies_across_pairs = sum(Shared_RefSpecies),
    .groups = "drop"
  )

metal_selection_matrix_96 <- candidate_support_96 %>%
  left_join(metal_support_96, by = c("Metal", "Results_96", "Tests_96", "References_96", "Species_96")) %>%
  left_join(ref_dominance_96, by = "Metal", suffix = c("", "_refdom")) %>%
  left_join(species_dominance_96, by = "Metal", suffix = c("", "_spdom")) %>%
  left_join(stage_missing_96, by = "Metal", suffix = c("", "_stage")) %>%
  left_join(env_support_96, by = "Metal", suffix = c("", "_env")) %>%
  left_join(ref_connectivity_summary_96, by = "Metal") %>%
  left_join(species_connectivity_summary_96, by = "Metal") %>%
  left_join(joint_connectivity_summary_96, by = "Metal") %>%
  mutate(
    Quantitative_96h_Status = case_when(
      Common_Molar_Ready ~ "Quantitatively harmonized at 96 h",
      Metal == "Cr" ~ "Not eligible for common-molar modelling",
      TRUE ~ "Review"
    ),
    Previous_Pooled_Panel = Metal %in% c("Cu", "Cd", "Zn")
  ) %>%
  arrange(factor(Metal, levels = candidate_metals))

write_csv(
  metal_selection_matrix_96,
  file.path(table_dir, "19_metal_selection_matrix_96h.csv")
)


# ============================================================
# 12. STRUCTURAL MULTI-METAL CONTEXTS FOR MATCHED FOLLOW-UP
# ============================================================
#
# IMPORTANT: these are candidate contexts only. They are NOT automatically
# accepted as source-verified matched comparisons. The within-reference metal script will still apply
# the tighter provenance / biological / protocol checks.

multi_metal_ref_species_duration <- common_pool %>%
  filter(
    !is.na(Reference_ID),
    !is.na(Species), Species != "",
    !is.na(Duration_days)
  ) %>%
  group_by(Reference_ID, Species, Duration_days) %>%
  summarise(
    n_metals = n_distinct(Metal),
    Metals = paste(sort(unique(Metal)), collapse = " | "),
    Results = n(),
    Tests = n_distinct(Test_ID),
    Source_Origins = paste(sort(unique(Source_Origin)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_metals >= 2) %>%
  arrange(desc(n_metals), Reference_ID, Species, Duration_days)

multi_metal_reference_duration <- common_pool %>%
  filter(!is.na(Reference_ID), !is.na(Duration_days)) %>%
  group_by(Reference_ID, Duration_days) %>%
  summarise(
    n_metals = n_distinct(Metal),
    Metals = paste(sort(unique(Metal)), collapse = " | "),
    Results = n(),
    Species = n_distinct(Species),
    Tests = n_distinct(Test_ID),
    Source_Origins = paste(sort(unique(Source_Origin)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_metals >= 2) %>%
  arrange(desc(n_metals), Reference_ID, Duration_days)

write_csv(
  multi_metal_ref_species_duration,
  file.path(table_dir, "20_candidate_matched_ref_species_duration_contexts.csv")
)

write_csv(
  multi_metal_reference_duration,
  file.path(table_dir, "21_candidate_multi_metal_reference_duration_contexts.csv")
)


# ============================================================
# 13. SOURCE-ORIGIN CONTRIBUTION TABLES
# ============================================================

source_origin_by_metal <- common_pool %>%
  count(Metal, Source_Origin, name = "Results") %>%
  group_by(Metal) %>%
  mutate(
    Metal_total = sum(Results),
    Share_pct = round(100 * Results / Metal_total, 1)
  ) %>%
  ungroup()

source_origin_by_metal_duration <- common_pool %>%
  count(Metal, Duration_days, Source_Origin, name = "Results") %>%
  group_by(Metal, Duration_days) %>%
  mutate(
    Cell_total = sum(Results),
    Share_pct = round(100 * Results / Cell_total, 1)
  ) %>%
  ungroup()

write_csv(
  source_origin_by_metal,
  file.path(table_dir, "22_source_origin_by_metal.csv")
)

write_csv(
  source_origin_by_metal_duration,
  file.path(table_dir, "23_source_origin_by_metal_duration.csv")
)


# ============================================================
# 14. COMPACT DECISION SNAPSHOT
# ============================================================

previous_panel_96 <- pool_96 %>%
  filter(Metal %in% c("Cu", "Cd", "Zn")) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species)
  )

support_snapshot <- metal_support_96 %>%
  select(
    Metal,
    Results_96,
    Tests_96,
    References_96,
    Species_96,
    ECOTOX_Results_96,
    WoS_Results_96,
    WoS_Share_96_pct
  ) %>%
  left_join(
    ref_dominance_96 %>%
      select(Metal, Top_reference_share_pct, Top_2_reference_share_pct),
    by = "Metal"
  ) %>%
  left_join(
    species_dominance_96 %>%
      select(Metal, Top_species_share_pct, Top_2_species_share_pct),
    by = "Metal"
  ) %>%
  left_join(ref_connectivity_summary_96, by = "Metal") %>%
  left_join(joint_connectivity_summary_96, by = "Metal") %>%
  mutate(
    Previously_selected = Metal %in% c("Cu", "Cd", "Zn")
  ) %>%
  arrange(desc(References_96), desc(Species_96), desc(Results_96))

write_csv(
  support_snapshot,
  file.path(table_dir, "24_DECISION_SNAPSHOT_96h.csv")
)

write_csv(
  previous_panel_96,
  file.path(table_dir, "25_previous_CuCdZn_panel_current_size.csv")
)


# ============================================================
# 15. AUDIT FIGURES
# ============================================================

# Figure A: Reference breadth across duration by metal.
fig_reference_duration <- ggplot(
  metal_duration_support,
  aes(x = Duration_h, y = References)
) +
  geom_point(size = 2.6) +
  facet_wrap(~ Metal) +
  scale_x_continuous(breaks = c(24, 48, 72, 96)) +
  labs(
    x = "Exposure duration (h)",
    y = "Distinct References",
    title = "Combined ECOTOX + WoS Reference support across acute durations",
    subtitle = "Points are support counts; they are not longitudinal trajectories."
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figure_dir, "01_reference_support_across_duration.png"),
  fig_reference_duration,
  width = 10,
  height = 6,
  dpi = 300
)

# Figure B: 96-h Results / References / Species side by side as separate facets.
fig_96_data <- metal_support_96 %>%
  select(Metal, Results_96, References_96, Species_96) %>%
  pivot_longer(
    cols = c(Results_96, References_96, Species_96),
    names_to = "Metric",
    values_to = "Count"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Results_96 = "Results",
      References_96 = "References",
      Species_96 = "Species"
    ),
    Metal = factor(Metal, levels = harmonized_metals)
  )

fig_96_support <- ggplot(
  fig_96_data,
  aes(x = Metal, y = Count)
) +
  geom_col() +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    x = NULL,
    y = "Count",
    title = "Combined 96-h evidence support by metal"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(figure_dir, "02_combined_96h_support_by_metal.png"),
  fig_96_support,
  width = 9,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 16. OPTIONAL EXCEL AUDIT WORKBOOK
# ============================================================

excel_file <- file.path(
  output_root,
  "Combined_support_audit.xlsx"
)

excel_sheets <- list(
  Metal_Support = metal_support_common_pool,
  Duration = duration_overall,
  Metal_x_Duration = metal_duration_support,
  Metal_96h = metal_support_96,
  Selection_Matrix_96h = metal_selection_matrix_96,
  Pair_Connectivity_96h = pair_connectivity_96,
  Reference_Dominance_96h = ref_dominance_96,
  Species_Dominance_96h = species_dominance_96,
  Stage_96h = stage_support_96,
  Order_96h = order_support_96,
  Environment_96h = env_support_96,
  Candidate_Matched = multi_metal_ref_species_duration,
  Source_Origin = source_origin_by_metal,
  Decision_Snapshot = support_snapshot
)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(excel_sheets, path = excel_file)
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(excel_sheets, file = excel_file, overwrite = TRUE)
} else {
  warning(
    "Neither writexl nor openxlsx is installed. CSV outputs were written, but the convenience Excel workbook was not created."
  )
}


# ============================================================
# 17. SAVE OBJECTS
# ============================================================

saveRDS(
  list(
    candidate_pool = candidate_pool,
    common_pool = common_pool,
    pool_96 = pool_96,
    metal_support_common_pool = metal_support_common_pool,
    duration_overall = duration_overall,
    metal_duration_support = metal_duration_support,
    pair_connectivity_all = pair_connectivity_all,
    pair_connectivity_96 = pair_connectivity_96,
    metal_support_96 = metal_support_96,
    metal_selection_matrix_96 = metal_selection_matrix_96,
    stage_support_96 = stage_support_96,
    order_support_96 = order_support_96,
    env_support_96 = env_support_96,
    candidate_matched_contexts = multi_metal_ref_species_duration
  ),
  file.path(object_dir, "metal_support_audit_objects.rds")
)


# ============================================================
# 18. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n============================================\n")
cat("METAL SUPPORT / COMPARABILITY AUDIT COMPLETE\n")
cat("============================================\n\n")

cat("Combined candidate rows:", nrow(candidate_pool), "\n")
cat("Combined harmonized rows:", nrow(common_pool), "\n")
cat("Combined References:", n_distinct(candidate_pool$Reference_ID), "\n")
cat("Harmonized References:", n_distinct(common_pool$Reference_ID), "\n\n")

cat("Harmonized support by metal:\n")
print(metal_support_common_pool)

cat("\nSupport by duration:\n")
print(duration_overall)

cat("\n96-h support by metal:\n")
print(metal_support_96)

cat("\n96-h pairwise connectivity:\n")
print(
  pair_connectivity_96 %>%
    select(
      Metal_A,
      Metal_B,
      Shared_References,
      Shared_Species,
      Shared_RefSpecies
    ) %>%
    arrange(desc(Shared_References), desc(Shared_RefSpecies))
)

cat("\nPreviously used Cu-Cd-Zn panel on the combined data:\n")
print(previous_panel_96)

cat("\nCandidate multi-metal Reference x Species x Duration contexts:",
    nrow(multi_metal_ref_species_duration), "\n")

cat("\nIMPORTANT: no pooled metal panel has been locked by this script.\n")
cat("Review 24_DECISION_SNAPSHOT_96h.csv and 19_metal_selection_matrix_96h.csv as the audit trail for the pooled metal domain.\n")
cat("For matched work, also inspect 20_candidate_matched_ref_species_duration_contexts.csv.\n")


# ============================================================
# 19. REPRODUCIBILITY RECORDS + FINAL STATUS
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  file.path(output_root, "sessionInfo.txt")
)

input_output_manifest <- tibble(
  Role = c(
    "Input",
    "Primary output",
    "Primary output",
    "Primary output",
    "Reusable R object"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "19_metal_selection_matrix_96h.csv"),
    file.path(table_dir, "24_DECISION_SNAPSHOT_96h.csv"),
    file.path(table_dir, "20_candidate_matched_ref_species_duration_contexts.csv"),
    file.path(object_dir, "metal_support_audit_objects.rds")
  )
)

write_csv(
  input_output_manifest,
  file.path(output_root, "input_output_manifest.csv")
)

stopifnot(
  nrow(candidate_pool) == 353,
  nrow(common_pool) == 304,
  nrow(pool_96) == 177,
  all(common_by_metal_qa$Match),
  all(metal_support_96$Count_QA)
)

cat("\nSTATUS: 03 AUDIT METAL SUPPORT = PASS\n")
