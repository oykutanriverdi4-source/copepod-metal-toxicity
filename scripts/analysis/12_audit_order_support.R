# ============================================================
# 12_audit_order_support.R
# ============================================================
# Purpose
# Audit the structural support for the final 96-h taxonomic-order
# comparison before fitting the Order model.
#
# Scientific question supported
# Does the Calanoida-versus-Harpacticoida LC50 contrast differ
# between Cu and Cd?
#
# Final model domain
# - exposure duration: 96 h
# - metals: Cu and Cd
# - Orders: Calanoida and Harpacticoida
#
# Input
# data/processed/combined_ecotox_wos_harmonized.csv
#
# Main outputs
# outputs/12_audit_order_support/
#
# Scientific workflow notes
# - This is a pre-model support audit; model outcomes are not used to
#   choose the comparison domain.
# - All harmonized metals are screened at 96 h before the final Cu-Cd
#   domain is documented.
# - Species naturally belong to taxonomic Orders; that biological
#   nesting is not treated as a data defect.
# - Structural support is assessed from Result, Reference and Species
#   breadth, cross-metal connectivity, dominance, developmental-stage
#   composition, and temperature/salinity coverage.
# - The scientific support calculations and frozen numerical QA gates
#   are preserved from the final thesis workflow.
# - Figures generated here are support/QA figures rather than the
#   polished thesis-facing figure suite.
#
# Reproducibility rule
# Run from the repository root with the RStudio project open.
# Script 02 must have produced the combined harmonized dataset.
# ============================================================


# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)

options(width = 220)

output_root <- file.path(
  "outputs",
  "12_audit_order_support"
)

table_dir <- file.path(output_root, "01_tables")
figure_dir <- file.path(output_root, "02_figures")
object_dir <- file.path(output_root, "03_objects")

completion_path <- file.path(output_root, "RUN_COMPLETE.txt")
if (file.exists(completion_path) && !file.remove(completion_path)) {
  stop("Could not remove the previous completion marker. Close open output files and retry.")
}

for (d in c(output_root, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 1. IMPORT + MASTER QA
# ============================================================

input_file <- file.path(
  "data",
  "processed",
  "combined_ecotox_wos_harmonized.csv"
)

stopifnot(file.exists(input_file))

dat <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "Reference_Number",
  "Metal",
  "Species",
  "Order",
  "Lifestage",
  "Duration_days",
  "Exposure_Type",
  "Conc_Type",
  "Temperature",
  "Salinity",
  "Harmonization_Status",
  "LC50_umol_L",
  "Source_Origin"
)

stopifnot(all(required_columns %in% names(dat)))

stopifnot(
  nrow(dat) == 353,
  n_distinct(dat$Reference_ID) == 84,
  sum(dat$Source_Origin == "ECOTOX") == 296,
  sum(dat$Source_Origin == "WoS_supplemental") == 57,
  sum(dat$Harmonization_Status == "HARMONIZED", na.rm = TRUE) == 304
)


# ============================================================
# 2. COMMON HARMONIZED QUANTITATIVE POOL
# ============================================================

architecture_data <- dat %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0,
    !is.na(Order),
    !is.na(Duration_days)
  ) %>%
  mutate(
    ln_LC50 = log(LC50_umol_L),
    Lifestage = na_if(trimws(as.character(Lifestage)), "")
  )

architecture_overall <- architecture_data %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    Metals = n_distinct(Metal),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental")
  )

stopifnot(architecture_overall$Results == 304)


# ============================================================
# 3. OVERALL ORDER / SPECIES ARCHITECTURE
# ============================================================

order_architecture <- architecture_data %>%
  group_by(Order) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    Metals = n_distinct(Metal),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(desc(References), desc(Species))

order_metal_support <- architecture_data %>%
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
  arrange(Metal, Order)

species_architecture <- architecture_data %>%
  group_by(Order, Species) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Metals = n_distinct(Metal),
    metal_set = paste(sort(unique(Metal)), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(Order, desc(References), desc(Results))


# ============================================================
# 4. ALL-METAL ORDER x DURATION SUPPORT LANDSCAPE
# ============================================================

order_two_main <- architecture_data %>%
  filter(Order %in% c("Calanoida", "Harpacticoida"))

order_duration_support <- order_two_main %>%
  group_by(Metal, Order, Duration_days) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Duration_days, Metal, Order)

order_duration_range <- order_two_main %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    min_duration = min(Duration_days),
    max_duration = max(Duration_days),
    durations = paste(sort(unique(Duration_days)), collapse = " | "),
    .groups = "drop"
  )

order_duration_reference_plot_data <- order_duration_support %>%
  select(Metal, Order, Duration_days, References) %>%
  tidyr::complete(
    Metal,
    Order,
    Duration_days = 1:4,
    fill = list(References = 0)
  )

order_duration_reference_plot <- ggplot(
  order_duration_reference_plot_data,
  aes(
    x = Duration_days,
    y = References,
    shape = Order
  )
) +
  geom_point(
    size = 2.6,
    position = position_dodge(width = 0.16)
  ) +
  facet_wrap(~ Metal) +
  scale_x_continuous(breaks = 1:4) +
  labs(
    x = "Exposure duration (days)",
    y = "Number of References",
    shape = "Order"
  ) +
  theme_classic()


# ============================================================
# 5. 96-h ALL-METAL CALANOIDA / HARPACTICOIDA SUPPORT
# ============================================================
# This is the main post-WoS domain re-check.

order_96_all <- order_two_main %>%
  filter(Duration_days == 4)

order_96_support <- order_96_all %>%
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
  arrange(Metal, Order)

# Regression checks against the already completed combined support audit.
expected_96_support <- tribble(
  ~Metal, ~Order, ~Results, ~References, ~Species, ~WoS_Results,
  "Cu", "Calanoida",      23L,  9L, 6L, 3L,
  "Cu", "Harpacticoida",  43L, 17L,10L, 0L,
  "Cd", "Calanoida",      24L,  9L, 5L, 6L,
  "Cd", "Harpacticoida",  11L,  7L, 5L, 1L,
  "Zn", "Calanoida",       8L,  2L, 4L, 0L,
  "Zn", "Harpacticoida",  19L, 11L, 6L, 4L,
  "Ag", "Calanoida",      12L,  4L, 3L, 1L,
  "Ag", "Harpacticoida",   6L,  3L, 2L, 0L,
  "Hg", "Calanoida",      11L,  5L, 7L, 4L,
  "Hg", "Harpacticoida",   1L,  1L, 1L, 0L,
  "Ni", "Calanoida",       9L,  5L, 4L, 4L,
  "Ni", "Harpacticoida",   2L,  2L, 2L, 0L
)

support_check <- order_96_support %>%
  inner_join(
    expected_96_support,
    by = c("Metal", "Order"),
    suffix = c("_Observed", "_Expected")
  )

stopifnot(
  nrow(support_check) == nrow(expected_96_support),
  all(support_check$Results_Observed == support_check$Results_Expected),
  all(support_check$References_Observed == support_check$References_Expected),
  all(support_check$Species_Observed == support_check$Species_Expected),
  all(support_check$WoS_Results_Observed == support_check$WoS_Results_Expected)
)

order_96_candidate_matrix <- order_96_support %>%
  select(
    Metal,
    Order,
    Results,
    References,
    Species,
    WoS_Results
  ) %>%
  pivot_wider(
    names_from = Order,
    values_from = c(Results, References, Species, WoS_Results),
    values_fill = 0
  ) %>%
  mutate(
    Both_Orders =
      Results_Calanoida > 0 &
      Results_Harpacticoida > 0,
    
    Min_Cell_Results =
      pmin(Results_Calanoida, Results_Harpacticoida),
    
    Min_Cell_References =
      pmin(References_Calanoida, References_Harpacticoida),
    
    Min_Cell_Species =
      pmin(Species_Calanoida, Species_Harpacticoida),
    
    Structural_note = case_when(
      !Both_Orders ~
        "No complete Calanoida-Harpacticoida support",
      
      Min_Cell_References >= 5 &
        Min_Cell_Species >= 5 ~
        "Broad two-Order support",
      
      Min_Cell_References >= 3 &
        Min_Cell_Species >= 2 ~
        "Two-Order support present but limited",
      
      TRUE ~
        "Two-Order support highly sparse"
    )
  ) %>%
  arrange(
    desc(Both_Orders),
    desc(Min_Cell_References),
    desc(Min_Cell_Species),
    Metal
  )


# ============================================================
# 6. 96-h CROSS-METAL CONNECTIVITY WITHIN EACH ORDER
# ============================================================
# For every metal pair, ask whether the same References or Species
# contribute evidence under both metals within a given Order.

available_metals <- sort(unique(order_96_all$Metal))

if (length(available_metals) >= 2) {
  
  metal_pairs <- combn(
    available_metals,
    2,
    simplify = FALSE
  )
  
  order_crossmetal_overlap_96h <- purrr::map_dfr(
    metal_pairs,
    function(pair) {
      
      metal_A <- pair[1]
      metal_B <- pair[2]
      
      order_96_all %>%
        filter(Metal %in% c(metal_A, metal_B)) %>%
        group_by(Order) %>%
        summarise(
          Refs_A = n_distinct(Reference_ID[Metal == metal_A]),
          Refs_B = n_distinct(Reference_ID[Metal == metal_B]),
          Shared_Refs = length(intersect(
            unique(Reference_ID[Metal == metal_A]),
            unique(Reference_ID[Metal == metal_B])
          )),
          Shared_Ref_IDs = paste(
            sort(intersect(
              unique(Reference_ID[Metal == metal_A]),
              unique(Reference_ID[Metal == metal_B])
            )),
            collapse = " | "
          ),
          
          Species_A = n_distinct(Species[Metal == metal_A]),
          Species_B = n_distinct(Species[Metal == metal_B]),
          Shared_Species = length(intersect(
            unique(Species[Metal == metal_A]),
            unique(Species[Metal == metal_B])
          )),
          Shared_Species_Names = paste(
            sort(intersect(
              unique(Species[Metal == metal_A]),
              unique(Species[Metal == metal_B])
            )),
            collapse = " | "
          ),
          .groups = "drop"
        ) %>%
        mutate(
          Metal_A = metal_A,
          Metal_B = metal_B,
          .before = 1
        ) %>%
        filter(Refs_A > 0, Refs_B > 0)
    }
  )
  
} else {
  
  order_crossmetal_overlap_96h <- tibble()
}


# ============================================================
# 7. REFERENCE-LEVEL ORDER OVERLAP WITHIN METAL
# ============================================================
# Diagnostic only. It is not required that the same Reference contain
# both Orders, because Orders consist of different Species.

order_reference_order_overlap <- order_96_all %>%
  distinct(
    Reference_ID,
    Metal,
    Order
  ) %>%
  group_by(
    Reference_ID,
    Metal
  ) %>%
  summarise(
    n_orders = n_distinct(Order),
    orders = paste(sort(unique(Order)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_orders == 2) %>%
  arrange(Metal, Reference_ID)


# ============================================================
# 8. PRIMARY CANDIDATE DOMAIN: 96-h Cu-Cd
# ============================================================

order_96h_cucd <- order_96_all %>%
  filter(Metal %in% c("Cu", "Cd"))

order_96h_support_cucd <- order_96h_cucd %>%
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
  arrange(Metal, Order)

stopifnot(
  nrow(order_96h_cucd) == 101,
  sum(order_96h_support_cucd$Results) == 101
)

order_model_overall <- order_96h_cucd %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental")
  )


# ============================================================
# 9. BEFORE / AFTER AUGMENTATION FOR Cu-Cd ORDER CELLS
# ============================================================

order_cucd_ecotox <- order_96h_cucd %>%
  filter(Source_Origin == "ECOTOX") %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(Evidence_Set = "ECOTOX only")

order_cucd_combined <- order_96h_cucd %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(Evidence_Set = "Combined ECOTOX + WoS")

order_cucd_before_after <- bind_rows(
  order_cucd_ecotox,
  order_cucd_combined
) %>%
  select(
    Evidence_Set,
    Metal,
    Order,
    Results,
    References,
    Species
  ) %>%
  arrange(Metal, Order, Evidence_Set)


# ============================================================
# 10. SPECIES CONNECTIVITY WITHIN EACH ORDER
# ============================================================

order_96h_species_connectivity <- order_96h_cucd %>%
  distinct(
    Order,
    Species,
    Metal,
    Reference_ID
  ) %>%
  group_by(
    Order,
    Species
  ) %>%
  summarise(
    n_metals = n_distinct(Metal),
    metals = paste(sort(unique(Metal)), collapse = " | "),
    n_references = n_distinct(Reference_ID),
    .groups = "drop"
  ) %>%
  arrange(
    Order,
    desc(n_metals),
    desc(n_references),
    Species
  )

order_96h_species_connectivity_summary <-
  order_96h_species_connectivity %>%
  group_by(Order) %>%
  summarise(
    n_species_total = n(),
    n_species_both_metals = sum(n_metals == 2),
    proportion_species_both_metals =
      n_species_both_metals / n_species_total,
    .groups = "drop"
  )


# ============================================================
# 11. REFERENCE CONNECTIVITY WITHIN EACH ORDER
# ============================================================

order_96h_reference_connectivity <- order_96h_cucd %>%
  distinct(
    Order,
    Reference_ID,
    Metal,
    Species
  ) %>%
  group_by(
    Order,
    Reference_ID
  ) %>%
  summarise(
    n_metals = n_distinct(Metal),
    metals = paste(sort(unique(Metal)), collapse = " | "),
    n_species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(
    Order,
    desc(n_metals),
    Reference_ID
  )

order_96h_reference_connectivity_summary <-
  order_96h_reference_connectivity %>%
  group_by(Order) %>%
  summarise(
    n_references_total = n(),
    n_references_both_metals = sum(n_metals == 2),
    proportion_references_both_metals =
      n_references_both_metals / n_references_total,
    .groups = "drop"
  )


# ============================================================
# 12. REFERENCE / SPECIES DOMINANCE
# ============================================================

order_96h_reference_dominance <- order_96h_cucd %>%
  count(
    Metal,
    Order,
    Reference_ID,
    name = "n_results"
  ) %>%
  group_by(Metal, Order) %>%
  mutate(
    total_results = sum(n_results),
    reference_share = n_results / total_results
  ) %>%
  arrange(Metal, Order, desc(reference_share)) %>%
  ungroup()

order_96h_reference_dominance_top <-
  order_96h_reference_dominance %>%
  group_by(Metal, Order) %>%
  slice_max(
    reference_share,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

order_96h_species_dominance <- order_96h_cucd %>%
  count(
    Metal,
    Order,
    Species,
    name = "n_results"
  ) %>%
  group_by(Metal, Order) %>%
  mutate(
    total_results = sum(n_results),
    species_share = n_results / total_results
  ) %>%
  arrange(Metal, Order, desc(species_share)) %>%
  ungroup()

order_96h_species_dominance_top <-
  order_96h_species_dominance %>%
  group_by(Metal, Order) %>%
  slice_max(
    species_share,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


# ============================================================
# 13. DEVELOPMENTAL-STAGE COMPOSITION
# ============================================================

order_96h_stage_support <- order_96h_cucd %>%
  group_by(
    Metal,
    Order,
    Lifestage
  ) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  arrange(Metal, Order, Lifestage)

order_96h_stage_completeness <- order_96h_cucd %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    stage_known = sum(!is.na(Lifestage)),
    stage_known_percent = 100 * mean(!is.na(Lifestage)),
    stage_labels = n_distinct(Lifestage, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 14. TEMPERATURE / SALINITY SUPPORT
# ============================================================
# This is a support audit for the Order model only.
# The separate environmental RQ will be revisited independently.
# pH is assessed separately in Script 08 and is not part of this Order-support audit.

order_96h_environment_support <- order_96h_cucd %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    
    temp_known = sum(!is.na(Temperature)),
    temp_percent = 100 * mean(!is.na(Temperature)),
    
    sal_known = sum(!is.na(Salinity)),
    sal_percent = 100 * mean(!is.na(Salinity)),
    
    temp_sal_both_known =
      sum(!is.na(Temperature) & !is.na(Salinity)),
    
    temp_sal_both_percent =
      100 * mean(!is.na(Temperature) & !is.na(Salinity)),
    
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    
    .groups = "drop"
  )

order_96h_joint_support <- order_96h_cucd %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

order_96h_temperature_range <- order_96h_cucd %>%
  filter(!is.na(Temperature)) %>%
  group_by(Metal, Order) %>%
  summarise(
    n = n(),
    min_temperature = min(Temperature),
    median_temperature = median(Temperature),
    max_temperature = max(Temperature),
    .groups = "drop"
  )

order_96h_salinity_range <- order_96h_cucd %>%
  filter(!is.na(Salinity)) %>%
  group_by(Metal, Order) %>%
  summarise(
    n = n(),
    min_salinity = min(Salinity),
    median_salinity = median(Salinity),
    max_salinity = max(Salinity),
    .groups = "drop"
  )

order_96h_environment_plot_data <- order_96h_cucd %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  distinct(
    Metal,
    Order,
    Reference_ID,
    Species,
    Temperature,
    Salinity,
    Source_Origin
  )

order_96h_environment_plot <- ggplot(
  order_96h_environment_plot_data,
  aes(
    x = Temperature,
    y = Salinity,
    shape = Order
  )
) +
  geom_point(
    size = 3,
    alpha = 0.70,
    position = position_jitter(
      width = 0.15,
      height = 0.15
    )
  ) +
  facet_wrap(~ Metal) +
  labs(
    x = "Temperature (C)",
    y = "Salinity",
    shape = "Order"
  ) +
  theme_classic()


# ============================================================
# 15. REFERENCE-LEVEL OUTCOME LANDSCAPE
# ============================================================
# Preserve the canonical support/landscape visual logic.

order_96h_ref <- order_96h_cucd %>%
  group_by(
    Reference_ID,
    Order,
    Metal
  ) %>%
  summarise(
    LC50_ref = exp(mean(log(LC50_umol_L))),
    n_results_in_ref = n(),
    Source_Origin = first(Source_Origin),
    .groups = "drop"
  )

order_96h_ref_summary <- order_96h_ref %>%
  group_by(Order, Metal) %>%
  summarise(
    References = n(),
    median_ref_LC50 = median(LC50_ref),
    .groups = "drop"
  )

order_96h_ref_plot <- ggplot(
  order_96h_ref,
  aes(
    x = Metal,
    y = log(LC50_ref)
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(shape = Source_Origin),
    width = 0.10,
    height = 0,
    alpha = 0.75,
    size = 2.2
  ) +
  facet_wrap(~ Order) +
  labs(
    x = "Metal",
    y = "ln(Reference-level LC50, umol/L)",
    shape = "Source"
  ) +
  theme_classic()


# ============================================================
# 16. FINAL PRE-MODEL SUPPORT DECISION
# ============================================================
# The audit documents why the final model domain is Cu-Cd at 96 h,
# comparing Calanoida with Harpacticoida. The decision is based on
# structural support, not on fitted model outcomes.

order_pre_model_decision <- tibble(
  analysis =
    "96-h Metal x Order LMM",
  
  final_model_domain =
    "96-h Cu-Cd; Calanoida vs Harpacticoida",
  
  decision_status =
    "FINAL DOMAIN SUPPORTED",
  
  rationale = paste(
    "Cu and Cd retain the broadest reusable two-Order structure after WoS augmentation.",
    "Zn has both Orders represented but Calanoida remains supported by only two References,",
    "and other metals retain still sparser two-Order support.",
    "Cross-metal Reference/Species connectivity, dominance, Stage composition and",
    "temperature/salinity common support are retained as structural support diagnostics.",
    "No metal is promoted solely because its Result count increased."
  )
)


# ============================================================
# 17. BUILD FINAL Cu-Cd ORDER MODEL DATASET
# ============================================================

order_model_data <- order_96h_cucd %>%
  mutate(
    Metal = factor(
      Metal,
      levels = c("Cu", "Cd")
    ),
    Order = factor(
      Order,
      levels = c("Calanoida", "Harpacticoida")
    ),
    Reference_ID = factor(Reference_ID),
    Species = factor(Species)
  ) %>%
  droplevels()


# Frozen provenance check for the final 101-result domain:
# 91 ECOTOX Results + 10 source-verified WoS supplemental Results.
stopifnot(
  nrow(order_model_data) == 101,
  n_distinct(order_model_data$Test_ID) == 101,
  n_distinct(order_model_data$Reference_ID) == 37,
  n_distinct(order_model_data$Species) == 18,
  sum(order_model_data$Source_Origin == "ECOTOX") == 91,
  sum(order_model_data$Source_Origin == "WoS_supplemental") == 10
)


# ============================================================
# 18. SAVE TABLES
# ============================================================

output_tables <- list(
  architecture_overall = architecture_overall,
  order_architecture = order_architecture,
  order_metal_support = order_metal_support,
  species_architecture = species_architecture,
  order_duration_support = order_duration_support,
  order_duration_range = order_duration_range,
  order_96_support = order_96_support,
  order_96_candidate_matrix = order_96_candidate_matrix,
  order_crossmetal_overlap_96h = order_crossmetal_overlap_96h,
  order_reference_order_overlap = order_reference_order_overlap,
  order_96h_cucd_support = order_96h_support_cucd,
  order_cucd_before_after = order_cucd_before_after,
  order_model_overall = order_model_overall,
  order_96h_species_connectivity = order_96h_species_connectivity,
  order_96h_species_connectivity_summary = order_96h_species_connectivity_summary,
  order_96h_reference_connectivity = order_96h_reference_connectivity,
  order_96h_reference_connectivity_summary = order_96h_reference_connectivity_summary,
  order_96h_reference_dominance = order_96h_reference_dominance,
  order_96h_reference_dominance_top = order_96h_reference_dominance_top,
  order_96h_species_dominance = order_96h_species_dominance,
  order_96h_species_dominance_top = order_96h_species_dominance_top,
  order_96h_stage_support = order_96h_stage_support,
  order_96h_stage_completeness = order_96h_stage_completeness,
  order_96h_environment_support = order_96h_environment_support,
  order_96h_joint_support = order_96h_joint_support,
  order_96h_temperature_range = order_96h_temperature_range,
  order_96h_salinity_range = order_96h_salinity_range,
  order_96h_ref_summary = order_96h_ref_summary,
  order_pre_model_decision = order_pre_model_decision
)

purrr::iwalk(
  output_tables,
  ~ write_csv(
    .x,
    file.path(
      table_dir,
      paste0(.y, ".csv")
    )
  )
)


# ============================================================
# 19. SAVE FIGURES
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "order_duration_reference_support.png"
  ),
  order_duration_reference_plot,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "order_96h_temperature_salinity_common_support.png"
  ),
  order_96h_environment_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "order_96h_reference_level_landscape.png"
  ),
  order_96h_ref_plot,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 20. SAVE OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    order_model_data =
      order_model_data,
    
    order_96_candidate_matrix =
      order_96_candidate_matrix,
    
    order_96h_cucd =
      order_96h_cucd,
    
    output_tables =
      output_tables
  ),
  file.path(
    object_dir,
    "order_support_audit_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_root,
    "sessionInfo.txt"
  )
)


input_output_manifest <- tibble::tibble(
  Role = c(
    "Input: combined harmonized dataset",
    "Output: 96-h candidate matrix",
    "Output: final Cu-Cd support table",
    "Output: final model-domain decision",
    "Output: reusable support-audit objects"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "order_96_candidate_matrix.csv"),
    file.path(table_dir, "order_96h_cucd_support.csv"),
    file.path(table_dir, "order_pre_model_decision.csv"),
    file.path(object_dir, "order_support_audit_objects.rds")
  )
)

write_csv(
  input_output_manifest,
  file.path(output_root, "input_output_manifest.csv")
)


# ============================================================
# 21. FINAL OUTPUT QA
# ============================================================

required_outputs <- c(
  file.path(
    table_dir,
    "order_96_candidate_matrix.csv"
  ),
  file.path(
    table_dir,
    "order_96h_cucd_support.csv"
  ),
  file.path(
    table_dir,
    "order_cucd_before_after.csv"
  ),
  file.path(
    table_dir,
    "order_crossmetal_overlap_96h.csv"
  ),
  file.path(
    table_dir,
    "order_96h_species_connectivity_summary.csv"
  ),
  file.path(
    table_dir,
    "order_96h_reference_connectivity_summary.csv"
  ),
  file.path(
    table_dir,
    "order_pre_model_decision.csv"
  ),
  file.path(
    figure_dir,
    "order_96h_reference_level_landscape.png"
  ),
  file.path(
    object_dir,
    "order_support_audit_objects.rds"
  )
)

stopifnot(all(file.exists(required_outputs)))


# ============================================================
# 22. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("12 TAXONOMIC-ORDER SUPPORT AUDIT\n")
cat("============================================================\n")

cat("\n--- 96-h ALL-METAL ORDER CANDIDATE MATRIX ---\n")
print(
  order_96_candidate_matrix,
  n = Inf,
  width = Inf
)

cat("\n--- 96-h Cu-Cd ORDER SUPPORT ---\n")
print(
  order_96h_support_cucd,
  n = Inf,
  width = Inf
)

cat("\n--- Cu-Cd BEFORE / AFTER AUGMENTATION ---\n")
print(
  order_cucd_before_after,
  n = Inf,
  width = Inf
)

cat("\n--- FINAL MODEL DATASET OVERALL ---\n")
print(
  order_model_overall,
  width = Inf
)

cat("\n--- SPECIES CONNECTIVITY WITHIN ORDER ---\n")
print(
  order_96h_species_connectivity_summary,
  n = Inf,
  width = Inf
)

cat("\n--- REFERENCE CONNECTIVITY WITHIN ORDER ---\n")
print(
  order_96h_reference_connectivity_summary,
  n = Inf,
  width = Inf
)

cat("\n--- TOP REFERENCE DOMINANCE ---\n")
print(
  order_96h_reference_dominance_top,
  n = Inf,
  width = Inf
)

cat("\n--- TOP SPECIES DOMINANCE ---\n")
print(
  order_96h_species_dominance_top,
  n = Inf,
  width = Inf
)

cat("\n--- ENVIRONMENTAL SUPPORT ---\n")
print(
  order_96h_environment_support,
  n = Inf,
  width = Inf
)

cat("\n--- FINAL PRE-MODEL SUPPORT DECISION ---\n")
print(
  order_pre_model_decision,
  width = Inf
)

writeLines(
  c(
    "STATUS: 12 AUDIT ORDER SUPPORT = PASS",
    "Final domain: 96-h Cu-Cd; Calanoida vs Harpacticoida",
    paste0("Results: ", nrow(order_model_data)),
    paste0("References: ", n_distinct(order_model_data$Reference_ID)),
    paste0("Species: ", n_distinct(order_model_data$Species))
  ),
  file.path(output_root, "RUN_COMPLETE.txt")
)

cat(
  "\nSTATUS: 12 AUDIT ORDER SUPPORT = PASS\n"
)
# ============================================================
# END
# ============================================================
