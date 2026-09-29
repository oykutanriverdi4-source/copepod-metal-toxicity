# ============================================================
# COPEPOD METAL TOXICITY THESIS REPRODUCIBILITY WORKFLOW
# 10_audit_stage_support.R
# ============================================================
# PURPOSE
# Audit structural support for the 96-h developmental-stage model
# after combining ECOTOX and source-verified supplementary records.
#
# This is a pre-model evidence audit. It does NOT fit the Stage model.
# It evaluates:
#   - Adult/Nauplii support by metal and duration;
#   - Reference and Species breadth;
#   - cross-metal connectivity;
#   - dominance by individual References or Species;
#   - environmental metadata support; and
#   - the final 96-h Cu-Cd Adult/Nauplii model domain used in the thesis.
#
# SCIENTIFIC SCOPE
# The analytical definitions and support calculations are preserved from
# the original stage-support audit. Repository cleanup changes only paths,
# file names, documentation, output organization, and final-status wording.
# No developmental-stage model is fitted in this script.
# ============================================================


# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)

options(width = 220)

output_root <- file.path(
  "outputs",
  "10_audit_stage_support"
)

table_dir <- file.path(output_root, "01_tables")
figure_dir <- file.path(output_root, "02_figures")
object_dir <- file.path(output_root, "03_objects")

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
  "Species",
  "Metal",
  "Duration_days",
  "Lifestage",
  "External_Stage_Sex",
  "Harmonization_Status",
  "LC50_umol_L",
  "Temperature",
  "Salinity",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Source_Origin"
)

stopifnot(all(required_columns %in% names(dat)))

stopifnot(
  nrow(dat) == 353,
  n_distinct(dat$Reference_ID) == 84,
  sum(dat$Source_Origin == "ECOTOX") == 296,
  sum(dat$Source_Origin == "WoS_supplemental") == 57
)


# ============================================================
# 2. GENERAL STAGE ARCHITECTURE
# ============================================================

stage_data <- dat %>%
  mutate(
    Lifestage = na_if(trimws(as.character(Lifestage)), ""),
    External_Stage_Sex = na_if(trimws(as.character(External_Stage_Sex)), ""),
    Stage_text = coalesce(External_Stage_Sex, Lifestage)
  ) %>%
  filter(
    Harmonization_Status == "HARMONIZED",
    !is.na(LC50_umol_L),
    LC50_umol_L > 0,
    !is.na(Stage_text)
  ) %>%
  mutate(
    Stage_group = case_when(
      str_detect(str_to_lower(Stage_text), "adult") ~ "Adult",
      str_detect(str_to_lower(Stage_text), "naupl") ~ "Nauplii",
      str_detect(str_to_lower(Stage_text), "copepodite|copepodid|cv|iv") ~ "Copepodite",
      str_detect(str_to_lower(Stage_text), "egg") ~ "Egg",
      str_detect(str_to_lower(Stage_text), "neonate") ~ "Neonate",
      str_detect(str_to_lower(Stage_text), "juvenile") ~ "Juvenile",
      TRUE ~ NA_character_
    ),
    Stage_role = case_when(
      Stage_group %in% c("Adult", "Nauplii") ~ "PRIMARY_STAGE_CONTRAST",
      Stage_group %in% c("Copepodite", "Egg") ~ "DEVELOPMENTAL_SECONDARY",
      Stage_group %in% c("Neonate", "Juvenile") ~ "AMBIGUOUS_OR_SPARSE",
      TRUE ~ "NOT_STAGE_ANALYSIS"
    ),
    ln_LC50 = log(LC50_umol_L)
  )

stage_label_audit <- stage_data %>%
  count(
    Source_Origin,
    Stage_text,
    Stage_group,
    Stage_role,
    Metal,
    sort = TRUE
  )

stage_support_overall <- stage_data %>%
  filter(!is.na(Stage_group)) %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Metal, Stage_group)


# ============================================================
# 3. PRIMARY ADULT / NAUPLII ARCHITECTURE
# ============================================================

stage_primary <- stage_data %>%
  filter(Stage_role == "PRIMARY_STAGE_CONTRAST")

stage_primary_support <- stage_primary %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Metal, Stage_group)

stage_duration_support <- stage_primary %>%
  group_by(Metal, Stage_group, Duration_days) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Metal, Stage_group, Duration_days)


# ============================================================
# 4. ADULT-NAUPLII WITHIN-METAL COMPOSITION OVERLAP
# ============================================================

stage_composition_overlap <- stage_primary %>%
  group_by(Metal, Duration_days) %>%
  summarise(
    Adult_Results = sum(Stage_group == "Adult"),
    Nauplii_Results = sum(Stage_group == "Nauplii"),
    
    Adult_References = n_distinct(
      Reference_ID[Stage_group == "Adult"]
    ),
    Nauplii_References = n_distinct(
      Reference_ID[Stage_group == "Nauplii"]
    ),
    
    Shared_References = length(intersect(
      unique(Reference_ID[Stage_group == "Adult"]),
      unique(Reference_ID[Stage_group == "Nauplii"])
    )),
    Shared_Reference_IDs = paste(
      sort(intersect(
        unique(Reference_ID[Stage_group == "Adult"]),
        unique(Reference_ID[Stage_group == "Nauplii"])
      )),
      collapse = " | "
    ),
    
    Adult_Species = n_distinct(
      Species[Stage_group == "Adult"]
    ),
    Nauplii_Species = n_distinct(
      Species[Stage_group == "Nauplii"]
    ),
    
    Shared_Species = length(intersect(
      unique(Species[Stage_group == "Adult"]),
      unique(Species[Stage_group == "Nauplii"])
    )),
    Shared_Species_Names = paste(
      sort(intersect(
        unique(Species[Stage_group == "Adult"]),
        unique(Species[Stage_group == "Nauplii"])
      )),
      collapse = " | "
    ),
    .groups = "drop"
  ) %>%
  filter(
    Adult_References > 0,
    Nauplii_References > 0
  ) %>%
  arrange(Metal, Duration_days)


# ============================================================
# 5. 96-H ALL-METAL STAGE SUPPORT MATRIX
# ============================================================
# This matrix documents the all-metal 96-h Adult/Nauplii evidence
# used to assess whether a broader model domain was structurally supported.

stage_96 <- stage_primary %>%
  filter(Duration_days == 4)

stage_96_support <- stage_96 %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Metal, Stage_group)

stage_96_candidate_matrix <- stage_96_support %>%
  select(
    Metal,
    Stage_group,
    Results,
    References,
    Species,
    WoS_Results
  ) %>%
  pivot_wider(
    names_from = Stage_group,
    values_from = c(
      Results,
      References,
      Species,
      WoS_Results
    ),
    values_fill = 0
  ) %>%
  mutate(
    Both_Adult_Nauplii =
      Results_Adult > 0 & Results_Nauplii > 0,
    
    Min_Cell_Results =
      pmin(Results_Adult, Results_Nauplii),
    
    Min_Cell_References =
      pmin(References_Adult, References_Nauplii),
    
    Min_Cell_Species =
      pmin(Species_Adult, Species_Nauplii),
    
    Structural_note = case_when(
      !Both_Adult_Nauplii ~
        "No complete Adult-Nauplii support",
      
      Min_Cell_References >= 4 &
        Min_Cell_Species >= 4 ~
        "Comparatively broad two-stage support",
      
      Min_Cell_References >= 2 ~
        "Two-stage support present but limited",
      
      TRUE ~
        "Two-stage support highly sparse"
    )
  ) %>%
  arrange(
    desc(Both_Adult_Nauplii),
    desc(Min_Cell_References),
    desc(Min_Cell_Species),
    Metal
  )


# ============================================================
# 6. STAGE-SPECIFIC CROSS-METAL CONNECTIVITY LANDSCAPE
# ============================================================

stage_crossmetal <- stage_primary %>%
  filter(Duration_days == 4)

available_metals <- sort(unique(stage_crossmetal$Metal))

if (length(available_metals) >= 2) {
  
  metal_pairs <- combn(
    available_metals,
    2,
    simplify = FALSE
  )
  
  stage_crossmetal_overlap <- purrr::map_dfr(
    metal_pairs,
    function(pair) {
      
      metal_A <- pair[1]
      metal_B <- pair[2]
      
      stage_crossmetal %>%
        filter(Metal %in% c(metal_A, metal_B)) %>%
        group_by(Stage_group) %>%
        summarise(
          Refs_A = n_distinct(
            Reference_ID[Metal == metal_A]
          ),
          Refs_B = n_distinct(
            Reference_ID[Metal == metal_B]
          ),
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
          
          Species_A = n_distinct(
            Species[Metal == metal_A]
          ),
          Species_B = n_distinct(
            Species[Metal == metal_B]
          ),
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
  
  stage_crossmetal_overlap <- tibble()
}


# ============================================================
# 7. FINAL MODEL DOMAIN: 96-H Cu-Cd ADULT/NAUPLII
# ============================================================
# The final thesis model domain is Cu-Cd at 96 h with Adult and Nauplii
# represented in every required Metal x Stage cell. The all-metal support
# matrix above is retained to document why other metals were not promoted
# to the same model.

stage_cucd_96 <- stage_96 %>%
  filter(Metal %in% c("Cu", "Cd"))

stage_cucd_96_support <- stage_cucd_96 %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(Metal, Stage_group)

stage_model_overall <- stage_cucd_96 %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental")
  )


# ============================================================
# 8. BEFORE / AFTER CELL COUNTS
# ============================================================

stage_cucd_96_ECOTOX <- stage_cucd_96 %>%
  filter(Source_Origin == "ECOTOX") %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(Evidence_Set = "ECOTOX only")

stage_cucd_96_combined <- stage_cucd_96 %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(Evidence_Set = "Combined ECOTOX + WoS")

stage_cucd_96_before_after <- bind_rows(
  stage_cucd_96_ECOTOX,
  stage_cucd_96_combined
) %>%
  select(
    Evidence_Set,
    Metal,
    Stage_group,
    Results,
    References,
    Species
  ) %>%
  arrange(Metal, Stage_group, Evidence_Set)


# ============================================================
# 9. PRE-MODEL ENVIRONMENTAL SUPPORT
# ============================================================

stage_environment_support <- stage_cucd_96 %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    
    Temp_known =
      sum(!is.na(Temperature)),
    Temp_percent =
      100 * mean(!is.na(Temperature)),
    
    Sal_known =
      sum(!is.na(Salinity)),
    Sal_percent =
      100 * mean(!is.na(Salinity)),
    
    Temp_sal_both_known =
      sum(!is.na(Temperature) & !is.na(Salinity)),
    Temp_sal_both_percent =
      100 * mean(!is.na(Temperature) & !is.na(Salinity)),
    
    References =
      n_distinct(Reference_ID),
    
    Species =
      n_distinct(Species),
    
    .groups = "drop"
  )

stage_joint_support <- stage_cucd_96 %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

stage_temperature_overlap <- stage_cucd_96 %>%
  filter(!is.na(Temperature)) %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    n = n(),
    min_temperature = min(Temperature),
    median_temperature = median(Temperature),
    max_temperature = max(Temperature),
    .groups = "drop"
  )

stage_salinity_overlap <- stage_cucd_96 %>%
  filter(!is.na(Salinity)) %>%
  group_by(Metal, Stage_group) %>%
  summarise(
    n = n(),
    min_salinity = min(Salinity),
    median_salinity = median(Salinity),
    max_salinity = max(Salinity),
    .groups = "drop"
  )

stage_environment_plot_data <- stage_cucd_96 %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  distinct(
    Metal,
    Stage_group,
    Reference_ID,
    Species,
    Temperature,
    Salinity,
    Source_Origin
  )

stage_environment_plot <- ggplot(
  stage_environment_plot_data,
  aes(
    x = Temperature,
    y = Salinity,
    shape = Stage_group
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
    shape = "Stage"
  ) +
  theme_classic()


# ============================================================
# 10. CROSS-METAL REFERENCE / SPECIES CONNECTIVITY
# ============================================================

stage_reference_connectivity_detail <- stage_cucd_96 %>%
  distinct(
    Stage_group,
    Reference_ID,
    Metal
  ) %>%
  group_by(
    Stage_group,
    Reference_ID
  ) %>%
  summarise(
    n_metals = n_distinct(Metal),
    metals = paste(
      sort(unique(Metal)),
      collapse = " | "
    ),
    .groups = "drop"
  )

stage_reference_connectivity <- stage_reference_connectivity_detail %>%
  group_by(Stage_group) %>%
  summarise(
    References_total = n(),
    References_both_CuCd =
      sum(n_metals == 2),
    .groups = "drop"
  )

stage_species_connectivity_detail <- stage_cucd_96 %>%
  distinct(
    Stage_group,
    Species,
    Metal
  ) %>%
  group_by(
    Stage_group,
    Species
  ) %>%
  summarise(
    n_metals = n_distinct(Metal),
    metals = paste(
      sort(unique(Metal)),
      collapse = " | "
    ),
    .groups = "drop"
  )

stage_species_connectivity <- stage_species_connectivity_detail %>%
  group_by(Stage_group) %>%
  summarise(
    Species_total = n(),
    Species_both_CuCd =
      sum(n_metals == 2),
    .groups = "drop"
  )


# ============================================================
# 11. REFERENCE / SPECIES DOMINANCE
# ============================================================

stage_reference_dominance <- stage_cucd_96 %>%
  count(
    Metal,
    Stage_group,
    Reference_ID,
    name = "n_results"
  ) %>%
  group_by(
    Metal,
    Stage_group
  ) %>%
  mutate(
    cell_total = sum(n_results),
    reference_share =
      n_results / cell_total
  ) %>%
  arrange(
    Metal,
    Stage_group,
    desc(reference_share)
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

stage_species_dominance <- stage_cucd_96 %>%
  count(
    Metal,
    Stage_group,
    Species,
    name = "n_results"
  ) %>%
  group_by(
    Metal,
    Stage_group
  ) %>%
  mutate(
    cell_total = sum(n_results),
    species_share =
      n_results / cell_total
  ) %>%
  arrange(
    Metal,
    Stage_group,
    desc(species_share)
  ) %>%
  slice_head(n = 1) %>%
  ungroup()


# ============================================================
# 12. FULLY-CROSSED SUPPORT
# ============================================================
# Does any Reference x Species contain all four cells?
# Adult-Cd, Adult-Cu, Nauplii-Cd, Nauplii-Cu.

stage_cucd_96_cross <- stage_cucd_96 %>%
  distinct(
    Reference_ID,
    Species,
    Stage_group,
    Metal
  ) %>%
  group_by(
    Reference_ID,
    Species
  ) %>%
  summarise(
    n_metals = n_distinct(Metal),
    n_stages = n_distinct(Stage_group),
    n_cells = n(),
    metals = paste(
      sort(unique(Metal)),
      collapse = ", "
    ),
    stages = paste(
      sort(unique(Stage_group)),
      collapse = ", "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells))

fully_crossed_stage_cucd <- stage_cucd_96_cross %>%
  filter(
    n_metals == 2,
    n_stages == 2,
    n_cells == 4
  )


# ============================================================
# 13. DESCRIPTIVE OUTCOME LANDSCAPE
# ============================================================
# Preserves the original reference-level landscape figure.

stage_cucd_96_result_summary <- stage_cucd_96 %>%
  group_by(Stage_group, Metal) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    median_LC50 = median(LC50_umol_L),
    median_ln_LC50 = median(ln_LC50),
    .groups = "drop"
  )

stage_cucd_96_ref <- stage_cucd_96 %>%
  group_by(
    Reference_ID,
    Stage_group,
    Metal
  ) %>%
  summarise(
    LC50_ref =
      exp(mean(log(LC50_umol_L))),
    n_results_in_ref = n(),
    Source_Origin = first(Source_Origin),
    .groups = "drop"
  )

stage_cucd_96_ref_summary <- stage_cucd_96_ref %>%
  group_by(Stage_group, Metal) %>%
  summarise(
    References = n(),
    median_ref_LC50 = median(LC50_ref),
    .groups = "drop"
  )

stage_cucd_96_ref_plot <- ggplot(
  stage_cucd_96_ref,
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
  facet_wrap(~ Stage_group) +
  labs(
    x = "Metal",
    y = "ln(Reference-level LC50, umol/L)",
    shape = "Source"
  ) +
  theme_classic()


# ============================================================
# 14. BROADER DURATION DIAGNOSTIC
# ============================================================

stage_cucd_broad <- stage_primary %>%
  filter(
    Metal %in% c("Cu", "Cd"),
    Duration_days %in% c(1, 2, 3, 4)
  )

stage_cucd_broad_refdur <- stage_cucd_broad %>%
  group_by(
    Reference_ID,
    Stage_group,
    Metal,
    Duration_days
  ) %>%
  summarise(
    LC50_refdur =
      exp(mean(log(LC50_umol_L))),
    n_results_in_cell = n(),
    Source_Origin = first(Source_Origin),
    .groups = "drop"
  )

stage_cucd_broad_summary <- stage_cucd_broad_refdur %>%
  group_by(
    Stage_group,
    Metal,
    Duration_days
  ) %>%
  summarise(
    Reference_duration_cells = n(),
    median_LC50 =
      median(LC50_refdur),
    .groups = "drop"
  ) %>%
  arrange(
    Stage_group,
    Duration_days,
    Metal
  )

stage_cucd_broad_plot <- ggplot(
  stage_cucd_broad_refdur,
  aes(
    x = Metal,
    y = log(LC50_refdur)
  )
) +
  geom_boxplot(
    aes(group = Metal),
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(
      shape = factor(Duration_days)
    ),
    width = 0.10,
    height = 0,
    alpha = 0.75,
    size = 2.4
  ) +
  facet_wrap(~ Stage_group) +
  labs(
    x = "Metal",
    y = "ln(Reference-duration-level LC50, umol/L)",
    shape = "Duration (days)"
  ) +
  theme_classic()


# ============================================================
# 15. SHARED-SPECIES COMPOSITION SENSITIVITY
# ============================================================

shared_species_96 <- stage_cucd_96 %>%
  distinct(
    Stage_group,
    Species,
    Metal
  ) %>%
  count(
    Stage_group,
    Species,
    name = "n_metals"
  ) %>%
  filter(n_metals == 2)

stage_cucd_96_sharedsp <- stage_cucd_96 %>%
  semi_join(
    shared_species_96,
    by = c(
      "Stage_group",
      "Species"
    )
  )

stage_shared_species_summary <- stage_cucd_96_sharedsp %>%
  group_by(
    Stage_group,
    Metal
  ) %>%
  summarise(
    Results = n(),
    References =
      n_distinct(Reference_ID),
    Species =
      n_distinct(Species),
    median_LC50 =
      median(LC50_umol_L),
    .groups = "drop"
  )

stage_cucd_96_sharedsp_refsp <- stage_cucd_96_sharedsp %>%
  group_by(
    Reference_ID,
    Species,
    Stage_group,
    Metal
  ) %>%
  summarise(
    LC50_refsp =
      exp(mean(log(LC50_umol_L))),
    n_results_in_cell = n(),
    .groups = "drop"
  )

stage_shared_species_refsp_summary <- stage_cucd_96_sharedsp_refsp %>%
  group_by(
    Stage_group,
    Metal
  ) %>%
  summarise(
    Ref_species_cells = n(),
    References =
      n_distinct(Reference_ID),
    Species =
      n_distinct(Species),
    median_LC50 =
      median(LC50_refsp),
    .groups = "drop"
  )


# ============================================================
# 16. INSPECT THE ORIGINAL FULLY-CROSSED ECOTOX CONTEXT
# ============================================================

ref14137_stage_2x2 <- stage_cucd_96 %>%
  filter(
    Reference_ID == "E_14137",
    Species == "Tisbe battagliai"
  ) %>%
  select(
    Reference_ID,
    Test_ID,
    Result_ID,
    Source_Origin,
    Species,
    Metal,
    Stage_group,
    Stage_text,
    Duration_days,
    Temperature,
    Salinity,
    Exposure_Type,
    Conc_Type,
    LC50_umol_L
  ) %>%
  arrange(Stage_group, Metal)


# ============================================================
# 17. FINAL MODEL-DOMAIN SNAPSHOT
# ============================================================
# No automatic model expansion is performed here.
# This table records the final thesis domain selected from the support audit.

stage_pre_model_decision <- tibble(
  analysis =
    "96-h Metal x Adult/Nauplii Stage LMM",
  
  final_model_domain =
    "Cu-Cd",
  
  decision_status =
    "FINAL THESIS DOMAIN CONFIRMED",
  
  rationale = paste(
    "Cu and Cd provide the strongest joint 96-h Adult/Nauplii support across",
    "Results, References, and Species groups, with both required stages represented",
    "for each metal. Zn contains Adult and Nauplii evidence but with narrower",
    "Reference and Species support, while the remaining metals do not provide",
    "complete Adult/Nauplii support. The selected model domain therefore remains",
    "Cu-Cd; no metal is promoted solely because supplementary records increased",
    "its Result count."
  )
)


# ============================================================
# 18. BUILD FINAL Cu-Cd MODEL DATASET
# ============================================================

stage_model_data <- stage_cucd_96 %>%
  mutate(
    Metal = factor(
      Metal,
      levels = c("Cu", "Cd")
    ),
    Stage = factor(
      Stage_group,
      levels = c("Adult", "Nauplii")
    ),
    Reference_ID = factor(Reference_ID),
    Species = factor(Species)
  ) %>%
  droplevels()


# ============================================================
# 19. SAVE TABLES
# ============================================================

output_tables <- list(
  stage_label_audit = stage_label_audit,
  stage_support_overall = stage_support_overall,
  stage_primary_support = stage_primary_support,
  stage_duration_support = stage_duration_support,
  stage_composition_overlap = stage_composition_overlap,
  stage_96_support = stage_96_support,
  stage_96_candidate_matrix = stage_96_candidate_matrix,
  stage_crossmetal_overlap_96h = stage_crossmetal_overlap,
  stage_cucd_96_support = stage_cucd_96_support,
  stage_cucd_96_before_after = stage_cucd_96_before_after,
  stage_model_overall = stage_model_overall,
  stage_environment_support = stage_environment_support,
  stage_joint_support = stage_joint_support,
  stage_temperature_overlap = stage_temperature_overlap,
  stage_salinity_overlap = stage_salinity_overlap,
  stage_reference_connectivity_detail = stage_reference_connectivity_detail,
  stage_reference_connectivity = stage_reference_connectivity,
  stage_species_connectivity_detail = stage_species_connectivity_detail,
  stage_species_connectivity = stage_species_connectivity,
  stage_reference_dominance = stage_reference_dominance,
  stage_species_dominance = stage_species_dominance,
  stage_cucd_96_cross = stage_cucd_96_cross,
  fully_crossed_stage_cucd = fully_crossed_stage_cucd,
  stage_pre_model_decision = stage_pre_model_decision,
  stage_cucd_96_result_summary = stage_cucd_96_result_summary,
  stage_cucd_96_ref_summary = stage_cucd_96_ref_summary,
  stage_cucd_broad_summary = stage_cucd_broad_summary,
  shared_species_96 = shared_species_96,
  stage_shared_species_summary = stage_shared_species_summary,
  stage_shared_species_refsp_summary = stage_shared_species_refsp_summary,
  stage_environment_plot_data = stage_environment_plot_data
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

write_csv(
  ref14137_stage_2x2,
  file.path(
    table_dir,
    "ref14137_stage_2x2_combined.csv"
  )
)

write_csv(
  stage_model_data %>%
    mutate(
      Metal = as.character(Metal),
      Stage = as.character(Stage),
      Reference_ID = as.character(Reference_ID),
      Species = as.character(Species)
    ),
  file.path(table_dir, "stage_model_data.csv")
)


# ============================================================
# 20. SAVE FIGURES
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "stage_96h_environment_common_support.png"
  ),
  stage_environment_plot,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "stage_96h_reference_level_landscape.png"
  ),
  stage_cucd_96_ref_plot,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "stage_24_96h_duration_diagnostic.png"
  ),
  stage_cucd_broad_plot,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 21. SAVE OBJECTS
# ============================================================

saveRDS(
  list(
    stage_model_data =
      stage_model_data,
    
    stage_96_candidate_matrix =
      stage_96_candidate_matrix,
    
    stage_cucd_96 =
      stage_cucd_96,
    
    output_tables =
      output_tables
  ),
  file.path(
    object_dir,
    "stage_support_audit_objects.rds"
  )
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_root,
    "sessionInfo.txt"
  )
)

# Reproducibility manifest
manifest <- tibble(
  Role = c(
    "Input: canonical combined dataset",
    "Output: all-metal 96-h stage candidate matrix",
    "Output: final Cu-Cd 96-h stage support",
    "Output: final stage model data object",
    "Output: R analysis objects"
  ),
  Path = c(
    input_file,
    file.path(table_dir, "stage_96_candidate_matrix.csv"),
    file.path(table_dir, "stage_cucd_96_support.csv"),
    file.path(table_dir, "stage_model_data.csv"),
    file.path(object_dir, "stage_support_audit_objects.rds")
  )
) %>%
  mutate(
    Exists = file.exists(Path),
    MD5 = if_else(Exists, as.character(tools::md5sum(Path)), NA_character_)
  )

write_csv(
  manifest,
  file.path(output_root, "input_output_manifest.csv")
)


# ============================================================
# 22. SCIENTIFIC BENCHMARK QA
# ============================================================
# These checks encode the final support architecture reported in the thesis.
# They are assertions only; they do not alter the data or analysis.

expected_stage_cells <- tribble(
  ~Metal, ~Stage_group, ~Results, ~Tests, ~References, ~Species, ~ECOTOX_Results, ~WoS_Results,
  "Cu", "Adult",   31L, 31L, 17L, 13L, 28L, 3L,
  "Cu", "Nauplii", 14L, 14L,  7L,  5L, 14L, 0L,
  "Cd", "Adult",   19L, 19L, 11L, 11L, 14L, 5L,
  "Cd", "Nauplii",  5L,  5L,  4L,  5L,  4L, 1L
)

observed_stage_cells <- stage_cucd_96_support %>%
  arrange(Metal, Stage_group) %>%
  select(Metal, Stage_group, Results, Tests, References, Species, ECOTOX_Results, WoS_Results)

expected_stage_cells <- expected_stage_cells %>%
  arrange(Metal, Stage_group)

stopifnot(identical(observed_stage_cells, expected_stage_cells))

stopifnot(
  nrow(stage_model_data) == 69,
  n_distinct(stage_model_data$Test_ID) == 69,
  n_distinct(stage_model_data$Reference_ID) == 31,
  n_distinct(stage_model_data$Species) == 19,
  sum(stage_model_data$Source_Origin == "ECOTOX") == 60,
  sum(stage_model_data$Source_Origin == "WoS_supplemental") == 9
)

stopifnot(
  nrow(fully_crossed_stage_cucd) == 1,
  as.character(fully_crossed_stage_cucd$Reference_ID[[1]]) == "E_14137",
  as.character(fully_crossed_stage_cucd$Species[[1]]) == "Tisbe battagliai"
)


# ============================================================
# 23. FINAL OUTPUT QA
# ============================================================

required_outputs <- c(
  file.path(
    table_dir,
    "stage_96_candidate_matrix.csv"
  ),
  file.path(
    table_dir,
    "stage_cucd_96_support.csv"
  ),
  file.path(
    table_dir,
    "stage_cucd_96_before_after.csv"
  ),
  file.path(
    table_dir,
    "stage_crossmetal_overlap_96h.csv"
  ),
  file.path(
    table_dir,
    "stage_pre_model_decision.csv"
  ),
  file.path(
    table_dir,
    "fully_crossed_stage_cucd.csv"
  ),
  file.path(
    figure_dir,
    "stage_96h_environment_common_support.png"
  ),
  file.path(
    figure_dir,
    "stage_96h_reference_level_landscape.png"
  ),
  file.path(
    object_dir,
    "stage_support_audit_objects.rds"
  ),
  file.path(table_dir, "stage_model_data.csv"),
  file.path(output_root, "input_output_manifest.csv")
)

stopifnot(all(file.exists(required_outputs)))

writeLines(
  c(
    "STATUS: 10 AUDIT STAGE SUPPORT = PASS",
    paste0("Stage model Results: ", nrow(stage_model_data)),
    paste0("References: ", n_distinct(stage_model_data$Reference_ID)),
    paste0("Species: ", n_distinct(stage_model_data$Species)),
    "Final model domain: 96-h Cu-Cd Adult/Nauplii"
  ),
  con = file.path(output_root, "RUN_COMPLETE.txt")
)


# ============================================================
# 24. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("10 STAGE SUPPORT AUDIT\n")
cat("============================================================\n")

cat("\n--- 96-h ALL-METAL ADULT/NAUPLII CANDIDATE MATRIX ---\n")
print(stage_96_candidate_matrix, n = Inf, width = Inf)

cat("\n--- Cu-Cd 96-h CELL SUPPORT ---\n")
print(stage_cucd_96_support, n = Inf, width = Inf)

cat("\n--- Cu-Cd BEFORE / AFTER AUGMENTATION ---\n")
print(stage_cucd_96_before_after, n = Inf, width = Inf)

cat("\n--- FINAL Cu-Cd 96-h DATASET OVERALL ---\n")
print(stage_model_overall, width = Inf)

cat("\n--- REFERENCE CONNECTIVITY ---\n")
print(stage_reference_connectivity, n = Inf, width = Inf)

cat("\n--- SPECIES CONNECTIVITY ---\n")
print(stage_species_connectivity, n = Inf, width = Inf)

cat("\n--- FULLY-CROSSED Cu-Cd SUPPORT ---\n")
print(fully_crossed_stage_cucd, n = Inf, width = Inf)

cat("\n--- ENVIRONMENTAL SUPPORT ---\n")
print(stage_environment_support, n = Inf, width = Inf)

cat("\n--- FINAL MODEL-DOMAIN DECISION ---\n")
print(stage_pre_model_decision, width = Inf)

cat("\nSTATUS: 10 AUDIT STAGE SUPPORT = PASS\n")
# ============================================================
# END
# ============================================================
