# ============================================================
# 11_fit_stage_model.R
# ============================================================
# Purpose
# Fit the final 96-h Cu-Cd developmental-stage linear mixed-effects
# model and run the diagnostics, robustness and sensitivity analyses
# used in the thesis.
#
# Scientific question
# Within the represented 96-h Cu-Cd copepod evidence, does the
# Nauplii-versus-Adult LC50 contrast differ between Cu and Cd?
#
# Primary model
# ln(LC50) ~ Metal * Stage +
#            (1 | Reference_ID) +
#            (1 | Species)
#
# Primary estimands
# - Cu: Nauplii / Adult LC50 ratio
# - Cd: Nauplii / Adult LC50 ratio
# - Cd-versus-Cu ratio of stage ratios
#
# Input
# outputs/10_audit_stage_support/03_objects/stage_support_audit_objects.rds
#
# Main outputs
# outputs/11_fit_stage_model/
#
# Scientific workflow notes
# - The final model domain (96-h Cu-Cd Adult/Nauplii) is inherited
#   from Script 10 and is not re-selected here.
# - The Metal x Stage interaction is retained because the difference
#   in the developmental-stage contrast between Cu and Cd is a
#   prespecified estimand.
# - Reference and Species are retained as random intercepts.
# - LORO and LOSO analyses are influence checks, not alternative
#   primary models.
# - The targeted negative-residual refit is a diagnostic sensitivity,
#   not an exclusion rule.
# - The ECOTOX-only fit is retained as a reproducibility/augmentation
#   check for the effect of adding source-verified WoS records.
# - Temperature/salinity sensitivity separates complete-case
#   restriction from covariate adjustment.
# - Holm adjustment is a post-model multiplicity sensitivity for the
#   three planned developmental-stage contrasts; it does not replace
#   the prespecified primary contrasts.
# - Figures created here are analysis/QA outputs. Thesis-facing
#   polished figures are generated later by the figure-suite script.
#
# Reproducibility rule
# Run from the repository root with the RStudio project open.
# Script 10 must be completed first.
# Scientific model definitions, contrasts, robustness checks and
# numerical QA gates are preserved from the final thesis workflow.
# ============================================================


# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)

# Reporting-only helper: does not fit/select models or suppress messages.
source(file.path("scripts", "helpers", "reporting_checks.R"), local = TRUE)

options(width = 220)
options(na.action = "na.fail")

output_root <- file.path(
  "outputs",
  "11_fit_stage_model"
)

model_dir <- file.path(output_root, "01_model")
diagnostic_dir <- file.path(output_root, "02_diagnostics")
robustness_dir <- file.path(output_root, "03_robustness")
sensitivity_dir <- file.path(output_root, "04_sensitivity")
synthesis_dir <- file.path(output_root, "05_synthesis")
object_dir <- file.path(output_root, "06_objects")

completion_path <- file.path(output_root, "RUN_COMPLETE.txt")
if (file.exists(completion_path) && !file.remove(completion_path)) {
  stop("Could not remove the previous completion marker. Close open output files and retry.")
}

for (d in c(
  output_root,
  model_dir,
  diagnostic_dir,
  robustness_dir,
  sensitivity_dir,
  synthesis_dir,
  object_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# ============================================================
# 1. LOAD FINAL STAGE-MODEL DATASET FROM SCRIPT 10
# ============================================================

stage_object_file <- file.path(
  "outputs",
  "10_audit_stage_support",
  "03_objects",
  "stage_support_audit_objects.rds"
)

stopifnot(file.exists(stage_object_file))

stage_objects <- readRDS(stage_object_file)

stopifnot(
  "stage_model_data" %in% names(stage_objects)
)

stage_model_data <- stage_objects$stage_model_data %>%
  mutate(
    Metal = factor(Metal, levels = c("Cu", "Cd")),
    Stage = factor(Stage, levels = c("Adult", "Nauplii")),
    Reference_ID = factor(Reference_ID),
    Species = factor(Species)
  ) %>%
  droplevels()

required_columns <- c(
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "Species",
  "Metal",
  "Stage",
  "LC50_umol_L",
  "ln_LC50",
  "Temperature",
  "Salinity",
  "Source_Origin"
)

stopifnot(all(required_columns %in% names(stage_model_data)))

stage_support <- stage_model_data %>%
  group_by(Metal, Stage) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  )

# Frozen final four-cell support inherited from Script 10.
stopifnot(
  nrow(stage_model_data) == 70,
  
  stage_support %>%
    filter(Metal == "Cu", Stage == "Adult") %>%
    pull(Results) == 32,
  
  stage_support %>%
    filter(Metal == "Cu", Stage == "Nauplii") %>%
    pull(Results) == 14,
  
  stage_support %>%
    filter(Metal == "Cd", Stage == "Adult") %>%
    pull(Results) == 19,
  
  stage_support %>%
    filter(Metal == "Cd", Stage == "Nauplii") %>%
    pull(Results) == 5
)

stage_overall_support <- stage_model_data %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental")
  )

write_csv(
  stage_support,
  file.path(model_dir, "stage_model_dataset_support.csv")
)

write_csv(
  stage_overall_support,
  file.path(model_dir, "stage_model_dataset_overall_support.csv")
)


# ============================================================
# 2. HELPERS
# ============================================================

prepare_stage_data <- function(dat) {
  dat %>%
    mutate(
      Metal = factor(Metal, levels = c("Cu", "Cd")),
      Stage = factor(Stage, levels = c("Adult", "Nauplii")),
      Reference_ID = droplevels(factor(Reference_ID)),
      Species = droplevels(factor(Species))
    ) %>%
    droplevels()
}

fit_stage_model <- function(dat) {
  dat <- prepare_stage_data(dat)
  
  lmer(
    ln_LC50 ~
      Metal * Stage +
      (1 | Reference_ID) +
      (1 | Species),
    data = dat,
    REML = TRUE
  )
}

extract_stage_metrics <- function(model) {
  
  emm_stage <- emmeans(
    model,
    ~ Stage | Metal,
    lmer.df = "kenward-roger"
  )
  
  stage_contrasts <- contrast(
    emm_stage,
    method = list(
      "Nauplii_vs_Adult" = c(-1, 1)
    )
  )
  
  stage_df <- summary(
    stage_contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
  
  interaction_contrast <- contrast(
    stage_contrasts,
    method = list(
      "Cd_vs_Cu_stage_contrast" = c(-1, 1)
    ),
    by = NULL
  )
  
  interaction_df <- summary(
    interaction_contrast,
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
  
  cu <- stage_df %>% filter(Metal == "Cu")
  cd <- stage_df %>% filter(Metal == "Cd")
  
  conv_msg <- model@optinfo$conv$lme4$messages
  
  if (is.null(conv_msg)) {
    conv_msg <- NA_character_
  } else {
    conv_msg <- paste(conv_msg, collapse = " | ")
  }
  
  tibble(
    Cu_log_contrast = cu$estimate,
    Cu_ratio = exp(cu$estimate),
    Cu_lower95 = exp(cu$lower.CL),
    Cu_upper95 = exp(cu$upper.CL),
    Cu_p = cu$p.value,
    
    Cd_log_contrast = cd$estimate,
    Cd_ratio = exp(cd$estimate),
    Cd_lower95 = exp(cd$lower.CL),
    Cd_upper95 = exp(cd$upper.CL),
    Cd_p = cd$p.value,
    
    interaction_log = interaction_df$estimate,
    ratio_of_ratios = exp(interaction_df$estimate),
    interaction_lower95 = exp(interaction_df$lower.CL),
    interaction_upper95 = exp(interaction_df$upper.CL),
    interaction_p = interaction_df$p.value,
    
    singular = isSingular(model, tol = 1e-4),
    convergence_message = conv_msg
  ) %>%
    bind_cols(
      report_fit_diagnostics(model) %>%
        dplyr::select(-singular, -lme4_messages)
    )
}


# ============================================================
# 3. PRIMARY COMBINED INTERACTION MODEL
# ============================================================

m_stage_int <- fit_stage_model(stage_model_data)

primary_singular <- isSingular(
  m_stage_int,
  tol = 1e-4
)

primary_convergence_messages <-
  m_stage_int@optinfo$conv$lme4$messages

primary_convergence_clean <-
  is.null(primary_convergence_messages)

primary_fit_checks <- tibble(
  Check = c(
    "Fitted observations",
    "References",
    "Species",
    "Singular fit",
    "Convergence clean"
  ),
  Value = c(
    as.character(nobs(m_stage_int)),
    as.character(n_distinct(stage_model_data$Reference_ID)),
    as.character(n_distinct(stage_model_data$Species)),
    as.character(primary_singular),
    as.character(primary_convergence_clean)
  )
)

writeLines(
  capture.output(summary(m_stage_int)),
  con = file.path(
    model_dir,
    "stage_interaction_model_summary.txt"
  )
)

write_csv(
  primary_fit_checks,
  file.path(
    model_dir,
    "stage_interaction_model_fit_checks.csv"
  )
)

write_csv(
  as.data.frame(VarCorr(m_stage_int)),
  file.path(
    model_dir,
    "stage_interaction_variance_components.csv"
  )
)

if (!primary_convergence_clean) {
  writeLines(
    primary_convergence_messages,
    con = file.path(
      model_dir,
      "stage_interaction_convergence_messages.txt"
    )
  )
}

stopifnot(
  nobs(m_stage_int) == 70,
  primary_convergence_clean
)


# ============================================================
# 4. FORMAL TYPE III TESTS
# ============================================================

stage_anova <- anova(
  m_stage_int,
  type = 3,
  ddf = "Kenward-Roger"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term")

write_csv(
  stage_anova,
  file.path(
    model_dir,
    "stage_type3_kenward_roger_tests.csv"
  )
)


# ============================================================
# 5. PRIMARY ESTIMANDS
# ============================================================

emm_stage <- emmeans(
  m_stage_int,
  ~ Stage | Metal,
  lmer.df = "kenward-roger"
)

stage_contrasts <- contrast(
  emm_stage,
  method = list(
    "Nauplii_vs_Adult" = c(-1, 1)
  )
)

stage_contrasts_summary <- summary(
  stage_contrasts,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame() %>%
  mutate(
    LC50_ratio_Nauplii_Adult = exp(estimate),
    ratio_lower_95 = exp(lower.CL),
    ratio_upper_95 = exp(upper.CL)
  )

interaction_contrast <- contrast(
  stage_contrasts,
  method = list(
    "Cd_vs_Cu_stage_contrast" = c(-1, 1)
  ),
  by = NULL
)

interaction_summary <- summary(
  interaction_contrast,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame() %>%
  mutate(
    ratio_of_ratios = exp(estimate),
    ratio_lower_95 = exp(lower.CL),
    ratio_upper_95 = exp(upper.CL)
  )

stage_emmeans <- emmeans(
  m_stage_int,
  ~ Metal * Stage,
  lmer.df = "kenward-roger"
)

stage_emmeans_table <- summary(
  stage_emmeans,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

baseline_metrics <- extract_stage_metrics(m_stage_int)

# Additional numerical-status report; existing primary QA gates are unchanged.
write_csv(
  report_fit_diagnostics(m_stage_int),
  file.path(model_dir, "stage_primary_numerical_diagnostics.csv")
)

write_csv(
  stage_contrasts_summary,
  file.path(
    model_dir,
    "stage_metal_specific_contrasts.csv"
  )
)

write_csv(
  interaction_summary,
  file.path(
    model_dir,
    "stage_interaction_contrast.csv"
  )
)

write_csv(
  stage_emmeans_table,
  file.path(
    model_dir,
    "stage_estimated_marginal_means.csv"
  )
)

write_csv(
  baseline_metrics,
  file.path(
    model_dir,
    "stage_baseline_metrics.csv"
  )
)


# ============================================================
# 5B. POST-MODEL MULTIPLICITY SENSITIVITY: HOLM ADJUSTMENT
# ============================================================

stage_holm_sensitivity <- bind_rows(
  stage_contrasts_summary %>%
    transmute(
      Model = "Developmental-stage",
      Contrast = paste0(
        as.character(Metal),
        ": Nauplii / Adult"
      ),
      Effect_type = "LC50 ratio",
      Effect_estimate = LC50_ratio_Nauplii_Adult,
      lower_95 = ratio_lower_95,
      upper_95 = ratio_upper_95,
      raw_p = p.value
    ),
  
  interaction_summary %>%
    transmute(
      Model = "Developmental-stage",
      Contrast = "Cd-vs-Cu ratio of ratios",
      Effect_type = "Ratio of ratios",
      Effect_estimate = ratio_of_ratios,
      lower_95 = ratio_lower_95,
      upper_95 = ratio_upper_95,
      raw_p = p.value
    )
) %>%
  mutate(
    holm_p = p.adjust(raw_p, method = "holm"),
    raw_below_0_05 = raw_p < 0.05,
    holm_below_0_05 = holm_p < 0.05,
    interpretation_changed =
      raw_below_0_05 != holm_below_0_05
  )

write_csv(
  stage_holm_sensitivity,
  file.path(
    sensitivity_dir,
    "stage_multiplicity_holm_sensitivity.csv"
  )
)


# ============================================================
# 6. ECOTOX-ONLY REPRODUCTION / AUGMENTATION CHECK
# ============================================================
# This is not a competing primary model. It checks reproduction of
# the pre-augmentation ECOTOX-only benchmark and quantifies the effect
# of adding source-verified WoS evidence.

stage_ecotox_only <- stage_model_data %>%
  filter(Source_Origin == "ECOTOX") %>%
  droplevels()

stopifnot(nrow(stage_ecotox_only) == 61)

m_stage_ecotox <- fit_stage_model(stage_ecotox_only)

ecotox_metrics <- extract_stage_metrics(m_stage_ecotox)

# Preserve the OLD 60-row numerical regression test on the OLD membership.
# It is not a target for the corrected 61-row ECOTOX-only model above.
# C6 was absent from the original Stage fit because of its former Copepodite label.
stage_ecotox_legacy <- stage_ecotox_only %>%
  filter(Result_ID != "E_R115201") %>% droplevels()
stopifnot(nrow(stage_ecotox_legacy) == 60)
legacy_ecotox_metrics <- extract_stage_metrics(fit_stage_model(stage_ecotox_legacy))
# Frozen historical values; they do not define the corrected combined result.
ecotox_regression_check <- tibble(
  Estimand = c(
    "Cu Nauplii / Adult",
    "Cd Nauplii / Adult",
    "Cd-vs-Cu ratio of ratios"
  ),
  Expected = c(
    0.3625,
    0.6745,
    1.8608
  ),
  ECOTOX_reproduction = c(
    legacy_ecotox_metrics$Cu_ratio,
    legacy_ecotox_metrics$Cd_ratio,
    legacy_ecotox_metrics$ratio_of_ratios
  )
) %>%
  mutate(
    Absolute_difference =
      abs(ECOTOX_reproduction - Expected),
    Match = Absolute_difference <
      c(0.02, 0.03, 0.08)
  )

stopifnot(all(ecotox_regression_check$Match))

ecotox_vs_combined <- tibble(
  Estimand = c(
    "Cu Nauplii / Adult",
    "Cd Nauplii / Adult",
    "Cd-vs-Cu ratio of ratios"
  ),
  Legacy_ECOTOX_benchmark = c(
    0.3625,
    0.6745,
    1.8608
  ),
  ECOTOX_reproduction = c(
    ecotox_metrics$Cu_ratio,
    ecotox_metrics$Cd_ratio,
    ecotox_metrics$ratio_of_ratios
  ),
  Combined = c(
    baseline_metrics$Cu_ratio,
    baseline_metrics$Cd_ratio,
    baseline_metrics$ratio_of_ratios
  )
) %>%
  mutate(
    Combined_vs_ECOTOX_percent_shift =
      100 * (Combined / ECOTOX_reproduction - 1),
    Direction_same_as_ECOTOX = case_when(
      Estimand == "Cd-vs-Cu ratio of ratios" ~
        (Combined > 1) == (ECOTOX_reproduction > 1),
      TRUE ~
        (Combined < 1) == (ECOTOX_reproduction < 1)
    )
  )

write_csv(
  ecotox_regression_check,
  file.path(
    synthesis_dir,
    "ecotox_only_stage_regression_check.csv"
  )
)

write_csv(
  ecotox_vs_combined,
  file.path(
    synthesis_dir,
    "ecotox_vs_combined_stage_shift.csv"
  )
)

saveRDS(
  m_stage_ecotox,
  file.path(
    sensitivity_dir,
    "stage_ECOTOX_only_model.rds"
  )
)


# ============================================================
# 7. INTERMEDIATE / QA STAGE EFFECT-RATIO FIGURE
# ============================================================

stage_effect_plot_data <- tibble(
  Estimand = factor(
    c(
      "Cu: Nauplii / Adult",
      "Cd: Nauplii / Adult",
      "Cd / Cu stage-contrast ratio"
    ),
    levels = rev(
      c(
        "Cu: Nauplii / Adult",
        "Cd: Nauplii / Adult",
        "Cd / Cu stage-contrast ratio"
      )
    )
  ),
  Ratio = c(
    baseline_metrics$Cu_ratio,
    baseline_metrics$Cd_ratio,
    baseline_metrics$ratio_of_ratios
  ),
  Lower = c(
    baseline_metrics$Cu_lower95,
    baseline_metrics$Cd_lower95,
    baseline_metrics$interaction_lower95
  ),
  Upper = c(
    baseline_metrics$Cu_upper95,
    baseline_metrics$Cd_upper95,
    baseline_metrics$interaction_upper95
  )
)

stage_effect_plot <- ggplot(
  stage_effect_plot_data,
  aes(
    x = Ratio,
    y = Estimand
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = Lower,
      xmax = Upper
    ),
    orientation = "y",
    width = 0.18
  ) +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(
    x = "LC50 ratio (log scale)",
    y = NULL
  ) +
  theme_classic()

ggsave(
  file.path(
    model_dir,
    "stage_interaction_effect_ratios.png"
  ),
  stage_effect_plot,
  width = 8.5,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    model_dir,
    "stage_interaction_effect_ratios.pdf"
  ),
  stage_effect_plot,
  width = 8.5,
  height = 4.8
)


# ============================================================
# 8. DIAGNOSTICS
# ============================================================

stage_diagnostics <- stage_model_data %>%
  mutate(
    fitted_value = fitted(m_stage_int),
    residual = residuals(m_stage_int),
    scaled_residual = residuals(
      m_stage_int,
      type = "pearson",
      scaled = TRUE
    ),
    abs_scaled_residual = abs(scaled_residual)
  )

top_residuals <- stage_diagnostics %>%
  arrange(desc(abs_scaled_residual)) %>%
  select(
    any_of(
      c(
        "Reference_ID",
        "Test_ID",
        "Result_ID",
        "Source_Origin",
        "Species",
        "Metal",
        "Stage",
        "LC50_umol_L",
        "ln_LC50",
        "fitted_value",
        "residual",
        "scaled_residual",
        "abs_scaled_residual"
      )
    )
  ) %>%
  slice_head(n = 10)

write_csv(
  stage_diagnostics,
  file.path(
    diagnostic_dir,
    "stage_observation_diagnostics.csv"
  )
)

write_csv(
  top_residuals,
  file.path(
    diagnostic_dir,
    "stage_top_10_scaled_residuals.csv"
  )
)

png(
  file.path(
    diagnostic_dir,
    "stage_residuals_vs_fitted.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  fitted(m_stage_int),
  resid(m_stage_int),
  xlab = "Fitted ln(LC50)",
  ylab = "Residual",
  main = "Residuals versus fitted values"
)

abline(h = 0, lty = 2)
dev.off()

png(
  file.path(
    diagnostic_dir,
    "stage_residuals_normal_qq.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

qqnorm(
  resid(m_stage_int),
  main = "Normal Q-Q plot: residuals"
)

qqline(resid(m_stage_int))
dev.off()

ref_re <- ranef(m_stage_int)$Reference_ID[, 1]

png(
  file.path(
    diagnostic_dir,
    "stage_reference_random_intercepts_qq.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

qqnorm(
  ref_re,
  main = "Normal Q-Q plot: Reference random intercepts"
)

qqline(ref_re)
dev.off()

species_re <- ranef(m_stage_int)$Species[, 1]

png(
  file.path(
    diagnostic_dir,
    "stage_species_random_intercepts_qq.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

qqnorm(
  species_re,
  main = "Normal Q-Q plot: Species random intercepts"
)

qqline(species_re)
dev.off()


# ============================================================
# 9. ROBUSTNESS: LEAVE-ONE-REFERENCE-OUT
# ============================================================

refs <- unique(
  as.character(stage_model_data$Reference_ID)
)

loro_results <- purrr::map_dfr(
  refs,
  function(ref_out) {
    
    dat_i <- stage_model_data %>%
      filter(
        as.character(Reference_ID) != ref_out
      )
    
    tryCatch({
      
      mod_i <- fit_stage_model(dat_i)
      metrics_i <- extract_stage_metrics(mod_i)
      
      metrics_i %>%
        mutate(
          removed_reference = ref_out,
          n_results = nrow(dat_i),
          n_references =
            n_distinct(dat_i$Reference_ID),
          n_species =
            n_distinct(dat_i$Species),
          error = NA_character_
        )
      
    }, error = function(e) {
      
      tibble(
        removed_reference = ref_out,
        error = conditionMessage(e)
      )
    })
  }
)

loro_summary <- loro_results %>%
  filter(!is.na(Cu_ratio)) %>%
  summarise(
    n_successful_refits = n(),
    
    Cu_ratio_min = min(Cu_ratio),
    Cu_ratio_median = median(Cu_ratio),
    Cu_ratio_max = max(Cu_ratio),
    
    Cd_ratio_min = min(Cd_ratio),
    Cd_ratio_median = median(Cd_ratio),
    Cd_ratio_max = max(Cd_ratio),
    
    RoR_min = min(ratio_of_ratios),
    RoR_median = median(ratio_of_ratios),
    RoR_max = max(ratio_of_ratios),
    
    Cu_Nauplii_lower_n =
      sum(Cu_ratio < 1),
    
    Cd_Nauplii_lower_n =
      sum(Cd_ratio < 1),
    
    interaction_RoR_above1_n =
      sum(ratio_of_ratios > 1),
    
    singular_refits =
      sum(singular, na.rm = TRUE),
    
    # Same raw message may report singularity; do not count it again
    # as a separate optimizer/convergence failure.
    any_lme4_message_refits =
      sum(!is.na(convergence_message)),
    singularity_message_refits =
      sum(!is.na(singularity_message)),
    other_numerical_alert_refits =
      sum(other_numerical_alert, na.rm = TRUE),
    optimizer_nonzero_code_refits =
      sum(optimizer_nonzero_code, na.rm = TRUE),
    optimizer_warning_refits =
      sum(!is.na(optimizer_warning)),
    optimizer_code_uninterpretable_refits =
      sum(optimizer_code_uninterpretable, na.rm = TRUE)
  )

write_csv(
  loro_results,
  file.path(
    robustness_dir,
    "stage_LORO_all_refits.csv"
  )
)

write_csv(
  loro_summary,
  file.path(
    robustness_dir,
    "stage_LORO_summary.csv"
  )
)

stopifnot(
  loro_summary$n_successful_refits ==
    n_distinct(stage_model_data$Reference_ID)
)


# ============================================================
# 10. ROBUSTNESS: LEAVE-ONE-SPECIES-OUT
# ============================================================

species_list <- unique(
  as.character(stage_model_data$Species)
)

loso_results <- purrr::map_dfr(
  species_list,
  function(sp_out) {
    
    dat_i <- stage_model_data %>%
      filter(
        as.character(Species) != sp_out
      )
    
    tryCatch({
      
      mod_i <- fit_stage_model(dat_i)
      metrics_i <- extract_stage_metrics(mod_i)
      
      metrics_i %>%
        mutate(
          removed_species = sp_out,
          n_results = nrow(dat_i),
          n_references =
            n_distinct(dat_i$Reference_ID),
          n_species =
            n_distinct(dat_i$Species),
          error = NA_character_
        )
      
    }, error = function(e) {
      
      tibble(
        removed_species = sp_out,
        error = conditionMessage(e)
      )
    })
  }
)

loso_summary <- loso_results %>%
  filter(!is.na(Cu_ratio)) %>%
  summarise(
    n_successful_refits = n(),
    
    Cu_ratio_min = min(Cu_ratio),
    Cu_ratio_median = median(Cu_ratio),
    Cu_ratio_max = max(Cu_ratio),
    
    Cd_ratio_min = min(Cd_ratio),
    Cd_ratio_median = median(Cd_ratio),
    Cd_ratio_max = max(Cd_ratio),
    
    RoR_min = min(ratio_of_ratios),
    RoR_median = median(ratio_of_ratios),
    RoR_max = max(ratio_of_ratios),
    
    Cu_Nauplii_lower_n =
      sum(Cu_ratio < 1),
    
    Cd_Nauplii_lower_n =
      sum(Cd_ratio < 1),
    
    interaction_RoR_above1_n =
      sum(ratio_of_ratios > 1),
    
    singular_refits =
      sum(singular, na.rm = TRUE),
    
    # Same raw message may report singularity; do not count it again
    # as a separate optimizer/convergence failure.
    any_lme4_message_refits =
      sum(!is.na(convergence_message)),
    singularity_message_refits =
      sum(!is.na(singularity_message)),
    other_numerical_alert_refits =
      sum(other_numerical_alert, na.rm = TRUE),
    optimizer_nonzero_code_refits =
      sum(optimizer_nonzero_code, na.rm = TRUE),
    optimizer_warning_refits =
      sum(!is.na(optimizer_warning)),
    optimizer_code_uninterpretable_refits =
      sum(optimizer_code_uninterpretable, na.rm = TRUE)
  )

write_csv(
  loso_results,
  file.path(
    robustness_dir,
    "stage_LOSO_all_refits.csv"
  )
)

write_csv(
  loso_summary,
  file.path(
    robustness_dir,
    "stage_LOSO_summary.csv"
  )
)

stopifnot(
  loso_summary$n_successful_refits ==
    n_distinct(stage_model_data$Species)
)


# ============================================================
# 11. TARGETED NEGATIVE-TAIL OBSERVATION SENSITIVITY
# ============================================================
# Preserved from the final developmental-stage workflow.
# This is a one-time diagnostic sensitivity, NOT an exclusion rule.

worst_negative_index <-
  which.min(stage_diagnostics$scaled_residual)

worst_negative_observation <-
  stage_diagnostics[worst_negative_index, ]

stage_without_negative <-
  stage_model_data[-worst_negative_index, ]

m_stage_no_negative <-
  fit_stage_model(stage_without_negative)

negative_tail_metrics <-
  extract_stage_metrics(m_stage_no_negative)

baseline_vs_negative <- bind_rows(
  baseline_metrics %>%
    mutate(model = "BASELINE"),
  negative_tail_metrics %>%
    mutate(
      model =
        "REMOVE_NEGATIVE_TAIL_OBSERVATION"
    )
) %>%
  select(
    model,
    Cu_ratio,
    Cu_lower95,
    Cu_upper95,
    Cd_ratio,
    Cd_lower95,
    Cd_upper95,
    ratio_of_ratios,
    interaction_lower95,
    interaction_upper95,
    singular,
    convergence_message
  )

write_csv(
  worst_negative_observation,
  file.path(
    robustness_dir,
    "stage_negative_tail_observation.csv"
  )
)

write_csv(
  baseline_vs_negative,
  file.path(
    robustness_dir,
    "stage_baseline_vs_negative_tail.csv"
  )
)

saveRDS(
  m_stage_no_negative,
  file.path(
    robustness_dir,
    "stage_model_without_negative_tail_observation.rds"
  )
)


# ============================================================
# 12. ENVIRONMENTAL COMPLETE-CASE SUPPORT
# ============================================================
# Temperature and salinity sensitivity is retained because those fields
# are available in the combined schema. This is a robustness analysis
# for the Stage model, not the separate environmental RQ audit.
# pH is handled separately in Script 08 and is not added to this sensitivity model.

stage_env_cc <- stage_model_data %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  droplevels()

stage_env_cc_overall <- stage_env_cc %>%
  summarise(
    Results = n(),
    References =
      n_distinct(Reference_ID),
    Species =
      n_distinct(Species)
  )

stage_env_cc_support <- stage_env_cc %>%
  group_by(Metal, Stage) %>%
  summarise(
    Results = n(),
    References =
      n_distinct(Reference_ID),
    Species =
      n_distinct(Species),
    .groups = "drop"
  )

write_csv(
  stage_env_cc_overall,
  file.path(
    sensitivity_dir,
    "stage_environment_complete_case_overall.csv"
  )
)

write_csv(
  stage_env_cc_support,
  file.path(
    sensitivity_dir,
    "stage_environment_complete_case_support.csv"
  )
)


# ============================================================
# 13. ENVIRONMENTAL SENSITIVITY MODELS
# ============================================================

environmental_models_status <- tibble(
  Analysis = c(
    "Complete-case same formula",
    "Complete-case + Temperature/Salinity"
  ),
  Status = "NOT_RUN"
)

stage_environment_comparison <- tibble()
cc_fit_checks <- tibble()
temp_sal_correlation <- tibble()

# Require all four cells and enough data before fitting.
cc_cell_support <- stage_env_cc %>%
  count(Metal, Stage, name = "n")

can_fit_env <- (
  nrow(cc_cell_support) == 4 &&
    all(cc_cell_support$n >= 3) &&
    n_distinct(stage_env_cc$Reference_ID) >= 8 &&
    n_distinct(stage_env_cc$Species) >= 6 &&
    sd(stage_env_cc$Temperature, na.rm = TRUE) > 0 &&
    sd(stage_env_cc$Salinity, na.rm = TRUE) > 0
)

if (can_fit_env) {
  
  m_stage_cc <- fit_stage_model(stage_env_cc)
  
  stage_env_cc_scaled <- stage_env_cc %>%
    mutate(
      Temperature_z =
        as.numeric(scale(Temperature)),
      Salinity_z =
        as.numeric(scale(Salinity))
    )
  
  m_stage_env <- lmer(
    ln_LC50 ~
      Metal * Stage +
      Temperature_z +
      Salinity_z +
      (1 | Reference_ID) +
      (1 | Species),
    data = stage_env_cc_scaled,
    REML = TRUE
  )
  
  cc_fit_checks <- tibble(
    Model = c(
      "Complete-case same formula",
      "Complete-case + Temperature/Salinity"
    ),
    Singular = c(
      isSingular(m_stage_cc, tol = 1e-4),
      isSingular(m_stage_env, tol = 1e-4)
    ),
    Convergence_clean = c(
      is.null(
        m_stage_cc@optinfo$conv$lme4$messages
      ),
      is.null(
        m_stage_env@optinfo$conv$lme4$messages
      )
    )
  )
  
  baseline_env_metrics <-
    extract_stage_metrics(m_stage_int) %>%
    mutate(model = "PRIMARY_FULL_DATA")
  
  cc_metrics <-
    extract_stage_metrics(m_stage_cc) %>%
    mutate(
      model =
        "TEMP_SAL_COMPLETE_CASE_SAME_FORMULA"
    )
  
  adjusted_metrics <-
    extract_stage_metrics(m_stage_env) %>%
    mutate(
      model =
        "TEMP_SAL_COMPLETE_CASE_ADJUSTED"
    )
  
  stage_environment_comparison <- bind_rows(
    baseline_env_metrics,
    cc_metrics,
    adjusted_metrics
  ) %>%
    select(
      model,
      Cu_ratio,
      Cu_lower95,
      Cu_upper95,
      Cd_ratio,
      Cd_lower95,
      Cd_upper95,
      ratio_of_ratios,
      interaction_lower95,
      interaction_upper95,
      singular,
      convergence_message
    )
  
  temp_sal_correlation <- tibble(
    Statistic = "Pearson correlation",
    Value = cor(
      stage_env_cc$Temperature,
      stage_env_cc$Salinity,
      use = "complete.obs"
    )
  )
  
  environmental_models_status$Status <- "RUN"
  
  writeLines(
    capture.output(summary(m_stage_cc)),
    con = file.path(
      sensitivity_dir,
      "stage_complete_case_same_formula_model_summary.txt"
    )
  )
  
  writeLines(
    capture.output(summary(m_stage_env)),
    con = file.path(
      sensitivity_dir,
      "stage_environment_adjusted_model_summary.txt"
    )
  )
  
  write_csv(
    cc_fit_checks,
    file.path(
      sensitivity_dir,
      "stage_environment_model_fit_checks.csv"
    )
  )
  
  write_csv(
    stage_environment_comparison,
    file.path(
      sensitivity_dir,
      "stage_environmental_sensitivity_comparison.csv"
    )
  )
  
  write_csv(
    temp_sal_correlation,
    file.path(
      sensitivity_dir,
      "stage_temperature_salinity_correlation.csv"
    )
  )
  
  saveRDS(
    m_stage_cc,
    file.path(
      sensitivity_dir,
      "stage_complete_case_same_formula_model.rds"
    )
  )
  
  saveRDS(
    m_stage_env,
    file.path(
      sensitivity_dir,
      "stage_complete_case_environment_adjusted_model.rds"
    )
  )
  
} else {
  
  environmental_models_status$Status <-
    "NOT RUN: insufficient four-cell/common-support structure"
  
  write_csv(
    environmental_models_status,
    file.path(
      sensitivity_dir,
      "stage_environmental_sensitivity_status.csv"
    )
  )
}

if (can_fit_env) {
  write_csv(
    environmental_models_status,
    file.path(
      sensitivity_dir,
      "stage_environmental_sensitivity_status.csv"
    )
  )
}


# ============================================================
# 14. ROBUSTNESS SYNTHESIS
# ============================================================

stage_robustness_synthesis <- tibble(
  Analysis = c(
    "Combined baseline",
    "ECOTOX-only reproduction",
    "LORO range",
    "LOSO range"
  ),
  Cu_Nauplii_Adult = c(
    baseline_metrics$Cu_ratio,
    ecotox_metrics$Cu_ratio,
    NA_real_,
    NA_real_
  ),
  Cu_range = c(
    paste0(
      round(baseline_metrics$Cu_lower95, 3),
      " - ",
      round(baseline_metrics$Cu_upper95, 3)
    ),
    paste0(
      round(ecotox_metrics$Cu_lower95, 3),
      " - ",
      round(ecotox_metrics$Cu_upper95, 3)
    ),
    paste0(
      round(loro_summary$Cu_ratio_min, 3),
      " - ",
      round(loro_summary$Cu_ratio_max, 3)
    ),
    paste0(
      round(loso_summary$Cu_ratio_min, 3),
      " - ",
      round(loso_summary$Cu_ratio_max, 3)
    )
  ),
  Cd_Nauplii_Adult = c(
    baseline_metrics$Cd_ratio,
    ecotox_metrics$Cd_ratio,
    NA_real_,
    NA_real_
  ),
  Cd_range = c(
    paste0(
      round(baseline_metrics$Cd_lower95, 3),
      " - ",
      round(baseline_metrics$Cd_upper95, 3)
    ),
    paste0(
      round(ecotox_metrics$Cd_lower95, 3),
      " - ",
      round(ecotox_metrics$Cd_upper95, 3)
    ),
    paste0(
      round(loro_summary$Cd_ratio_min, 3),
      " - ",
      round(loro_summary$Cd_ratio_max, 3)
    ),
    paste0(
      round(loso_summary$Cd_ratio_min, 3),
      " - ",
      round(loso_summary$Cd_ratio_max, 3)
    )
  ),
  Cd_vs_Cu_RoR = c(
    baseline_metrics$ratio_of_ratios,
    ecotox_metrics$ratio_of_ratios,
    NA_real_,
    NA_real_
  ),
  RoR_range = c(
    paste0(
      round(baseline_metrics$interaction_lower95, 3),
      " - ",
      round(baseline_metrics$interaction_upper95, 3)
    ),
    paste0(
      round(ecotox_metrics$interaction_lower95, 3),
      " - ",
      round(ecotox_metrics$interaction_upper95, 3)
    ),
    paste0(
      round(loro_summary$RoR_min, 3),
      " - ",
      round(loro_summary$RoR_max, 3)
    ),
    paste0(
      round(loso_summary$RoR_min, 3),
      " - ",
      round(loso_summary$RoR_max, 3)
    )
  )
)

if (can_fit_env) {
  stage_robustness_synthesis <- bind_rows(
    stage_robustness_synthesis,
    tibble(
      Analysis = c(
        "Temperature/salinity complete-case",
        "Temperature/salinity adjusted"
      ),
      Cu_Nauplii_Adult = c(
        cc_metrics$Cu_ratio,
        adjusted_metrics$Cu_ratio
      ),
      Cu_range = c(
        paste0(
          round(cc_metrics$Cu_lower95, 3),
          " - ",
          round(cc_metrics$Cu_upper95, 3)
        ),
        paste0(
          round(adjusted_metrics$Cu_lower95, 3),
          " - ",
          round(adjusted_metrics$Cu_upper95, 3)
        )
      ),
      Cd_Nauplii_Adult = c(
        cc_metrics$Cd_ratio,
        adjusted_metrics$Cd_ratio
      ),
      Cd_range = c(
        paste0(
          round(cc_metrics$Cd_lower95, 3),
          " - ",
          round(cc_metrics$Cd_upper95, 3)
        ),
        paste0(
          round(adjusted_metrics$Cd_lower95, 3),
          " - ",
          round(adjusted_metrics$Cd_upper95, 3)
        )
      ),
      Cd_vs_Cu_RoR = c(
        cc_metrics$ratio_of_ratios,
        adjusted_metrics$ratio_of_ratios
      ),
      RoR_range = c(
        paste0(
          round(cc_metrics$interaction_lower95, 3),
          " - ",
          round(cc_metrics$interaction_upper95, 3)
        ),
        paste0(
          round(adjusted_metrics$interaction_lower95, 3),
          " - ",
          round(adjusted_metrics$interaction_upper95, 3)
        )
      )
    )
  )
}

write_csv(
  stage_robustness_synthesis,
  file.path(
    robustness_dir,
    "stage_robustness_synthesis.csv"
  )
)


# ============================================================
# 15. CONCLUSION-STABILITY SNAPSHOT
# ============================================================

conclusion_stability <- tibble(
  Question = c(
    "Cu: are Nauplii point estimates lower than Adults?",
    "Cd: are Nauplii point estimates lower than Adults?",
    "Is the unadjusted Cd-vs-Cu interaction p-value below 0.05?"
  ),
  Criterion = c(
    "POINT_ESTIMATE_DIRECTION_ONLY",
    "POINT_ESTIMATE_DIRECTION_ONLY",
    "UNADJUSTED_INTERACTION_P_LT_0_05"
  ),
  Interpretation_note = c(
    "TRUE describes point direction, not statistical support.",
    "TRUE describes point direction, not statistical support.",
    "Unadjusted planned interaction test; Holm sensitivity is reported separately."
  ),
  Combined_primary = c(
    baseline_metrics$Cu_ratio < 1,
    baseline_metrics$Cd_ratio < 1,
    baseline_metrics$interaction_p < 0.05
  ),
  ECOTOX_only = c(
    ecotox_metrics$Cu_ratio < 1,
    ecotox_metrics$Cd_ratio < 1,
    ecotox_metrics$interaction_p < 0.05
  ),
  LORO_direction_consistency = c(
    loro_summary$Cu_Nauplii_lower_n ==
      loro_summary$n_successful_refits,
    loro_summary$Cd_Nauplii_lower_n ==
      loro_summary$n_successful_refits,
    NA
  ),
  LOSO_direction_consistency = c(
    loso_summary$Cu_Nauplii_lower_n ==
      loso_summary$n_successful_refits,
    loso_summary$Cd_Nauplii_lower_n ==
      loso_summary$n_successful_refits,
    NA
  )
)

write_csv(
  conclusion_stability,
  file.path(
    synthesis_dir,
    "stage_conclusion_stability_summary.csv"
  )
)


# ============================================================
# 16. SAVE OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    stage_model_data = stage_model_data,
    primary_model = m_stage_int,
    stage_anova = stage_anova,
    stage_emmeans = stage_emmeans,
    stage_contrasts = stage_contrasts_summary,
    interaction_contrast = interaction_summary,
    baseline_metrics = baseline_metrics,
    ECOTOX_only_model = m_stage_ecotox,
    ECOTOX_only_metrics = ecotox_metrics,
    ecotox_vs_combined = ecotox_vs_combined,
    stage_effect_plot = stage_effect_plot,
    diagnostics = stage_diagnostics,
    loro_results = loro_results,
    loso_results = loso_results,
    negative_tail_model = m_stage_no_negative,
    environmental_models_run = can_fit_env,
    environmental_comparison = stage_environment_comparison
  ),
  file.path(
    object_dir,
    "stage_model_analysis_objects.rds"
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
    "Input: stage support audit objects",
    "Output: primary model contrasts",
    "Output: LORO summary",
    "Output: LOSO summary",
    "Output: environmental sensitivity",
    "Output: reusable analysis objects"
  ),
  Path = c(
    stage_object_file,
    file.path(model_dir, "stage_metal_specific_contrasts.csv"),
    file.path(robustness_dir, "stage_LORO_summary.csv"),
    file.path(robustness_dir, "stage_LOSO_summary.csv"),
    file.path(sensitivity_dir, "stage_environmental_sensitivity_comparison.csv"),
    file.path(object_dir, "stage_model_analysis_objects.rds")
  )
)

write_csv(
  input_output_manifest,
  file.path(output_root, "input_output_manifest.csv")
)


# ============================================================
# 17. FINAL OUTPUT QA
# ============================================================

required_outputs <- c(
  file.path(model_dir, "stage_primary_numerical_diagnostics.csv"),
  file.path(
    object_dir,
    "stage_model_analysis_objects.rds"
  ),
  file.path(
    model_dir,
    "stage_interaction_effect_ratios.png"
  ),
  file.path(
    model_dir,
    "stage_baseline_metrics.csv"
  ),
  file.path(
    model_dir,
    "stage_type3_kenward_roger_tests.csv"
  ),
  file.path(
    diagnostic_dir,
    "stage_residuals_vs_fitted.png"
  ),
  file.path(
    diagnostic_dir,
    "stage_residuals_normal_qq.png"
  ),
  file.path(
    robustness_dir,
    "stage_LORO_summary.csv"
  ),
  file.path(
    robustness_dir,
    "stage_LOSO_summary.csv"
  ),
  file.path(
    sensitivity_dir,
    "stage_multiplicity_holm_sensitivity.csv"
  ),
  file.path(
    synthesis_dir,
    "ecotox_vs_combined_stage_shift.csv"
  ),
  file.path(
    synthesis_dir,
    "stage_conclusion_stability_summary.csv"
  )
)

stopifnot(all(file.exists(required_outputs)))


# ============================================================
# 18. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("11 DEVELOPMENTAL-STAGE MODEL - FINAL SUMMARY\n")
cat("============================================================\n")

cat("\n--- COMBINED MODEL SUPPORT ---\n")
print(stage_support, n = Inf, width = Inf)
print(stage_overall_support, width = Inf)

cat("\n--- PRIMARY MODEL FIT ---\n")
print(primary_fit_checks, n = Inf)

cat("\n--- TYPE III TESTS ---\n")
print(
  tibble::as_tibble(stage_anova),
  n = Inf,
  width = Inf
)

cat("\n--- PRIMARY STAGE RATIOS ---\n")
print(baseline_metrics, width = Inf)

cat("\n--- ECOTOX-ONLY / COMBINED SHIFT ---\n")
print(ecotox_vs_combined, n = Inf, width = Inf)

cat("\n--- LORO SUMMARY ---\n")
print(loro_summary, width = Inf)

cat("\n--- LOSO SUMMARY ---\n")
print(loso_summary, width = Inf)

cat("\n--- TEMPERATURE/SALINITY SENSITIVITY STATUS ---\n")
print(environmental_models_status, n = Inf, width = Inf)

if (can_fit_env) {
  cat("\n--- TEMPERATURE/SALINITY SENSITIVITY ---\n")
  print(stage_environment_comparison, n = Inf, width = Inf)
}

cat("\n--- CONCLUSION STABILITY ---\n")
print(conclusion_stability, n = Inf, width = Inf)

cat(
  "\nInterpretive rule:\n",
  "The combined 96-h Cu-Cd Adult/Nauplii model is the primary developmental-stage analysis.\n",
  "The ECOTOX-only model is a reproducibility/augmentation check, not a competing primary model.\n",
  "The model domain is inherited from Script 10 and is not re-selected in this script.\n",
  sep = ""
)

writeLines(
  c(
    "STATUS: 11 FIT STAGE MODEL = PASS",
    "Primary domain: 96-h Cu-Cd Adult/Nauplii",
    paste0("Results: ", nrow(stage_model_data)),
    paste0("References: ", n_distinct(stage_model_data$Reference_ID)),
    paste0("Species: ", n_distinct(stage_model_data$Species))
  ),
  file.path(output_root, "RUN_COMPLETE.txt")
)

cat(
  "\nSTATUS: 11 FIT STAGE MODEL = PASS\n"
)
# ============================================================
# END
# ============================================================
