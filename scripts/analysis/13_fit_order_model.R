# ============================================================
# 13_fit_order_model.R
# ============================================================
# Purpose
# Fit the final 96-h Cu-Cd taxonomic-order linear mixed-effects
# model and run the associated diagnostics, robustness and
# sensitivity analyses used in the thesis.
#
# Scientific question
# Within the represented 96-h Cu-Cd copepod evidence, does the
# Harpacticoida-versus-Calanoida LC50 contrast differ between Cu and Cd?
#
# Primary model
# ln(LC50) ~ Metal * Order +
#            (1 | Reference_ID) +
#            (1 | Species)
#
# Primary estimands
# - Cu: Harpacticoida / Calanoida LC50 ratio
# - Cd: Harpacticoida / Calanoida LC50 ratio
# - Cd-versus-Cu ratio of Order ratios
#
# Input
# outputs/12_audit_order_support/03_objects/order_support_audit_objects.rds
#
# Main outputs
# outputs/13_fit_order_model/
#
# Scientific workflow notes
# - The final model domain (96-h Cu-Cd; Calanoida/Harpacticoida)
#   is inherited from Script 12 and is not re-selected here.
# - The Metal x Order interaction is retained because the difference
#   in the Order contrast between Cu and Cd is a prespecified estimand.
# - Reference and Species are retained as random intercepts.
# - LORO and LOSO analyses are influence checks, not alternative
#   primary models.
# - The ECOTOX-only fit is retained as a reproducibility/augmentation
#   check for the effect of adding source-verified WoS records.
# - Temperature/salinity sensitivity separates complete-case
#   restriction from covariate adjustment.
# - Holm adjustment is a post-model multiplicity sensitivity for the
#   three planned Order contrasts; it does not replace the
#   prespecified primary contrasts.
# - Figures created here are analysis/QA outputs. Thesis-facing
#   polished figures are generated later by the figure-suite script.
#
# Reproducibility rule
# Run from the repository root with the RStudio project open.
# Script 12 must be completed first.
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

options(width = 220)
options(na.action = "na.fail")

output_root <- file.path(
  "outputs",
  "13_fit_order_model"
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
# 1. LOAD FINAL ORDER-MODEL DATASET FROM SCRIPT 12
# ============================================================

audit_rds <- file.path(
  "outputs",
  "12_audit_order_support",
  "03_objects",
  "order_support_audit_objects.rds"
)

stopifnot(file.exists(audit_rds))

order_objects <- readRDS(audit_rds)

stopifnot(
  "order_model_data" %in% names(order_objects)
)

order_model_data <- order_objects$order_model_data %>%
  mutate(
    Metal = factor(Metal, levels = c("Cu", "Cd")),
    Order = factor(Order, levels = c("Calanoida", "Harpacticoida")),
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
  "Order",
  "LC50_umol_L",
  "ln_LC50",
  "Temperature",
  "Salinity",
  "Source_Origin"
)

stopifnot(all(required_columns %in% names(order_model_data)))

order_support <- order_model_data %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  )

# Frozen final four-cell support inherited from Script 12.
stopifnot(
  nrow(order_model_data) == 101,
  
  order_support %>%
    filter(Metal == "Cu", Order == "Calanoida") %>%
    pull(Results) == 23,
  
  order_support %>%
    filter(Metal == "Cu", Order == "Harpacticoida") %>%
    pull(Results) == 43,
  
  order_support %>%
    filter(Metal == "Cd", Order == "Calanoida") %>%
    pull(Results) == 24,
  
  order_support %>%
    filter(Metal == "Cd", Order == "Harpacticoida") %>%
    pull(Results) == 11
)

order_overall_support <- order_model_data %>%
  summarise(
    Results = n(),
    Tests = n_distinct(Test_ID),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental")
  )

write_csv(
  order_support,
  file.path(model_dir, "order_model_dataset_support.csv")
)

write_csv(
  order_overall_support,
  file.path(model_dir, "order_model_dataset_overall_support.csv")
)


# ============================================================
# 2. HELPERS
# ============================================================

prepare_order_data <- function(dat) {
  dat %>%
    mutate(
      Metal = factor(Metal, levels = c("Cu", "Cd")),
      Order = factor(Order, levels = c("Calanoida", "Harpacticoida")),
      Reference_ID = droplevels(factor(Reference_ID)),
      Species = droplevels(factor(Species))
    ) %>%
    droplevels()
}

fit_order_model <- function(dat) {
  
  dat <- prepare_order_data(dat)
  
  cell_check <- table(dat$Metal, dat$Order)
  
  if (any(cell_check == 0)) {
    stop("At least one Metal x Order cell became empty.")
  }
  
  lmer(
    ln_LC50 ~
      Metal * Order +
      (1 | Reference_ID) +
      (1 | Species),
    data = dat,
    REML = TRUE
  )
}

extract_order_metrics <- function(model) {
  
  emm_order <- emmeans(
    model,
    ~ Order | Metal,
    lmer.df = "kenward-roger"
  )
  
  order_contrasts <- contrast(
    emm_order,
    method = list(
      "Harpacticoida_vs_Calanoida" = c(-1, 1)
    )
  )
  
  order_df <- summary(
    order_contrasts,
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
  
  interaction_contrast <- contrast(
    order_contrasts,
    method = list(
      "Cd_vs_Cu_order_contrast" = c(-1, 1)
    ),
    by = NULL
  )
  
  interaction_df <- summary(
    interaction_contrast,
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
  
  cu <- order_df %>% filter(Metal == "Cu")
  cd <- order_df %>% filter(Metal == "Cd")
  
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
  )
}


# ============================================================
# 3. PRIMARY COMBINED INTERACTION MODEL
# ============================================================

m_order_int <- fit_order_model(order_model_data)

primary_singular <- isSingular(
  m_order_int,
  tol = 1e-4
)

primary_convergence_messages <-
  m_order_int@optinfo$conv$lme4$messages

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
    as.character(nobs(m_order_int)),
    as.character(n_distinct(order_model_data$Reference_ID)),
    as.character(n_distinct(order_model_data$Species)),
    as.character(primary_singular),
    as.character(primary_convergence_clean)
  )
)

writeLines(
  capture.output(summary(m_order_int)),
  con = file.path(
    model_dir,
    "order_interaction_model_summary.txt"
  )
)

write_csv(
  primary_fit_checks,
  file.path(
    model_dir,
    "order_interaction_model_fit_checks.csv"
  )
)

write_csv(
  as.data.frame(VarCorr(m_order_int)),
  file.path(
    model_dir,
    "order_interaction_variance_components.csv"
  )
)

if (!primary_convergence_clean) {
  writeLines(
    primary_convergence_messages,
    con = file.path(
      model_dir,
      "order_interaction_convergence_messages.txt"
    )
  )
}

# A singular primary fit is reported rather than silently hidden.
# Convergence failure remains a hard stop.
stopifnot(
  nobs(m_order_int) == 101,
  primary_convergence_clean
)


# ============================================================
# 4. FORMAL TYPE III TESTS
# ============================================================

order_anova <- anova(
  m_order_int,
  type = 3,
  ddf = "Kenward-Roger"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term")

write_csv(
  order_anova,
  file.path(
    model_dir,
    "order_type3_kenward_roger_tests.csv"
  )
)


# ============================================================
# 5. PRIMARY ORDER CONTRASTS + INTERACTION
# ============================================================

emm_order <- emmeans(
  m_order_int,
  ~ Order | Metal,
  lmer.df = "kenward-roger"
)

order_contrasts <- contrast(
  emm_order,
  method = list(
    "Harpacticoida_vs_Calanoida" = c(-1, 1)
  )
)

order_contrasts_summary <- summary(
  order_contrasts,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame() %>%
  mutate(
    LC50_ratio_Harp_Cal = exp(estimate),
    ratio_lower_95 = exp(lower.CL),
    ratio_upper_95 = exp(upper.CL)
  )

interaction_contrast <- contrast(
  order_contrasts,
  method = list(
    "Cd_vs_Cu_order_contrast" = c(-1, 1)
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

order_emmeans <- emmeans(
  m_order_int,
  ~ Metal * Order,
  lmer.df = "kenward-roger"
)

order_emmeans_table <- summary(
  order_emmeans,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

baseline_metrics <- extract_order_metrics(m_order_int)

write_csv(
  order_contrasts_summary,
  file.path(
    model_dir,
    "order_metal_specific_contrasts.csv"
  )
)

write_csv(
  interaction_summary,
  file.path(
    model_dir,
    "order_interaction_contrast.csv"
  )
)

write_csv(
  order_emmeans_table,
  file.path(
    model_dir,
    "order_estimated_marginal_means.csv"
  )
)

write_csv(
  baseline_metrics,
  file.path(
    model_dir,
    "order_baseline_metrics.csv"
  )
)


# ============================================================
# 5B. POST-MODEL MULTIPLICITY SENSITIVITY: HOLM ADJUSTMENT
# ============================================================

order_holm_sensitivity <- bind_rows(
  order_contrasts_summary %>%
    transmute(
      Model = "Taxonomic-order",
      Contrast = paste0(
        as.character(Metal),
        ": Harpacticoida / Calanoida"
      ),
      Effect_type = "LC50 ratio",
      Effect_estimate = LC50_ratio_Harp_Cal,
      lower_95 = ratio_lower_95,
      upper_95 = ratio_upper_95,
      raw_p = p.value
    ),
  
  interaction_summary %>%
    transmute(
      Model = "Taxonomic-order",
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
  order_holm_sensitivity,
  file.path(
    sensitivity_dir,
    "order_multiplicity_holm_sensitivity.csv"
  )
)


# ============================================================
# 6. ECOTOX-ONLY REPRODUCTION / AUGMENTATION CHECK
# ============================================================
# The pre-augmentation ECOTOX-only analysis used 91 Results.
# This section verifies that frozen benchmark and then quantifies the
# effect of adding 10 source-verified WoS 96-h Cu-Cd Order Results.

order_ecotox_only <- order_model_data %>%
  filter(Source_Origin == "ECOTOX") %>%
  droplevels()

stopifnot(nrow(order_ecotox_only) == 91)

m_order_ecotox <- fit_order_model(order_ecotox_only)

ecotox_metrics <- extract_order_metrics(m_order_ecotox)

# Frozen ECOTOX-only benchmark values retained as a regression test.
# They do not define the combined analysis result.
ecotox_regression_check <- tibble(
  Estimand = c(
    "Cu Harpacticoida / Calanoida",
    "Cd Harpacticoida / Calanoida",
    "Cd-vs-Cu ratio of ratios"
  ),
  Expected = c(
    9.8301,
    5.4786,
    0.5573
  ),
  ECOTOX_reproduction = c(
    ecotox_metrics$Cu_ratio,
    ecotox_metrics$Cd_ratio,
    ecotox_metrics$ratio_of_ratios
  )
) %>%
  mutate(
    Absolute_difference =
      abs(ECOTOX_reproduction - Expected),
    Match = Absolute_difference <
      c(0.15, 0.12, 0.03)
  )

stopifnot(all(ecotox_regression_check$Match))

ecotox_vs_combined <- tibble(
  Estimand = c(
    "Cu Harpacticoida / Calanoida",
    "Cd Harpacticoida / Calanoida",
    "Cd-vs-Cu ratio of ratios"
  ),
  Legacy_ECOTOX_benchmark = c(
    9.8301,
    5.4786,
    0.5573
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
    
    Direction_same_as_ECOTOX = c(
      (Combined[1] > 1) == (ECOTOX_reproduction[1] > 1),
      (Combined[2] > 1) == (ECOTOX_reproduction[2] > 1),
      (Combined[3] < 1) == (ECOTOX_reproduction[3] < 1)
    )
  )

write_csv(
  ecotox_regression_check,
  file.path(
    synthesis_dir,
    "ecotox_only_order_regression_check.csv"
  )
)

write_csv(
  ecotox_vs_combined,
  file.path(
    synthesis_dir,
    "ecotox_vs_combined_order_shift.csv"
  )
)

saveRDS(
  m_order_ecotox,
  file.path(
    sensitivity_dir,
    "order_ECOTOX_only_model.rds"
  )
)


# ============================================================
# 7. INTERMEDIATE / QA ORDER EFFECT-RATIO FIGURE
# ============================================================

order_effect_plot_data <- tibble(
  Estimand = factor(
    c(
      "Cu: Harpacticoida / Calanoida",
      "Cd: Harpacticoida / Calanoida",
      "Cd / Cu order-contrast ratio"
    ),
    levels = rev(
      c(
        "Cu: Harpacticoida / Calanoida",
        "Cd: Harpacticoida / Calanoida",
        "Cd / Cu order-contrast ratio"
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

order_effect_plot <- ggplot(
  order_effect_plot_data,
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
    "order_interaction_effect_ratios.png"
  ),
  order_effect_plot,
  width = 8.8,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    model_dir,
    "order_interaction_effect_ratios.pdf"
  ),
  order_effect_plot,
  width = 8.8,
  height = 4.8
)


# ============================================================
# 8. DIAGNOSTICS
# ============================================================

order_diagnostics <- order_model_data %>%
  mutate(
    fitted_value = fitted(m_order_int),
    residual = residuals(m_order_int),
    scaled_residual = residuals(
      m_order_int,
      type = "pearson",
      scaled = TRUE
    ),
    abs_scaled_residual = abs(scaled_residual)
  )

top_residuals <- order_diagnostics %>%
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
        "Order",
        "Lifestage",
        "LC50_umol_L",
        "ln_LC50",
        "Temperature",
        "Salinity",
        "fitted_value",
        "residual",
        "scaled_residual",
        "abs_scaled_residual"
      )
    )
  ) %>%
  slice_head(n = 10)

write_csv(
  order_diagnostics,
  file.path(
    diagnostic_dir,
    "order_observation_diagnostics.csv"
  )
)

write_csv(
  top_residuals,
  file.path(
    diagnostic_dir,
    "order_top_10_scaled_residuals.csv"
  )
)

png(
  file.path(
    diagnostic_dir,
    "order_residuals_vs_fitted.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

plot(
  fitted(m_order_int),
  resid(m_order_int),
  xlab = "Fitted ln(LC50)",
  ylab = "Residual",
  main = "Residuals versus fitted values"
)

abline(h = 0, lty = 2)
dev.off()

png(
  file.path(
    diagnostic_dir,
    "order_residuals_normal_qq.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

qqnorm(
  resid(m_order_int),
  main = "Normal Q-Q plot: residuals"
)

qqline(resid(m_order_int))
dev.off()

ref_re <- ranef(m_order_int)$Reference_ID[, 1]

png(
  file.path(
    diagnostic_dir,
    "order_reference_random_intercepts_qq.png"
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

species_re <- ranef(m_order_int)$Species[, 1]

png(
  file.path(
    diagnostic_dir,
    "order_species_random_intercepts_qq.png"
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
  as.character(order_model_data$Reference_ID)
)

loro_results <- purrr::map_dfr(
  refs,
  function(ref_out) {
    
    dat_i <- order_model_data %>%
      filter(
        as.character(Reference_ID) != ref_out
      )
    
    tryCatch({
      
      mod_i <- fit_order_model(dat_i)
      metrics_i <- extract_order_metrics(mod_i)
      
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
    
    Cu_Harp_higher_n =
      sum(Cu_ratio > 1),
    
    Cd_Harp_higher_n =
      sum(Cd_ratio > 1),
    
    RoR_below_1_n =
      sum(ratio_of_ratios < 1),
    
    singular_refits =
      sum(singular, na.rm = TRUE),
    
    convergence_warning_refits =
      sum(!is.na(convergence_message))
  )

write_csv(
  loro_results,
  file.path(
    robustness_dir,
    "order_LORO_all_refits.csv"
  )
)

write_csv(
  loro_summary,
  file.path(
    robustness_dir,
    "order_LORO_summary.csv"
  )
)

stopifnot(
  loro_summary$n_successful_refits ==
    n_distinct(order_model_data$Reference_ID)
)


# ============================================================
# 10. ROBUSTNESS: LEAVE-ONE-SPECIES-OUT
# ============================================================

species_list <- unique(
  as.character(order_model_data$Species)
)

loso_results <- purrr::map_dfr(
  species_list,
  function(sp_out) {
    
    dat_i <- order_model_data %>%
      filter(
        as.character(Species) != sp_out
      )
    
    tryCatch({
      
      mod_i <- fit_order_model(dat_i)
      metrics_i <- extract_order_metrics(mod_i)
      
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
    
    Cu_Harp_higher_n =
      sum(Cu_ratio > 1),
    
    Cd_Harp_higher_n =
      sum(Cd_ratio > 1),
    
    RoR_below_1_n =
      sum(ratio_of_ratios < 1),
    
    singular_refits =
      sum(singular, na.rm = TRUE),
    
    convergence_warning_refits =
      sum(!is.na(convergence_message))
  )

write_csv(
  loso_results,
  file.path(
    robustness_dir,
    "order_LOSO_all_refits.csv"
  )
)

write_csv(
  loso_summary,
  file.path(
    robustness_dir,
    "order_LOSO_summary.csv"
  )
)

stopifnot(
  loso_summary$n_successful_refits ==
    n_distinct(order_model_data$Species)
)


# ============================================================
# 11. DIRECTION-CHANGE CHECKS + INTERACTION EXTREMES
# ============================================================

baseline_Cu_above1 <-
  baseline_metrics$Cu_ratio > 1

baseline_Cd_above1 <-
  baseline_metrics$Cd_ratio > 1

baseline_RoR_below1 <-
  baseline_metrics$ratio_of_ratios < 1

loro_direction_changes <- loro_results %>%
  filter(!is.na(Cu_ratio)) %>%
  mutate(
    Cu_direction_changed =
      (Cu_ratio > 1) != baseline_Cu_above1,
    Cd_direction_changed =
      (Cd_ratio > 1) != baseline_Cd_above1,
    RoR_direction_changed =
      (ratio_of_ratios < 1) != baseline_RoR_below1
  ) %>%
  filter(
    Cu_direction_changed |
      Cd_direction_changed |
      RoR_direction_changed
  )

loso_direction_changes <- loso_results %>%
  filter(!is.na(Cu_ratio)) %>%
  mutate(
    Cu_direction_changed =
      (Cu_ratio > 1) != baseline_Cu_above1,
    Cd_direction_changed =
      (Cd_ratio > 1) != baseline_Cd_above1,
    RoR_direction_changed =
      (ratio_of_ratios < 1) != baseline_RoR_below1
  ) %>%
  filter(
    Cu_direction_changed |
      Cd_direction_changed |
      RoR_direction_changed
  )

loro_interaction_extremes <- loro_results %>%
  filter(!is.na(ratio_of_ratios)) %>%
  arrange(ratio_of_ratios) %>%
  bind_rows(
    loro_results %>%
      filter(!is.na(ratio_of_ratios)) %>%
      arrange(desc(ratio_of_ratios)) %>%
      slice_head(n = 1)
  ) %>%
  distinct(removed_reference, .keep_all = TRUE) %>%
  slice_head(n = 2)

loso_interaction_extremes <- loso_results %>%
  filter(!is.na(ratio_of_ratios)) %>%
  arrange(ratio_of_ratios) %>%
  bind_rows(
    loso_results %>%
      filter(!is.na(ratio_of_ratios)) %>%
      arrange(desc(ratio_of_ratios)) %>%
      slice_head(n = 1)
  ) %>%
  distinct(removed_species, .keep_all = TRUE) %>%
  slice_head(n = 2)

write_csv(
  loro_direction_changes,
  file.path(
    robustness_dir,
    "order_LORO_direction_changes.csv"
  )
)

write_csv(
  loso_direction_changes,
  file.path(
    robustness_dir,
    "order_LOSO_direction_changes.csv"
  )
)

write_csv(
  loro_interaction_extremes,
  file.path(
    robustness_dir,
    "order_LORO_interaction_extremes.csv"
  )
)

write_csv(
  loso_interaction_extremes,
  file.path(
    robustness_dir,
    "order_LOSO_interaction_extremes.csv"
  )
)


# ============================================================
# 12. ENVIRONMENTAL COMPLETE-CASE SUPPORT
# ============================================================
# Temperature/salinity adjustment remains a sensitivity of the Order
# model, not the separate RQ3 environmental analysis.
# pH is assessed separately in Script 08 and is not added to this sensitivity model.

order_env_cc <- order_model_data %>%
  filter(
    !is.na(Temperature),
    !is.na(Salinity)
  ) %>%
  droplevels()

order_env_cc_overall <- order_env_cc %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species)
  )

order_env_cc_support <- order_env_cc %>%
  group_by(Metal, Order) %>%
  summarise(
    Results = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    .groups = "drop"
  )

write_csv(
  order_env_cc_overall,
  file.path(
    sensitivity_dir,
    "order_environment_complete_case_overall.csv"
  )
)

write_csv(
  order_env_cc_support,
  file.path(
    sensitivity_dir,
    "order_environment_complete_case_support.csv"
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

order_environment_comparison <- tibble()
environment_fit_checks <- tibble()
temp_sal_correlation <- tibble()

cc_cell_support <- order_env_cc %>%
  count(Metal, Order, name = "n")

can_fit_env <- (
  nrow(cc_cell_support) == 4 &&
    all(cc_cell_support$n >= 3) &&
    n_distinct(order_env_cc$Reference_ID) >= 8 &&
    n_distinct(order_env_cc$Species) >= 6 &&
    sd(order_env_cc$Temperature, na.rm = TRUE) > 0 &&
    sd(order_env_cc$Salinity, na.rm = TRUE) > 0
)

if (can_fit_env) {
  
  m_order_cc <- fit_order_model(order_env_cc)
  
  order_env_cc_scaled <- order_env_cc %>%
    mutate(
      Temperature_z =
        as.numeric(scale(Temperature)),
      Salinity_z =
        as.numeric(scale(Salinity))
    )
  
  m_order_env <- lmer(
    ln_LC50 ~
      Metal * Order +
      Temperature_z +
      Salinity_z +
      (1 | Reference_ID) +
      (1 | Species),
    data = order_env_cc_scaled,
    REML = TRUE
  )
  
  environment_fit_checks <- tibble(
    Model = c(
      "Complete-case same formula",
      "Complete-case + Temperature/Salinity"
    ),
    Singular = c(
      isSingular(m_order_cc, tol = 1e-4),
      isSingular(m_order_env, tol = 1e-4)
    ),
    Convergence_clean = c(
      is.null(
        m_order_cc@optinfo$conv$lme4$messages
      ),
      is.null(
        m_order_env@optinfo$conv$lme4$messages
      )
    )
  )
  
  baseline_env_metrics <-
    extract_order_metrics(m_order_int) %>%
    mutate(model = "PRIMARY_FULL_DATA")
  
  cc_metrics <-
    extract_order_metrics(m_order_cc) %>%
    mutate(
      model =
        "TEMP_SAL_COMPLETE_CASE_SAME_FORMULA"
    )
  
  adjusted_metrics <-
    extract_order_metrics(m_order_env) %>%
    mutate(
      model =
        "TEMP_SAL_COMPLETE_CASE_ADJUSTED"
    )
  
  order_environment_comparison <- bind_rows(
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
      order_env_cc$Temperature,
      order_env_cc$Salinity,
      use = "complete.obs"
    )
  )
  
  environmental_models_status$Status <- "RUN"
  
  writeLines(
    capture.output(summary(m_order_cc)),
    con = file.path(
      sensitivity_dir,
      "order_complete_case_same_formula_model_summary.txt"
    )
  )
  
  writeLines(
    capture.output(summary(m_order_env)),
    con = file.path(
      sensitivity_dir,
      "order_environment_adjusted_model_summary.txt"
    )
  )
  
  write_csv(
    environment_fit_checks,
    file.path(
      sensitivity_dir,
      "order_environment_model_fit_checks.csv"
    )
  )
  
  write_csv(
    order_environment_comparison,
    file.path(
      sensitivity_dir,
      "order_environmental_sensitivity_comparison.csv"
    )
  )
  
  write_csv(
    temp_sal_correlation,
    file.path(
      sensitivity_dir,
      "order_temperature_salinity_correlation.csv"
    )
  )
  
  saveRDS(
    m_order_cc,
    file.path(
      sensitivity_dir,
      "order_complete_case_same_formula_model.rds"
    )
  )
  
  saveRDS(
    m_order_env,
    file.path(
      sensitivity_dir,
      "order_complete_case_environment_adjusted_model.rds"
    )
  )
  
} else {
  
  environmental_models_status$Status <-
    "NOT RUN: insufficient four-cell/common-support structure"
}

write_csv(
  environmental_models_status,
  file.path(
    sensitivity_dir,
    "order_environmental_sensitivity_status.csv"
  )
)


# ============================================================
# 14. ROBUSTNESS SYNTHESIS
# ============================================================

order_robustness_synthesis <- tibble(
  Analysis = c(
    "Combined baseline",
    "ECOTOX-only reproduction",
    "LORO range",
    "LOSO range"
  ),
  
  Cu_Harp_Cal = c(
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
  
  Cd_Harp_Cal = c(
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
  
  order_robustness_synthesis <- bind_rows(
    order_robustness_synthesis,
    tibble(
      Analysis = c(
        "Temperature/salinity complete-case",
        "Temperature/salinity adjusted"
      ),
      
      Cu_Harp_Cal = c(
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
      
      Cd_Harp_Cal = c(
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
  order_robustness_synthesis,
  file.path(
    robustness_dir,
    "order_robustness_synthesis.csv"
  )
)


# ============================================================
# 15. CONCLUSION-STABILITY SNAPSHOT
# ============================================================

conclusion_stability <- tibble(
  Question = c(
    "Cu: is Harpacticoida LC50 higher than Calanoida?",
    "Cd: is Harpacticoida LC50 higher than Calanoida?",
    "Is the Cd order contrast smaller than the Cu order contrast?"
  ),
  
  Combined_primary = c(
    baseline_metrics$Cu_ratio > 1,
    baseline_metrics$Cd_ratio > 1,
    baseline_metrics$ratio_of_ratios < 1
  ),
  
  ECOTOX_only = c(
    ecotox_metrics$Cu_ratio > 1,
    ecotox_metrics$Cd_ratio > 1,
    ecotox_metrics$ratio_of_ratios < 1
  ),
  
  LORO_direction_consistency = c(
    loro_summary$Cu_Harp_higher_n ==
      loro_summary$n_successful_refits,
    
    loro_summary$Cd_Harp_higher_n ==
      loro_summary$n_successful_refits,
    
    loro_summary$RoR_below_1_n ==
      loro_summary$n_successful_refits
  ),
  
  LOSO_direction_consistency = c(
    loso_summary$Cu_Harp_higher_n ==
      loso_summary$n_successful_refits,
    
    loso_summary$Cd_Harp_higher_n ==
      loso_summary$n_successful_refits,
    
    loso_summary$RoR_below_1_n ==
      loso_summary$n_successful_refits
  )
)

write_csv(
  conclusion_stability,
  file.path(
    synthesis_dir,
    "order_conclusion_stability_summary.csv"
  )
)


# ============================================================
# 16. SAVE OBJECTS + SESSION INFO
# ============================================================

saveRDS(
  list(
    order_model_data = order_model_data,
    primary_model = m_order_int,
    order_anova = order_anova,
    order_emmeans = order_emmeans,
    order_contrasts = order_contrasts_summary,
    interaction_contrast = interaction_summary,
    baseline_metrics = baseline_metrics,
    ECOTOX_only_model = m_order_ecotox,
    ECOTOX_only_metrics = ecotox_metrics,
    ecotox_vs_combined = ecotox_vs_combined,
    order_effect_plot = order_effect_plot,
    diagnostics = order_diagnostics,
    loro_results = loro_results,
    loso_results = loso_results,
    loro_direction_changes = loro_direction_changes,
    loso_direction_changes = loso_direction_changes,
    environmental_models_run = can_fit_env,
    environmental_comparison = order_environment_comparison
  ),
  file.path(
    object_dir,
    "order_model_analysis_objects.rds"
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
    "Input: Order support audit objects",
    "Output: primary model contrasts",
    "Output: LORO summary",
    "Output: LOSO summary",
    "Output: environmental sensitivity",
    "Output: reusable analysis objects"
  ),
  Path = c(
    audit_rds,
    file.path(model_dir, "order_metal_specific_contrasts.csv"),
    file.path(robustness_dir, "order_LORO_summary.csv"),
    file.path(robustness_dir, "order_LOSO_summary.csv"),
    file.path(sensitivity_dir, "order_environmental_sensitivity_comparison.csv"),
    file.path(object_dir, "order_model_analysis_objects.rds")
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
  file.path(
    object_dir,
    "order_model_analysis_objects.rds"
  ),
  file.path(
    model_dir,
    "order_interaction_effect_ratios.png"
  ),
  file.path(
    model_dir,
    "order_baseline_metrics.csv"
  ),
  file.path(
    model_dir,
    "order_type3_kenward_roger_tests.csv"
  ),
  file.path(
    diagnostic_dir,
    "order_residuals_vs_fitted.png"
  ),
  file.path(
    diagnostic_dir,
    "order_residuals_normal_qq.png"
  ),
  file.path(
    robustness_dir,
    "order_LORO_summary.csv"
  ),
  file.path(
    robustness_dir,
    "order_LOSO_summary.csv"
  ),
  file.path(
    sensitivity_dir,
    "order_multiplicity_holm_sensitivity.csv"
  ),
  file.path(
    synthesis_dir,
    "ecotox_vs_combined_order_shift.csv"
  ),
  file.path(
    synthesis_dir,
    "order_conclusion_stability_summary.csv"
  )
)

stopifnot(all(file.exists(required_outputs)))


# ============================================================
# 18. FINAL CONSOLE SUMMARY
# ============================================================

cat("\n\n")
cat("============================================================\n")
cat("13 TAXONOMIC-ORDER MODEL - FINAL SUMMARY\n")
cat("============================================================\n")

cat("\n--- COMBINED MODEL SUPPORT ---\n")
print(order_support, n = Inf, width = Inf)
print(order_overall_support, width = Inf)

cat("\n--- PRIMARY MODEL FIT ---\n")
print(primary_fit_checks, n = Inf)

cat("\n--- TYPE III TESTS ---\n")
print(
  tibble::as_tibble(order_anova),
  n = Inf,
  width = Inf
)

cat("\n--- PRIMARY ORDER RATIOS ---\n")
print(baseline_metrics, width = Inf)

cat("\n--- ECOTOX-ONLY / COMBINED SHIFT ---\n")
print(ecotox_vs_combined, n = Inf, width = Inf)

cat("\n--- LORO SUMMARY ---\n")
print(loro_summary, width = Inf)

cat("\n--- LOSO SUMMARY ---\n")
print(loso_summary, width = Inf)

cat("\n--- LORO DIRECTION CHANGES ---\n")
print(loro_direction_changes, n = Inf, width = Inf)

cat("\n--- LOSO DIRECTION CHANGES ---\n")
print(loso_direction_changes, n = Inf, width = Inf)

cat("\n--- TEMPERATURE/SALINITY SENSITIVITY STATUS ---\n")
print(environmental_models_status, n = Inf, width = Inf)

if (can_fit_env) {
  cat("\n--- TEMPERATURE/SALINITY SENSITIVITY ---\n")
  print(
    order_environment_comparison,
    n = Inf,
    width = Inf
  )
}

cat("\n--- CONCLUSION STABILITY ---\n")
print(conclusion_stability, n = Inf, width = Inf)

cat(
  "\nInterpretive rule:\n",
  "The combined 96-h Cu-Cd Calanoida/Harpacticoida interaction model is the primary Order analysis.\n",
  "The ECOTOX-only model is retained only to quantify the effect of WoS augmentation.\n",
  "Alternative metals remain in descriptive support outputs but are not comparably supported for the\n",
  "primary Metal x Order estimand: Zn has only two Calanoida References; Ag remains limited; Ni/Hg\n",
  "have sparse Harpacticoida support; and Pb has no Harpacticoida cell.\n",
  sep = ""
)

writeLines(
  c(
    "STATUS: 13 FIT ORDER MODEL = PASS",
    "Primary domain: 96-h Cu-Cd; Calanoida vs Harpacticoida",
    paste0("Results: ", nrow(order_model_data)),
    paste0("References: ", n_distinct(order_model_data$Reference_ID)),
    paste0("Species: ", n_distinct(order_model_data$Species))
  ),
  file.path(output_root, "RUN_COMPLETE.txt")
)

cat(
  "\nSTATUS: 13 FIT ORDER MODEL = PASS\n"
)
# ============================================================
# END
# ============================================================
