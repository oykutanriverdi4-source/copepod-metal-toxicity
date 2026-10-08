# ============================================================
# 04_fit_pooled_metal_model.R
# ============================================================
# Purpose
# Fit the primary pooled 96-h Cu-Cd-Zn linear mixed-effects model
# and run the associated diagnostics, robustness and sensitivity
# analyses used in the thesis.
#
# Primary estimand
# Pairwise mean differences in ln(LC50) among Cu, Cd and Zn,
# reported on the original scale as multiplicative LC50 ratios.
#
# Primary model
# ln(LC50) ~ Metal + (1 | Reference_ID) + (1 | Species)
#
# Input
# outputs/03_audit_metal_support/03_objects/metal_support_audit_objects.rds
#
# Main outputs
# outputs/04_fit_pooled_metal_model/
#
# Scientific workflow notes
# - The 96-h Cu-Cd-Zn domain is inherited from the pre-model support
#   audit; this script does not re-select the duration or metal panel.
# - Leave-one-Reference-out and leave-one-Species-out analyses are
#   influence checks, not alternative primary models.
# - Reference fixed-effects, Reference-only random-intercept and
#   Reference x Species x Metal cell-collapse analyses assess model-
#   structure and weighting sensitivity.
# - The ECOTOX-only analysis evaluates sensitivity to the addition of
#   source-verified WoS records.
# - Exploratory analyses that were not retained in the final thesis
#   inferential framework are excluded from this public reproducibility script.
# - Environmental sensitivity is evaluated on complete-case records
#   so that data restriction can be distinguished from covariate
#   adjustment.
# - No model is selected or discarded solely on the basis of p-values.
#
# Metadata limitation
# Exposure-design and chemical-analysis metadata are not encoded in
# fully equivalent fields across ECOTOX and the supplementary WoS
# records. Protocol-specific sensitivities that require harmonized
# cross-source metadata are therefore not silently imposed.
#
# Reproducibility rule
# Run from the repository root. Script 03 must be completed first.
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
  "04_fit_pooled_metal_model"
)

model_dir <- file.path(output_root, "01_model")
diagnostic_dir <- file.path(output_root, "02_diagnostics")
robustness_dir <- file.path(output_root, "03_robustness")
sensitivity_dir <- file.path(output_root, "04_sensitivity")
synthesis_dir <- file.path(output_root, "05_synthesis")
object_dir <- file.path(output_root, "06_objects")

# Invalidate the previous completion marker before overwriting this run's outputs.
# No output directories are renamed or archived.
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
# 1. LOAD COMBINED SUPPORT-AUDIT OBJECTS FROM SCRIPT 03
# ============================================================

audit_object_file <- file.path(
  "outputs",
  "03_audit_metal_support",
  "03_objects",
  "metal_support_audit_objects.rds"
)

if (!file.exists(audit_object_file)) {
  stop(
    paste0(
      "Could not find the combined support-audit object:\n  ",
      audit_object_file,
      "\n\nRun scripts/analysis/03_audit_metal_support.R first."
    )
  )
}

audit_objects <- readRDS(audit_object_file)

required_objects <- c(
  "common_pool",
  "pool_96",
  "metal_support_96",
  "pair_connectivity_96"
)

missing_objects <- setdiff(required_objects, names(audit_objects))

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Script 03 support-audit object is missing: ",
      paste(missing_objects, collapse = ", ")
    )
  )
}

pool_96_all <- audit_objects$pool_96

required_columns <- c(
  "Metal",
  "Species",
  "Order",
  "Lifestage",
  "Duration_days",
  "LC50_umol_L",
  "ln_LC50_umol_L",
  "Reference_ID",
  "Test_ID",
  "Result_ID",
  "Source_Origin",
  "Exposure_Type",
  "Conc_Type",
  "Basis_Interpretation",
  "Chemical_Name"
)

missing_columns <- setdiff(required_columns, names(pool_96_all))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The 96-h combined pool is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Script-03 benchmark for the complete harmonized 96-h domain.
stopifnot(
  nrow(pool_96_all) == 177,
  all(pool_96_all$Duration_days == 4),
  dplyr::n_distinct(pool_96_all$Result_ID) == nrow(pool_96_all)
)


# ============================================================
# 2. LOCK PRIMARY 96-h Cu-Cd-Zn DOMAIN
# ============================================================

primary_metals <- c("Cu", "Cd", "Zn")

primary96 <- pool_96_all %>%
  dplyr::filter(Metal %in% primary_metals) %>%
  dplyr::mutate(
    Metal = factor(Metal, levels = primary_metals),
    Reference_ID = factor(Reference_ID),
    Species = factor(Species)
  ) %>%
  droplevels()

stopifnot(nrow(primary96) == 131L,
          dplyr::n_distinct(primary96$Reference_ID) == 47L,
          dplyr::n_distinct(primary96$Species) == 22L,
          !anyNA(primary96$Species), !anyNA(primary96$Reference_ID),
          all(nzchar(trimws(as.character(primary96$Species)))),
          all(is.finite(primary96$ln_LC50_umol_L)))

primary_support <- primary96 %>%
  dplyr::group_by(Metal) %>%
  dplyr::summarise(
    Results = n(),
    Tests = dplyr::n_distinct(Test_ID),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    ECOTOX_Results = sum(Source_Origin == "ECOTOX"),
    WoS_Results = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  )

expected_primary_results <- tribble(
  ~Metal, ~Expected_Results,
  "Cu", 67L,
  "Cd", 37L,
  "Zn", 27L
)

primary_support_qa <- primary_support %>%
  dplyr::mutate(Metal = as.character(Metal)) %>%
  dplyr::left_join(expected_primary_results, by = "Metal") %>%
  dplyr::mutate(Match = Results == Expected_Results)

stopifnot(
  nrow(primary96) == 131,
  dplyr::n_distinct(primary96$Reference_ID) == 47,
  sum(primary96$Result_ID == "W_017_T01_R01") == 1L,
  all(primary_support_qa$Match),
  dplyr::n_distinct(primary96$Result_ID) == 131,
  all(levels(primary96$Metal) == primary_metals)
)

write_csv(
  primary_support,
  file.path(model_dir, "primary_model_dataset_support.csv")
)

write_csv(
  primary_support_qa,
  file.path(model_dir, "primary_model_dataset_support_QA.csv")
)

source_origin_support <- primary96 %>%
  dplyr::count(Metal, Source_Origin, name = "Results") %>%
  dplyr::group_by(Metal) %>%
  dplyr::mutate(Percent_results = 100 * Results / sum(Results)) %>%
  dplyr::ungroup()

write_csv(
  source_origin_support,
  file.path(model_dir, "primary_model_source_origin_support.csv")
)


# ============================================================
# 3. GROUPING-STRUCTURE AUDIT
# ============================================================

reference_structure <- primary96 %>%
  dplyr::group_by(Reference_ID) %>%
  dplyr::summarise(
    Results = n(),
    Tests = dplyr::n_distinct(Test_ID),
    Species = dplyr::n_distinct(Species),
    Metals = dplyr::n_distinct(Metal),
    Source_Origin = paste(sort(unique(Source_Origin)), collapse = " | "),
    Metal_list = paste(sort(unique(as.character(Metal))), collapse = " | "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(Results), Reference_ID)

species_structure <- primary96 %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    Results = n(),
    Tests = dplyr::n_distinct(Test_ID),
    References = dplyr::n_distinct(Reference_ID),
    Metals = dplyr::n_distinct(Metal),
    Metal_list = paste(sort(unique(as.character(Metal))), collapse = " | "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(Results), Species)

grouping_structure_summary <- tibble(
  Grouping_factor = c("Reference", "Species"),
  Levels = c(
    dplyr::n_distinct(primary96$Reference_ID),
    dplyr::n_distinct(primary96$Species)
  ),
  Levels_with_multiple_results = c(
    sum(reference_structure$Results > 1),
    sum(species_structure$Results > 1)
  ),
  Levels_with_multiple_metals = c(
    sum(reference_structure$Metals > 1),
    sum(species_structure$Metals > 1)
  ),
  Median_results_per_level = c(
    median(reference_structure$Results),
    median(species_structure$Results)
  ),
  Max_results_per_level = c(
    max(reference_structure$Results),
    max(species_structure$Results)
  )
)

write_csv(reference_structure, file.path(model_dir, "reference_grouping_structure.csv"))
write_csv(species_structure, file.path(model_dir, "species_grouping_structure.csv"))
write_csv(grouping_structure_summary, file.path(model_dir, "grouping_structure_summary.csv"))


# Joint support audit: review alongside the separate grouping summaries.
write_csv(primary96 %>% dplyr::count(Reference_ID, Species, Metal, name = "Records"),
          file.path(model_dir, "reference_species_metal_support.csv"))
write_csv(primary96 %>% dplyr::count(Species, Metal, name = "Records"),
          file.path(model_dir, "species_metal_support.csv"))

# ============================================================
# 4. HELPERS
# ============================================================

model_convergence_clean <- function(fit) {
  is.null(fit@optinfo$conv$lme4$messages) &&
    (is.null(fit@optinfo$conv$opt) || all(fit@optinfo$conv$opt == 0))
}

extract_pairwise_ratios <- function(fit, df_method = "kenward-roger") {
  emm <- emmeans(fit, ~ Metal, lmer.df = df_method)
  out <- summary(
    pairs(emm, adjust = "tukey"),
    infer = c(TRUE, TRUE)
  ) %>%
    as.data.frame()
  
  out %>%
    dplyr::mutate(
      LC50_ratio = exp(estimate),
      ratio_lower = exp(lower.CL),
      ratio_upper = exp(upper.CL)
    )
}

make_variant_table <- function(table, variant, family) {
  table %>%
    dplyr::transmute(
      Variant = variant,
      Variant_family = family,
      contrast,
      LC50_ratio,
      ratio_lower,
      ratio_upper,
      p.value = if ("p.value" %in% names(table)) p.value else NA_real_
    )
}

extract_variance <- function(vc_table, group_name) {
  out <- vc_table %>%
    dplyr::filter(grp == group_name) %>%
    pull(vcov)
  
  if (length(out) == 0) return(NA_real_)
  out[1]
}

print_full <- function(x) {
  print(tibble::as_tibble(x), n = Inf, width = Inf)
}


# ============================================================
# 5. POOLED REFERENCE + SPECIES RANDOM-INTERCEPT LMM
# ============================================================

m_metal_ref_species <- lmer(
  ln_LC50_umol_L ~ Metal + (1 | Reference_ID) + (1 | Species),
  data = primary96,
  REML = TRUE
)

primary_singular <- isSingular(m_metal_ref_species, tol = 1e-4)
primary_convergence_messages <- m_metal_ref_species@optinfo$conv$lme4$messages
primary_convergence_clean <- model_convergence_clean(m_metal_ref_species)

primary_fit_checks <- tibble(
  Check = c(
    "Fitted observations",
    "Reference levels",
    "Species random-intercept levels",
    "ECOTOX Results",
    "WoS Results",
    "Singular fit",
    "Convergence clean"
  ),
  Value = c(
    as.character(nobs(m_metal_ref_species)),
    as.character(dplyr::n_distinct(primary96$Reference_ID)),
    as.character(dplyr::n_distinct(primary96$Species)),
    as.character(sum(primary96$Source_Origin == "ECOTOX")),
    as.character(sum(primary96$Source_Origin == "WoS_supplemental")),
    as.character(primary_singular),
    as.character(primary_convergence_clean)
  )
)

primary_variance_components <- as.data.frame(VarCorr(m_metal_ref_species))

write_csv(primary_fit_checks, file.path(model_dir, "primary_model_fit_checks.csv"))
write_csv(primary_variance_components, file.path(model_dir, "primary_model_variance_components.csv"))

writeLines(
  capture.output(summary(m_metal_ref_species, ddf = "Kenward-Roger")),
  con = file.path(model_dir, "primary_model_summary.txt")
)

if (!primary_convergence_clean) {
  writeLines(
    primary_convergence_messages,
    con = file.path(model_dir, "primary_model_convergence_messages.txt")
  )
  warning("Primary model produced a convergence message. Review the saved diagnostics.")
}

if (primary_singular) {
  warning("Primary model is singular. Review variance components before interpretation.")
}

stopifnot(nobs(m_metal_ref_species) == 131)


if (primary_singular || !primary_convergence_clean || nobs(m_metal_ref_species) != nrow(primary96)) {
  stop("Reference + Species fit failed the fit gate. Review 01_model fit checks; no model simplification is performed automatically.")
}

# ============================================================
# 6. GLOBAL METAL TEST, EMMs AND PAIRWISE LC50 RATIOS
# ============================================================

global_metal_test <- anova(
  m_metal_ref_species,
  ddf = "Kenward-Roger"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term")

emm_metal <- emmeans(
  m_metal_ref_species,
  ~ Metal,
  lmer.df = "kenward-roger"
)

emm_metal_table <- summary(
  emm_metal,
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

emm_metal_original_scale <- emm_metal_table %>%
  dplyr::transmute(
    Metal,
    emmean_log = emmean,
    SE_log = SE,
    df,
    lower_log = lower.CL,
    upper_log = upper.CL,
    LC50_model_estimate_umol_L = exp(emmean),
    LC50_lower_umol_L = exp(lower.CL),
    LC50_upper_umol_L = exp(upper.CL)
  )

metal_pairs_log <- summary(
  pairs(emm_metal, adjust = "tukey"),
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

metal_pairs_ratio <- metal_pairs_log %>%
  dplyr::mutate(
    LC50_ratio = exp(estimate),
    ratio_lower = exp(lower.CL),
    ratio_upper = exp(upper.CL)
  ) %>%
  dplyr::select(
    contrast,
    estimate,
    lower.CL,
    upper.CL,
    LC50_ratio,
    ratio_lower,
    ratio_upper,
    SE,
    df,
    t.ratio,
    p.value
  )

write_csv(global_metal_test, file.path(model_dir, "global_metal_test.csv"))
write_csv(emm_metal_table, file.path(model_dir, "metal_estimated_marginal_means_log_scale.csv"))
write_csv(emm_metal_original_scale, file.path(model_dir, "metal_model_estimates_original_scale.csv"))
write_csv(metal_pairs_log, file.path(model_dir, "metal_pairwise_contrasts_log_scale.csv"))
write_csv(metal_pairs_ratio, file.path(model_dir, "metal_pairwise_lc50_ratios.csv"))


# ============================================================
# 7. PRE-AUGMENTATION VS COMBINED RATIO-SHIFT AUDIT
# ============================================================
# Historical tracking table retained for audit compatibility; it is not a model-selection benchmark.

pre_augmentation_ratio_benchmark <- tibble(
  contrast = c("Cu - Cd", "Cu - Zn", "Cd - Zn"),
  pre_augmentation_ratio = c(0.7837, 0.2797, 0.3569)
)

pre_augmentation_vs_combined_ratio_shift <- metal_pairs_ratio %>%
  dplyr::select(
    contrast,
    combined_ratio = LC50_ratio,
    combined_lower = ratio_lower,
    combined_upper = ratio_upper,
    combined_p_value = p.value
  ) %>%
  dplyr::left_join(pre_augmentation_ratio_benchmark, by = "contrast") %>%
  dplyr::mutate(
    ratio_multiplier = combined_ratio / pre_augmentation_ratio,
    percent_change = 100 * (combined_ratio - pre_augmentation_ratio) / pre_augmentation_ratio,
    pre_augmentation_direction = case_when(
      pre_augmentation_ratio < 1 ~ "<1",
      pre_augmentation_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    combined_direction = case_when(
      combined_ratio < 1 ~ "<1",
      combined_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    direction_changed = pre_augmentation_direction != combined_direction
  )

write_csv(
  pre_augmentation_vs_combined_ratio_shift,
  file.path(model_dir, "pre_augmentation_vs_combined_ratio_shift.csv")
)


# ============================================================
# 8. PRIMARY PAIRWISE RATIO FIGURE
# ============================================================

primary_ratio_plot_data <- metal_pairs_ratio %>%
  dplyr::mutate(comparison = stringr::str_replace(contrast, " - ", " / "))

primary_ratio_plot <- ggplot(
  primary_ratio_plot_data,
  aes(
    x = LC50_ratio,
    y = forcats::fct_rev(comparison)
  )
) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(
    aes(
      x = ratio_lower,
      xend = ratio_upper,
      yend = forcats::fct_rev(comparison)
    ),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(
    title = "Combined pooled 96-h metal contrasts",
    subtitle = "ECOTOX + source-verified WoS; Reference random-intercept LMM",
    x = "LC50 ratio (first metal / second metal; log scale)",
    y = NULL
  ) +
  theme_classic()

ggsave(
  file.path(model_dir, "primary_pooled_pairwise_lc50_ratios.png"),
  primary_ratio_plot,
  width = 8,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 9. RAW DATA + MODEL ESTIMATE FIGURE
# ============================================================
# The final thesis styling will still be handled in Script 10.
# This plot provides an immediate model/data QA view from Script 02.

raw_model_plot_data <- primary96 %>%
  dplyr::mutate(Metal = factor(Metal, levels = primary_metals))

model_plot_data <- emm_metal_original_scale %>%
  dplyr::mutate(Metal = factor(Metal, levels = primary_metals))

raw_model_plot <- ggplot(raw_model_plot_data, aes(x = Metal, y = LC50_umol_L)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.55, size = 1.8) +
  geom_errorbar(
    data = model_plot_data,
    aes(
      x = Metal,
      ymin = LC50_lower_umol_L,
      ymax = LC50_upper_umol_L
    ),
    inherit.aes = FALSE,
    width = 0.08,
    linewidth = 0.8
  ) +
  geom_point(
    data = model_plot_data,
    aes(x = Metal, y = LC50_model_estimate_umol_L),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  scale_y_log10() +
  labs(
    title = "Observed and model-estimated 96-h LC50",
    subtitle = "Combined ECOTOX + source-verified WoS evidence",
    x = NULL,
    y = expression(LC[50]~(mu*mol~L^{-1}))
  ) +
  theme_classic()

ggsave(
  file.path(model_dir, "primary_pooled_raw_plus_model_estimates.png"),
  raw_model_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 10. POST-MODEL DIAGNOSTICS
# ============================================================

diag_table <- primary96 %>%
  dplyr::mutate(
    fitted = fitted(m_metal_ref_species),
    residual = resid(m_metal_ref_species),
    abs_residual = abs(residual)
  ) %>%
  dplyr::arrange(desc(abs_residual)) %>%
  dplyr::select(
    Source_Origin,
    Reference_ID,
    Test_ID,
    Result_ID,
    Species,
    Metal,
    LC50_umol_L,
    ln_LC50_umol_L,
    fitted,
    residual,
    abs_residual
  )

species_diag <- ranef(m_metal_ref_species)$Species %>%
  tibble::rownames_to_column("Species") %>%
  dplyr::rename(species_intercept = `(Intercept)`)
write_csv(species_diag, file.path(diagnostic_dir, "species_random_intercept_diagnostics.csv"))
png(file.path(diagnostic_dir, "species_random_intercepts_normal_qq.png"),
    width = 1800, height = 1400, res = 220)
qqnorm(species_diag$species_intercept, main = "Normal Q-Q plot: Species random intercepts")
qqline(species_diag$species_intercept)
dev.off()

ref_diag <- ranef(m_metal_ref_species)$Reference_ID %>%
  tibble::rownames_to_column("Reference_ID") %>%
  dplyr::rename(ref_intercept = `(Intercept)`) %>%
  dplyr::mutate(abs_ref_intercept = abs(ref_intercept)) %>%
  dplyr::arrange(desc(abs_ref_intercept))

write_csv(diag_table, file.path(diagnostic_dir, "observation_residual_diagnostics.csv"))
write_csv(head(diag_table, 10), file.path(diagnostic_dir, "top_10_absolute_residuals.csv"))
write_csv(ref_diag, file.path(diagnostic_dir, "reference_random_intercept_diagnostics.csv"))

png(
  file.path(diagnostic_dir, "residuals_vs_fitted.png"),
  width = 1800,
  height = 1400,
  res = 220
)
plot(
  fitted(m_metal_ref_species),
  resid(m_metal_ref_species),
  xlab = "Fitted ln(LC50)",
  ylab = "Residual",
  main = "Residuals versus fitted values"
)
abline(h = 0, lty = 2)
dev.off()

png(
  file.path(diagnostic_dir, "residuals_normal_qq.png"),
  width = 1800,
  height = 1400,
  res = 220
)
qqnorm(resid(m_metal_ref_species), main = "Normal Q-Q plot: residuals")
qqline(resid(m_metal_ref_species))
dev.off()

ref_re <- ranef(m_metal_ref_species)$Reference_ID[, 1]

png(
  file.path(diagnostic_dir, "reference_random_intercepts_normal_qq.png"),
  width = 1800,
  height = 1400,
  res = 220
)
qqnorm(ref_re, main = "Normal Q-Q plot: Reference random intercepts")
qqline(ref_re)
dev.off()


# ============================================================
# 11. ROBUSTNESS 1: LEAVE-ONE-REFERENCE-OUT (LORO)
# ============================================================

refs <- levels(primary96$Reference_ID)

loro_results <- purrr::map_dfr(
  refs,
  function(ref_out) {
    dat_loo <- primary96 %>%
      dplyr::filter(as.character(Reference_ID) != ref_out) %>%
      droplevels()
    
    if (dplyr::n_distinct(dat_loo$Metal) != 3) return(tibble())
    
    fit_loo <- lmer(
      ln_LC50_umol_L ~ Metal + (1 | Reference_ID) + (1 | Species),
      data = dat_loo,
      REML = TRUE
    )
    
    extract_pairwise_ratios(fit_loo) %>%
      dplyr::mutate(
        removed_reference = ref_out,
        n_results = nrow(dat_loo),
        n_references = dplyr::n_distinct(dat_loo$Reference_ID),
        n_species = dplyr::n_distinct(dat_loo$Species),
        singular = isSingular(fit_loo, tol = 1e-4),
        convergence_clean = model_convergence_clean(fit_loo)
      )
  }
)

loro_summary <- loro_results %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    n_refits = n(),
    min_ratio = min(LC50_ratio),
    max_ratio = max(LC50_ratio),
    min_lower = min(ratio_lower),
    max_upper = max(ratio_upper),
    ratios_below_1 = sum(LC50_ratio < 1),
    ratios_above_1 = sum(LC50_ratio > 1),
    ratios_equal_1 = sum(LC50_ratio == 1),
    all_point_estimates_below_1 = all(LC50_ratio < 1),
    all_point_estimates_above_1 = all(LC50_ratio > 1),
    singular_models = sum(singular),
    convergence_failures = sum(!convergence_clean),
    .groups = "drop"
  )

write_csv(loro_results, file.path(robustness_dir, "leave_one_reference_out_all_refits.csv"))
write_csv(loro_summary, file.path(robustness_dir, "leave_one_reference_out_summary.csv"))

stopifnot(
  dplyr::n_distinct(loro_results$removed_reference) == dplyr::n_distinct(primary96$Reference_ID),
  all(loro_summary$n_refits == dplyr::n_distinct(primary96$Reference_ID))
)


# ============================================================
# 12. ROBUSTNESS 2: LEAVE-ONE-SPECIES-OUT (LOSO)
# ============================================================

species_levels <- levels(primary96$Species)

loso_results <- purrr::map_dfr(
  species_levels,
  function(species_out) {
    dat_loso <- primary96 %>%
      dplyr::filter(as.character(Species) != species_out) %>%
      droplevels()
    
    if (dplyr::n_distinct(dat_loso$Metal) != 3) return(tibble())
    
    fit_loso <- lmer(
      ln_LC50_umol_L ~ Metal + (1 | Reference_ID) + (1 | Species),
      data = dat_loso,
      REML = TRUE
    )
    
    extract_pairwise_ratios(fit_loso) %>%
      dplyr::mutate(
        removed_species = species_out,
        n_results = nrow(dat_loso),
        n_references = dplyr::n_distinct(dat_loso$Reference_ID),
        n_species = dplyr::n_distinct(dat_loso$Species),
        singular = isSingular(fit_loso, tol = 1e-4),
        convergence_clean = model_convergence_clean(fit_loso)
      )
  }
)

loso_summary <- loso_results %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    n_refits = n(),
    min_ratio = min(LC50_ratio),
    max_ratio = max(LC50_ratio),
    min_lower = min(ratio_lower),
    max_upper = max(ratio_upper),
    ratios_below_1 = sum(LC50_ratio < 1),
    ratios_above_1 = sum(LC50_ratio > 1),
    ratios_equal_1 = sum(LC50_ratio == 1),
    all_point_estimates_below_1 = all(LC50_ratio < 1),
    all_point_estimates_above_1 = all(LC50_ratio > 1),
    singular_models = sum(singular),
    convergence_failures = sum(!convergence_clean),
    .groups = "drop"
  )

write_csv(loso_results, file.path(robustness_dir, "leave_one_species_out_all_refits.csv"))
write_csv(loso_summary, file.path(robustness_dir, "leave_one_species_out_summary.csv"))

stopifnot(
  dplyr::n_distinct(loso_results$removed_species) == dplyr::n_distinct(primary96$Species),
  all(loso_summary$n_refits == dplyr::n_distinct(primary96$Species))
)


# Reporting-only counts: same saved Tukey intervals, no new tests/refits.
delete_interval_report <- function(dat, check_name, deleted_id, expected_n) {
  bind_rows(lapply(split(dat, as.character(dat$contrast)), function(one_contrast) {
    bind_cols(
      tibble(Check = check_name, contrast = as.character(one_contrast$contrast[1])),
      report_ratio_support(
        one_contrast, ratio = "LC50_ratio", lower = "ratio_lower",
        upper = "ratio_upper", id = deleted_id, expected_refits = expected_n,
        interval_adjustment = "Tukey: three metal contrasts within each refit"
      )
    )
  }))
}
deletion_direction_and_interval_summary <- bind_rows(
  delete_interval_report(
    loro_results, "LORO", "removed_reference", n_distinct(primary96$Reference_ID)
  ),
  delete_interval_report(
    loso_results, "LOSO", "removed_species", n_distinct(primary96$Species)
  )
)
write_csv(
  deletion_direction_and_interval_summary,
  file.path(robustness_dir, "deletion_direction_and_interval_summary.csv")
)

# ============================================================
# 13. ROBUSTNESS 3: REFERENCE FIXED-EFFECTS SENSITIVITY
# ============================================================

m_metal_ref_FE <- lmer(
  ln_LC50_umol_L ~ Reference_ID + Metal + (1 | Species),
  data = primary96, REML = TRUE
)

reference_fixed_effects_fit_checks <- bind_cols(
  tibble(
    Model = "Reference fixed effects + Species random intercept",
    Results = nobs(m_metal_ref_FE),
    References = n_distinct(primary96$Reference_ID),
    Species = n_distinct(primary96$Species),
    Legacy_convergence_clean = model_convergence_clean(m_metal_ref_FE),
    Diagnostic_scope = "Numerical status only; not a test of model adequacy"
  ),
  report_fit_diagnostics(m_metal_ref_FE)
)
write_csv(
  reference_fixed_effects_fit_checks,
  file.path(robustness_dir, "reference_fixed_effects_fit_checks.csv")
)

emm_ref_FE <- emmeans(m_metal_ref_FE, ~ Metal, lmer.df = "kenward-roger")

ref_FE_pairs <- summary(
  pairs(emm_ref_FE, adjust = "tukey"),
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

ref_FE_ratios <- ref_FE_pairs %>%
  dplyr::mutate(
    LC50_ratio = exp(estimate),
    ratio_lower = exp(lower.CL),
    ratio_upper = exp(upper.CL)
  ) %>%
  dplyr::select(
    contrast,
    estimate,
    lower.CL,
    upper.CL,
    LC50_ratio,
    ratio_lower,
    ratio_upper,
    SE,
    df,
    t.ratio,
    p.value
  )

writeLines(
  capture.output(summary(m_metal_ref_FE, ddf = "Kenward-Roger")),
  con = file.path(robustness_dir, "reference_fixed_effects_model_summary.txt")
)

write_csv(
  ref_FE_ratios,
  file.path(robustness_dir, "reference_fixed_effects_lc50_ratios.csv")
)


# ============================================================
# 14. ROBUSTNESS 4: REFERENCE-ONLY RANDOM INTERCEPT
# ============================================================

m_metal_ref_only <- lmer(
  ln_LC50_umol_L ~ Metal + (1 | Reference_ID),
  data = primary96,
  REML = TRUE
)

ref_only_singular <- isSingular(m_metal_ref_only, tol = 1e-4)
ref_only_convergence_messages <- m_metal_ref_only@optinfo$conv$lme4$messages
ref_only_convergence_clean <- model_convergence_clean(m_metal_ref_only)

emm_ref_only <- emmeans(
  m_metal_ref_only,
  ~ Metal,
  lmer.df = "kenward-roger"
)

ref_only_pairs <- summary(
  pairs(emm_ref_only, adjust = "tukey"),
  infer = c(TRUE, TRUE)
) %>%
  as.data.frame()

ref_only_ratios <- ref_only_pairs %>%
  dplyr::mutate(
    LC50_ratio = exp(estimate),
    ratio_lower = exp(lower.CL),
    ratio_upper = exp(upper.CL)
  ) %>%
  dplyr::select(
    contrast,
    estimate,
    lower.CL,
    upper.CL,
    LC50_ratio,
    ratio_lower,
    ratio_upper,
    SE,
    df,
    t.ratio,
    p.value
  )

ref_only_fit_checks <- tibble(
  Check = c(
    "Fitted observations",
    "Reference levels",
    "Species levels",
    "Singular fit",
    "Convergence clean"
  ),
  Value = c(
    as.character(nobs(m_metal_ref_only)),
    as.character(dplyr::n_distinct(primary96$Reference_ID)),
    as.character(dplyr::n_distinct(primary96$Species)),
    as.character(ref_only_singular),
    as.character(ref_only_convergence_clean)
  )
)

ref_only_variance_components <- as.data.frame(VarCorr(m_metal_ref_only))

writeLines("Comparison columns: primary_* = Reference + Species; ref_only_* = Reference only. Percent change is (Reference-only / Reference+Species - 1) * 100.",
           file.path(robustness_dir, "model_comparison_column_definitions.txt"))
primary_vs_ref_only <- metal_pairs_ratio %>%
  dplyr::select(
    contrast,
    primary_ratio = LC50_ratio,
    primary_lower = ratio_lower,
    primary_upper = ratio_upper
  ) %>%
  dplyr::left_join(
    ref_only_ratios %>%
      dplyr::select(
        contrast,
        ref_only_ratio = LC50_ratio,
        ref_only_lower = ratio_lower,
        ref_only_upper = ratio_upper
      ),
    by = "contrast"
  ) %>%
  dplyr::mutate(
    ratio_multiplier = ref_only_ratio / primary_ratio,
    percent_change = 100 * (ref_only_ratio - primary_ratio) / primary_ratio,
    primary_direction = case_when(
      primary_ratio < 1 ~ "<1",
      primary_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    ref_only_direction = case_when(
      ref_only_ratio < 1 ~ "<1",
      ref_only_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    direction_changed = primary_direction != ref_only_direction
  )

writeLines(
  capture.output(summary(m_metal_ref_only, ddf = "Kenward-Roger")),
  con = file.path(robustness_dir, "reference_only_random_intercept_model_summary.txt")
)

if (!ref_only_convergence_clean) {
  writeLines(
    ref_only_convergence_messages,
    con = file.path(robustness_dir, "reference_only_random_intercept_convergence_messages.txt")
  )
}

write_csv(ref_only_fit_checks, file.path(robustness_dir, "reference_only_random_intercept_fit_checks.csv"))
write_csv(ref_only_variance_components, file.path(robustness_dir, "reference_only_random_intercept_variance_components.csv"))
write_csv(ref_only_ratios, file.path(robustness_dir, "reference_only_random_intercept_lc50_ratios.csv"))
write_csv(primary_vs_ref_only, file.path(robustness_dir, "primary_vs_reference_only_random_intercepts.csv"))


# ============================================================
# 15. ROBUSTNESS 5: REFERENCE x SPECIES x METAL CELL COLLAPSE
# ============================================================

cluster_sizes_96 <- primary96 %>%
  dplyr::count(Reference_ID, Species, Metal, name = "Results_in_cell") %>%
  dplyr::arrange(desc(Results_in_cell))

cluster_size_distribution <- cluster_sizes_96 %>%
  dplyr::count(Results_in_cell, name = "Number_of_cells") %>%
  dplyr::arrange(Results_in_cell)

cluster_size_summary <- cluster_sizes_96 %>%
  dplyr::summarise(
    Total_cells = n(),
    Single_result_cells = sum(Results_in_cell == 1),
    Multi_result_cells = sum(Results_in_cell > 1),
    Max_results_in_cell = max(Results_in_cell)
  )

primary96_collapsed <- primary96 %>%
  dplyr::group_by(Reference_ID, Species, Metal) %>%
  dplyr::summarise(
    Results_in_cell = n(),
    ln_LC50_cell = mean(ln_LC50_umol_L),
    LC50_cell_umol_L = exp(mean(ln_LC50_umol_L)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Metal = factor(Metal, levels = primary_metals),
    Reference_ID = factor(Reference_ID),
    Species = factor(Species)
  ) %>%
  droplevels()

collapsed_support <- primary96_collapsed %>%
  dplyr::group_by(Metal) %>%
  dplyr::summarise(
    Cells = n(),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  )

m_metal_collapsed <- lmer(
  ln_LC50_cell ~ Metal + (1 | Reference_ID) + (1 | Species),
  data = primary96_collapsed,
  REML = TRUE
)

collapsed_singular <- isSingular(m_metal_collapsed, tol = 1e-4)
collapsed_convergence_clean <- model_convergence_clean(m_metal_collapsed)

collapsed_ratios <- extract_pairwise_ratios(m_metal_collapsed) %>%
  dplyr::select(
    contrast,
    estimate,
    lower.CL,
    upper.CL,
    LC50_ratio,
    ratio_lower,
    ratio_upper,
    SE,
    df,
    t.ratio,
    p.value
  )

weighting_comparison <- metal_pairs_ratio %>%
  dplyr::select(
    contrast,
    primary_ratio = LC50_ratio,
    primary_lower = ratio_lower,
    primary_upper = ratio_upper
  ) %>%
  dplyr::left_join(
    collapsed_ratios %>%
      dplyr::select(
        contrast,
        collapsed_ratio = LC50_ratio,
        collapsed_lower = ratio_lower,
        collapsed_upper = ratio_upper
      ),
    by = "contrast"
  ) %>%
  dplyr::mutate(ratio_change = collapsed_ratio / primary_ratio)

collapsed_fit_checks <- tibble(
  Check = c("Fitted cells", "Singular fit", "Convergence clean"),
  Value = c(
    as.character(nobs(m_metal_collapsed)),
    as.character(collapsed_singular),
    as.character(collapsed_convergence_clean)
  )
)

write_csv(cluster_sizes_96, file.path(robustness_dir, "reference_only_metal_cell_sizes.csv"))
write_csv(cluster_size_distribution, file.path(robustness_dir, "cell_size_distribution.csv"))
write_csv(cluster_size_summary, file.path(robustness_dir, "cell_size_summary.csv"))
write_csv(collapsed_support, file.path(robustness_dir, "collapsed_dataset_support.csv"))
writeLines(capture.output(summary(m_metal_collapsed, ddf = "Kenward-Roger")), con = file.path(robustness_dir, "collapsed_model_summary.txt"))
write_csv(collapsed_ratios, file.path(robustness_dir, "collapsed_model_lc50_ratios.csv"))
write_csv(weighting_comparison, file.path(robustness_dir, "primary_vs_cell_collapsed_comparison.csv"))
write_csv(collapsed_fit_checks, file.path(robustness_dir, "collapsed_model_fit_checks.csv"))


# ============================================================
# 16. AUGMENTATION SENSITIVITY: COMBINED VS ECOTOX-ONLY
# ============================================================
# This is the most direct check of whether the WoS augmentation changes
# the pooled Cu-Cd-Zn estimates relative to the original sampling frame.

primary96_ecotox <- primary96 %>%
  dplyr::filter(Source_Origin == "ECOTOX") %>%
  droplevels()

stopifnot(
  nrow(primary96_ecotox) == 117,
  dplyr::n_distinct(primary96_ecotox$Metal) == 3
)

ecotox_only_support <- primary96_ecotox %>%
  dplyr::group_by(Metal) %>%
  dplyr::summarise(
    Results = n(),
    Tests = dplyr::n_distinct(Test_ID),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  )

m_metal_ecotox_only <- lmer(
  ln_LC50_umol_L ~ Metal + (1 | Reference_ID) + (1 | Species),
  data = primary96_ecotox,
  REML = TRUE
)

ecotox_only_ratios <- extract_pairwise_ratios(m_metal_ecotox_only) %>%
  dplyr::select(
    contrast,
    LC50_ratio,
    ratio_lower,
    ratio_upper,
    SE,
    df,
    p.value
  )

combined_vs_ecotox_only <- metal_pairs_ratio %>%
  dplyr::select(
    contrast,
    Combined_ratio = LC50_ratio,
    Combined_lower = ratio_lower,
    Combined_upper = ratio_upper,
    Combined_p = p.value
  ) %>%
  dplyr::left_join(
    ecotox_only_ratios %>%
      dplyr::select(
        contrast,
        ECOTOX_only_ratio = LC50_ratio,
        ECOTOX_only_lower = ratio_lower,
        ECOTOX_only_upper = ratio_upper,
        ECOTOX_only_p = p.value
      ),
    by = "contrast"
  ) %>%
  dplyr::mutate(
    Ratio_multiplier_combined_vs_ECOTOX = Combined_ratio / ECOTOX_only_ratio,
    Percent_change = 100 * (Combined_ratio - ECOTOX_only_ratio) / ECOTOX_only_ratio,
    ECOTOX_direction = case_when(
      ECOTOX_only_ratio < 1 ~ "<1",
      ECOTOX_only_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    Combined_direction = case_when(
      Combined_ratio < 1 ~ "<1",
      Combined_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    Direction_changed = ECOTOX_direction != Combined_direction
  )

write_csv(ecotox_only_support, file.path(sensitivity_dir, "ecotox_only_primary_support.csv"))
writeLines(capture.output(summary(m_metal_ecotox_only, ddf = "Kenward-Roger")), con = file.path(sensitivity_dir, "ecotox_only_primary_model_summary.txt"))
write_csv(ecotox_only_ratios, file.path(sensitivity_dir, "ecotox_only_primary_lc50_ratios.csv"))
write_csv(combined_vs_ecotox_only, file.path(sensitivity_dir, "combined_vs_ecotox_only_lc50_ratios.csv"))


# ============================================================
# 17. PROTOCOL / CHEMISTRY COMPARABILITY AUDIT
# ============================================================
# Exposure_Type is tabulated by source because its semantics are not yet
# equivalent across ECOTOX and WoS. Other fields are retained as descriptive
# metadata audits rather than automatic inclusion/exclusion criteria.

primary_exposure_audit_by_source <- primary96 %>%
  dplyr::group_by(Source_Origin, Metal, Exposure_Type) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Source_Origin, Metal) %>%
  dplyr::mutate(Percent_results = 100 * Results / sum(Results)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Source_Origin, Metal, desc(Results))

primary_conc_type_audit <- primary96 %>%
  dplyr::group_by(Source_Origin, Metal, Conc_Type) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Source_Origin, Metal, desc(Results))

primary_chemical_name_audit <- primary96 %>%
  dplyr::group_by(Source_Origin, Metal, Chemical_Name) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Source_Origin, Metal, desc(Results))

primary_basis_audit <- primary96 %>%
  dplyr::group_by(Source_Origin, Metal, Basis_Interpretation) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Source_Origin, Metal, desc(Results))

# Protocol-specific sensitivities requiring cross-source equivalent
# exposure-design or chemical-analysis fields are not fitted. The available
# ECOTOX and WoS fields are not harmonized for those purposes, so imposing
# those filters would create an unsupported comparison.

write_csv(primary_exposure_audit_by_source, file.path(sensitivity_dir, "protocol_audit_exposure_type_by_source.csv"))
write_csv(primary_conc_type_audit, file.path(sensitivity_dir, "protocol_audit_concentration_type.csv"))
write_csv(primary_chemical_name_audit, file.path(sensitivity_dir, "protocol_audit_chemical_name.csv"))
write_csv(primary_basis_audit, file.path(sensitivity_dir, "protocol_audit_basis_interpretation.csv"))


# ============================================================
# 18. ROBUSTNESS SYNTHESIS: POINT-ESTIMATE VARIANTS
# ============================================================
# Protocol-only sensitivities are deliberately absent until their cross-source
# metadata are harmonized. ECOTOX-only is included as augmentation sensitivity.

robustness_point_estimates <- bind_rows(
  make_variant_table(
    metal_pairs_ratio,
    "Combined Reference + Species random intercepts",
    "Primary model"
  ),
  make_variant_table(
    ecotox_only_ratios,
    "ECOTOX-only Reference + Species random intercepts",
    "Augmentation sensitivity"
  ),
  make_variant_table(
    ref_FE_ratios,
    "Reference fixed effects + Species random intercept",
    "Source-structure sensitivity"
  ),
  make_variant_table(
    ref_only_ratios,
    "Reference-only random intercept",
    "Model-structure sensitivity"
  ),
  make_variant_table(
    collapsed_ratios,
    "Reference x Species x Metal cell collapse",
    "Weighting sensitivity"
  )
) %>%
  dplyr::mutate(
    Direction = case_when(
      LC50_ratio < 1 ~ "<1",
      LC50_ratio > 1 ~ ">1",
      TRUE ~ "=1"
    ),
    CI_excludes_1 = ratio_upper < 1 | ratio_lower > 1,
    Interval_level = 0.95,
    Interval_adjustment = "Tukey: three metal contrasts within each model",
    Across_variants_adjustment = "None; models are sensitivity checks, not independent confirmations"
  )

write_csv(
  robustness_point_estimates,
  file.path(synthesis_dir, "robustness_point_estimate_synthesis.csv")
)


# ============================================================
# 19. ROBUSTNESS SYNTHESIS: LORO + LOSO RANGES
# ============================================================

leave_one_out_range_synthesis <- bind_rows(
  loro_summary %>%
    dplyr::transmute(
      Check = "LORO",
      contrast,
      n_refits,
      min_ratio,
      max_ratio,
      min_lower,
      max_upper,
      ratios_below_1,
      ratios_above_1,
      all_point_estimates_below_1,
      all_point_estimates_above_1,
      singular_models,
      convergence_failures
    ),
  loso_summary %>%
    dplyr::transmute(
      Check = "LOSO",
      contrast,
      n_refits,
      min_ratio,
      max_ratio,
      min_lower,
      max_upper,
      ratios_below_1,
      ratios_above_1,
      all_point_estimates_below_1,
      all_point_estimates_above_1,
      singular_models,
      convergence_failures
    )
)

write_csv(
  leave_one_out_range_synthesis,
  file.path(synthesis_dir, "leave_one_out_range_synthesis.csv")
)


# ============================================================
# 20. MODEL-STRUCTURE SYNTHESIS
# ============================================================

primary_ref_var <- extract_variance(primary_variance_components, "Reference_ID")
primary_residual_var <- extract_variance(primary_variance_components, "Residual")
ref_only_ref_var <- extract_variance(ref_only_variance_components, "Reference_ID")
ref_only_species_var <- extract_variance(ref_only_variance_components, "Species")
ref_only_residual_var <- extract_variance(ref_only_variance_components, "Residual")

model_structure_variance_summary <- tibble(
  Model = c("Reference + Species RI", "Reference RI"),
  Reference_variance = c(primary_ref_var, ref_only_ref_var),
  Species_variance = c(extract_variance(primary_variance_components, "Species"), ref_only_species_var),
  Residual_variance = c(primary_residual_var, ref_only_residual_var),
  Singular = c(primary_singular, ref_only_singular),
  Convergence_clean = c(primary_convergence_clean, ref_only_convergence_clean)
) %>%
  dplyr::mutate(
    Total_variance = rowSums(
      cbind(
        replace_na(Reference_variance, 0),
        replace_na(Species_variance, 0),
        replace_na(Residual_variance, 0)
      )
    ),
    Reference_variance_share = Reference_variance / Total_variance,
    Species_variance_share = Species_variance / Total_variance,
    Residual_variance_share = Residual_variance / Total_variance
  )

write_csv(
  model_structure_variance_summary,
  file.path(synthesis_dir, "model_structure_variance_summary.csv")
)

write_csv(
  primary_vs_ref_only,
  file.path(synthesis_dir, "primary_vs_reference_only_contrast_summary.csv")
)


# ============================================================
# 21. CONCLUSION-STABILITY SUMMARY
# ============================================================

point_direction_summary <- robustness_point_estimates %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    Variants = n(),
    Min_ratio = min(LC50_ratio),
    Max_ratio = max(LC50_ratio),
    All_below_1 = all(LC50_ratio < 1),
    All_above_1 = all(LC50_ratio > 1),
    Any_direction_change_across_variants = !(All_below_1 | All_above_1),
    Variants_with_CI_excluding_1 = sum(CI_excludes_1),
    .groups = "drop"
  )

conclusion_stability_summary <- metal_pairs_ratio %>%
  dplyr::select(
    contrast,
    Primary_ratio = LC50_ratio,
    Primary_lower = ratio_lower,
    Primary_upper = ratio_upper
  ) %>%
  dplyr::left_join(point_direction_summary, by = "contrast") %>%
  dplyr::left_join(
    loro_summary %>%
      dplyr::select(
        contrast,
        LORO_min = min_ratio,
        LORO_max = max_ratio,
        LORO_all_below_1 = all_point_estimates_below_1,
        LORO_all_above_1 = all_point_estimates_above_1
      ),
    by = "contrast"
  ) %>%
  dplyr::left_join(
    loso_summary %>%
      dplyr::select(
        contrast,
        LOSO_min = min_ratio,
        LOSO_max = max_ratio,
        LOSO_all_below_1 = all_point_estimates_below_1,
        LOSO_all_above_1 = all_point_estimates_above_1
      ),
    by = "contrast"
  )

write_csv(
  conclusion_stability_summary,
  file.path(synthesis_dir, "conclusion_stability_summary.csv")
)


# ============================================================
# 22. TEMPERATURE / SALINITY SENSITIVITY ON IDENTICAL RECORDS
# ============================================================
# Compare FULL -> COMPLETE SAME FORMULA -> COMPLETE ADJUSTED.
# This separates record selection from adding environmental covariates.
# Keep the primary Reference random intercept and Tukey contrasts unchanged.
# Linear common slopes; no causal claim or metal-specific slope claim.
# pH is outside this specific control. No missing values are imputed.
# Reuses lmer and emmeans settings already used in this script.

if (!all(c("Temperature", "Salinity") %in% names(primary96))) {
  stop("Temperature/Salinity missing from script 03 support-audit objects; check the input version.")
}
stopifnot(is.numeric(primary96$Temperature), is.numeric(primary96$Salinity))
metal_env_membership <- primary96 %>% dplyr::mutate(
  Environment_complete = is.finite(Temperature) & is.finite(Salinity))
write_csv(metal_env_membership,
          file.path(sensitivity_dir, "metal_environment_record_membership.csv"))
metal_env_support <- metal_env_membership %>%
  dplyr::group_by(Environment_complete, Metal) %>% dplyr::summarise(
    Results = n(), References = dplyr::n_distinct(Reference_ID),
    Species = dplyr::n_distinct(Species), .groups = "drop")
write_csv(metal_env_support,
          file.path(sensitivity_dir, "metal_environment_complete_case_support.csv"))
metal_env_cc <- metal_env_membership %>% dplyr::filter(Environment_complete) %>% droplevels()
metal_env_models <- list(PRIMARY_FULL_DATA = m_metal_ref_species)
metal_env_checks <- list()
metal_env_comparisons <- list()
metal_env_status <- "INSUFFICIENT_SUPPORT_OR_NO_COVARIATE_VARIATION"
metal_env_gate <- nrow(metal_env_cc) > 5 &&
  dplyr::n_distinct(metal_env_cc$Metal) == 3 &&
  dplyr::n_distinct(metal_env_cc$Reference_ID) >= 2 &&
  dplyr::n_distinct(metal_env_cc$Reference_ID) < nrow(metal_env_cc) &&
  sd(metal_env_cc$Temperature) > 0 && sd(metal_env_cc$Salinity) > 0
if (metal_env_gate) {
  metal_env_scaling <- tibble(Variable = c("Temperature", "Salinity"),
                              Center = c(mean(metal_env_cc$Temperature), mean(metal_env_cc$Salinity)),
                              Scale = c(sd(metal_env_cc$Temperature), sd(metal_env_cc$Salinity)))
  write_csv(metal_env_scaling, file.path(sensitivity_dir, "metal_environment_scaling.csv"))
  metal_env_cc <- metal_env_cc %>% dplyr::mutate(
    Temperature_z = (Temperature-mean(Temperature))/sd(Temperature),
    Salinity_z = (Salinity-mean(Salinity))/sd(Salinity))
  env_x <- model.matrix(~ Metal + Temperature_z + Salinity_z, data = metal_env_cc)
  metal_env_gate <- qr(env_x)$rank == ncol(env_x)
  if (!metal_env_gate) metal_env_status <- "RANK_DEFICIENT_FIXED_EFFECTS"
}
if (metal_env_gate) {
  env_formulas <- list(
    TEMP_SAL_COMPLETE_CASE_SAME_FORMULA = ln_LC50_umol_L ~ Metal + (1 | Reference_ID) + (1 | Species),
    TEMP_SAL_COMPLETE_CASE_ADJUSTED = ln_LC50_umol_L ~ Metal + Temperature_z + Salinity_z + (1 | Reference_ID) + (1 | Species))
  for (label in names(env_formulas)) {
    result <- tryCatch(lmer(env_formulas[[label]], data = metal_env_cc,
                            REML = TRUE, na.action = na.fail), error = function(e) e)
    if (inherits(result, "error")) {
      metal_env_checks[[length(metal_env_checks)+1L]] <- tibble(
        Model = label, Results = nrow(metal_env_cc), References = dplyr::n_distinct(metal_env_cc$Reference_ID),
        Species = dplyr::n_distinct(metal_env_cc$Species), Singular = NA, Convergence_clean = FALSE,
        Status = "FIT_ERROR", Message = conditionMessage(result))
    } else {
      stopifnot(nobs(result) == nrow(metal_env_cc))
      metal_env_models[[label]] <- result
      writeLines(capture.output(summary(result, ddf = "Kenward-Roger")),
                 file.path(sensitivity_dir, paste0("metal_environment_", label, "_summary.txt")))
    }
  }
  metal_env_status <- if (length(metal_env_models) == 3) "FITTED_REVIEW_DIAGNOSTICS" else "PARTIAL_FIT_FAILURE"
}
for (label in names(metal_env_models)) {
  fit <- metal_env_models[[label]]
  dat_env <- if (label == "PRIMARY_FULL_DATA") primary96 else metal_env_cc
  msgs <- fit@optinfo$conv$lme4$messages
  sing <- isSingular(fit, tol = 1e-4)
  metal_env_checks[[length(metal_env_checks)+1L]] <- tibble(Model = label,
                                                            Results = nobs(fit), References = dplyr::n_distinct(dat_env$Reference_ID),
                                                            Species = dplyr::n_distinct(dat_env$Species), Singular = sing,
                                                            Convergence_clean = model_convergence_clean(fit), Status = "FITTED",
                                                            Message = if (is.null(msgs)) NA_character_ else paste(msgs, collapse = " | "))
  ratios_env <- tryCatch(extract_pairwise_ratios(fit), error = function(e) e)
  if (inherits(ratios_env, "error")) {
    warning(paste("Environmental contrast extraction failed:", label, conditionMessage(ratios_env)))
    ratios_env <- tibble(contrast = c("Cu - Cd", "Cu - Zn", "Cd - Zn"),
                         LC50_ratio = NA_real_, ratio_lower = NA_real_, ratio_upper = NA_real_, p.value = NA_real_)
    contrast_status <- "CONTRAST_ERROR"
  } else contrast_status <- "COMPUTED"
  metal_env_comparisons[[length(metal_env_comparisons)+1L]] <- ratios_env %>%
    dplyr::transmute(Model = label, Results = nobs(fit), contrast, LC50_ratio,
                     ratio_lower, ratio_upper, p.value, Singular = sing,
                     Convergence_clean = model_convergence_clean(fit), Contrast_status = contrast_status)
  write_csv(tibble(Record = as.character(dat_env$Result_ID),
                   Fitted_log_LC50 = fitted(fit), Residual_log_LC50 = residuals(fit)),
            file.path(sensitivity_dir, paste0("metal_environment_", label, "_residuals.csv")))
}
metal_environmental_sensitivity_comparison <- bind_rows(metal_env_comparisons)
metal_environment_fit_checks <- bind_rows(metal_env_checks)
write_csv(metal_environmental_sensitivity_comparison,
          file.path(sensitivity_dir, "metal_environmental_sensitivity_comparison.csv"))
write_csv(metal_environment_fit_checks,
          file.path(sensitivity_dir, "metal_environment_fit_checks.csv"))
write_csv(tibble(Status = metal_env_status),
          file.path(sensitivity_dir, "metal_environmental_sensitivity_status.csv"))
cat("\n--- METAL TEMPERATURE / SALINITY CONTROL ---\n")
print(metal_environmental_sensitivity_comparison)
print(metal_environment_fit_checks)

# ============================================================
# 23. SAVE ANALYSIS OBJECTS
# ============================================================

analysis_objects <- list(
  model_revision = "repo_04_pooled_metal_model_v1",
  reported_formula = paste(deparse(formula(m_metal_ref_species)), collapse = " "),
  validation_note = "Numerical fit gate passed; residual and grouping audits require review.",
  primary_model_data = primary96,
  metal_environment_models = metal_env_models,
  metal_environment_complete_data = metal_env_cc,
  metal_environment_comparison = metal_environmental_sensitivity_comparison,
  metal_environment_fit_checks = metal_environment_fit_checks,
  primary_support = primary_support,
  grouping_structure_summary = grouping_structure_summary,
  reference_structure = reference_structure,
  species_structure = species_structure,
  primary_model = m_metal_ref_species,
  global_metal_test = global_metal_test,
  metal_emmeans = emm_metal,
  metal_emmeans_original_scale = emm_metal_original_scale,
  metal_pairwise_ratios = metal_pairs_ratio,
  pre_augmentation_vs_combined_ratio_shift = pre_augmentation_vs_combined_ratio_shift,
  leave_one_reference_out = loro_results,
  leave_one_species_out = loso_results,
  reference_fixed_effects_model = m_metal_ref_FE,
  reference_fixed_effects_fit_checks = reference_fixed_effects_fit_checks,
  deletion_direction_and_interval_summary = deletion_direction_and_interval_summary,
  reference_only_model = m_metal_ref_only,
  primary_vs_reference_only = primary_vs_ref_only,
  collapsed_model_data = primary96_collapsed,
  collapsed_model = m_metal_collapsed,
  ecotox_only_model = m_metal_ecotox_only,
  combined_vs_ecotox_only = combined_vs_ecotox_only,
  robustness_point_estimates = robustness_point_estimates,
  leave_one_out_range_synthesis = leave_one_out_range_synthesis,
  model_structure_variance_summary = model_structure_variance_summary,
  conclusion_stability_summary = conclusion_stability_summary
)

saveRDS(
  analysis_objects,
  file.path(object_dir, "combined_pooled_analysis_objects.rds")
)

# Compatibility copy for later revised figure scripts.
saveRDS(
  analysis_objects,
  file.path(model_dir, "primary_pooled_analysis_objects.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)

input_output_manifest <- tibble::tibble(
  Type = c(
    "Input",
    "Primary output directory",
    "Reusable analysis object",
    "Primary ratio table",
    "Primary model estimates"
  ),
  Path = c(
    audit_object_file,
    output_root,
    file.path(object_dir, "combined_pooled_analysis_objects.rds"),
    file.path(model_dir, "metal_pairwise_lc50_ratios.csv"),
    file.path(model_dir, "metal_model_estimates_original_scale.csv")
  )
)

readr::write_csv(
  input_output_manifest,
  file.path(output_root, "input_output_manifest.csv")
)


# ============================================================
# 24. FINAL QA + CONSOLE SUMMARY
# ============================================================

cat("\n\n")
cat("============================================\n")
cat("COMBINED PRIMARY POOLED METAL ANALYSIS SUMMARY\n")
cat("============================================\n")

cat("\n--- PRIMARY MODEL DATASET ---\n")
print_full(primary_support)

cat("\n--- GROUPING STRUCTURE ---\n")
print_full(grouping_structure_summary)

cat("\n--- PRIMARY MODEL FIT CHECKS ---\n")
print_full(primary_fit_checks)

cat("\n--- GLOBAL METAL TEST ---\n")
print_full(global_metal_test)

cat("\n--- PRIMARY PAIRWISE LC50 RATIOS ---\n")
print_full(metal_pairs_ratio)

cat("\n--- SUBMITTED VS COMBINED RATIO SHIFT ---\n")
print_full(pre_augmentation_vs_combined_ratio_shift)

cat("\n--- COMBINED VS ECOTOX-ONLY ---\n")
print_full(combined_vs_ecotox_only)

cat("\n--- LORO SUMMARY ---\n")
print_full(loro_summary)

cat("\n--- LOSO SUMMARY ---\n")
print_full(loso_summary)

cat("\n--- REFERENCE-ONLY RANDOM-INTERCEPT FIT ---\n")
print_full(ref_only_fit_checks)

cat("\n--- REPORTED VS REFERENCE-ONLY CONTRASTS ---\n")
print_full(primary_vs_ref_only)

cat("\n--- PRIMARY VS CELL-COLLAPSED ---\n")
print_full(weighting_comparison)

cat("\n--- POINT-ESTIMATE ROBUSTNESS SYNTHESIS ---\n")
print_full(robustness_point_estimates)

cat("\n--- LORO / LOSO RANGE SYNTHESIS ---\n")
print_full(leave_one_out_range_synthesis)

cat("\n--- MODEL-STRUCTURE VARIANCE SUMMARY ---\n")
print_full(model_structure_variance_summary)

cat("\n--- CONCLUSION-STABILITY SUMMARY ---\n")
print_full(conclusion_stability_summary)

required_outputs <- c(
  file.path(robustness_dir, "reference_fixed_effects_fit_checks.csv"),
  file.path(robustness_dir, "deletion_direction_and_interval_summary.csv"),
  file.path(model_dir, "metal_pairwise_lc50_ratios.csv"),
  file.path(model_dir, "metal_model_estimates_original_scale.csv"),
  file.path(model_dir, "pre_augmentation_vs_combined_ratio_shift.csv"),
  file.path(model_dir, "primary_pooled_pairwise_lc50_ratios.png"),
  file.path(model_dir, "primary_pooled_raw_plus_model_estimates.png"),
  file.path(diagnostic_dir, "residuals_vs_fitted.png"),
  file.path(diagnostic_dir, "residuals_normal_qq.png"),
  file.path(diagnostic_dir, "reference_random_intercepts_normal_qq.png"),
  file.path(robustness_dir, "leave_one_reference_out_summary.csv"),
  file.path(robustness_dir, "leave_one_species_out_summary.csv"),
  file.path(robustness_dir, "primary_vs_reference_only_random_intercepts.csv"),
  file.path(robustness_dir, "primary_vs_cell_collapsed_comparison.csv"),
  file.path(sensitivity_dir, "combined_vs_ecotox_only_lc50_ratios.csv"),
  file.path(synthesis_dir, "robustness_point_estimate_synthesis.csv"),
  file.path(synthesis_dir, "leave_one_out_range_synthesis.csv"),
  file.path(synthesis_dir, "model_structure_variance_summary.csv"),
  file.path(synthesis_dir, "conclusion_stability_summary.csv"),
  file.path(object_dir, "combined_pooled_analysis_objects.rds")
)

stopifnot(all(file.exists(required_outputs)))

writeLines("04_fit_pooled_metal_model_PASS", file.path(output_root, "RUN_COMPLETE.txt"))
cat("\nSTATUS: 04 FIT POOLED METAL MODEL = PASS\n")
cat("Review residual, grouping, robustness and sensitivity outputs before interpretation.\n")

# ============================================================
# END
# ============================================================
