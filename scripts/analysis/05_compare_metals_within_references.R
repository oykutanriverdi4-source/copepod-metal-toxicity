# ============================================================
# ESS THESIS REPRODUCIBILITY WORKFLOW
# 05_compare_metals_within_references.R
# ============================================================
# PURPOSE
# Reproduce the source-verified within-Reference metal comparisons
# reported in Section 3.7.1 / 4.2.2 of the thesis.
#
# SCIENTIFIC DESIGN
# - Comparisons are constructed only from manually source-verified
#   rows retained in the frozen adjudication ledgers.
# - Matching is within Reference and source-defined context.
# - Exposure duration is held constant within each comparison, but
#   eligible comparisons can come from the full acute 24-96 h range.
# - Repeated Results within a context x metal cell are collapsed by
#   geometric mean before pairwise comparison.
# - Context-level log LC50 ratios are averaged within Reference;
#   References are then given equal weight.
# - Pointwise 95% Student-t intervals describe uncertainty in the
#   mean Reference-level log ratio where estimable.
# - Exact two-sided sign tests describe directional consistency
#   across References; Holm adjustment uses the full family of
#   21 planned metal-pair hypotheses separately within each scope.
#
# PRIMARY / SUPPORTING INTERPRETATION
# - Core: Cu-Cd, Cu-Zn, Cd-Zn
# - Supporting: Cu-Hg, Cd-Hg
# - Other verified pairs remain in the full 21-pair accounting and
#   reproducibility outputs but are not promoted to separate thesis
#   conclusions when evidence is sparse.
#
# IMPORTANT
# This script reads the reviewed source-verification ledgers as
# immutable inputs. It never rewrites or updates those decisions.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

options(width = 220)

# ============================================================
# 0. PATHS
# ============================================================

combined_file <- file.path(
  "data", "processed", "combined_ecotox_wos_harmonized.csv"
)

adjudication_file <- file.path(
  "data", "curated", "source_verification",
  "combined_source_adjudication_rows.csv"
)

verified_file <- file.path(
  "data", "curated", "source_verification",
  "combined_source_verified_rows_used.csv"
)

output_root <- file.path("outputs", "05_compare_metals_within_references")
table_dir   <- file.path(output_root, "01_tables")
figure_dir  <- file.path(output_root, "02_figures")
object_dir  <- file.path(output_root, "03_objects")

for (d in c(output_root, table_dir, figure_dir, object_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

required_inputs <- c(combined_file, adjudication_file, verified_file)
if (!all(file.exists(required_inputs))) {
  stop(
    "Missing required input file(s):\n",
    paste(required_inputs[!file.exists(required_inputs)], collapse = "\n")
  )
}

ledger_files <- c(adjudication_file, verified_file)
ledger_hashes_before <- unname(tools::md5sum(ledger_files))

# ============================================================
# 1. LOAD FINAL COMBINED DATASET
# ============================================================

combined <- read_csv(combined_file, show_col_types = FALSE)

required_columns <- c(
  "Reference_Number", "Test_Number", "Result_Number",
  "Reference_ID", "Test_ID", "Result_ID", "Source_Origin",
  "Metal", "Species", "Order", "Lifestage", "Duration_days",
  "LC50_umol_L", "ln_LC50_umol_L", "Harmonization_Status",
  "External_Stage_Sex", "External_Source_Location", "External_Source_Note"
)

missing_columns <- setdiff(required_columns, names(combined))
if (length(missing_columns) > 0) {
  stop(
    "Combined input is missing required columns: ",
    paste(missing_columns, collapse = " | ")
  )
}

stopifnot(
  nrow(combined) == 353,
  n_distinct(combined$Reference_ID) == 84,
  n_distinct(combined$Result_ID) == 353,
  sum(combined$Harmonization_Status == "HARMONIZED", na.rm = TRUE) == 304
)

# ============================================================
# 2. LOAD AND VERIFY FROZEN SOURCE-DECISION LEDGERS
# ============================================================

audited <- read_csv(adjudication_file, show_col_types = FALSE)
verified_rows <- read_csv(verified_file, show_col_types = FALSE)

decision_columns <- c(
  "Result_ID", "Source_Context_ID", "Verification_Status",
  "Final_Pairwise_Use", "Decision_Note"
)

for (x in list(audited, verified_rows)) {
  if (!all(c(required_columns, decision_columns) %in% names(x))) {
    stop("A source-verification ledger is missing required columns.")
  }
  if (anyNA(x$Result_ID) || anyDuplicated(x$Result_ID)) {
    stop("Missing or duplicated Result_ID in a source-verification ledger.")
  }
}

accepted <- audited %>% filter(Final_Pairwise_Use == "YES")

if (!setequal(accepted$Result_ID, verified_rows$Result_ID)) {
  stop("Adjudication and accepted-row ledgers disagree on retained Result_ID values.")
}

same_field <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    same <- (is.na(a) & is.na(b)) |
      (!is.na(a) & !is.na(b) &
         abs(a - b) <= 1e-10 * pmax(1, abs(a), abs(b)))
  } else {
    a <- as.character(a)
    b <- as.character(b)
    same <- (is.na(a) & is.na(b)) |
      (!is.na(a) & !is.na(b) & a == b)
  }
  isTRUE(all(same))
}

if (!setequal(names(accepted), names(verified_rows))) {
  stop("Adjudication and accepted-row ledger schemas differ.")
}

accepted <- accepted[match(verified_rows$Result_ID, accepted$Result_ID), ]

for (f in names(verified_rows)) {
  if (!same_field(accepted[[f]], verified_rows[[f]])) {
    stop("Accepted-row ledger differs from adjudication in field: ", f)
  }
}

# Verify retained analysis rows against the canonical combined dataset.
ix <- match(verified_rows$Result_ID, combined$Result_ID)
if (anyNA(ix)) {
  stop("At least one accepted Result_ID is absent from the combined dataset.")
}

for (f in names(combined)) {
  if (!f %in% names(verified_rows) ||
      !same_field(verified_rows[[f]], combined[[f]][ix])) {
    stop("Accepted rows disagree with the combined dataset in field: ", f)
  }
}

stopifnot(
  nrow(verified_rows) == 70,
  n_distinct(verified_rows$Reference_ID) == 14,
  nrow(distinct(verified_rows, Reference_ID, Source_Context_ID)) == 29,
  sum(verified_rows$Source_Origin == "ECOTOX") == 61,
  sum(verified_rows$Source_Origin == "WoS_supplemental") == 9,
  all(verified_rows$Harmonization_Status == "HARMONIZED"),
  all(is.finite(verified_rows$LC50_umol_L)),
  all(verified_rows$LC50_umol_L > 0),
  all(!is.na(verified_rows$Source_Context_ID)),
  all(nzchar(verified_rows$Source_Context_ID))
)

# Compact copies retained in the saved RDS for provenance/inspection.
source_audit_key <- audited %>%
  filter(Source_Origin == "ECOTOX") %>%
  select(
    Reference_Number, Test_Number, Result_Number, Source_Context_ID,
    Verification_Status, Final_Pairwise_Use, Decision_Note
  )

wos_source_audit_key <- audited %>%
  filter(Source_Origin == "WoS_supplemental") %>%
  select(
    Reference_ID, Test_ID, Metal, Duration_days, Expected_Stage_Sex,
    Source_Context_ID, Verification_Status, Final_Pairwise_Use, Decision_Note
  )

# ============================================================
# 3. COLLAPSE VERIFIED CONTEXT x METAL CELLS
# ============================================================

context_metal_verified <- verified_rows %>%
  group_by(Reference_ID, Source_Context_ID, Metal) %>%
  summarise(
    Source_Origin = first(Source_Origin),
    LC50_context_umol_L = exp(mean(log(LC50_umol_L))),
    n_results = n(),
    n_tests = n_distinct(Test_ID),
    .groups = "drop"
  ) %>%
  mutate(
    ECOTOX_Reference_Number = if_else(
      Source_Origin == "ECOTOX",
      as.numeric(str_remove(
        if_else(Source_Origin == "ECOTOX", Reference_ID, NA_character_),
        "^E_"
      )),
      NA_real_
    )
  )

multiplicity_QA <- context_metal_verified %>%
  filter(n_results > 1)

context_metadata <- verified_rows %>%
  mutate(
    Biological_Context = case_when(
      Source_Origin == "WoS_supplemental" & !is.na(External_Stage_Sex) ~ External_Stage_Sex,
      TRUE ~ Lifestage
    )
  ) %>%
  group_by(Reference_ID, Source_Context_ID) %>%
  summarise(
    Source_Origin = first(Source_Origin),
    n_species = n_distinct(Species),
    n_durations = n_distinct(Duration_days),
    Species = first(Species),
    Duration_days = first(Duration_days),
    Lifestage = paste(sort(unique(na.omit(Lifestage))), collapse = " | "),
    Biological_Context = paste(sort(unique(na.omit(Biological_Context))), collapse = " | "),
    .groups = "drop"
  )

stopifnot(
  all(context_metadata$n_species == 1),
  all(context_metadata$n_durations == 1)
)

context_inventory <- context_metal_verified %>%
  group_by(Reference_ID, Source_Context_ID, Source_Origin) %>%
  summarise(
    n_metals = n_distinct(Metal),
    Metals = paste(sort(unique(Metal)), collapse = " | "),
    Results = sum(n_results),
    .groups = "drop"
  ) %>%
  left_join(
    context_metadata %>%
      select(
        Reference_ID, Source_Context_ID, Species,
        Duration_days, Biological_Context
      ),
    by = c("Reference_ID", "Source_Context_ID")
  ) %>%
  arrange(desc(n_metals), Reference_ID, Source_Context_ID)

write_csv(
  context_metal_verified,
  file.path(table_dir, "verified_context_metal_values.csv")
)
write_csv(
  multiplicity_QA,
  file.path(table_dir, "verified_context_multiplicity_QA.csv")
)
write_csv(
  context_inventory,
  file.path(table_dir, "verified_context_inventory.csv")
)
write_csv(
  context_inventory %>% filter(Source_Origin == "WoS_supplemental"),
  file.path(table_dir, "wos_verified_matched_contexts.csv")
)

# ============================================================
# 4. BUILD WIDE VERIFIED-CONTEXT TABLE
# ============================================================

wide_verified <- context_metal_verified %>%
  select(
    Reference_ID, Source_Context_ID, Source_Origin,
    Metal, LC50_context_umol_L
  ) %>%
  pivot_wider(
    names_from = Metal,
    values_from = LC50_context_umol_L
  ) %>%
  left_join(
    context_metadata %>%
      select(
        Reference_ID, Source_Context_ID, Species,
        Duration_days, Lifestage, Biological_Context
      ),
    by = c("Reference_ID", "Source_Context_ID")
  )

write_csv(
  wide_verified,
  file.path(table_dir, "verified_contexts_wide.csv")
)

# ============================================================
# 5. CONSTRUCT ALL VERIFIED METAL PAIRS
# ============================================================

make_pair <- function(wide_data, metal_A, metal_B) {
  if (!(metal_A %in% names(wide_data)) || !(metal_B %in% names(wide_data))) {
    return(tibble())
  }

  wide_data %>%
    filter(!is.na(.data[[metal_A]]), !is.na(.data[[metal_B]])) %>%
    transmute(
      Reference_ID,
      Source_Context_ID,
      Source_Origin,
      Species,
      Duration_days,
      Lifestage,
      Biological_Context,
      Comparison = paste0(metal_A, "-", metal_B),
      Metal_A = metal_A,
      Metal_B = metal_B,
      LC50_A = .data[[metal_A]],
      LC50_B = .data[[metal_B]],
      LC50_ratio = LC50_A / LC50_B,
      ln_ratio = log(LC50_ratio)
    )
}

metal_order <- c("Cu", "Cd", "Zn", "Ni", "Hg", "Ag", "Pb")
available_metals <- metal_order[metal_order %in% names(wide_verified)]

pair_specs <- combn(available_metals, 2, simplify = FALSE)

all_pair_data <- map_dfr(
  pair_specs,
  ~ make_pair(wide_verified, .x[1], .x[2])
) %>%
  mutate(
    Analysis_Role = case_when(
      Comparison %in% c("Cu-Cd", "Cu-Zn", "Cd-Zn") ~ "CORE",
      Comparison %in% c("Cu-Hg", "Cd-Hg") ~ "SUPPORTING",
      TRUE ~ "OTHER_VERIFIED"
    )
  )

pair_support_inventory <- all_pair_data %>%
  group_by(Comparison, Analysis_Role) %>%
  summarise(
    Matched_contexts = n(),
    References = n_distinct(Reference_ID),
    Species = n_distinct(Species),
    ECOTOX_contexts = sum(Source_Origin == "ECOTOX"),
    WoS_contexts = sum(Source_Origin == "WoS_supplemental"),
    .groups = "drop"
  ) %>%
  arrange(desc(Matched_contexts), Comparison)

write_csv(
  all_pair_data,
  file.path(table_dir, "all_verified_pair_contexts_long.csv")
)
write_csv(
  pair_support_inventory,
  file.path(table_dir, "all_verified_pair_support_inventory.csv")
)

get_pair <- function(label) {
  all_pair_data %>% filter(Comparison == label)
}

# ============================================================
# 6. BENCHMARK QA: PRESERVE PRE-AUGMENTATION ECOTOX EVIDENCE
# ============================================================

pair_summary <- function(dat, comparison, scope) {
  dat %>%
    filter(Comparison == comparison) %>%
    summarise(
      Comparison = comparison,
      Scope = scope,
      Matched_contexts = n(),
      References = n_distinct(Reference_ID),
      Species = n_distinct(Species),
      .groups = "drop"
    )
}

ecotox_only_pairs <- all_pair_data %>%
  filter(Source_Origin == "ECOTOX")

ecotox_pairwise_QA <- bind_rows(
  pair_summary(ecotox_only_pairs, "Cu-Cd", "ECOTOX_FROZEN"),
  pair_summary(ecotox_only_pairs, "Cd-Zn", "ECOTOX_FROZEN"),
  pair_summary(ecotox_only_pairs, "Cu-Zn", "ECOTOX_FROZEN"),
  pair_summary(ecotox_only_pairs, "Cu-Hg", "ECOTOX_FROZEN"),
  pair_summary(ecotox_only_pairs, "Cd-Hg", "ECOTOX_FROZEN")
)


stopifnot(
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Cd") %>% pull(Matched_contexts) == 19,
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Cd") %>% pull(References) == 9,
  ecotox_pairwise_QA %>% filter(Comparison == "Cd-Zn") %>% pull(Matched_contexts) == 7,
  ecotox_pairwise_QA %>% filter(Comparison == "Cd-Zn") %>% pull(References) == 4,
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Zn") %>% pull(Matched_contexts) == 6,
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Zn") %>% pull(References) == 4,
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Hg") %>% pull(Matched_contexts) == 6,
  ecotox_pairwise_QA %>% filter(Comparison == "Cu-Hg") %>% pull(References) == 4,
  ecotox_pairwise_QA %>% filter(Comparison == "Cd-Hg") %>% pull(Matched_contexts) == 5,
  ecotox_pairwise_QA %>% filter(Comparison == "Cd-Hg") %>% pull(References) == 3
)

combined_key_pair_QA <- bind_rows(
  pair_summary(all_pair_data, "Cu-Cd", "COMBINED"),
  pair_summary(all_pair_data, "Cd-Zn", "COMBINED"),
  pair_summary(all_pair_data, "Cu-Zn", "COMBINED"),
  pair_summary(all_pair_data, "Cu-Hg", "COMBINED"),
  pair_summary(all_pair_data, "Cd-Hg", "COMBINED")
)

stopifnot(
  combined_key_pair_QA %>% filter(Comparison == "Cu-Cd") %>% pull(Matched_contexts) == 21,
  combined_key_pair_QA %>% filter(Comparison == "Cu-Cd") %>% pull(References) == 10,
  combined_key_pair_QA %>% filter(Comparison == "Cd-Zn") %>% pull(Matched_contexts) == 7,
  combined_key_pair_QA %>% filter(Comparison == "Cd-Zn") %>% pull(References) == 4,
  combined_key_pair_QA %>% filter(Comparison == "Cu-Zn") %>% pull(Matched_contexts) == 6,
  combined_key_pair_QA %>% filter(Comparison == "Cu-Zn") %>% pull(References) == 4,
  combined_key_pair_QA %>% filter(Comparison == "Cu-Hg") %>% pull(Matched_contexts) == 6,
  combined_key_pair_QA %>% filter(Comparison == "Cu-Hg") %>% pull(References) == 4,
  combined_key_pair_QA %>% filter(Comparison == "Cd-Hg") %>% pull(Matched_contexts) == 5,
  combined_key_pair_QA %>% filter(Comparison == "Cd-Hg") %>% pull(References) == 3
)

write_csv(
  ecotox_pairwise_QA,
  file.path(table_dir, "matched_metal_pairwise_QA_ecotox_frozen.csv")
)
write_csv(
  combined_key_pair_QA,
  file.path(table_dir, "matched_metal_pairwise_QA_combined.csv")
)

# ============================================================
# 7. DIRECTION SUMMARIES
# ============================================================

summarise_context_direction <- function(dat) {
  dat %>%
    group_by(Comparison, Analysis_Role) %>%
    summarise(
      Matched_contexts = n(),
      References = n_distinct(Reference_ID),
      First_metal_lower_contexts = sum(ln_ratio < 0),
      Second_metal_lower_contexts = sum(ln_ratio > 0),
      Equal_contexts = sum(ln_ratio == 0),
      .groups = "drop"
    )
}

context_direction_summary <- summarise_context_direction(all_pair_data)

context_direction_ecotox <- summarise_context_direction(ecotox_only_pairs) %>%
  mutate(Scope = "ECOTOX_FROZEN")

context_direction_combined <- context_direction_summary %>%
  mutate(Scope = "COMBINED")

reference_direction <- all_pair_data %>%
  group_by(Comparison, Analysis_Role, Reference_ID, Source_Origin) %>%
  summarise(
    Contexts = n(),
    mean_ln_ratio = mean(ln_ratio),
    geometric_mean_ratio = exp(mean_ln_ratio),
    .groups = "drop"
  )

reference_direction_summary <- reference_direction %>%
  group_by(Comparison, Analysis_Role) %>%
  summarise(
    References = n(),
    First_metal_lower_references = sum(mean_ln_ratio < 0),
    Second_metal_lower_references = sum(mean_ln_ratio > 0),
    Equal_references = sum(mean_ln_ratio == 0),
    .groups = "drop"
  )

reference_direction_ecotox <- reference_direction %>%
  filter(Source_Origin == "ECOTOX") %>%
  group_by(Comparison, Analysis_Role) %>%
  summarise(
    References = n(),
    First_metal_lower_references = sum(mean_ln_ratio < 0),
    Second_metal_lower_references = sum(mean_ln_ratio > 0),
    Equal_references = sum(mean_ln_ratio == 0),
    .groups = "drop"
  ) %>%
  mutate(Scope = "ECOTOX_FROZEN")

reference_direction_combined <- reference_direction_summary %>%
  mutate(Scope = "COMBINED")

core_direction_shift <- bind_rows(
  context_direction_ecotox,
  context_direction_combined
) %>%
  filter(Comparison %in% c("Cu-Cd", "Cd-Zn", "Cu-Zn")) %>%
  arrange(Comparison, Scope)

core_reference_direction_shift <- bind_rows(
  reference_direction_ecotox,
  reference_direction_combined
) %>%
  filter(Comparison %in% c("Cu-Cd", "Cd-Zn", "Cu-Zn")) %>%
  arrange(Comparison, Scope)

write_csv(
  context_direction_summary,
  file.path(table_dir, "pairwise_context_direction_summary.csv")
)
write_csv(
  reference_direction,
  file.path(table_dir, "pairwise_reference_level_direction.csv")
)
write_csv(
  reference_direction_summary,
  file.path(table_dir, "pairwise_reference_direction_summary.csv")
)
write_csv(
  core_direction_shift,
  file.path(table_dir, "core_pair_context_direction_shift_ecotox_vs_combined.csv")
)
write_csv(
  core_reference_direction_shift,
  file.path(table_dir, "core_pair_reference_direction_shift_ecotox_vs_combined.csv")
)

# Frozen pre-augmentation regression checks.
stopifnot(
  context_direction_ecotox %>%
    filter(Comparison == "Cu-Cd") %>%
    pull(First_metal_lower_contexts) == 15,
  reference_direction_ecotox %>%
    filter(Comparison == "Cu-Cd") %>%
    pull(First_metal_lower_references) == 7
)

# ============================================================
# 8. REFERENCE-EQUAL METAL EFFECT SUMMARIES
# ============================================================
# Ratio A/B is an LC50 ratio.
# - ratio < 1: first-listed metal has lower LC50
# - ratio > 1: second-listed metal has lower LC50
#
# ALL_DURATIONS and 96H are separate descriptive scopes.
# The 96H scope holds exposure duration constant at 96 h.
#
# Exact sign tests assess direction among non-tied References.
# They do not test effect magnitude. Holm adjustment is applied
# over the full family of 21 planned metal pairs within each scope.

reference_effect_rows <- list()
effect_summary_rows <- list()
loro_rows <- list()

all_metal_specs <- combn(metal_order, 2, simplify = FALSE)

for (scope in c("ALL_DURATIONS", "96H")) {

  scope_data <- all_pair_data
  if (scope == "96H") {
    scope_data <- filter(scope_data, Duration_days == 4)
  }

  for (ab in all_metal_specs) {

    label <- paste(ab, collapse = "-")
    d <- filter(scope_data, Comparison == label)

    r <- d %>%
      group_by(Reference_ID) %>%
      summarise(
        Contexts = n(),
        Mean_log_ratio = mean(ln_ratio),
        Ratio_A_over_B = exp(Mean_log_ratio),
        .groups = "drop"
      ) %>%
      mutate(
        Scope = scope,
        Comparison = label,
        .before = 1
      )

    reference_effect_rows[[length(reference_effect_rows) + 1L]] <- r

    x <- r$Mean_log_ratio
    k <- length(x)

    estimate <- if (k > 0L) exp(mean(x)) else NA_real_
    lo <- NA_real_
    hi <- NA_real_

    interval_status <- if (k < 2L) {
      "FEWER_THAN_TWO_REFERENCES"
    } else {
      "ZERO_VARIANCE"
    }

    if (k >= 2L && is.finite(sd(x)) && sd(x) > 0) {
      margin <- qt(0.975, df = k - 1L) * sd(x) / sqrt(k)
      lo <- exp(mean(x) - margin)
      hi <- exp(mean(x) + margin)
      interval_status <- "POINTWISE_T_INTERVAL_ASSUMPTION_DEPENDENT"
    }

    # Numerical equality tolerance is used only for sign classification.
    negative <- sum(x < -1e-12)
    positive <- sum(x > 1e-12)
    tied <- k - negative - positive
    n_direction <- negative + positive

    sign_p <- if (k >= 2L && n_direction > 0L) {
      binom.test(negative, n_direction, p = 0.5)$p.value
    } else {
      NA_real_
    }

    sign_status <- if (k < 2L) {
      "FEWER_THAN_TWO_REFERENCES"
    } else if (n_direction == 0L) {
      "ALL_TIED"
    } else {
      "EXACT_TWO_SIDED_SIGN_TEST"
    }

    loro_min <- NA_real_
    loro_max <- NA_real_

    if (k >= 2L) {
      omitted_ratios <- vapply(
        seq_len(k),
        function(j) exp(mean(x[-j])),
        numeric(1)
      )

      loro_min <- min(omitted_ratios)
      loro_max <- max(omitted_ratios)

      loro_rows[[length(loro_rows) + 1L]] <- tibble(
        Scope = scope,
        Comparison = label,
        Omitted_Reference = r$Reference_ID,
        Remaining_References = k - 1L,
        Ratio_A_over_B = omitted_ratios
      )
    }

    effect_summary_rows[[length(effect_summary_rows) + 1L]] <- tibble(
      Scope = scope,
      Comparison = label,
      Metal_A = ab[1],
      Metal_B = ab[2],
      Matched_contexts = nrow(d),
      References = k,
      Ratio_A_over_B = estimate,
      CI95_low = lo,
      CI95_high = hi,
      Interval_status = interval_status,
      A_lower_references = negative,
      B_lower_references = positive,
      Tied_references = tied,
      Sign_p = sign_p,
      Sign_status = sign_status,
      LORO_min = loro_min,
      LORO_max = loro_max,
      Support_note = case_when(
        k == 0L ~ "NO_MATCHED_CONTEXT",
        k == 1L ~ "ONE_REFERENCE_DESCRIPTIVE_ONLY",
        TRUE ~ "REPORT_REFERENCE_COUNT_AND_ASSUMPTIONS"
      )
    )
  }
}

metal_reference_effects <- bind_rows(reference_effect_rows)

metal_effect_summary <- bind_rows(effect_summary_rows) %>%
  group_by(Scope) %>%
  mutate(
    Sign_p_Holm = p.adjust(Sign_p, method = "holm", n = 21L)
  ) %>%
  ungroup()

metal_effect_loro <- bind_rows(loro_rows)

stopifnot(
  nrow(metal_effect_summary) == 42L,
  all(
    is.na(metal_effect_summary$Sign_p_Holm) |
      metal_effect_summary$Sign_p_Holm >= metal_effect_summary$Sign_p
  )
)

write_csv(
  metal_reference_effects,
  file.path(table_dir, "metal_reference_effects.csv")
)
write_csv(
  metal_effect_summary,
  file.path(table_dir, "metal_effect_summary.csv")
)
write_csv(
  metal_effect_loro,
  file.path(table_dir, "metal_effect_loro.csv")
)

# ============================================================
# 9. CORE-PAIR SPECIES DECOMPOSITION
# ============================================================

core_species_decomposition <- all_pair_data %>%
  filter(Comparison %in% c("Cu-Cd", "Cd-Zn", "Cu-Zn")) %>%
  group_by(Comparison, Species) %>%
  summarise(
    Matched_contexts = n(),
    References = n_distinct(Reference_ID),
    mean_ln_ratio = mean(ln_ratio),
    geometric_mean_ratio = exp(mean_ln_ratio),
    First_metal_lower = mean_ln_ratio < 0,
    .groups = "drop"
  ) %>%
  arrange(Comparison, Species)

core_species_summary <- core_species_decomposition %>%
  group_by(Comparison) %>%
  summarise(
    Species = n(),
    First_metal_lower_species = sum(First_metal_lower),
    Second_metal_lower_species = sum(!First_metal_lower),
    .groups = "drop"
  )

write_csv(
  core_species_decomposition,
  file.path(table_dir, "core_pair_species_decomposition.csv")
)
write_csv(
  core_species_summary,
  file.path(table_dir, "core_pair_species_summary.csv")
)

# ============================================================
# 10. SAME-CONTEXT MULTI-METAL EVIDENCE
# ============================================================

cu_cd_zn_contexts <- wide_verified %>%
  filter(!is.na(Cu), !is.na(Cd), !is.na(Zn)) %>%
  mutate(
    Cu_Cd_ratio = Cu / Cd,
    Cu_Zn_ratio = Cu / Zn,
    Cd_Zn_ratio = Cd / Zn
  )

cu_cd_zn_96h <- cu_cd_zn_contexts %>%
  filter(Duration_days == 4)

stopifnot(
  nrow(cu_cd_zn_contexts) == 5,
  n_distinct(cu_cd_zn_contexts$Reference_ID) == 3,
  nrow(cu_cd_zn_96h) == 2,
  n_distinct(cu_cd_zn_96h$Reference_ID) == 2
)

write_csv(
  cu_cd_zn_contexts,
  file.path(table_dir, "Cu_Cd_Zn_same_context_all_durations.csv")
)
write_csv(
  cu_cd_zn_96h,
  file.path(table_dir, "Cu_Cd_Zn_same_context_96h.csv")
)

# ============================================================
# 11. THESIS-INTERPRETED SUMMARY TABLE
# ============================================================
# The thesis interprets the three core pairs and the two smaller
# Hg-linked supporting comparisons. Other verified pairs remain in
# the complete 21-pair outputs for reproducibility and Holm accounting.

thesis_pairs <- c("Cu-Cd", "Cu-Zn", "Cu-Hg", "Cd-Zn", "Cd-Hg")

thesis_pair_summary <- metal_effect_summary %>%
  filter(Scope == "ALL_DURATIONS", Comparison %in% thesis_pairs) %>%
  select(
    Comparison,
    Matched_contexts,
    References,
    Ratio_A_over_B,
    CI95_low,
    CI95_high,
    A_lower_references,
    B_lower_references,
    Sign_p,
    Sign_p_Holm
  ) %>%
  left_join(
    context_direction_summary %>%
      select(
        Comparison,
        First_metal_lower_contexts,
        Second_metal_lower_contexts,
        Equal_contexts
      ),
    by = "Comparison"
  ) %>%
  mutate(
    Analytical_role = case_when(
      Comparison %in% c("Cu-Cd", "Cu-Zn", "Cd-Zn") ~ "Core",
      TRUE ~ "Supporting"
    )
  ) %>%
  arrange(match(Comparison, thesis_pairs))

write_csv(
  thesis_pair_summary,
  file.path(table_dir, "thesis_within_reference_metal_summary.csv")
)

# Published/working-thesis benchmark values from the frozen verified dataset.
benchmark <- thesis_pair_summary %>%
  select(Comparison, Matched_contexts, References, Ratio_A_over_B)

stopifnot(
  benchmark %>% filter(Comparison == "Cu-Cd") %>% pull(Matched_contexts) == 21,
  benchmark %>% filter(Comparison == "Cu-Cd") %>% pull(References) == 10,
  benchmark %>% filter(Comparison == "Cu-Zn") %>% pull(Matched_contexts) == 6,
  benchmark %>% filter(Comparison == "Cu-Zn") %>% pull(References) == 4,
  benchmark %>% filter(Comparison == "Cd-Zn") %>% pull(Matched_contexts) == 7,
  benchmark %>% filter(Comparison == "Cd-Zn") %>% pull(References) == 4,
  benchmark %>% filter(Comparison == "Cu-Hg") %>% pull(Matched_contexts) == 6,
  benchmark %>% filter(Comparison == "Cu-Hg") %>% pull(References) == 4,
  benchmark %>% filter(Comparison == "Cd-Hg") %>% pull(Matched_contexts) == 5,
  benchmark %>% filter(Comparison == "Cd-Hg") %>% pull(References) == 3
)

# ============================================================
# 12. ANALYSIS FIGURE
# ============================================================
# This plot is an analysis/reproducibility figure. The polished thesis
# figure is generated later by the dedicated final figure script.

make_composite <- function(dat, comparisons, title, subtitle) {

  plot_data <- dat %>%
    filter(Comparison %in% comparisons) %>%
    mutate(
      Comparison = factor(Comparison, levels = comparisons),
      Reference_label = Reference_ID
    )

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  ranges <- plot_data %>%
    group_by(Comparison, Reference_label) %>%
    summarise(
      n_contexts = n(),
      min_ln_ratio = min(ln_ratio),
      max_ln_ratio = max(ln_ratio),
      .groups = "drop"
    ) %>%
    filter(n_contexts > 1)

  ggplot(plot_data, aes(x = ln_ratio, y = Reference_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
    geom_segment(
      data = ranges,
      aes(
        x = min_ln_ratio,
        xend = max_ln_ratio,
        y = Reference_label,
        yend = Reference_label
      ),
      inherit.aes = FALSE,
      linewidth = 0.6,
      alpha = 0.25
    ) +
    geom_point(aes(shape = Source_Origin), size = 2.8, alpha = 0.84) +
    facet_wrap(~ Comparison, ncol = 1, scales = "free_y") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "ln(LC50 first metal / LC50 second metal)",
      y = "Reference_ID",
      shape = "Source",
      caption = paste0(
        "Negative values indicate lower LC50 for the first-listed metal; ",
        "positive values indicate lower LC50 for the second-listed metal."
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0)
    )
}

core_composite_plot <- make_composite(
  all_pair_data,
  c("Cu-Cd", "Cu-Zn", "Cd-Zn"),
  "Within-Reference core metal comparisons",
  "Source-verified combined ECOTOX + supplementary WoS evidence"
)

if (!is.null(core_composite_plot)) {
  ggsave(
    file.path(figure_dir, "core_within_reference_metal_comparisons.png"),
    core_composite_plot,
    width = 9,
    height = 11,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, "core_within_reference_metal_comparisons.pdf"),
    core_composite_plot,
    width = 9,
    height = 11,
    units = "in"
  )
}

# ============================================================
# 13. SAVE PAIR-SPECIFIC DATASETS
# ============================================================

supported_pair_labels <- pair_support_inventory %>%
  filter(Matched_contexts > 0) %>%
  pull(Comparison)

for (lab in supported_pair_labels) {
  dat <- get_pair(lab)
  if (nrow(dat) == 0) next

  write_csv(
    dat,
    file.path(
      table_dir,
      paste0(str_replace_all(lab, "-", "_"), "_verified_contexts.csv")
    )
  )
}

# ============================================================
# 14. SAVE R OBJECTS
# ============================================================

analysis_objects <- list(
  ecotox_source_audit_key = source_audit_key,
  wos_source_audit_key = wos_source_audit_key,
  audited = audited,
  verified_rows = verified_rows,
  context_metal_verified = context_metal_verified,
  context_inventory = context_inventory,
  wide_verified = wide_verified,
  all_pair_data = all_pair_data,
  pair_support_inventory = pair_support_inventory,
  metal_effect_summary = metal_effect_summary,
  metal_reference_effects = metal_reference_effects,
  metal_effect_loro = metal_effect_loro,
  context_direction_summary = context_direction_summary,
  reference_direction = reference_direction,
  core_species_decomposition = core_species_decomposition,
  cu_cd_zn_contexts = cu_cd_zn_contexts,
  cu_cd_zn_96h = cu_cd_zn_96h,
  thesis_pair_summary = thesis_pair_summary,
  core_composite_plot = core_composite_plot
)

object_file <- file.path(
  object_dir,
  "within_reference_metal_analysis_objects.rds"
)

saveRDS(analysis_objects, object_file)

# ============================================================
# 15. MANIFEST + SESSION INFO
# ============================================================

manifest <- tibble(
  Role = c(
    "Input: canonical combined dataset",
    "Input: source adjudication ledger",
    "Input: accepted source-verified rows",
    "Output: main metal effect summary",
    "Output: thesis-interpreted summary",
    "Output: R analysis objects"
  ),
  Path = c(
    combined_file,
    adjudication_file,
    verified_file,
    file.path(table_dir, "metal_effect_summary.csv"),
    file.path(table_dir, "thesis_within_reference_metal_summary.csv"),
    object_file
  )
) %>%
  mutate(
    Exists = file.exists(Path),
    MD5 = if_else(
      Exists,
      unname(tools::md5sum(Path)),
      NA_character_
    )
  )

write_csv(
  manifest,
  file.path(output_root, "input_output_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)

# Confirm that the frozen source ledgers were not modified.
ledger_hashes_after <- unname(tools::md5sum(ledger_files))
stopifnot(identical(ledger_hashes_before, ledger_hashes_after))

# ============================================================
# 16. FINAL OUTPUT CHECKS + CONSOLE SUMMARY
# ============================================================

required_outputs <- c(
  file.path(table_dir, "metal_effect_summary.csv"),
  file.path(table_dir, "metal_reference_effects.csv"),
  file.path(table_dir, "metal_effect_loro.csv"),
  file.path(table_dir, "all_verified_pair_support_inventory.csv"),
  file.path(table_dir, "thesis_within_reference_metal_summary.csv"),
  file.path(table_dir, "Cu_Cd_Zn_same_context_96h.csv"),
  file.path(figure_dir, "core_within_reference_metal_comparisons.png"),
  object_file,
  file.path(output_root, "input_output_manifest.csv"),
  file.path(output_root, "sessionInfo.txt")
)

stopifnot(all(file.exists(required_outputs)))

cat("\n")
cat("=====================================================\n")
cat("05 WITHIN-REFERENCE METAL COMPARISONS\n")
cat("=====================================================\n")

cat("\n--- SOURCE-VERIFIED INPUT ---\n")
cat("ECOTOX accepted rows:", sum(verified_rows$Source_Origin == "ECOTOX"), "\n")
cat("WoS accepted rows:", sum(verified_rows$Source_Origin == "WoS_supplemental"), "\n")
cat("Total accepted rows:", nrow(verified_rows), "\n")
cat("References:", n_distinct(verified_rows$Reference_ID), "\n")
cat("Source-defined contexts:", nrow(distinct(verified_rows, Reference_ID, Source_Context_ID)), "\n")

cat("\n--- THESIS-INTERPRETED ALL-DURATION PAIRS ---\n")
print(thesis_pair_summary)

cat("\n--- ALL VERIFIED PAIR SUPPORT ---\n")
print(pair_support_inventory)

cat("\n--- SAME-CONTEXT Cu-Cd-Zn SUPPORT ---\n")
cat(
  "All durations:", nrow(cu_cd_zn_contexts), "contexts /",
  n_distinct(cu_cd_zn_contexts$Reference_ID), "References\n"
)
cat(
  "96 h:", nrow(cu_cd_zn_96h), "contexts /",
  n_distinct(cu_cd_zn_96h$Reference_ID), "References\n"
)

writeLines(
  c(
    "STATUS: 05 COMPARE METALS WITHIN REFERENCES = PASS",
    paste0("Accepted source-verified rows: ", nrow(verified_rows)),
    "Frozen source-decision ledgers were read but not modified."
  ),
  con = file.path(output_root, "RUN_COMPLETE.txt")
)

cat("\nSTATUS: 05 COMPARE METALS WITHIN REFERENCES = PASS\n")
cat("Frozen source-decision ledgers preserved.\n")
cat("Core, supporting, and full 21-pair reproducibility outputs created.\n")

# ============================================================
# END
# ============================================================
