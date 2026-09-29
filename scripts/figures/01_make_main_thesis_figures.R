# ============================================================
# 01_make_main_thesis_figures.R
# ============================================================
# Purpose
# Regenerate the reader-facing thesis figures from the frozen outputs
# of the cleaned analysis pipeline. This script does NOT refit models,
# alter eligibility decisions, change contrasts, or recompute the
# scientific analyses.
#
# Analysis sources
# - Script 04: pooled 96-h Cu-Cd-Zn model
# - Script 05: within-Reference metal comparisons
# - Script 06: exposure-duration comparisons
# - Script 07: salinity and temperature evidence
# - Script 08: pH/confounding evidence
# - Script 09: within-Reference developmental-stage comparisons
# - Script 11: developmental-stage model
# - Script 13: taxonomic-order model
# - Script 14: matched-metal covariance/correlation
#
# Reader-facing figures generated
# Main Results:
#   Figure 4.1  Evidence landscape for the 96-h Cu-Cd-Zn comparison
#   Figure 4.2  Pooled observed LC50 values + model estimates
#   Figure 4.3  Pooled pairwise LC50 contrasts
#   Figure 4.4  Within-Reference metal comparisons
#   Figure 4.5  Matched-metal associations for Cu-Cd-Zn
#   Figure 4.6  Between- and within-Reference covariance components
#   Figure 4.7  Developmental-stage observed values + model estimates
#   Figure 4.8  Developmental-stage contrasts + interaction
#   Figure 4.9  Within-Reference developmental-stage comparisons
#   Figure 4.10 Taxonomic-order observed values + model estimates
#   Figure 4.11 Taxonomic-order contrasts + interaction
#   Figure 4.12 Exposure-duration comparisons
#   Figure 4.13 Salinity comparisons
#   Figure 4.14 Temperature comparisons
#   Figure 4.15 Pooled-model robustness/sensitivity
#   Figure 4.16 Developmental-stage robustness/sensitivity
#   Figure 4.17 Taxonomic-order robustness/sensitivity
#
# Appendix:
#   Figure C.1  pH-salinity confounding audit
#   Figure D.7  Additional matched-metal associations (Cu-Hg, Cd-Hg)
#
# Figure policy
# - Figures are presentation-only reconstructions of already-frozen
#   numerical results.
# - No significance-based figure selection is performed.
# - No source-restriction reproduction is shown as a reader-facing
#   result.
# - LORO = leave-one-Reference-out; LOSO = leave-one-Species-out.
# - Covariance/correlation figures describe association, not relative
#   toxicity, shared mechanism, mixture toxicity, or causation.
#
# Outputs
# Reader-facing PNG/PDF files:
#   results/figures/
# Plot-data/provenance/object files:
#   outputs/16_make_main_thesis_figures/
#
# Reproducibility rule
# Run from the repository root with copepod-metal-toxicity.Rproj open.
# Scripts 04-14 must have completed successfully.
# ============================================================


# ============================================================
# 0. PACKAGES + OUTPUT DIRECTORIES
# ============================================================

library(tidyverse)

options(width = 200)

# Required published summaries: fail before creating any figures if absent.
# Run the updated analysis scripts in your existing project first.
new_summary_paths <- c(
  metal_environment = "outputs/04_fit_pooled_metal_model/04_sensitivity/metal_environmental_sensitivity_comparison.csv",
  stage_environment = "outputs/11_fit_stage_model/04_sensitivity/stage_environmental_sensitivity_comparison.csv",
  order_environment = "outputs/13_fit_order_model/04_sensitivity/order_environmental_sensitivity_comparison.csv",
  metal_pairs = "outputs/05_compare_metals_within_references/01_tables/metal_effect_summary.csv",
  stage_pairs = "outputs/09_compare_stages_within_references/02_tables/stage_pair_effect_summary.csv",
  duration_pairs = "outputs/06_compare_exposure_duration/02_tables/duration_pair_effect_summary.csv",
  duration_repeated = "outputs/06_compare_exposure_duration/02_tables/duration_repeated_measures_summary.csv"
)
missing_summaries <- new_summary_paths[!file.exists(new_summary_paths)]
if (length(missing_summaries)) stop(paste(
  "Missing updated analysis outputs; no models are fitted by this figure-suite script:",
  paste(missing_summaries, collapse = "\n"), sep = "\n"))



output_root <- file.path(
  "outputs",
  "16_make_main_thesis_figures"
)

plot_data_dir <- file.path(
  output_root,
  "01_plot_data"
)

figure_dir <- file.path(
  "results",
  "figures"
)

object_dir <- file.path(
  output_root,
  "02_objects"
)

for (d in c(
  output_root,
  plot_data_dir,
  figure_dir,
  object_dir
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 1. HELPERS
# ============================================================

blue <- "#1677C8"
raw_grey <- "grey42"
band_grey <- "grey87"
grid_major <- "grey84"
grid_minor <- "grey92"
dark_grey <- "grey30"
mid_grey <- "grey55"
light_blue <- "#AFCFEA"
interaction_colour <- "#7A5AA6"
interaction_fill <- "#F4F0F8"
wos_colour <- "#7A5AA6"


# Export layout: wrap text against the full figure width, not the panel width.
# Only presentation is changed; data, model objects and intervals are untouched.
wrap_figure_text <- function(x, width_in, size_pt) {
  if (is.null(x) || !is.character(x) || length(x) != 1L || is.na(x)) return(x)
  nch <- max(24L, floor((width_in - 0.8) * 72 / (0.65 * size_pt)))
  paste(vapply(strsplit(x, "\n", fixed = TRUE)[[1]], function(line) {
    paste(strwrap(line, width = nch), collapse = "\n")
  }, character(1)), collapse = "\n")
}

prepare_figure_layout <- function(plot, width) {
  sizes <- c(title = 15, subtitle = 11, caption = 9)
  added_height <- 0
  labels <- list()
  for (key in names(sizes)) {
    old <- plot$labels[[key]]
    new <- wrap_figure_text(old, width, sizes[[key]])
    if (is.character(old) && length(old) == 1L && !is.na(old)) {
      old_n <- length(strsplit(old, "\n", fixed = TRUE)[[1]])
      new_n <- length(strsplit(new, "\n", fixed = TRUE)[[1]])
      added_height <- added_height + max(0, new_n - old_n) * sizes[[key]] * 1.3 / 72
      labels[[key]] <- new
    }
  }
  plot <- plot + do.call(ggplot2::labs, labels) + ggplot2::theme(
    plot.title.position = "plot", plot.caption.position = "plot",
    plot.title = ggplot2::element_text(size = 15, face = "bold", hjust = 0, lineheight = 1.1),
    plot.subtitle = ggplot2::element_text(size = 11, hjust = 0, lineheight = 1.1),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0, lineheight = 1.15,
                                         margin = ggplot2::margin(t = 10)),
    plot.margin = ggplot2::margin(t = 14, r = 24, b = 18, l = 18)
  )
  list(plot = plot, extra_height = added_height)
}

save_final_figure <- function(
    plot,
    filename,
    width,
    height
) {
  layout <- prepare_figure_layout(plot, width)
  plot <- layout$plot
  height <- height + layout$extra_height
  
  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename,
        ".png"
      )
    ),
    plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
  
  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename,
        ".pdf"
      )
    ),
    plot,
    width = width,
    height = height
  )
}


save_vertical_figure <- function(
    upper_plot,
    lower_plot,
    filename,
    width,
    height,
    upper_fraction = 0.68
) {
  upper <- prepare_figure_layout(upper_plot, width)
  lower <- prepare_figure_layout(lower_plot, width)
  upper_height <- height * upper_fraction + upper$extra_height
  # Five rows include multiline method labels; reserve room independently
  # of the upper panels and of the wrapped caption.
  lower_height <- max(4.6, height * (1 - upper_fraction)) + lower$extra_height
  height <- upper_height + lower_height
  upper_fraction <- upper_height / height
  upper_plot <- upper$plot
  lower_plot <- lower$plot
  
  png_file <- file.path(
    figure_dir,
    paste0(
      filename,
      ".png"
    )
  )
  
  pdf_file <- file.path(
    figure_dir,
    paste0(
      filename,
      ".pdf"
    )
  )
  
  png(
    filename = png_file,
    width = width,
    height = height,
    units = "in",
    res = 320,
    bg = "white"
  )
  
  grid::grid.newpage()
  
  print(
    upper_plot,
    vp = grid::viewport(
      x = 0.5,
      y = 1 - upper_fraction / 2,
      width = 1,
      height = upper_fraction
    )
  )
  
  print(
    lower_plot,
    vp = grid::viewport(
      x = 0.5,
      y = (1 - upper_fraction) / 2,
      width = 1,
      height = 1 - upper_fraction
    )
  )
  
  dev.off()
  
  pdf(
    file = pdf_file,
    width = width,
    height = height
  )
  
  grid::grid.newpage()
  
  print(
    upper_plot,
    vp = grid::viewport(
      x = 0.5,
      y = 1 - upper_fraction / 2,
      width = 1,
      height = upper_fraction
    )
  )
  
  print(
    lower_plot,
    vp = grid::viewport(
      x = 0.5,
      y = (1 - upper_fraction) / 2,
      width = 1,
      height = 1 - upper_fraction
    )
  )
  
  dev.off()
}


format_lc50 <- function(x) {
  
  ifelse(
    x < 10,
    sprintf(
      "%.2f",
      x
    ),
    sprintf(
      "%.1f",
      x
    )
  )
}


# Standard 1-2-5 sequence on a log axis.
# This makes positions such as 3 and 11 interpretable from the
# axis itself instead of leaving only decade marks (1, 10, 100).
make_125_breaks <- function(
    lower,
    upper
) {
  
  stopifnot(
    lower > 0,
    upper > lower
  )
  
  exponents <- seq(
    floor(
      log10(
        lower
      )
    ) - 1,
    ceiling(
      log10(
        upper
      )
    ) + 1
  )
  
  values <- sort(
    unique(
      as.vector(
        outer(
          c(
            1,
            2,
            5
          ),
          10 ^ exponents
        )
      )
    )
  )
  
  values[
    values >=
      lower *
      0.95 &
      values <=
      upper *
      1.05
  ]
}


log_tick_labels <- function(x) {
  
  vapply(
    x,
    function(z) {
      
      # ggplot2 may pass NA/NaN placeholders to a label function when
      # panel-specific scales are trained. Return a valid character label
      # instead of evaluating an indeterminate logical condition.
      if (
        is.na(z)
      ) {
        return(
          NA_character_
        )
      }
      
      if (
        z < 1
      ) {
        format(
          z,
          scientific = FALSE,
          trim = TRUE,
          digits = 3
        )
      } else {
        format(
          z,
          scientific = FALSE,
          trim = TRUE,
          big.mark = ",",
          digits = 6
        )
      }
    },
    character(
      1
    )
  )
}


raw_model_limits <- function(
    raw_y,
    lower_ci,
    upper_ci
) {
  
  all_y <- c(
    raw_y,
    lower_ci,
    upper_ci
  )
  
  all_y <- all_y[
    is.finite(
      all_y
    ) &
      all_y >
      0
  ]
  
  lower <- min(
    all_y
  ) *
    0.88
  
  upper <- max(
    all_y
  ) *
    1.08
  
  c(
    lower,
    upper
  )
}


theme_final <- theme_minimal(
  base_size = 13
) +
  theme(
    panel.grid.major.x =
      element_blank(),
    
    panel.grid.major.y =
      element_line(
        colour = grid_major,
        linewidth = 0.42
      ),
    
    panel.grid.minor.x =
      element_blank(),
    
    panel.grid.minor.y =
      element_line(
        colour = grid_minor,
        linewidth = 0.30
      ),
    
    axis.line.x =
      element_line(
        colour = "grey30",
        linewidth = 0.5
      ),
    
    axis.line.y =
      element_line(
        colour = "grey30",
        linewidth = 0.5
      ),
    
    axis.ticks =
      element_line(
        colour = "grey35",
        linewidth = 0.4
      ),
    
    axis.title =
      element_text(
        size = 12.5,
        colour = "black"
      ),
    
    axis.text =
      element_text(
        size = 11.5,
        colour = "black"
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    
    plot.subtitle =
      element_text(
        size = 11.5,
        margin = margin(
          b = 9
        )
      ),
    
    plot.caption =
      element_text(
        size = 9.3,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(
          t = 9
        )
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11.5
      ),
    
    strip.background =
      element_blank(),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        t = 12,
        r = 18,
        b = 15,
        l = 12
      )
  )


theme_ratio <- theme_classic(
  base_size = 13
) +
  theme(
    axis.line =
      element_line(
        colour = "grey25",
        linewidth = 0.55
      ),
    
    axis.ticks =
      element_line(
        colour = "grey25",
        linewidth = 0.5
      ),
    
    axis.title =
      element_text(
        size = 12.5,
        colour = "black"
      ),
    
    axis.text =
      element_text(
        size = 11.5,
        colour = "black"
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    
    plot.subtitle =
      element_text(
        size = 11.5,
        margin = margin(
          b = 10
        )
      ),
    
    plot.caption =
      element_text(
        size = 10,
        hjust = 0.5,
        lineheight = 1.08,
        margin = margin(
          t = 9
        )
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11.5
      ),
    
    strip.background =
      element_blank(),
    
    legend.position =
      "none",
    
    plot.margin =
      margin(
        t = 12,
        r = 18,
        b = 16,
        l = 12
      )
  )


make_raw_model_plot <- function(
    raw_data,
    model_data,
    x_breaks,
    x_labels,
    title,
    subtitle,
    caption,
    y_title,
    width = 8.6,
    height = 7.5
) {
  
  y_limits <- raw_model_limits(
    raw_y =
      raw_data$LC50_umol_L,
    lower_ci =
      model_data$lower,
    upper_ci =
      model_data$upper
  )
  
  y_breaks <- make_125_breaks(
    y_limits[
      1
    ],
    y_limits[
      2
    ]
  )
  
  model_for_plot <- model_data %>%
    dplyr::mutate(
      xmin =
        x -
        0.20,
      
      xmax =
        x +
        0.20,
      
      estimate_label =
        format_lc50(
          estimate
        ),
      
      label_x =
        x +
        0.22,
      
      label_y =
        estimate *
        1.04
    )
  
  p <- ggplot() +
    
    # CI first: raw points remain visible on top.
    geom_rect(
      data =
        model_for_plot,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = lower,
        ymax = upper
      ),
      fill = band_grey,
      colour = NA,
      alpha = 0.72
    ) +
    
    # One point per retained database Result.
    geom_point(
      data =
        raw_data,
      aes(
        x = x,
        y = LC50_umol_L
      ),
      position =
        position_jitter(
          width = 0.13,
          height = 0,
          seed = 20260828
        ),
      colour = raw_grey,
      alpha = 0.76,
      size = 2.0
    ) +
    
    # Model estimate.
    geom_segment(
      data =
        model_for_plot,
      aes(
        x = xmin,
        xend = xmax,
        y = estimate,
        yend = estimate
      ),
      colour = blue,
      linewidth = 1.35,
      lineend = "round"
    ) +
    
    # Exact estimate beside the line, not in the middle of the CI block.
    geom_text(
      data =
        model_for_plot,
      aes(
        x = label_x,
        y = label_y,
        label = estimate_label
      ),
      colour = blue,
      fontface = "bold",
      size = 3.7,
      hjust = 0,
      vjust = 0
    ) +
    
    scale_x_continuous(
      breaks =
        x_breaks,
      labels =
        x_labels,
      expand =
        expansion(
          add = 0.55
        )
    ) +
    
    scale_y_log10(
      breaks =
        y_breaks,
      labels =
        log_tick_labels(
          y_breaks
        ),
      minor_breaks =
        NULL
    ) +
    
    coord_cartesian(
      ylim =
        y_limits,
      clip =
        "off"
    ) +
    
    labs(
      title =
        title,
      subtitle =
        subtitle,
      x =
        NULL,
      y =
        y_title,
      caption =
        caption
    ) +
    
    theme_final
  
  attr(
    p,
    "final_width"
  ) <- width
  
  attr(
    p,
    "final_height"
  ) <- height
  
  p
}


make_leave_one_ranges <- function(
    points
) {
  
  points %>%
    dplyr::group_by(
      Estimand,
      Method
    ) %>%
    dplyr::summarise(
      min_ratio =
        min(
          Ratio,
          na.rm = TRUE
        ),
      max_ratio =
        max(
          Ratio,
          na.rm = TRUE
        ),
      .groups = "drop"
    )
}


# ============================================================
# 2. LOAD FROZEN ANALYSIS OBJECTS
# ============================================================

pooled_object_file <- file.path(
  "outputs",
  "04_fit_pooled_metal_model",
  "01_model",
  "primary_pooled_analysis_objects.rds"
)

stage_object_file <- file.path(
  "outputs",
  "11_fit_stage_model",
  "06_objects",
  "stage_model_analysis_objects.rds"
)

order_object_file <- file.path(
  "outputs",
  "13_fit_order_model",
  "06_objects",
  "order_model_analysis_objects.rds"
)

within_ref_dir <- file.path(
  "outputs",
  "05_compare_metals_within_references",
  "01_tables"
)

environmental_table_dir <- file.path(
  "outputs",
  "07_assess_salinity_temperature",
  "02_tables"
)

ph_table_dir <- file.path(
  "outputs",
  "08_assess_ph_evidence",
  "02_tables"
)

required_files <- c(
  pooled_object_file,
  stage_object_file,
  order_object_file,
  file.path("outputs", "04_fit_pooled_metal_model", "RUN_COMPLETE.txt"),
  file.path("outputs", "11_fit_stage_model", "RUN_COMPLETE.txt"),
  file.path("outputs", "13_fit_order_model", "RUN_COMPLETE.txt"),
  file.path("outputs", "14_analyze_metal_covariance", "RUN_COMPLETE.txt"),
  
  file.path(
    "outputs",
    "11_fit_stage_model",
    "01_model",
    "stage_estimated_marginal_means.csv"
  ),
  
  file.path(
    "outputs",
    "11_fit_stage_model",
    "01_model",
    "stage_metal_specific_contrasts.csv"
  ),
  
  file.path(
    "outputs",
    "11_fit_stage_model",
    "03_robustness",
    "stage_robustness_synthesis.csv"
  ),
  
  file.path(
    "outputs",
    "13_fit_order_model",
    "01_model",
    "order_estimated_marginal_means.csv"
  ),
  
  file.path(
    "outputs",
    "13_fit_order_model",
    "01_model",
    "order_metal_specific_contrasts.csv"
  ),
  
  file.path(
    "outputs",
    "13_fit_order_model",
    "03_robustness",
    "order_robustness_synthesis.csv"
  ),
  
  file.path(
    within_ref_dir,
    "Cu_Cd_verified_contexts.csv"
  ),
  
  file.path(
    within_ref_dir,
    "Cd_Zn_verified_contexts.csv"
  ),
  
  file.path(
    within_ref_dir,
    "Cu_Zn_verified_contexts.csv"
  ),
  
  file.path(
    environmental_table_dir,
    "salinity_main_direct_values.csv"
  ),
  
  file.path(
    environmental_table_dir,
    "salinity_extended_E2332_values.csv"
  ),
  
  file.path(
    environmental_table_dir,
    "temperature_final_direct_values.csv"
  ),
  
  file.path(
    ph_table_dir,
    "ph_within_reference_candidate_rows.csv"
  ),
  
  file.path(
    ph_table_dir,
    "ph_metadata_coverage_overall.csv"
  ),
  
  file.path(
    ph_table_dir,
    "wei2021_direct_acidification_values.csv"
  )
)

stopifnot(
  all(
    file.exists(
      required_files
    )
  )
)


pooled_obj <- readRDS(
  pooled_object_file
)

# Confirm that the pooled object contains the final Reference + Species
# random-intercept model used by the cleaned repository workflow.
pooled_groups <- names(lme4::getME(pooled_obj$primary_model, "flist"))
stopifnot(setequal(pooled_groups, c("Reference_ID", "Species")))

stage_obj <- readRDS(
  stage_object_file
)

order_obj <- readRDS(
  order_object_file
)


stopifnot(
  all(
    c(
      "primary_model_data",
      "metal_emmeans_original_scale",
      "metal_pairwise_ratios",
      "leave_one_reference_out",
      "leave_one_species_out",
      "robustness_point_estimates"
    ) %in%
      names(
        pooled_obj
      )
  ),
  
  all(
    c(
      "stage_model_data",
      "baseline_metrics",
      "loro_results",
      "loso_results"
    ) %in%
      names(
        stage_obj
      )
  ),
  
  all(
    c(
      "order_model_data",
      "baseline_metrics",
      "loro_results",
      "loso_results"
    ) %in%
      names(
        order_obj
      )
  )
)


# ============================================================
# 3. POOLED DATA PREPARATION
# ============================================================

pooled_raw <- pooled_obj$primary_model_data %>%
  as_tibble() %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels = c(
          "Cu",
          "Cd",
          "Zn"
        )
      ),
    
    x =
      as.numeric(
        Metal
      )
  )

pooled_model <- pooled_obj$metal_emmeans_original_scale %>%
  as_tibble() %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels = c(
          "Cu",
          "Cd",
          "Zn"
        )
      ),
    
    x =
      as.numeric(
        Metal
      )
  ) %>%
  dplyr::transmute(
    Metal,
    x,
    
    estimate =
      LC50_model_estimate_umol_L,
    
    lower =
      LC50_lower_umol_L,
    
    upper =
      LC50_upper_umol_L
  )

pooled_counts <- pooled_raw %>%
  dplyr::group_by(
    Metal
  ) %>%
  dplyr::summarise(
    Results =
      n(),
    References =
      dplyr::n_distinct(
        Reference_ID
      ),
    .groups =
      "drop"
  )

pooled_x_labels <- pooled_counts %>%
  dplyr::arrange(
    Metal
  ) %>%
  dplyr::transmute(
    label =
      paste0(
        as.character(
          Metal
        ),
        "\n",
        Results,
        " records | ",
        References,
        " publications"
      )
  ) %>%
  pull(
    label
  )

stopifnot(
  nrow(
    pooled_raw
  ) == 131,
  nrow(
    pooled_model
  ) == 3
)

write_csv(
  pooled_raw,
  file.path(
    plot_data_dir,
    "F01_pooled_raw_results.csv"
  )
)

write_csv(
  pooled_model,
  file.path(
    plot_data_dir,
    "Figure_4_02_model_estimates.csv"
  )
)


# ============================================================
# 3B. F00 96-h Cu-Cd-Zn EVIDENCE LANDSCAPE
# ============================================================
# Dedicated support/composition figure retained because the pooled model
# is only interpretable in relation to the biological evidence structure.
# Same molar LC50 scale is used in every panel.

landscape_data <- pooled_raw %>%
  dplyr::filter(
    !is.na(Order),
    !is.na(Species),
    Metal %in% c("Cu", "Cd", "Zn")
  ) %>%
  dplyr::mutate(
    Order = factor(
      Order,
      levels = c("Calanoida", "Cyclopoida", "Harpacticoida")
    ),
    Species = as.character(Species)
  )

landscape_species_levels <- landscape_data %>%
  dplyr::distinct(Order, Species) %>%
  dplyr::arrange(Order, Species) %>%
  pull(Species) %>%
  as.character()

landscape_data <- landscape_data %>%
  dplyr::mutate(
    Species = factor(
      as.character(Species),
      levels = rev(unique(landscape_species_levels))
    ),
    Source_label = recode(
      Source_Origin,
      "ECOTOX" = "ECOTOX",
      "WoS_supplemental" = "WoS supplemental",
      .default = as.character(Source_Origin)
    )
  )

landscape_limits <- c(
  min(landscape_data$LC50_umol_L, na.rm = TRUE) * 0.8,
  max(landscape_data$LC50_umol_L, na.rm = TRUE) * 1.2
)

landscape_breaks <- c(0.1, 1, 10, 100, 500)
landscape_breaks <- landscape_breaks[
  landscape_breaks >= landscape_limits[1] &
    landscape_breaks <= landscape_limits[2]
]

write_csv(
  landscape_data,
  file.path(
    plot_data_dir,
    "Figure_4_01_plot_data.csv"
  )
)

F00 <- ggplot(
  landscape_data,
  aes(
    x = LC50_umol_L,
    y = Species
  )
) +
  geom_point(
    aes(shape = Source_label, colour = Metal),
    alpha = 0.78,
    size = 2.7,
    position = position_jitter(height = 0.05, width = 0)
  ) +
  facet_grid(
    rows = vars(Order),
    cols = vars(Metal),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_log10(
    limits = landscape_limits,
    breaks = landscape_breaks,
    labels = log_tick_labels
  ) +
  scale_colour_manual(
    values = c(
      "Cu" = "#F26B5E",
      "Cd" = "#24B86A",
      "Zn" = "#5794E8"
    ),
    guide = "none"
  ) +
  labs(
    title =
      "Available 96-h LC50 data for Cu, Cd and Zn",
    subtitle =
      "Data available for the pooled metal comparison",
    x =
      expression(
        LC[50]~(mu*mol~metal~L^{-1})~"(log scale)"
      ),
    y =
      NULL,
    shape =
      "Data source",
    caption =
      paste0(
        "Each point is one retained 96-h LC50 record. ",
        "Panels use the same LC50 scale; empty metal-by-order cells contain no retained records.\n",
        "The figure shows which metals and taxonomic orders are represented in the data."
      )
  ) +
  theme_final +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    strip.text.y.right = element_text(
      face = "bold",
      angle = 0,
      margin = margin(l = 7, r = 7)
    ),
    axis.text.y = element_text(size = 8.3, face = "italic"),
    panel.grid.major.x = element_line(
      colour = "grey86",
      linewidth = 0.38
    ),
    panel.grid.minor.x = element_line(
      colour = "grey94",
      linewidth = 0.25
    ),
    panel.spacing = grid::unit(0.70, "lines"),
    plot.margin = margin(
      t = 12,
      r = 58,
      b = 18,
      l = 12
    )
  )

save_final_figure(
  F00,
  "Figure_4_01_evidence_landscape",
  width = 15.2,
  height = 9.4
)


# ============================================================
# 4. F01 POOLED RESULTS + MODEL ESTIMATES
# ============================================================

F01 <- make_raw_model_plot(
  raw_data = pooled_raw,
  model_data = pooled_model,
  x_breaks = c(1, 2, 3),
  x_labels = pooled_x_labels,
  title = "96-h LC50 for Cu, Cd, and Zn",
  subtitle = paste0(
    "96-h mixed-effects model with publication and species random intercepts\n",
    "Grey points = LC50 records | grey band = 95% confidence interval | blue line = model estimate"
  ),
  caption = paste0(
    "Records were compiled from ECOTOX and the supplementary literature search; one publication can contribute several records.\n",
    "Lower LC50 indicates greater acute toxicity."
  ),
  y_title = expression(
    LC[50]~(mu*"mol metal L"^{-1}*"; log scale")
  ),
  width = 8.8,
  height = 7.8
)

save_final_figure(
  F01,
  "Figure_4_02_pooled_results_and_model_estimates",
  width = 8.8,
  height = 7.8
)


# ============================================================
# 5. F02 POOLED PAIRWISE METAL CONTRASTS
# ============================================================
#
# Representation intentionally follows the existing thesis figure.

pooled_pairwise <- pooled_obj$metal_pairwise_ratios %>%
  as_tibble() %>%
  dplyr::mutate(
    Comparison =
      stringr::str_replace(
        contrast,
        " - ",
        " / "
      ),
    
    Comparison =
      factor(
        Comparison,
        levels =
          rev(
            c(
              "Cu / Cd",
              "Cu / Zn",
              "Cd / Zn"
            )
          )
      ),
    
    ratio_label =
      sprintf(
        "%.2f",
        LC50_ratio
      )
  )

write_csv(
  pooled_pairwise,
  file.path(
    plot_data_dir,
    "Figure_4_03_pairwise_contrasts.csv"
  )
)


F02 <- ggplot(
  pooled_pairwise,
  aes(
    x =
      LC50_ratio,
    y =
      Comparison
  )
) +
  
  geom_vline(
    xintercept =
      1,
    linetype =
      "dashed",
    linewidth =
      0.72,
    colour =
      "grey35"
  ) +
  
  geom_errorbar(
    aes(
      xmin =
        ratio_lower,
      xmax =
        ratio_upper
    ),
    orientation =
      "y",
    width =
      0.16,
    linewidth =
      0.95,
    colour =
      blue
  ) +
  
  geom_point(
    size =
      3.7,
    shape =
      21,
    stroke =
      0.9,
    fill =
      "white",
    colour =
      blue
  ) +
  
  geom_text(
    aes(
      label =
        ratio_label
    ),
    nudge_y =
      0.16,
    colour =
      blue,
    fontface =
      "bold",
    size =
      3.7
  ) +
  
  scale_x_log10(
    breaks =
      c(
        0.1,
        0.2,
        0.3,
        0.5,
        1,
        2
      ),
    labels =
      c(
        "0.1",
        "0.2",
        "0.3",
        "0.5",
        "1",
        "2"
      ),
    limits =
      c(
        0.085,
        2.05
      )
  ) +
  
  annotate(
    "text",
    x =
      0.94,
    y =
      0.62,
    label =
      "Equal LC50",
    hjust =
      1,
    size =
      3.2,
    colour =
      "grey30"
  ) +
  
  labs(
    title =
      "96-h LC50 Ratios Among Cu, Cd, and Zn",
    
    subtitle =
      "96-h mixed-effects model with publication and species random intercepts",
    
    x =
      "LC50 ratio (first-listed / second-listed metal; log scale)",
    
    y =
      NULL,
    
    caption =
      paste0(
        "< 1: first-listed metal has lower LC50   |   > 1: second-listed metal has lower LC50\n",
        "Bars: Tukey-adjusted 95% confidence intervals. Lower LC50 indicates greater acute toxicity."
      )
  ) +
  
  coord_cartesian(
    clip =
      "off"
  ) +
  
  theme_ratio +
  
  theme(
    plot.margin =
      margin(
        t = 12,
        r = 18,
        b = 30,
        l = 12
      )
  )

save_final_figure(
  F02,
  "Figure_4_03_pooled_pairwise_contrasts",
  width =
    8.4,
  height =
    5.8
)


# ============================================================
# 6. F03 POOLED POST-MODEL ROBUSTNESS
# ============================================================

pooled_variant_map <- c(
  "Combined Reference + Species random intercepts" =
    "Combined pooled
Cu-Cd-Zn LMM",
  
  "Reference fixed effects + Species random intercept" =
    "Publication FE + Species RI",
  
  "Reference-only random intercept" =
    "Publication-only RI",
  
  "Reference x Species x Metal cell collapse" =
    "Cell collapse"
)

pooled_method_order <- c(
  "Combined pooled
Cu-Cd-Zn LMM",
  "LORO",
  "LOSO",
  "Publication FE + Species RI",
  "Publication-only RI",
  "Cell collapse",
  "Complete-case",
  "Temperature +
salinity adjusted"
)

pooled_estimand_order <- c(
  "Cu / Cd",
  "Cu / Zn",
  "Cd / Zn"
)

pooled_variant_points <- pooled_obj$robustness_point_estimates %>%
  as_tibble() %>%
  dplyr::filter(
    Variant %in%
      names(
        pooled_variant_map
      )
  ) %>%
  dplyr::mutate(
    Estimand =
      stringr::str_replace(
        contrast,
        " - ",
        " / "
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          pooled_estimand_order
      ),
    
    Method =
      recode(
        Variant,
        !!!pooled_variant_map
      ),
    
    Method =
      factor(
        Method,
        levels =
          rev(
            pooled_method_order
          )
      ),
    
    Baseline =
      Variant ==
      "Combined Reference + Species random intercepts"
  ) %>%
  dplyr::transmute(
    Estimand,
    Method,
    Baseline,
    Ratio =
      LC50_ratio,
    Lower =
      ratio_lower,
    Upper =
      ratio_upper
  )


# Environmental-data sensitivity is shown in the same reader-facing
# robustness figure as the other pooled-model sensitivity checks.
# The full-data model is already represented by the baseline row.
pooled_environment_points <- readr::read_csv(
  new_summary_paths[["metal_environment"]],
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  dplyr::filter(
    Model %in% c(
      "TEMP_SAL_COMPLETE_CASE_SAME_FORMULA",
      "TEMP_SAL_COMPLETE_CASE_ADJUSTED"
    )
  ) %>%
  dplyr::mutate(
    Estimand =
      stringr::str_replace(
        contrast,
        " - ",
        " / "
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          pooled_estimand_order
      ),
    
    Method =
      dplyr::recode(
        Model,
        "TEMP_SAL_COMPLETE_CASE_SAME_FORMULA" =
          "Complete-case",
        "TEMP_SAL_COMPLETE_CASE_ADJUSTED" =
          "Temperature +\nsalinity adjusted"
      ),
    
    Method =
      factor(
        Method,
        levels =
          rev(
            pooled_method_order
          )
      ),
    
    Baseline =
      FALSE
  ) %>%
  dplyr::transmute(
    Estimand,
    Method,
    Baseline,
    Ratio =
      LC50_ratio,
    Lower =
      ratio_lower,
    Upper =
      ratio_upper
  )

stopifnot(
  nrow(pooled_environment_points) == 6,
  all(is.finite(pooled_environment_points$Ratio)),
  all(pooled_environment_points$Ratio > 0)
)

pooled_variant_points <- dplyr::bind_rows(
  pooled_variant_points,
  pooled_environment_points
)

# Final-reader robustness outputs must not contain a source-restriction row.
stopifnot(
  !any(
    stringr::str_detect(
      as.character(pooled_variant_points$Method),
      stringr::fixed("ECOTOX")
    ),
    na.rm = TRUE
  )
)


pooled_loro_points <- pooled_obj$leave_one_reference_out %>%
  as_tibble() %>%
  dplyr::mutate(
    Estimand =
      stringr::str_replace(
        contrast,
        " - ",
        " / "
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          pooled_estimand_order
      ),
    
    Method =
      factor(
        "LORO",
        levels =
          rev(
            pooled_method_order
          )
      )
  ) %>%
  dplyr::transmute(
    Estimand,
    Method,
    Ratio =
      LC50_ratio
  )


pooled_loso_points <- pooled_obj$leave_one_species_out %>%
  as_tibble() %>%
  dplyr::mutate(
    Estimand =
      stringr::str_replace(
        contrast,
        " - ",
        " / "
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          pooled_estimand_order
      ),
    
    Method =
      factor(
        "LOSO",
        levels =
          rev(
            pooled_method_order
          )
      )
  ) %>%
  dplyr::transmute(
    Estimand,
    Method,
    Ratio =
      LC50_ratio
  )


pooled_loro_ranges <- bind_rows(
  pooled_loro_points,
  pooled_loso_points
) %>%
  make_leave_one_ranges()


write_csv(
  pooled_variant_points,
  file.path(
    plot_data_dir,
    "Figure_4_15_robustness_variants.csv"
  )
)

write_csv(
  pooled_loro_points,
  file.path(
    plot_data_dir,
    "Figure_4_15_LORO_points.csv"
  )
)

write_csv(
  pooled_loso_points,
  file.path(
    plot_data_dir,
    "Figure_4_15_LOSO_points.csv"
  )
)


F03 <- ggplot() +
  
  geom_vline(
    xintercept =
      1,
    linetype =
      "dashed",
    colour =
      "grey35",
    linewidth =
      0.70
  ) +
  
  geom_segment(
    data =
      pooled_loro_ranges,
    aes(
      x =
        min_ratio,
      xend =
        max_ratio,
      y =
        Method,
      yend =
        Method
    ),
    colour =
      "grey78",
    linewidth =
      3.0,
    lineend =
      "round"
  ) +
  
  geom_point(
    data =
      pooled_loro_points,
    aes(
      x =
        Ratio,
      y =
        Method
    ),
    position =
      position_jitter(
        width = 0,
        height = 0.065,
        seed = 20260828
      ),
    colour =
      "grey52",
    alpha =
      0.70,
    size =
      1.7
  ) +
  
  geom_point(
    data =
      pooled_loso_points,
    aes(
      x =
        Ratio,
      y =
        Method
    ),
    position =
      position_jitter(
        width = 0,
        height = 0.065,
        seed = 20260829
      ),
    colour =
      dark_grey,
    alpha =
      0.82,
    size =
      1.8,
    shape =
      17
  ) +
  
  geom_errorbar(
    data =
      pooled_variant_points %>%
      dplyr::filter(
        !Baseline
      ),
    aes(
      x =
        Ratio,
      xmin =
        Lower,
      xmax =
        Upper,
      y =
        Method
    ),
    orientation =
      "y",
    width =
      0.15,
    colour =
      dark_grey,
    linewidth =
      0.72
  ) +
  
  geom_point(
    data =
      pooled_variant_points %>%
      dplyr::filter(
        !Baseline
      ),
    aes(
      x =
        Ratio,
      y =
        Method
    ),
    shape =
      21,
    fill =
      "white",
    colour =
      dark_grey,
    stroke =
      0.7,
    size =
      2.5
  ) +
  
  geom_errorbar(
    data =
      pooled_variant_points %>%
      dplyr::filter(
        Baseline
      ),
    aes(
      x =
        Ratio,
      xmin =
        Lower,
      xmax =
        Upper,
      y =
        Method
    ),
    orientation =
      "y",
    width =
      0.16,
    colour =
      blue,
    linewidth =
      1.05
  ) +
  
  geom_point(
    data =
      pooled_variant_points %>%
      dplyr::filter(
        Baseline
      ),
    aes(
      x =
        Ratio,
      y =
        Method
    ),
    shape =
      21,
    fill =
      "white",
    colour =
      blue,
    stroke =
      1.0,
    size =
      3.0
  ) +
  
  scale_x_log10(
    breaks =
      c(
        0.1,
        0.2,
        0.5,
        1,
        2,
        5
      ),
    labels =
      c(
        "0.1",
        "0.2",
        "0.5",
        "1",
        "2",
        "5"
      )
  ) +
  
  scale_y_discrete(
    limits =
      rev(
        pooled_method_order
      ),
    drop =
      FALSE
  ) +
  
  facet_grid(
    rows =
      vars(
        Estimand
      )
  ) +
  
  labs(
    title =
      "Sensitivity of the pooled metal comparisons",
    
    subtitle =
      paste0(
        "96-h Cu-Cd-Zn LMM\n",
        "Primary LC50 ratios, omission checks, alternative model specifications and environmental sensitivity"
      ),
    
    x =
      "LC50 ratio (log scale)",
    
    y =
      NULL,
    
    caption =
      paste0(
        "LORO = one publication omitted at a time; LOSO = one species group omitted at a time. Grey bands show the resulting point-estimate ranges.\n",
        "RI = random intercept; FE = fixed effect. Cell collapse = one mean ln(LC50) per publication x species x metal group.\n",
        "CC = records with both temperature and salinity reported. CC + T + S = the same subset with temperature and salinity added to the model.\n",
        "Open circles and bars: model ratios and Tukey-adjusted 95% confidence intervals. Grey bands show point-estimate ranges, not confidence intervals."
      )
  ) +
  
  theme_ratio +
  
  theme(
    panel.grid.major.x =
      element_line(
        colour =
          grid_minor,
        linewidth =
          0.35
      ),
    
    panel.grid.major.y =
      element_blank(),
    
    strip.text.y =
      element_text(
        angle =
          0,
        hjust =
          0
      ),
    
    strip.placement =
      "outside",
    
    panel.spacing.y =
      unit(
        0.9,
        "lines"
      ),
    plot.caption = element_text(
      size = 8.6,
      hjust = 0,
      lineheight = 1.10,
      margin = margin(t = 10)
    ),
    plot.margin = margin(
      t = 12,
      r = 22,
      b = 24,
      l = 16
    )
  )

save_final_figure(
  F03,
  "Figure_4_15_pooled_model_robustness",
  width =
    10.2,
  height =
    10.9
)


# ============================================================
# 7. F04 WITHIN-REFERENCE METAL COMPARISONS
# ============================================================
#
# LOCKED DESIGN:
# - one figure, three metal pairs
# - pair order matches pooled pairwise contrasts:
#   Cu-Cd -> Cu-Zn -> Cd-Zn
# - landscape layout
# - no right-side facet boxes
# - pair support appears in compact panel headers
# - same Reference may appear in more than one pair
# - subtle colour separation aids pair recognition

within_ref_data <- bind_rows(
  read_csv(
    file.path(
      within_ref_dir,
      "Cu_Cd_verified_contexts.csv"
    ),
    show_col_types = FALSE
  ),
  read_csv(
    file.path(
      within_ref_dir,
      "Cu_Zn_verified_contexts.csv"
    ),
    show_col_types = FALSE
  ),
  read_csv(
    file.path(
      within_ref_dir,
      "Cd_Zn_verified_contexts.csv"
    ),
    show_col_types = FALSE
  )
) %>%
  dplyr::mutate(
    Comparison = factor(
      Comparison,
      levels = c(
        "Cu-Cd",
        "Cu-Zn",
        "Cd-Zn"
      )
    ),
    Ref_label = as.character(Reference_ID)
  )

within_ref_support <- within_ref_data %>%
  dplyr::group_by(
    Comparison
  ) %>%
  dplyr::summarise(
    Contexts = n(),
    References = dplyr::n_distinct(
      Reference_ID
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Header = paste0(
      as.character(
        Comparison
      ),
      " | ",
      Contexts,
      " matched comparisons | ",
      References,
      " publications"
    )
  )

pair_colours <- c(
  # These colours identify comparison blocks, not individual metals.
  # They are deliberately distinct from the landscape metal palette.
  "Cu-Cd" = "#6C5B7B",
  "Cu-Zn" = "#A06A3B",
  "Cd-Zn" = "#5D6A6E"
)

comparison_order <- c(
  # Numerical y positions increase upward. Build bottom-to-top so
  # the displayed order is Cu-Cd -> Cu-Zn -> Cd-Zn.
  "Cd-Zn",
  "Cu-Zn",
  "Cu-Cd"
)

row_key_list <- list()
current_top <- 0

for (
  comp in comparison_order
) {
  
  refs <- within_ref_data %>%
    dplyr::filter(
      as.character(
        Comparison
      ) == comp
    ) %>%
    dplyr::distinct(
      Reference_ID,
      Ref_label
    ) %>%
    dplyr::arrange(
      Reference_ID
    )
  
  n_ref <- nrow(
    refs
  )
  
  refs <- refs %>%
    dplyr::mutate(
      y = current_top +
        rev(
          seq_len(
            n_ref
          )
        )
    )
  
  row_key_list[[comp]] <- refs %>%
    dplyr::mutate(
      Comparison = comp
    )
  
  current_top <- current_top +
    n_ref +
    2.2
}

within_ref_row_key <- bind_rows(
  row_key_list
)

within_ref_plot <- within_ref_data %>%
  dplyr::mutate(
    Comparison_chr = as.character(
      Comparison
    )
  ) %>%
  dplyr::left_join(
    within_ref_row_key %>%
      dplyr::transmute(
        Comparison_chr = Comparison,
        Reference_ID,
        y
      ),
    by = c(
      "Comparison_chr",
      "Reference_ID"
    )
  )

within_ref_segments <- within_ref_plot %>%
  dplyr::group_by(
    Comparison,
    Reference_ID,
    y
  ) %>%
  dplyr::summarise(
    xmin = min(
      ln_ratio
    ),
    xmax = max(
      ln_ratio
    ),
    n_contexts = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_contexts > 1
  )

within_ref_headers <- within_ref_support %>%
  dplyr::mutate(
    Comparison_chr = as.character(
      Comparison
    )
  ) %>%
  dplyr::left_join(
    within_ref_row_key %>%
      dplyr::group_by(
        Comparison
      ) %>%
      dplyr::summarise(
        header_y = max(
          y
        ) + 0.75,
        .groups = "drop"
      ) %>%
      dplyr::rename(
        Comparison_chr = Comparison
      ),
    by = "Comparison_chr"
  )

within_ref_separators <- within_ref_row_key %>%
  dplyr::group_by(
    Comparison
  ) %>%
  dplyr::summarise(
    # Two separators: one between each adjacent comparison block.
    separator_y = max(
      y
    ) + 1.55,
    .groups = "drop"
  ) %>%
  dplyr::filter(
    as.character(
      Comparison
    ) != tail(
      comparison_order,
      1
    )
  )

x_min <- min(
  within_ref_plot$ln_ratio,
  na.rm = TRUE
)

x_max <- max(
  within_ref_plot$ln_ratio,
  na.rm = TRUE
)

x_pad <- max(
  0.25,
  0.07 * (
    x_max - x_min
  )
)

write_csv(
  within_ref_support,
  file.path(
    plot_data_dir,
    "Figure_4_04_support.csv"
  )
)

write_csv(
  within_ref_plot,
  file.path(
    plot_data_dir,
    "Figure_4_04_points.csv"
  )
)

write_csv(
  within_ref_segments,
  file.path(
    plot_data_dir,
    "Figure_4_04_segments.csv"
  )
)

F04 <- ggplot() +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +
  
  geom_hline(
    data = within_ref_separators,
    aes(
      yintercept = separator_y
    ),
    colour = "grey80",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  
  geom_segment(
    data = within_ref_segments,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y,
      colour = Comparison
    ),
    linewidth = 1.10,
    alpha = 0.40
  ) +
  
  geom_point(
    data = within_ref_plot,
    aes(
      x = ln_ratio,
      y = y,
      colour = Comparison
    ),
    size = 2.65,
    alpha = 0.95
  ) +
  
  geom_text(
    data = within_ref_headers,
    aes(
      x = x_min - x_pad,
      y = header_y,
      label = Header,
      colour = Comparison
    ),
    hjust = 0,
    vjust = 0.5,
    fontface = "bold",
    size = 3.65,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = pair_colours
  ) +
  
  scale_y_continuous(
    breaks = within_ref_row_key$y,
    labels = within_ref_row_key$Ref_label,
    expand = expansion(
      mult = c(
        0.03,
        0.07
      )
    )
  ) +
  
  coord_cartesian(
    xlim = c(
      x_min - x_pad,
      x_max + x_pad
    ),
    clip = "off"
  ) +
  
  labs(
    title =
      "LC50 comparisons between metals",
    
    subtitle =
      paste0(
        "Each point compares two metals under matched test conditions within the same publication.\n",
        "The dashed line at 0 indicates equal LC50 for the two metals."
      ),
    
    x =
      "ln(LC50 ratio: first-listed / second-listed metal)",
    
    y =
      "Publication",
    
    caption =
      paste0(
        "< 0: first-listed metal has lower LC50   |   > 0: second-listed metal has lower LC50.\n",
        "Horizontal segments connect comparisons from the same publication; a publication can contribute to more than one metal pair."
      )
  ) +
  
  theme_ratio +
  
  theme(
    panel.grid.major.x =
      element_line(
        colour = grid_minor,
        linewidth = 0.35
      ),
    
    panel.grid.major.y =
      element_blank(),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9.6
      )
  )

save_final_figure(
  F04,
  "Figure_4_04_within_reference_metal_comparisons",
  width = 11.8,
  height = 8.4
)


# ============================================================
# 8. STAGE DATA PREPARATION
# ============================================================

stage_raw <- stage_obj$stage_model_data %>%
  as_tibble() %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Stage =
      factor(
        Stage,
        levels =
          c(
            "Adult",
            "Nauplii"
          )
      )
  )


stage_x_map <- tribble(
  ~Metal, ~Stage, ~x,
  "Cu", "Adult", 1,
  "Cu", "Nauplii", 3,
  "Cd", "Adult", 6,
  "Cd", "Nauplii", 8
) %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Stage =
      factor(
        Stage,
        levels =
          c(
            "Adult",
            "Nauplii"
          )
      )
  )


stage_raw <- stage_raw %>%
  dplyr::left_join(
    stage_x_map,
    by =
      c(
        "Metal",
        "Stage"
      )
  )

stopifnot(
  !any(
    is.na(
      stage_raw$x
    )
  )
)


stage_emm <- read_csv(
  file.path(
    "outputs",
    "11_fit_stage_model",
    "01_model",
    "stage_estimated_marginal_means.csv"
  ),
  show_col_types =
    FALSE
) %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Stage =
      factor(
        Stage,
        levels =
          c(
            "Adult",
            "Nauplii"
          )
      ),
    
    estimate =
      exp(
        emmean
      ),
    
    lower =
      exp(
        lower.CL
      ),
    
    upper =
      exp(
        upper.CL
      )
  ) %>%
  dplyr::left_join(
    stage_x_map,
    by =
      c(
        "Metal",
        "Stage"
      )
  )


stage_counts <- stage_raw %>%
  dplyr::group_by(Metal, Stage) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    .groups = "drop"
  ) %>%
  dplyr::left_join(stage_x_map, by = c("Metal", "Stage")) %>%
  dplyr::arrange(x)

stage_x_labels <- paste0(
  as.character(stage_counts$Metal),
  " | ",
  as.character(stage_counts$Stage),
  "\n",
  stage_counts$Results,
  " records | ",
  stage_counts$References,
  " publications"
)


write_csv(
  stage_raw,
  file.path(
    plot_data_dir,
    "Figure_4_07_raw_results.csv"
  )
)

write_csv(
  stage_emm,
  file.path(
    plot_data_dir,
    "Figure_4_07_model_estimates.csv"
  )
)


# ============================================================
# 9. F05 STAGE RESULTS + MODEL ESTIMATES
# ============================================================

F05 <- make_raw_model_plot(
  raw_data = stage_raw,
  model_data = stage_emm %>%
    dplyr::transmute(Metal, Stage, x, estimate, lower, upper),
  x_breaks = stage_counts$x,
  x_labels = stage_x_labels,
  title = "96-h LC50 by Developmental Stage",
  subtitle = paste0(
    "96-h Cu-Cd developmental-stage mixed-effects model\n",
    "Grey points = LC50 records | grey band = 95% confidence interval | blue line = model estimate"
  ),
  caption = paste0(
    "Adult and Nauplii model estimates are shown separately for Cu and Cd.\n",
    "Lower LC50 indicates greater acute toxicity."
  ),
  y_title = expression(
    LC[50]~(mu*"mol metal L"^{-1}*"; log scale")
  ),
  width = 9.2,
  height = 7.6
)

save_final_figure(
  F05,
  "Figure_4_07_stage_results_and_model_estimates",
  width = 10.4,
  height = 7.8
)


# ============================================================
# 10. F06A DEVELOPMENTAL-STAGE MODEL CONTRASTS
# ============================================================
#
# LOCKED DESIGN:
# - upper panel: two metal-specific Nauplii / Adult contrasts
# - lower mini-panel: formal Metal x Stage interaction
# - interaction is a distinct purple diamond estimand, not a
#   third peer group
# - RoR is defined explicitly in the caption

stage_contrasts <- read_csv(
  file.path(
    "outputs",
    "11_fit_stage_model",
    "01_model",
    "stage_metal_specific_contrasts.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Metal = factor(
      Metal,
      levels = rev(
        c(
          "Cu",
          "Cd"
        )
      )
    ),
    ratio_label = sprintf(
      "%.2f",
      LC50_ratio_Nauplii_Adult
    )
  )

stage_base_for_figure <- stage_obj$baseline_metrics %>%
  as_tibble()

stopifnot(
  nrow(
    stage_base_for_figure
  ) == 1
)

stage_interaction <- tibble(
  Estimate =
    stage_base_for_figure$ratio_of_ratios,
  Lower =
    stage_base_for_figure$interaction_lower95,
  Upper =
    stage_base_for_figure$interaction_upper95,
  Label =
    "Interaction: Cd/Cu RoR"
)

write_csv(
  stage_contrasts,
  file.path(
    plot_data_dir,
    "Figure_4_08_stage_contrasts.csv"
  )
)

write_csv(
  stage_interaction,
  file.path(
    plot_data_dir,
    "Figure_4_08_stage_interaction.csv"
  )
)

F06_main <- ggplot(
  stage_contrasts,
  aes(
    x = LC50_ratio_Nauplii_Adult,
    y = Metal
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.72
  ) +
  
  geom_errorbar(
    aes(
      xmin = ratio_lower_95,
      xmax = ratio_upper_95
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.90,
    colour = blue
  ) +
  
  geom_point(
    size = 3.5,
    shape = 21,
    fill = "white",
    colour = blue,
    stroke = 0.9
  ) +
  
  geom_text(
    aes(
      label = ratio_label
    ),
    nudge_y = 0.13,
    colour = blue,
    fontface = "bold",
    size = 3.6
  ) +
  
  scale_x_log10(
    breaks = c(
      0.1,
      0.2,
      0.5,
      1,
      2,
      5
    ),
    labels = c(
      "0.1",
      "0.2",
      "0.5",
      "1",
      "2",
      "5"
    ),
    limits = c(
      0.08,
      5.0
    )
  ) +
  
  annotate(
    "text",
    x = 1,
    y = 2.38,
    label = "Equal LC50",
    hjust = 0.5,
    vjust = 0,
    size = 3.0,
    colour = "grey30"
  ) +
  
  labs(
    title =
      "Model-estimated LC50 ratios between developmental stages",
    
    subtitle =
      "Model-estimated Nauplii / Adult LC50 ratios with 95% confidence intervals",
    
    x =
      "Nauplii / Adult LC50 ratio (log scale)",
    
    y =
      NULL,
    
    caption =
      "Ratio < 1: Nauplii have lower LC50   |   Ratio > 1: Adults have lower LC50"
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_ratio +
  
  theme(
    plot.margin = margin(
      t = 12,
      r = 18,
      b = 8,
      l = 12
    )
  )

F06_interaction <- ggplot(
  stage_interaction,
  aes(
    x = Estimate,
    y = 1
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.62
  ) +
  
  geom_errorbar(
    aes(
      xmin = Lower,
      xmax = Upper
    ),
    orientation = "y",
    width = 0.11,
    linewidth = 0.82,
    colour = interaction_colour
  ) +
  
  geom_point(
    shape = 23,
    size = 3.5,
    fill = interaction_colour,
    colour = interaction_colour
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.2f", Estimate)
    ),
    nudge_y = 0.13,
    hjust = 0.5,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.15
  ) +
  
  scale_x_log10(
    breaks = c(
      0.2,
      0.5,
      1,
      2,
      5,
      10,
      20
    ),
    labels = c(
      "0.2",
      "0.5",
      "1",
      "2",
      "5",
      "10",
      "20"
    ),
    limits = c(
      0.18,
      22
    )
  ) +
  
  scale_y_continuous(
    breaks = NULL
  ) +
  
  labs(
    x =
      "Cd / Cu ratio of developmental-stage contrasts (log scale)",
    
    y =
      NULL,
    
    caption =
      paste0(
        "RoR = ratio of ratios; it compares the Cd Nauplii / Adult contrast with the corresponding Cu contrast. ",
        "RoR = 1 means the developmental-stage contrast is the same in Cu and Cd."
      )
  ) +
  
  theme_ratio +
  
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_text(
      hjust = 0,
      size = 8.8
    ),
    plot.margin = margin(
      t = 2,
      r = 18,
      b = 12,
      l = 12
    )
  )

# Present the two metal-specific contrasts and the formal interaction in one
# compact forest plot, matching the visual logic of pooled Figure F02.

stage_model_comparisons <- bind_rows(
  stage_contrasts %>%
    dplyr::transmute(
      Comparison = paste0(as.character(Metal), ": Nauplii / Adult"),
      Estimate = LC50_ratio_Nauplii_Adult,
      Lower = ratio_lower_95,
      Upper = ratio_upper_95,
      Estimand_type = "Metal-specific contrast"
    ),
  stage_interaction %>%
    dplyr::transmute(
      Comparison = "Interaction: Cd / Cu RoR",
      Estimate,
      Lower,
      Upper,
      Estimand_type = "Interaction"
    )
) %>%
  dplyr::mutate(
    Comparison = factor(
      Comparison,
      levels = rev(c(
        "Cu: Nauplii / Adult",
        "Cd: Nauplii / Adult",
        "Interaction: Cd / Cu RoR"
      ))
    ),
    ratio_label = sprintf("%.2f", Estimate),
    y_position = case_when(
      as.character(Comparison) == "Cu: Nauplii / Adult" ~ 3.0,
      as.character(Comparison) == "Cd: Nauplii / Adult" ~ 2.0,
      TRUE ~ 0.65
    )
  )

F06A <- ggplot(
  stage_model_comparisons,
  aes(x = Estimate, y = y_position)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.72,
    colour = "grey35"
  ) +
  geom_errorbar(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    aes(xmin = Lower, xmax = Upper),
    orientation = "y",
    width = 0.16,
    linewidth = 0.95,
    colour = blue
  ) +
  geom_point(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    size = 3.7,
    shape = 21,
    stroke = 0.9,
    fill = "white",
    colour = blue
  ) +
  geom_text(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    aes(label = ratio_label),
    nudge_y = 0.16,
    colour = blue,
    fontface = "bold",
    size = 3.7
  ) +
  geom_errorbar(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    aes(xmin = Lower, xmax = Upper),
    orientation = "y",
    width = 0.16,
    linewidth = 0.95,
    colour = interaction_colour
  ) +
  geom_point(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    size = 3.9,
    shape = 23,
    stroke = 0.9,
    fill = "white",
    colour = interaction_colour
  ) +
  geom_text(
    data = stage_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    aes(label = ratio_label),
    nudge_y = 0.16,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.7
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20"),
    limits = c(0.08, 22)
  ) +
  scale_y_continuous(
    breaks = c(3.0, 2.0, 0.65),
    labels = c(
      "Cu: Nauplii / Adult",
      "Cd: Nauplii / Adult",
      "Interaction: Cd / Cu RoR"
    )
  ) +
  labs(
    title = "Developmental-stage comparisons and their variation between metals",
    subtitle = "Ratios and 95% confidence intervals from the 96-h Cu-Cd mixed-effects model",
    x = "Ratio estimate (log scale)",
    y = NULL,
    caption = paste0(
      "Blue: Nauplii / Adult LC50 ratio within each metal. Purple: Cd / Cu ratio of those contrasts (RoR).\n",
      "The dashed line at 1 indicates equal LC50 for the blue contrasts and no interaction for the purple RoR."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    plot.margin = margin(t = 12, r = 18, b = 30, l = 12)
  )

save_final_figure(
  F06A,
  "Figure_4_08_stage_contrasts_and_interaction",
  width = 8.4,
  height = 5.8
)


# ============================================================
# 11. F06 STAGE POST-MODEL ROBUSTNESS
# ============================================================

stage_method_order <- c(
  "96-h Cu-Cd\nStage LMM",
  "LORO",
  "LOSO",
  "Complete-case",
  "Temperature +\nsalinity adjusted"
)

stage_method_labels <- c(
  "96-h Cu-Cd\nStage LMM" = "96-h Cu-Cd\nStage LMM",
  "LORO" = "LORO",
  "LOSO" = "LOSO",
  "Complete-case" = "CC\n(T + S reported)",
  "Temperature +\nsalinity adjusted" = "CC + T + S\n(T + S in model)"
)

stage_estimand_order <- c(
  "Cu: Nauplii / Adult",
  "Cd: Nauplii / Adult",
  "Interaction: Cd / Cu ratio of Stage contrasts"
)


stage_base <- stage_obj$baseline_metrics %>%
  as_tibble()

stopifnot(
  nrow(
    stage_base
  ) ==
    1
)


stage_baseline_points <- tibble(
  Estimand =
    factor(
      stage_estimand_order,
      levels =
        stage_estimand_order
    ),
  
  Method =
    factor(
      rep(
        "96-h Cu-Cd\nStage LMM",
        3
      ),
      levels =
        rev(
          stage_method_order
        )
    ),
  
  Baseline =
    TRUE,
  
  Ratio =
    c(
      stage_base$Cu_ratio,
      stage_base$Cd_ratio,
      stage_base$ratio_of_ratios
    ),
  
  Lower =
    c(
      stage_base$Cu_lower95,
      stage_base$Cd_lower95,
      stage_base$interaction_lower95
    ),
  
  Upper =
    c(
      stage_base$Cu_upper95,
      stage_base$Cd_upper95,
      stage_base$interaction_upper95
    )
)


stage_loro <- stage_obj$loro_results %>%
  as_tibble() %>%
  dplyr::filter(
    !is.na(
      Cu_ratio
    )
  )

stage_loso <- stage_obj$loso_results %>%
  as_tibble() %>%
  dplyr::filter(
    !is.na(
      Cu_ratio
    )
  )


stage_loro_points <- bind_rows(
  
  stage_loro %>%
    dplyr::transmute(
      Estimand =
        "Cu: Nauplii / Adult",
      Ratio =
        Cu_ratio
    ),
  
  stage_loro %>%
    dplyr::transmute(
      Estimand =
        "Cd: Nauplii / Adult",
      Ratio =
        Cd_ratio
    ),
  
  stage_loro %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Stage contrasts",
      Ratio =
        ratio_of_ratios
    )
) %>%
  dplyr::mutate(
    Estimand =
      factor(
        Estimand,
        levels =
          stage_estimand_order
      ),
    
    Method =
      factor(
        "LORO",
        levels =
          rev(
            stage_method_order
          )
      )
  )


stage_loso_points <- bind_rows(
  
  stage_loso %>%
    dplyr::transmute(
      Estimand =
        "Cu: Nauplii / Adult",
      Ratio =
        Cu_ratio
    ),
  
  stage_loso %>%
    dplyr::transmute(
      Estimand =
        "Cd: Nauplii / Adult",
      Ratio =
        Cd_ratio
    ),
  
  stage_loso %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Stage contrasts",
      Ratio =
        ratio_of_ratios
    )
) %>%
  dplyr::mutate(
    Estimand =
      factor(
        Estimand,
        levels =
          stage_estimand_order
      ),
    
    Method =
      factor(
        "LOSO",
        levels =
          rev(
            stage_method_order
          )
      )
  )


stage_syn <- read_csv(
  file.path(
    "outputs",
    "11_fit_stage_model",
    "03_robustness",
    "stage_robustness_synthesis.csv"
  ),
  show_col_types =
    FALSE
)


extract_range <- function(
    x,
    which =
      c(
        "lower",
        "upper"
      )
) {
  
  which <- match.arg(
    which
  )
  
  split <- stringr::str_split_fixed(
    x,
    "\\s+-\\s+",
    n =
      2
  )
  
  if (
    which ==
    "lower"
  ) {
    as.numeric(
      split[
        ,
        1
      ]
    )
  } else {
    as.numeric(
      split[
        ,
        2
      ]
    )
  }
}


stage_alt_points <- bind_rows(
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cu: Nauplii / Adult",
      Method =
        "Complete-case",
      Ratio =
        Cu_Nauplii_Adult,
      Range =
        Cu_range
    ),
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cd: Nauplii / Adult",
      Method =
        "Complete-case",
      Ratio =
        Cd_Nauplii_Adult,
      Range =
        Cd_range
    ),
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Stage contrasts",
      Method =
        "Complete-case",
      Ratio =
        Cd_vs_Cu_RoR,
      Range =
        RoR_range
    ),
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cu: Nauplii / Adult",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cu_Nauplii_Adult,
      Range =
        Cu_range
    ),
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cd: Nauplii / Adult",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cd_Nauplii_Adult,
      Range =
        Cd_range
    ),
  
  stage_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Stage contrasts",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cd_vs_Cu_RoR,
      Range =
        RoR_range
    )
) %>%
  dplyr::mutate(
    Lower =
      extract_range(
        Range,
        "lower"
      ),
    
    Upper =
      extract_range(
        Range,
        "upper"
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          stage_estimand_order
      ),
    
    Method =
      factor(
        Method,
        levels =
          rev(
            stage_method_order
          )
      ),
    
    Baseline =
      FALSE
  ) %>%
  dplyr::select(
    -Range
  )


stage_loro_ranges <- bind_rows(
  stage_loro_points,
  stage_loso_points
) %>%
  make_leave_one_ranges()


stage_variant_points <- bind_rows(
  stage_baseline_points,
  stage_alt_points
)


write_csv(
  stage_variant_points,
  file.path(
    plot_data_dir,
    "Figure_4_16_robustness_variants.csv"
  )
)

write_csv(
  stage_loro_points,
  file.path(
    plot_data_dir,
    "Figure_4_16_LORO_points.csv"
  )
)

write_csv(
  stage_loso_points,
  file.path(
    plot_data_dir,
    "Figure_4_16_LOSO_points.csv"
  )
)


stage_interaction_name <- "Interaction: Cd / Cu ratio of Stage contrasts"

F07_upper <- ggplot() +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +
  
  geom_segment(
    data = stage_loro_ranges %>%
      dplyr::filter(as.character(Estimand) != stage_interaction_name),
    aes(x = min_ratio, xend = max_ratio, y = Method, yend = Method),
    colour = "grey78",
    linewidth = 3.0,
    lineend = "round"
  ) +
  
  geom_point(
    data = stage_loro_points %>%
      dplyr::filter(as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260828),
    colour = "grey52",
    alpha = 0.70,
    size = 1.7
  ) +
  
  geom_point(
    data = stage_loso_points %>%
      dplyr::filter(as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260829),
    colour = dark_grey,
    alpha = 0.82,
    size = 1.8,
    shape = 17
  ) +
  
  geom_errorbar(
    data = stage_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.15,
    colour = dark_grey,
    linewidth = 0.72
  ) +
  
  geom_point(
    data = stage_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = dark_grey,
    stroke = 0.7,
    size = 2.5
  ) +
  
  geom_errorbar(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.16,
    colour = blue,
    linewidth = 1.05
  ) +
  
  geom_point(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != stage_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = blue,
    stroke = 1.0,
    size = 3.0
  ) +
  
  geom_text(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != stage_interaction_name) %>%
      dplyr::mutate(Label = sprintf("%.2f", Ratio), LabelX = Ratio * 1.10),
    aes(x = LabelX, y = Method, label = Label),
    hjust = 0,
    vjust = -0.8,
    colour = blue,
    fontface = "bold",
    size = 3.4,
    show.legend = FALSE
  ) +
  
  scale_x_log10(
    breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
    labels = c("0.05", "0.1", "0.2", "0.5", "1", "2", "5", "10")
  ) +
  
  scale_y_discrete(
    limits = rev(stage_method_order),
    labels = stage_method_labels,
    drop = FALSE
  ) +
  
  facet_wrap(
    ~ Estimand,
    ncol = 1,
    scales = "fixed",
    strip.position = "top"
  ) +
  
  labs(
    title = "Sensitivity of the developmental-stage comparisons",
    subtitle = "96-h Cu-Cd developmental-stage mixed-effects model",
    x = "Nauplii / Adult LC50 ratio (log scale)",
    y = NULL,
    caption = "1 = equal LC50 between stages."
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    panel.grid.major.x = element_line(colour = grid_minor, linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5, colour = "black"),
    strip.placement = "outside",
    panel.spacing.y = unit(0.85, "lines"),
    plot.caption = element_text(hjust = 0, size = 9.0),
    plot.margin = margin(t = 12, r = 24, b = 8, l = 12)
  )

F07_interaction <- ggplot() +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +
  
  geom_segment(
    data = stage_loro_ranges %>%
      dplyr::filter(as.character(Estimand) == stage_interaction_name),
    aes(x = min_ratio, xend = max_ratio, y = Method, yend = Method),
    colour = "grey78",
    linewidth = 3.0,
    lineend = "round"
  ) +
  
  geom_point(
    data = stage_loro_points %>%
      dplyr::filter(as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260828),
    colour = "grey52",
    alpha = 0.70,
    size = 1.7
  ) +
  
  geom_point(
    data = stage_loso_points %>%
      dplyr::filter(as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260829),
    colour = dark_grey,
    alpha = 0.82,
    size = 1.8,
    shape = 17
  ) +
  
  geom_errorbar(
    data = stage_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.15,
    colour = dark_grey,
    linewidth = 0.72
  ) +
  
  geom_point(
    data = stage_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = dark_grey,
    stroke = 0.7,
    size = 2.5
  ) +
  
  geom_errorbar(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.16,
    colour = interaction_colour,
    linewidth = 1.05
  ) +
  
  geom_point(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == stage_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 23,
    fill = "white",
    colour = interaction_colour,
    stroke = 1.0,
    size = 3.2
  ) +
  
  geom_text(
    data = stage_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == stage_interaction_name) %>%
      dplyr::mutate(Label = sprintf("%.2f", Ratio), LabelX = Ratio * 1.10),
    aes(x = LabelX, y = Method, label = Label),
    hjust = 0,
    vjust = -0.8,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.4,
    show.legend = FALSE
  ) +
  
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20")
  ) +
  
  scale_y_discrete(
    limits = rev(stage_method_order),
    labels = stage_method_labels,
    drop = FALSE
  ) +
  
  facet_wrap(
    ~ Estimand,
    ncol = 1,
    scales = "fixed",
    strip.position = "top"
  ) +
  
  labs(
    x = "Ratio of ratios (log scale)",
    y = NULL,
    caption = paste0(
      "RoR = ratio of ratios; 1 = the developmental-stage contrast is the same in Cu and Cd.\n",
      "Grey bands show point-estimate ranges when one publication (LORO) or one species group (LOSO) is omitted at a time.\n",
      "CC = records with both temperature and salinity reported.\n",
      "CC + T + S = the same subset with temperature and salinity added to the model."
    )
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.grid.major.x = element_line(colour = grid_minor, linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5, colour = interaction_colour),
    strip.placement = "outside",
    plot.caption = element_text(hjust = 0, size = 9.0),
    plot.margin = margin(t = 2, r = 24, b = 12, l = 12)
  )

save_vertical_figure(
  upper_plot = F07_upper,
  lower_plot = F07_interaction,
  filename = "Figure_4_16_stage_model_robustness",
  width = 9.2,
  height = 10.8,
  upper_fraction = 0.67
)

F06 <- list(
  contrasts = F07_upper,
  interaction = F07_interaction
)


# ============================================================
# 12. ORDER DATA PREPARATION
# ============================================================

order_raw <- order_obj$order_model_data %>%
  as_tibble() %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Order =
      factor(
        Order,
        levels =
          c(
            "Calanoida",
            "Harpacticoida"
          )
      )
  )


order_x_map <- tribble(
  ~Metal, ~Order, ~x,
  "Cu", "Calanoida", 1,
  "Cu", "Harpacticoida", 3,
  "Cd", "Calanoida", 6,
  "Cd", "Harpacticoida", 8
) %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Order =
      factor(
        Order,
        levels =
          c(
            "Calanoida",
            "Harpacticoida"
          )
      )
  )


order_raw <- order_raw %>%
  dplyr::left_join(
    order_x_map,
    by =
      c(
        "Metal",
        "Order"
      )
  )

stopifnot(
  !any(
    is.na(
      order_raw$x
    )
  )
)


order_emm <- read_csv(
  file.path(
    "outputs",
    "13_fit_order_model",
    "01_model",
    "order_estimated_marginal_means.csv"
  ),
  show_col_types =
    FALSE
) %>%
  dplyr::mutate(
    Metal =
      factor(
        Metal,
        levels =
          c(
            "Cu",
            "Cd"
          )
      ),
    
    Order =
      factor(
        Order,
        levels =
          c(
            "Calanoida",
            "Harpacticoida"
          )
      ),
    
    estimate =
      exp(
        emmean
      ),
    
    lower =
      exp(
        lower.CL
      ),
    
    upper =
      exp(
        upper.CL
      )
  ) %>%
  dplyr::left_join(
    order_x_map,
    by =
      c(
        "Metal",
        "Order"
      )
  )


order_counts <- order_raw %>%
  dplyr::group_by(Metal, Order) %>%
  dplyr::summarise(
    Results = n(),
    References = dplyr::n_distinct(Reference_ID),
    .groups = "drop"
  ) %>%
  dplyr::left_join(order_x_map, by = c("Metal", "Order")) %>%
  dplyr::arrange(x)

order_x_labels <- paste0(
  as.character(order_counts$Metal),
  " | ",
  as.character(order_counts$Order),
  "\n",
  order_counts$Results,
  " records | ",
  order_counts$References,
  " publications"
)


write_csv(
  order_raw,
  file.path(
    plot_data_dir,
    "Figure_4_10_raw_results.csv"
  )
)

write_csv(
  order_emm,
  file.path(
    plot_data_dir,
    "Figure_4_10_model_estimates.csv"
  )
)


# ============================================================
# 13. F07 ORDER RESULTS + MODEL ESTIMATES
# ============================================================

F07 <- make_raw_model_plot(
  raw_data = order_raw,
  model_data = order_emm %>%
    dplyr::transmute(Metal, Order, x, estimate, lower, upper),
  x_breaks = order_counts$x,
  x_labels = order_x_labels,
  title = "96-h LC50 by Taxonomic Order",
  subtitle = paste0(
    "96-h Cu-Cd taxonomic-order mixed-effects model\n",
    "Grey points = LC50 records | grey band = 95% confidence interval | blue line = model estimate"
  ),
  caption = paste0(
    "Calanoida and Harpacticoida model estimates are shown separately for Cu and Cd.\n",
    "Lower LC50 indicates greater acute toxicity."
  ),
  y_title = expression(
    LC[50]~(mu*"mol metal L"^{-1}*"; log scale")
  ),
  width = 9.5,
  height = 7.6
)

save_final_figure(
  F07,
  "Figure_4_10_order_results_and_model_estimates",
  width = 10.8,
  height = 7.8
)


# ============================================================
# 14. F08A TAXONOMIC-ORDER MODEL CONTRASTS
# ============================================================
#
# LOCKED DESIGN:
# - upper panel: metal-specific Harpacticoida / Calanoida contrasts
# - lower mini-panel: formal Metal x Order interaction
# - interaction is purple + diamond and explicitly labelled
# - RoR is defined in the caption

order_contrasts <- read_csv(
  file.path(
    "outputs",
    "13_fit_order_model",
    "01_model",
    "order_metal_specific_contrasts.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Metal = factor(
      Metal,
      levels = rev(
        c(
          "Cu",
          "Cd"
        )
      )
    ),
    ratio_label = sprintf(
      "%.2f",
      LC50_ratio_Harp_Cal
    )
  )

order_base_for_figure <- order_obj$baseline_metrics %>%
  as_tibble()

stopifnot(
  nrow(
    order_base_for_figure
  ) == 1
)

order_interaction <- tibble(
  Estimate =
    order_base_for_figure$ratio_of_ratios,
  Lower =
    order_base_for_figure$interaction_lower95,
  Upper =
    order_base_for_figure$interaction_upper95,
  Label =
    "Interaction: Cd/Cu RoR"
)

write_csv(
  order_contrasts,
  file.path(
    plot_data_dir,
    "Figure_4_11_order_contrasts.csv"
  )
)

write_csv(
  order_interaction,
  file.path(
    plot_data_dir,
    "Figure_4_11_order_interaction.csv"
  )
)

F09_main <- ggplot(
  order_contrasts,
  aes(
    x = LC50_ratio_Harp_Cal,
    y = Metal
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.72
  ) +
  
  geom_errorbar(
    aes(
      xmin = ratio_lower_95,
      xmax = ratio_upper_95
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.90,
    colour = blue
  ) +
  
  geom_point(
    size = 3.5,
    shape = 21,
    fill = "white",
    colour = blue,
    stroke = 0.9
  ) +
  
  geom_text(
    aes(
      label = ratio_label
    ),
    nudge_y = 0.13,
    colour = blue,
    fontface = "bold",
    size = 3.6
  ) +
  
  scale_x_log10(
    breaks = c(
      0.5,
      1,
      2,
      5,
      10,
      20,
      50
    ),
    labels = c(
      "0.5",
      "1",
      "2",
      "5",
      "10",
      "20",
      "50"
    ),
    limits = c(
      0.45,
      55
    )
  ) +
  
  annotate(
    "text",
    x = 1,
    y = 2.38,
    label = "Equal LC50",
    hjust = 0.5,
    vjust = 0,
    size = 3.0,
    colour = "grey30"
  ) +
  
  labs(
    title =
      "Model-estimated LC50 ratios between taxonomic orders",
    
    subtitle =
      "Model-estimated Harpacticoida / Calanoida LC50 ratios with 95% confidence intervals",
    
    x =
      "Harpacticoida / Calanoida LC50 ratio (log scale)",
    
    y =
      NULL,
    
    caption =
      "Ratio < 1: Harpacticoida have lower LC50   |   Ratio > 1: Harpacticoida have higher LC50"
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_ratio +
  
  theme(
    plot.margin = margin(
      t = 12,
      r = 18,
      b = 8,
      l = 12
    )
  )

F09_interaction <- ggplot(
  order_interaction,
  aes(
    x = Estimate,
    y = 1
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.62
  ) +
  
  geom_errorbar(
    aes(
      xmin = Lower,
      xmax = Upper
    ),
    orientation = "y",
    width = 0.11,
    linewidth = 0.82,
    colour = interaction_colour
  ) +
  
  geom_point(
    shape = 23,
    size = 3.5,
    fill = interaction_colour,
    colour = interaction_colour
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.2f", Estimate)
    ),
    nudge_y = 0.13,
    hjust = 0.5,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.15
  ) +
  
  scale_x_log10(
    breaks = c(
      0.1,
      0.2,
      0.5,
      1,
      2
    ),
    labels = c(
      "0.1",
      "0.2",
      "0.5",
      "1",
      "2"
    ),
    limits = c(
      0.12,
      2.3
    )
  ) +
  
  scale_y_continuous(
    breaks = NULL
  ) +
  
  labs(
    x =
      "Cd / Cu ratio of taxonomic-order contrasts (log scale)",
    
    y =
      NULL,
    
    caption =
      paste0(
        "RoR = ratio of ratios; it compares the Cd Harpacticoida / Calanoida contrast with the corresponding Cu contrast. ",
        "RoR = 1 means the taxonomic-order contrast is the same in Cu and Cd."
      )
  ) +
  
  theme_ratio +
  
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_text(
      hjust = 0,
      size = 8.8
    ),
    plot.margin = margin(
      t = 2,
      r = 18,
      b = 12,
      l = 12
    )
  )

# Present the two metal-specific contrasts and the formal interaction in one
# compact forest plot, matching the visual logic of pooled Figure F02.

order_model_comparisons <- bind_rows(
  order_contrasts %>%
    dplyr::transmute(
      Comparison = paste0(as.character(Metal), ": Harpacticoida / Calanoida"),
      Estimate = LC50_ratio_Harp_Cal,
      Lower = ratio_lower_95,
      Upper = ratio_upper_95,
      Estimand_type = "Metal-specific contrast"
    ),
  order_interaction %>%
    dplyr::transmute(
      Comparison = "Interaction: Cd / Cu RoR",
      Estimate,
      Lower,
      Upper,
      Estimand_type = "Interaction"
    )
) %>%
  dplyr::mutate(
    Comparison = factor(
      Comparison,
      levels = rev(c(
        "Cu: Harpacticoida / Calanoida",
        "Cd: Harpacticoida / Calanoida",
        "Interaction: Cd / Cu RoR"
      ))
    ),
    ratio_label = sprintf("%.2f", Estimate),
    y_position = case_when(
      as.character(Comparison) == "Cu: Harpacticoida / Calanoida" ~ 3.0,
      as.character(Comparison) == "Cd: Harpacticoida / Calanoida" ~ 2.0,
      TRUE ~ 0.65
    )
  )

F08A <- ggplot(
  order_model_comparisons,
  aes(x = Estimate, y = y_position)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.72,
    colour = "grey35"
  ) +
  geom_errorbar(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    aes(xmin = Lower, xmax = Upper),
    orientation = "y",
    width = 0.16,
    linewidth = 0.95,
    colour = blue
  ) +
  geom_point(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    size = 3.7,
    shape = 21,
    stroke = 0.9,
    fill = "white",
    colour = blue
  ) +
  geom_text(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Metal-specific contrast"),
    aes(label = ratio_label),
    nudge_y = 0.16,
    colour = blue,
    fontface = "bold",
    size = 3.7
  ) +
  geom_errorbar(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    aes(xmin = Lower, xmax = Upper),
    orientation = "y",
    width = 0.16,
    linewidth = 0.95,
    colour = interaction_colour
  ) +
  geom_point(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    size = 3.9,
    shape = 23,
    stroke = 0.9,
    fill = "white",
    colour = interaction_colour
  ) +
  geom_text(
    data = order_model_comparisons %>%
      dplyr::filter(Estimand_type == "Interaction"),
    aes(label = ratio_label),
    nudge_y = 0.16,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.7
  ) +
  scale_x_log10(
    breaks = c(0.2, 0.5, 1, 2, 5, 10, 20, 50),
    labels = c("0.2", "0.5", "1", "2", "5", "10", "20", "50"),
    limits = c(0.12, 55)
  ) +
  scale_y_continuous(
    breaks = c(3.0, 2.0, 0.65),
    labels = c(
      "Cu: Harpacticoida / Calanoida",
      "Cd: Harpacticoida / Calanoida",
      "Interaction: Cd / Cu RoR"
    )
  ) +
  labs(
    title = "Taxonomic-order comparisons and their variation between metals",
    subtitle = "Ratios and 95% confidence intervals from the 96-h Cu-Cd mixed-effects model",
    x = "Ratio estimate (log scale)",
    y = NULL,
    caption = paste0(
      "Blue: Harpacticoida / Calanoida LC50 ratio within each metal. Purple: Cd / Cu ratio of those contrasts (RoR).\n",
      "The dashed line at 1 indicates equal LC50 for the blue contrasts and no interaction for the purple RoR."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    plot.margin = margin(t = 12, r = 18, b = 30, l = 12)
  )

save_final_figure(
  F08A,
  "Figure_4_11_order_contrasts_and_interaction",
  width = 8.4,
  height = 5.8
)


# ============================================================
# 15. F08 ORDER POST-MODEL ROBUSTNESS
# ============================================================

order_method_order <- c(
  "96-h Cu-Cd\nOrder LMM",
  "LORO",
  "LOSO",
  "Complete-case",
  "Temperature +\nsalinity adjusted"
)

order_method_labels <- c(
  "96-h Cu-Cd\nOrder LMM" = "96-h Cu-Cd\nOrder LMM",
  "LORO" = "LORO",
  "LOSO" = "LOSO",
  "Complete-case" = "CC\n(T + S reported)",
  "Temperature +\nsalinity adjusted" = "CC + T + S\n(T + S in model)"
)

order_estimand_order <- c(
  "Cu: Harpacticoida / Calanoida",
  "Cd: Harpacticoida / Calanoida",
  "Interaction: Cd / Cu ratio of Order contrasts"
)


order_base <- order_obj$baseline_metrics %>%
  as_tibble()

stopifnot(
  nrow(
    order_base
  ) ==
    1
)


order_baseline_points <- tibble(
  Estimand =
    factor(
      order_estimand_order,
      levels =
        order_estimand_order
    ),
  
  Method =
    factor(
      rep(
        "96-h Cu-Cd\nOrder LMM",
        3
      ),
      levels =
        rev(
          order_method_order
        )
    ),
  
  Baseline =
    TRUE,
  
  Ratio =
    c(
      order_base$Cu_ratio,
      order_base$Cd_ratio,
      order_base$ratio_of_ratios
    ),
  
  Lower =
    c(
      order_base$Cu_lower95,
      order_base$Cd_lower95,
      order_base$interaction_lower95
    ),
  
  Upper =
    c(
      order_base$Cu_upper95,
      order_base$Cd_upper95,
      order_base$interaction_upper95
    )
)


order_loro <- order_obj$loro_results %>%
  as_tibble() %>%
  dplyr::filter(
    !is.na(
      Cu_ratio
    )
  )

order_loso <- order_obj$loso_results %>%
  as_tibble() %>%
  dplyr::filter(
    !is.na(
      Cu_ratio
    )
  )


order_loro_points <- bind_rows(
  
  order_loro %>%
    dplyr::transmute(
      Estimand =
        "Cu: Harpacticoida / Calanoida",
      Ratio =
        Cu_ratio
    ),
  
  order_loro %>%
    dplyr::transmute(
      Estimand =
        "Cd: Harpacticoida / Calanoida",
      Ratio =
        Cd_ratio
    ),
  
  order_loro %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Order contrasts",
      Ratio =
        ratio_of_ratios
    )
) %>%
  dplyr::mutate(
    Estimand =
      factor(
        Estimand,
        levels =
          order_estimand_order
      ),
    
    Method =
      factor(
        "LORO",
        levels =
          rev(
            order_method_order
          )
      )
  )


order_loso_points <- bind_rows(
  
  order_loso %>%
    dplyr::transmute(
      Estimand =
        "Cu: Harpacticoida / Calanoida",
      Ratio =
        Cu_ratio
    ),
  
  order_loso %>%
    dplyr::transmute(
      Estimand =
        "Cd: Harpacticoida / Calanoida",
      Ratio =
        Cd_ratio
    ),
  
  order_loso %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Order contrasts",
      Ratio =
        ratio_of_ratios
    )
) %>%
  dplyr::mutate(
    Estimand =
      factor(
        Estimand,
        levels =
          order_estimand_order
      ),
    
    Method =
      factor(
        "LOSO",
        levels =
          rev(
            order_method_order
          )
      )
  )


order_syn <- read_csv(
  file.path(
    "outputs",
    "13_fit_order_model",
    "03_robustness",
    "order_robustness_synthesis.csv"
  ),
  show_col_types =
    FALSE
)


order_alt_points <- bind_rows(
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cu: Harpacticoida / Calanoida",
      Method =
        "Complete-case",
      Ratio =
        Cu_Harp_Cal,
      Range =
        Cu_range
    ),
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cd: Harpacticoida / Calanoida",
      Method =
        "Complete-case",
      Ratio =
        Cd_Harp_Cal,
      Range =
        Cd_range
    ),
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity complete-case"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Order contrasts",
      Method =
        "Complete-case",
      Ratio =
        Cd_vs_Cu_RoR,
      Range =
        RoR_range
    ),
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cu: Harpacticoida / Calanoida",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cu_Harp_Cal,
      Range =
        Cu_range
    ),
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Cd: Harpacticoida / Calanoida",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cd_Harp_Cal,
      Range =
        Cd_range
    ),
  
  order_syn %>%
    dplyr::filter(
      Analysis ==
        "Temperature/salinity adjusted"
    ) %>%
    dplyr::transmute(
      Estimand =
        "Interaction: Cd / Cu ratio of Order contrasts",
      Method =
        "Temperature +\nsalinity adjusted",
      Ratio =
        Cd_vs_Cu_RoR,
      Range =
        RoR_range
    )
) %>%
  dplyr::mutate(
    Lower =
      extract_range(
        Range,
        "lower"
      ),
    
    Upper =
      extract_range(
        Range,
        "upper"
      ),
    
    Estimand =
      factor(
        Estimand,
        levels =
          order_estimand_order
      ),
    
    Method =
      factor(
        Method,
        levels =
          rev(
            order_method_order
          )
      ),
    
    Baseline =
      FALSE
  ) %>%
  dplyr::select(
    -Range
  )


order_loro_ranges <- bind_rows(
  order_loro_points,
  order_loso_points
) %>%
  make_leave_one_ranges()


order_variant_points <- bind_rows(
  order_baseline_points,
  order_alt_points
)


write_csv(
  order_variant_points,
  file.path(
    plot_data_dir,
    "Figure_4_17_robustness_variants.csv"
  )
)

write_csv(
  order_loro_points,
  file.path(
    plot_data_dir,
    "Figure_4_17_LORO_points.csv"
  )
)

write_csv(
  order_loso_points,
  file.path(
    plot_data_dir,
    "Figure_4_17_LOSO_points.csv"
  )
)


order_interaction_name <- "Interaction: Cd / Cu ratio of Order contrasts"

F10_upper <- ggplot() +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +
  
  geom_segment(
    data = order_loro_ranges %>%
      dplyr::filter(as.character(Estimand) != order_interaction_name),
    aes(x = min_ratio, xend = max_ratio, y = Method, yend = Method),
    colour = "grey78",
    linewidth = 3.0,
    lineend = "round"
  ) +
  
  geom_point(
    data = order_loro_points %>%
      dplyr::filter(as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260828),
    colour = "grey52",
    alpha = 0.70,
    size = 1.7
  ) +
  
  geom_point(
    data = order_loso_points %>%
      dplyr::filter(as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260829),
    colour = dark_grey,
    alpha = 0.82,
    size = 1.8,
    shape = 17
  ) +
  
  geom_errorbar(
    data = order_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.15,
    colour = dark_grey,
    linewidth = 0.72
  ) +
  
  geom_point(
    data = order_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = dark_grey,
    stroke = 0.7,
    size = 2.5
  ) +
  
  geom_errorbar(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.16,
    colour = blue,
    linewidth = 1.05
  ) +
  
  geom_point(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != order_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = blue,
    stroke = 1.0,
    size = 3.0
  ) +
  
  geom_text(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) != order_interaction_name) %>%
      dplyr::mutate(Label = sprintf("%.2f", Ratio), LabelX = Ratio * 1.10),
    aes(x = LabelX, y = Method, label = Label),
    hjust = 0,
    vjust = -0.8,
    colour = blue,
    fontface = "bold",
    size = 3.4,
    show.legend = FALSE
  ) +
  
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
    labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")
  ) +
  
  scale_y_discrete(
    limits = rev(order_method_order),
    labels = order_method_labels,
    drop = FALSE
  ) +
  
  facet_wrap(
    ~ Estimand,
    ncol = 1,
    scales = "fixed",
    strip.position = "top"
  ) +
  
  labs(
    title = "Sensitivity of the taxonomic-order comparisons",
    subtitle = "96-h Cu-Cd taxonomic-order mixed-effects model",
    x = "Harpacticoida / Calanoida LC50 ratio (log scale)",
    y = NULL,
    caption = "1 = equal LC50 between orders."
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    panel.grid.major.x = element_line(colour = grid_minor, linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5, colour = "black"),
    strip.placement = "outside",
    panel.spacing.y = unit(0.85, "lines"),
    plot.caption = element_text(hjust = 0, size = 9.0),
    plot.margin = margin(t = 12, r = 24, b = 8, l = 12)
  )

F10_interaction <- ggplot() +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.70
  ) +
  
  geom_segment(
    data = order_loro_ranges %>%
      dplyr::filter(as.character(Estimand) == order_interaction_name),
    aes(x = min_ratio, xend = max_ratio, y = Method, yend = Method),
    colour = "grey78",
    linewidth = 3.0,
    lineend = "round"
  ) +
  
  geom_point(
    data = order_loro_points %>%
      dplyr::filter(as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260828),
    colour = "grey52",
    alpha = 0.70,
    size = 1.7
  ) +
  
  geom_point(
    data = order_loso_points %>%
      dplyr::filter(as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, y = Method),
    position = position_jitter(width = 0, height = 0.065, seed = 20260829),
    colour = dark_grey,
    alpha = 0.82,
    size = 1.8,
    shape = 17
  ) +
  
  geom_errorbar(
    data = order_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.15,
    colour = dark_grey,
    linewidth = 0.72
  ) +
  
  geom_point(
    data = order_variant_points %>%
      dplyr::filter(!Baseline, as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 21,
    fill = "white",
    colour = dark_grey,
    stroke = 0.7,
    size = 2.5
  ) +
  
  geom_errorbar(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, xmin = Lower, xmax = Upper, y = Method),
    orientation = "y",
    width = 0.16,
    colour = interaction_colour,
    linewidth = 1.05
  ) +
  
  geom_point(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == order_interaction_name),
    aes(x = Ratio, y = Method),
    shape = 23,
    fill = "white",
    colour = interaction_colour,
    stroke = 1.0,
    size = 3.2
  ) +
  
  geom_text(
    data = order_variant_points %>%
      dplyr::filter(Baseline, as.character(Estimand) == order_interaction_name) %>%
      dplyr::mutate(Label = sprintf("%.2f", Ratio), LabelX = Ratio * 1.10),
    aes(x = LabelX, y = Method, label = Label),
    hjust = 0,
    vjust = -0.8,
    colour = interaction_colour,
    fontface = "bold",
    size = 3.4,
    show.legend = FALSE
  ) +
  
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5),
    labels = c("0.1", "0.2", "0.5", "1", "2", "5")
  ) +
  
  scale_y_discrete(
    limits = rev(order_method_order),
    labels = order_method_labels,
    drop = FALSE
  ) +
  
  facet_wrap(
    ~ Estimand,
    ncol = 1,
    scales = "fixed",
    strip.position = "top"
  ) +
  
  labs(
    x = "Ratio of ratios (log scale)",
    y = NULL,
    caption = paste0(
      "RoR = ratio of ratios; 1 = the taxonomic-order contrast is the same in Cu and Cd.\n",
      "Grey bands show point-estimate ranges when one publication (LORO) or one species group (LOSO) is omitted at a time.\n",
      "CC = records with both temperature and salinity reported.\n",
      "CC + T + S = the same subset with temperature and salinity added to the model."
    )
  ) +
  
  coord_cartesian(clip = "off") +
  theme_ratio +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.grid.major.x = element_line(colour = grid_minor, linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0, size = 10.5, colour = interaction_colour),
    strip.placement = "outside",
    plot.caption = element_text(hjust = 0, size = 9.0),
    plot.margin = margin(t = 2, r = 24, b = 12, l = 12)
  )

save_vertical_figure(
  upper_plot = F10_upper,
  lower_plot = F10_interaction,
  filename = "Figure_4_17_order_model_robustness",
  width = 9.2,
  height = 10.8,
  upper_fraction = 0.67
)

F08 <- list(
  contrasts = F10_upper,
  interaction = F10_interaction
)


# ============================================================
# 16. F09 WITHIN-REFERENCE DEVELOPMENTAL-STAGE COMPARISONS
# ============================================================
#
# Updated combined representation:
# - original exact-context ECOTOX series are preserved;
# - one new source-verified WoS series (Kadiene et al. 2019) is added;
# - one panel = one exact comparable developmental-stage series;
# - no connecting lines and no pooled regression;
# - common log-LC50 y-scale.

stage_within_file <- file.path(
  "outputs",
  "09_compare_stages_within_references",
  "02_tables",
  "stage_within_combined_table_data.csv"
)

stopifnot(file.exists(stage_within_file))

stage_within_summary <- read_csv(
  stage_within_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Source_label = recode(
      Source_Origin,
      "ECOTOX" = "ECOTOX",
      "WoS_supplemental" = "Supplementary literature"
    ),
    Panel =
      paste0(
        Reference_ID,
        " | ",
        Metal,
        "\n",
        Species
      )
  )

stage_within_points <- stage_within_summary %>%
  dplyr::select(
    Panel,
    Series_ID,
    Reference_ID,
    Species,
    Metal,
    Source_Origin,
    Source_label,
    Stage_profile
  ) %>%
  separate_rows(
    Stage_profile,
    sep = " \\| "
  ) %>%
  separate(
    Stage_profile,
    into = c("Stage_raw", "LC50_umol_L"),
    sep = ": ",
    convert = TRUE,
    extra = "merge"
  ) %>%
  dplyr::mutate(
    Stage_display =
      recode(
        Stage_raw,
        "Nauplius 1 d" = "Nauplius\n1 d",
        "Nauplius 5 d" = "Nauplius\n5 d",
        "Copepodid 10 d" = "Copepodid\n10 d",
        "Adult ovigerous band" = "Adult\nband",
        "Adult ovigerous sac" = "Adult\nsac"
      ),
    
    Stage_display =
      factor(
        Stage_display,
        levels = c(
          "Egg",
          "Nauplius\n1 d",
          "Nauplius\n5 d",
          "Nauplii",
          "Copepodid\n10 d",
          "Copepodid",
          "Copepodite",
          "Adult\nband",
          "Adult\nsac",
          "Adult"
        )
      ),
    
    Panel =
      factor(
        Panel,
        levels = stage_within_summary$Panel
      )
  )

stopifnot(
  dplyr::n_distinct(stage_within_points$Series_ID) == 10,
  dplyr::n_distinct(stage_within_points$Reference_ID) == 7,
  dplyr::n_distinct(stage_within_points$Species) == 6,
  all(stage_within_points$LC50_umol_L > 0),
  all(!is.na(stage_within_points$Stage_display))
)

stage_within_limits <- c(
  min(stage_within_points$LC50_umol_L) * 0.82,
  max(stage_within_points$LC50_umol_L) * 1.18
)

stage_within_breaks <- c(0.2, 1, 5, 20, 100, 500)
stage_within_breaks <- stage_within_breaks[stage_within_breaks >= stage_within_limits[1] &
                                             stage_within_breaks <= stage_within_limits[2]]

write_csv(
  stage_within_points,
  file.path(
    plot_data_dir,
    "Figure_4_09_points.csv"
  )
)

F09 <- ggplot(
  stage_within_points,
  aes(
    x = Stage_display,
    y = LC50_umol_L
  )
) +
  geom_point(
    aes(
      shape = Source_label,
      colour = Source_label
    ),
    size = 2.6,
    alpha = 0.95
  ) +
  facet_wrap(
    ~ Panel,
    ncol = 2,
    scales = "free_x"
  ) +
  scale_y_log10(
    limits = stage_within_limits,
    breaks = stage_within_breaks,
    labels = log_tick_labels
  ) +
  scale_shape_manual(
    values = c(
      "ECOTOX" = 16,
      "Supplementary literature" = 17
    )
  ) +
  scale_colour_manual(
    values = c(
      "ECOTOX" = blue,
      "Supplementary literature" = wos_colour
    )
  ) +
  labs(
    title =
      "LC50 across developmental stages",
    subtitle =
      "10 stage series | 7 publications | 6 species groups",
    x =
      "Developmental stage",
    y =
      expression(
        LC[50]~(mu*mol~metal~L^{-1})~"(log scale)"
      ),
    shape = "Data source",
    colour = "Data source",
    caption =
      paste0(
        "Each panel compares stages within one publication under matched test conditions. ",
        "In nine of ten series, the first stage shown had lower LC50 than the last stage shown.\n",
        "Lower LC50 indicates greater acute sensitivity."
      )
  ) +
  theme_final +
  theme(
    legend.position = "top",
    axis.text.x =
      element_text(
        size = 8.2,
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        lineheight = 0.95
      ),
    strip.text =
      element_text(
        face = "bold",
        size = 10.0,
        lineheight = 1.00
      ),
    panel.spacing =
      grid::unit(1.10, "lines"),
    plot.caption =
      element_text(
        size = 9.2,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(t = 10)
      )
  )

save_final_figure(
  F09,
  "Figure_4_09_within_reference_stage_comparisons",
  width = 10.8,
  height = 12.5
)


# ============================================================
# 17. F10 WITHIN-REFERENCE SALINITY COMPARISONS
# ============================================================
#
# Direct environmental evidence, not a common-support diagnostic.
# No connecting lines: salinity levels are separate LC50 estimates.
# All retained comparisons are displayed together.

salinity_main <- read_csv(
  file.path(
    environmental_table_dir,
    "salinity_main_direct_values.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::mutate(Data_set = "Previously Retained")

salinity_additional <- read_csv(
  file.path(
    environmental_table_dir,
    "salinity_extended_E2332_values.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::mutate(Data_set = "Additional Screened")

salinity_direct <- bind_rows(
  salinity_main %>%
    dplyr::select(
      Comparison_ID,
      Reference_ID,
      Reference_Number,
      Species,
      Metal,
      Salinity,
      LC50_umol_L,
      Data_set
    ),
  salinity_additional %>%
    dplyr::select(
      Comparison_ID,
      Reference_ID,
      Reference_Number,
      Species,
      Metal,
      Salinity,
      LC50_umol_L,
      Data_set
    )
) %>%
  dplyr::mutate(
    Panel =
      paste0(
        "Publication ",
        Reference_Number,
        " | ",
        Metal,
        "\n",
        Species
      ),
    Panel =
      factor(
        Panel,
        levels = unique(Panel)
      )
  )

stopifnot(
  dplyr::n_distinct(salinity_direct$Comparison_ID) == 7,
  dplyr::n_distinct(salinity_direct$Reference_Number) == 5,
  dplyr::n_distinct(salinity_direct$Species) == 5,
  all(salinity_direct$LC50_umol_L > 0)
)

salinity_limits <- c(
  min(salinity_direct$LC50_umol_L) * 0.85,
  max(salinity_direct$LC50_umol_L) * 1.18
)

salinity_breaks <- make_125_breaks(
  salinity_limits[1],
  salinity_limits[2]
)

write_csv(
  salinity_direct,
  file.path(
    plot_data_dir,
    "Figure_4_13_points.csv"
  )
)

F10 <- ggplot(
  salinity_direct,
  aes(
    x = Salinity,
    y = LC50_umol_L
  )
) +
  geom_point(
    colour = blue,
    size = 3.0,
    alpha = 0.95
  ) +
  facet_wrap(
    ~ Panel,
    ncol = 3,
    axes = "all_x",
    axis.labels = "all_x"
  ) +
  scale_x_continuous(
    breaks = c(5, 10, 15, 20, 25, 30)
  ) +
  scale_y_log10(
    limits = salinity_limits,
    breaks = salinity_breaks,
    labels = log_tick_labels
  ) +
  labs(
    title =
      "LC50 across salinity levels",
    subtitle =
      "7 comparisons | 5 publications | 5 species groups",
    x =
      "Salinity",
    y =
      expression(
        LC[50]~(mu*mol~metal~L^{-1})~"(log scale)"
      ),
    caption =
      paste0(
        "Each panel shows LC50 at different salinities within the same publication. ",
        "LC50 increased, decreased or changed direction across salinity levels; no trend line was fitted.\n",
        "For publication 2332, temperature, pH and life stage were not reported. Its Zn comparison includes salinities 15 and 25, where the same zinc salt was used."
      )
  ) +
  theme_final +
  theme(
    strip.text =
      element_text(
        face = "bold",
        size = 10.0,
        lineheight = 1.00
      ),
    panel.spacing =
      grid::unit(0.85, "lines"),
    axis.text.x =
      element_text(size = 8.2),
    plot.caption =
      element_text(
        size = 9.2,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(t = 10)
      )
  )

save_final_figure(
  F10,
  "Figure_4_13_salinity_comparisons",
  width = 10.6,
  height = 7.6
)


# ============================================================
# 18. F11 WITHIN-REFERENCE TEMPERATURE COMPARISONS
# ============================================================
#
# Sparse direct evidence only. No connecting lines and no common slope.

temperature_direct <- read_csv(
  file.path(
    environmental_table_dir,
    "temperature_final_direct_values.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::filter(Include_clean_direct == "YES") %>%
  dplyr::mutate(
    Source_label = recode(
      Source_Origin,
      "ECOTOX" = "ECOTOX",
      "SOURCE_REVIEWED_PUBLICATION" = "Supplementary literature"
    ),
    Source_label = factor(
      Source_label,
      levels = c(
        "ECOTOX",
        "Supplementary literature"
      )
    ),
    Panel = paste0(
      Reference_label,
      " | ",
      Metal,
      "\n",
      Species
    ),
    Panel = factor(Panel, levels = unique(Panel))
  )

temperature_direction_check <- temperature_direct %>%
  dplyr::arrange(Temperature_C) %>%
  dplyr::group_by(Temperature_Comparison_ID) %>%
  dplyr::summarise(
    LC50_change = last(LC50_umol_L) - first(LC50_umol_L),
    .groups = "drop"
  )

stopifnot(
  nrow(temperature_direct) == 8,
  dplyr::n_distinct(temperature_direct$Temperature_Comparison_ID) == 3,
  dplyr::n_distinct(temperature_direct$Reference_ID) == 3,
  dplyr::n_distinct(temperature_direct$Species) == 2,
  dplyr::n_distinct(temperature_direct$Metal) == 3,
  all(temperature_direct$LC50_umol_L > 0),
  all(temperature_direction_check$LC50_change < 0)
)

write_csv(
  temperature_direct,
  file.path(
    plot_data_dir,
    "Figure_4_14_points.csv"
  )
)

F11 <- ggplot(
  temperature_direct,
  aes(
    x = Temperature_C,
    y = LC50_umol_L,
    shape = Source_label,
    colour = Source_label
  )
) +
  geom_point(
    size = 3.0,
    alpha = 0.95
  ) +
  facet_wrap(
    ~ Panel,
    ncol = 3,
    scales = "free",
    axes = "all_x",
    axis.labels = "all_x"
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4)
  ) +
  scale_y_log10(
    labels = log_tick_labels
  ) +
  scale_shape_manual(
    values = c(
      "ECOTOX" = 16,
      "Supplementary literature" = 17
    )
  ) +
  scale_colour_manual(
    values = c(
      "ECOTOX" = blue,
      "Supplementary literature" = wos_colour
    )
  ) +
  labs(
    title =
      "LC50 across temperatures",
    subtitle =
      "3 comparisons | 3 publications | 2 species groups | 3 metals",
    x =
      expression("Temperature ("*degree*C*")"),
    y =
      expression(
        LC[50]~(mu*mol~metal~L^{-1})~"(log scale)"
      ),
    shape =
      "Data source",
    colour =
      "Data source",
    caption =
      paste0(
        "Each panel compares temperatures within one publication. Vertical axes use different log scales; no overall temperature slope was fitted.\n",
        "LC50 was lower at the higher temperature in all three comparisons. The comparisons represent few publications and different test conditions.",
        ""
      )
  ) +
  theme_final +
  theme(
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(
      face = "bold",
      size = 10.0,
      lineheight = 1.00
    ),
    panel.spacing = grid::unit(0.85, "lines"),
    axis.text.x = element_text(size = 8.2),
    plot.caption = element_text(
      size = 9.2,
      hjust = 0,
      lineheight = 1.08,
      margin = margin(t = 10)
    )
  )

save_final_figure(
  F11,
  "Figure_4_14_temperature_comparisons",
  width = 10.8,
  height = 6.9
)


# ============================================================
# 19. F12 WITHIN-ASSAY EXPOSURE-DURATION COMPARISONS
# ============================================================
#
# Updated combined source-verified duration evidence.
# Each series is normalized to its shortest reported exposure.
# These are repeated LC50 observation times within source-verified
# acute assay trajectories, not independent randomized duration treatments.

duration_verified_file <- file.path(
  "outputs",
  "06_compare_exposure_duration",
  "02_tables",
  "duration_verified_plot_data_combined.csv"
)

stopifnot(file.exists(duration_verified_file))

duration_verified <- read_csv(
  duration_verified_file,
  show_col_types = FALSE
) %>%
  dplyr::rename(
    Comparison_ID = series_id
  ) %>%
  dplyr::mutate(
    Source_label = recode(
      Source_Origin,
      "ECOTOX" = "ECOTOX",
      "WoS_supplemental" = "Supplementary literature"
    ),
    Source_label = factor(
      Source_label,
      levels = c("ECOTOX", "Supplementary literature")
    ),
    Metal = factor(
      Metal,
      levels = c("Cu", "Cd", "Zn", "Hg", "Ni")
    )
  )

duration_final_check <- duration_verified %>%
  dplyr::group_by(Comparison_ID) %>%
  slice_max(
    Duration_days,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

stopifnot(
  dplyr::n_distinct(duration_verified$Comparison_ID) == 19,
  dplyr::n_distinct(duration_verified$Reference_ID) == 10,
  dplyr::n_distinct(duration_verified$Species) == 9,
  sum(duration_final_check$ln_relative_LC50 < 0) == 18,
  sum(duration_final_check$ln_relative_LC50 > 0) == 1
)

write_csv(
  duration_verified,
  file.path(
    plot_data_dir,
    "Figure_4_12_points.csv"
  )
)

duration_y_min <- floor(
  min(duration_verified$ln_relative_LC50, na.rm = TRUE) * 10
) / 10

duration_y_max <- max(
  0.08,
  ceiling(max(duration_verified$ln_relative_LC50, na.rm = TRUE) * 10) / 10
)

duration_panel_labels <- duration_verified %>%
  dplyr::distinct(Comparison_ID, Metal) %>%
  dplyr::count(Metal, name = "n_series") %>%
  dplyr::mutate(
    facet_label = paste0(
      Metal,
      " (",
      n_series,
      if_else(n_series == 1, " series)", " series)")
    )
  )

duration_labeller <- setNames(
  duration_panel_labels$facet_label,
  duration_panel_labels$Metal
)

F12 <- ggplot(
  duration_verified,
  aes(
    x = Duration_days,
    y = ln_relative_LC50,
    group = Comparison_ID,
    shape = Source_label,
    colour = Source_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.42,
    colour = dark_grey
  ) +
  geom_line(
    linewidth = 0.48,
    alpha = 0.34
  ) +
  geom_point(
    aes(
      shape = Source_label,
      colour = Source_label
    ),
    size = 2.25,
    alpha = 0.90
  ) +
  facet_wrap(
    ~ Metal,
    ncol = 2,
    labeller = as_labeller(duration_labeller)
  ) +
  scale_shape_manual(
    values = c(
      "ECOTOX" = 16,
      "Supplementary literature" = 17
    )
  ) +
  scale_colour_manual(
    values = c(
      "ECOTOX" = blue,
      "Supplementary literature" = wos_colour
    )
  ) +
  scale_x_continuous(
    breaks = 1:4
  ) +
  scale_y_continuous(
    limits = c(duration_y_min, duration_y_max)
  ) +
  labs(
    title =
      "LC50 across exposure durations",
    subtitle =
      "19 test series | 10 publications | 9 species groups",
    x =
      "Exposure duration (days)",
    y =
      "ln(LC50 / LC50 at shortest exposure)",
    shape =
      "Data source",
    colour =
      "Data source",
    caption =
      paste0(
        "Each line follows LC50 over observation times within the same toxicity test. Zero marks the shortest exposure; ",
        "values below 0 indicate lower LC50 at later observation times.\n",
        "Eighteen of 19 series showed lower LC50 at the longer exposure; ",
        "one Ni series showed a small change in the opposite direction."
      )
  ) +
  theme_final +
  theme(
    legend.position = "top",
    legend.justification = "left",
    panel.spacing = grid::unit(0.85, "lines"),
    strip.text =
      element_text(
        face = "bold",
        size = 11.0
      ),
    plot.caption =
      element_text(
        size = 9.2,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(t = 10)
      )
  )

save_final_figure(
  F12,
  "Figure_4_12_exposure_duration_comparisons",
  width = 10.8,
  height = 9.2
)


# ============================================================
# 20. APP_A1 pH CANDIDATE CONFOUNDING AUDIT
# ============================================================
#
# This is an evidence-structure figure, not a pH effect estimate.
# Each candidate context changes salinity together with pH.

ph_candidates <- read_csv(
  file.path(
    ph_table_dir,
    "ph_within_reference_candidate_rows.csv"
  ),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Panel = paste0(
      Reference_ID,
      " | ",
      Metal,
      "\n",
      Species
    ),
    Panel = factor(Panel, levels = unique(Panel)),
    Salinity_label = paste0("S = ", format(Salinity, trim = TRUE))
  )

ph_coverage <- read_csv(
  file.path(
    ph_table_dir,
    "ph_metadata_coverage_overall.csv"
  ),
  show_col_types = FALSE
)

wei_acidification <- read_csv(
  file.path(
    ph_table_dir,
    "wei2021_direct_acidification_values.csv"
  ),
  show_col_types = FALSE
)

wei_lc50_ratio <- wei_acidification %>%
  dplyr::arrange(pCO2_uatm) %>%
  dplyr::summarise(
    ratio = last(LC50_reported_mg_L) / first(LC50_reported_mg_L)
  ) %>%
  pull(ratio)

stopifnot(
  nrow(ph_candidates) == 9,
  dplyr::n_distinct(ph_candidates$Reference_ID) == 3,
  dplyr::n_distinct(ph_candidates$Species) == 2,
  all(ph_candidates$LC50_umol_L > 0),
  nrow(ph_coverage) == 1,
  ph_coverage$Results == 304,
  ph_coverage$pH_known == 89,
  nrow(wei_acidification) == 2,
  abs(wei_lc50_ratio - 0.755) < 0.002
)

write_csv(
  ph_candidates,
  file.path(
    plot_data_dir,
    "Figure_C_01_candidate_points.csv"
  )
)

write_csv(
  ph_coverage,
  file.path(
    plot_data_dir,
    "Figure_C_01_metadata_coverage.csv"
  )
)

write_csv(
  wei_acidification,
  file.path(
    plot_data_dir,
    "Figure_C_01_Wei2021_acidification_context.csv"
  )
)

APP_A1 <- ggplot(
  ph_candidates,
  aes(
    x = pH_final,
    y = LC50_umol_L
  )
) +
  geom_point(
    colour = blue,
    size = 3.0,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = Salinity_label),
    colour = dark_grey,
    size = 3.15,
    vjust = -0.85,
    check_overlap = TRUE
  ) +
  facet_wrap(
    ~ Panel,
    ncol = 3,
    scales = "free_x",
    axes = "all_x",
    axis.labels = "all_x"
  ) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  scale_y_log10(
    labels = log_tick_labels,
    expand = expansion(mult = c(0.12, 0.28))
  ) +
  labs(
    title =
      "LC50 across pH and salinity conditions",
    subtitle =
      "3 comparisons with simultaneous changes in pH and salinity",
    x =
      "Reported pH",
    y =
      expression(
        LC[50]~(mu*mol~metal~L^{-1})~"(log scale)"
      ),
    caption =
      paste0(
        "Labels show the salinity reported at each pH level.\n",
        "Salinity changed with pH in all three comparisons, so an independent pH effect could not be isolated. ",
        "The pCO2 experiment reported by Wei et al. (2021) was considered separately because verified exposure-water pH values were unavailable."
      )
  ) +
  theme_final +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10.0,
      lineheight = 1.00
    ),
    panel.spacing = grid::unit(0.85, "lines"),
    axis.text.x = element_text(size = 8.2),
    plot.caption = element_text(
      size = 9.2,
      hjust = 0,
      lineheight = 1.08,
      margin = margin(t = 10)
    )
  )

save_final_figure(
  APP_A1,
  "Figure_C_01_pH_salinity_confounding",
  width = 10.8,
  height = 6.8
)



# ============================================================
# 20B. FIGURES 4.5, 4.6 AND D.7: MATCHED-METAL COVARIANCE
# ============================================================
#
# IMPORTANT:
# The covariance figures intentionally retain the mature display logic
# from the final covariance workflow used in the thesis. Script 14
# supplies the frozen numerical results; this section only redraws
# those results with thesis-facing filenames.
#
# Groups reflect the thesis metal panel and structural support,
# never p-values.

covariance_dir <- file.path(
  "outputs",
  "14_analyze_metal_covariance"
)

covariance_overview_file <- file.path(
  covariance_dir,
  "02_tables",
  "metal_covariance_main_scope_overview.csv"
)

covariance_pair_values_file <- file.path(
  covariance_dir,
  "02_tables",
  "metal_covariance_pair_context_values.csv"
)

covariance_components_file <- file.path(
  covariance_dir,
  "02_tables",
  "metal_covariance_covariance_all_pairs_and_scopes.csv"
)

stopifnot(
  file.exists(covariance_overview_file),
  file.exists(covariance_pair_values_file),
  file.exists(covariance_components_file)
)

covariance_overview <- readr::read_csv(
  covariance_overview_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

covariance_pair_values <- readr::read_csv(
  covariance_pair_values_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

covariance_components <- readr::read_csv(
  covariance_components_file,
  show_col_types = FALSE,
  na = c("", "NA")
)

# Keep only the approved-ledger scope for reader-facing covariance figures.
covariance_approved <- covariance_components %>%
  dplyr::filter(
    Subset == "approved_ledger"
  )

stopifnot(
  !anyDuplicated(covariance_overview$Pair),
  !anyDuplicated(
    covariance_pair_values[
      c(
        "Pair",
        "Context_Key"
      )
    ]
  )
)

# Check that saved context values agree with the support table.
for (i in seq_len(nrow(covariance_overview))) {

  pair_i <- covariance_overview$Pair[i]

  x <- covariance_pair_values[
    covariance_pair_values$Pair == pair_i,
    ,
    drop = FALSE
  ]

  if (
    nrow(x) != covariance_overview$n_contexts[i] ||
      length(
        unique(
          x$Reference_ID
        )
      ) != covariance_overview$n_references[i]
  ) {
    stop(
      "Covariance support tables disagree for pair: ",
      pair_i
    )
  }
}


# ----------------------------
# Base-graphics save helper
# ----------------------------

save_covariance_figure <- function(
  stem,
  width_px,
  height_px,
  width_in,
  height_in,
  draw
) {

  png(
    file.path(
      figure_dir,
      paste0(
        stem,
        ".png"
      )
    ),
    width = width_px,
    height = height_px,
    res = 300,
    pointsize = 12,
    family = "sans"
  )

  draw()

  dev.off()

  pdf(
    file.path(
      figure_dir,
      paste0(
        stem,
        ".pdf"
      )
    ),
    width = width_in,
    height = height_in,
    family = "sans",
    useDingbats = FALSE
  )

  draw()

  dev.off()
}


# ----------------------------
# Publication colours
# ----------------------------

covariance_references <- sort(
  unique(
    covariance_pair_values$Reference_ID
  )
)

covariance_colours <- setNames(
  grDevices::hcl.colors(
    length(
      covariance_references
    ),
    "Dark 3"
  ),
  covariance_references
)


draw_covariance_pair <- function(
  pair,
  tag
) {

  x <- covariance_pair_values[
    covariance_pair_values$Pair == pair,
    ,
    drop = FALSE
  ]

  z <- covariance_overview[
    covariance_overview$Pair == pair,
    ,
    drop = FALSE
  ]

  metals <- strsplit(
    pair,
    "-",
    fixed = TRUE
  )[[1]]

  plot(
    log(
      x$A_value
    ),
    log(
      x$B_value
    ),
    pch = 21,
    bg = covariance_colours[
      x$Reference_ID
    ],
    col = "grey25",
    cex = 1.4,
    xlab = paste(
      metals[1],
      "log-LC50"
    ),
    ylab = paste(
      metals[2],
      "log-LC50"
    ),
    main = paste0(
      tag,
      "  ",
      paste(
        metals,
        collapse = " - "
      )
    ),
    sub = sprintf(
      "%d matched comparisons | %d publications",
      z$n_contexts,
      z$n_references
    )
  )

  limited <- z$n_references < 3 ||
    z$n_contexts <= 2

  label <- if (limited) {

    z$Interpretation

  } else if (
    is.finite(
      z$r_equal_reference
    )
  ) {

    sprintf(
      "Equal-publication r = %.3f",
      z$r_equal_reference
    )

  } else {

    "No correlation available"
  }

  mtext(
    label,
    side = 3,
    line = 0.3,
    cex = 0.9,
    col = if (
      limited
    ) {
      "#8C3B17"
    } else {
      "grey20"
    }
  )
}


# ----------------------------
# Figure 4.5: focal Cu-Cd-Zn associations
# ----------------------------

focal_covariance_pairs <- c(
  "Cu-Cd",
  "Cu-Zn",
  "Cd-Zn"
)

covariance_focal_support <- covariance_overview %>%
  dplyr::filter(
    Pair %in% focal_covariance_pairs
  ) %>%
  dplyr::mutate(
    Pair = factor(
      Pair,
      levels = focal_covariance_pairs
    )
  ) %>%
  dplyr::arrange(
    Pair
  )

covariance_focal_points <- covariance_pair_values %>%
  dplyr::filter(
    Pair %in% focal_covariance_pairs
  )

stopifnot(
  nrow(covariance_focal_support) == 3,
  identical(
    as.integer(
      covariance_focal_support$n_contexts
    ),
    c(
      21L,
      6L,
      7L
    )
  ),
  identical(
    as.integer(
      covariance_focal_support$n_references
    ),
    c(
      10L,
      4L,
      4L
    )
  )
)

write_csv(
  covariance_focal_points,
  file.path(
    plot_data_dir,
    "Figure_4_05_points.csv"
  )
)

write_csv(
  covariance_focal_support,
  file.path(
    plot_data_dir,
    "Figure_4_05_support.csv"
  )
)

draw_figure_4_05 <- function() {

  par(
    mfrow = c(
      2,
      2
    ),
    mar = c(
      4.8,
      4.8,
      3.7,
      1.2
    ),
    oma = c(
      3.4,
      0,
      2.8,
      0
    ),
    cex = 1,
    las = 1
  )

  for (
    j in seq_along(
      focal_covariance_pairs
    )
  ) {

    draw_covariance_pair(
      focal_covariance_pairs[j],
      LETTERS[j]
    )
  }

  plot.new()

  shown <- covariance_references[
    covariance_references %in%
      covariance_pair_values$Reference_ID[
        covariance_pair_values$Pair %in%
          focal_covariance_pairs
      ]
  ]

  legend(
    "center",
    legend = shown,
    pt.bg = covariance_colours[
      shown
    ],
    pch = 21,
    pt.cex = 1.4,
    cex = 0.95,
    bty = "n",
    title = "Publications",
    ncol = if (
      length(
        shown
      ) > 8
    ) {
      2
    } else {
      1
    }
  )

  mtext(
    "Core metal comparisons",
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.2
  )

  mtext(
    "Natural-log LC50; concentrations expressed in micromol/L before transformation.",
    side = 1,
    outer = TRUE,
    line = 0.7,
    cex = 0.85
  )

  mtext(
    "Each point is a matched comparison within one publication. Colours identify publications.",
    side = 1,
    outer = TRUE,
    line = 1.8,
    cex = 0.85
  )
}

save_covariance_figure(
  "Figure_4_05_matched_metal_associations",
  width_px = 3600,
  height_px = 3000,
  width_in = 12,
  height_in = 10,
  draw = draw_figure_4_05
)


# ----------------------------
# Figure 4.6: covariance components
# ----------------------------
#
# Preserve the thesis display scope from the final covariance script:
# all approved-ledger pairs with at least two matched contexts.
# This is deliberately broader than only the three focal Cu-Cd-Zn
# pairs and is not a significance-based selection.

covariance_component_source <- covariance_approved %>%
  dplyr::filter(
    n_contexts >= 2
  )

covariance_component_plot_data <- dplyr::bind_rows(
  covariance_component_source %>%
    dplyr::transmute(
      Pair,
      n_contexts,
      n_references,
      Panel = "A  Total and between publications",
      Component = "Total",
      Covariance = total_log_covariance
    ),
  covariance_component_source %>%
    dplyr::transmute(
      Pair,
      n_contexts,
      n_references,
      Panel = "A  Total and between publications",
      Component = "Between publications",
      Covariance = between_reference_covariance
    ),
  covariance_component_source %>%
    dplyr::transmute(
      Pair,
      n_contexts,
      n_references,
      Panel = "B  Within publications",
      Component = "Within publications",
      Covariance = within_reference_covariance
    )
)

write_csv(
  covariance_component_plot_data,
  file.path(
    plot_data_dir,
    "Figure_4_06_covariance_components.csv"
  )
)

draw_figure_4_06 <- function() {

  z <- covariance_component_source

  par(
    mfrow = c(
      1,
      2
    ),
    oma = c(
      3.5,
      0,
      2,
      0
    ),
    las = 1
  )

  y <- rev(
    seq_len(
      nrow(
        z
      )
    )
  )

  labels <- paste0(
    z$Pair,
    ifelse(
      z$n_references < 3,
      " *",
      ""
    ),
    "\n",
    z$n_contexts,
    " comparisons / ",
    z$n_references,
    " publications"
  )

  label_lines <- unlist(
    strsplit(
      labels,
      "\n",
      fixed = TRUE
    )
  )

  label_margin <- max(
    strwidth(
      label_lines,
      units = "inches",
      cex = 0.85
    )
  ) + 0.35

  bounds <- function(x) {

    r <- range(
      c(
        0,
        x
      ),
      finite = TRUE
    )

    pad <- if (
      diff(
        r
      ) > 0
    ) {
      diff(
        r
      ) * 0.12
    } else {
      0.05
    }

    r + c(
      -pad,
      pad
    )
  }

  par(
    mai = c(
      0.85,
      label_margin,
      0.75,
      0.25
    )
  )

  plot(
    z$total_log_covariance,
    y,
    type = "n",
    xlim = bounds(
      c(
        z$total_log_covariance,
        z$between_reference_covariance
      )
    ),
    ylim = c(
      0.5,
      nrow(
        z
      ) + 1.6
    ),
    yaxt = "n",
    ylab = "",
    xlab = "Weighted log-LC50 covariance",
    main = "A  Total and between publications"
  )

  axis(
    2,
    at = y,
    labels = labels,
    tick = FALSE,
    cex.axis = 0.85
  )

  abline(
    v = 0,
    col = "grey70"
  )

  abline(
    h = y,
    col = "grey93"
  )

  points(
    z$total_log_covariance,
    y + 0.12,
    pch = 16,
    col = "#234E70",
    cex = 1.2
  )

  points(
    z$between_reference_covariance,
    y - 0.12,
    pch = 15,
    col = "#269DA3",
    cex = 1.1
  )

  legend(
    "top",
    legend = c(
      "Total",
      "Between publications"
    ),
    pch = c(
      16,
      15
    ),
    col = c(
      "#234E70",
      "#269DA3"
    ),
    bty = "n",
    horiz = TRUE,
    cex = 0.9
  )

  plot(
    z$within_reference_covariance,
    y,
    pch = 16,
    col = "#CE881D",
    cex = 1.2,
    xlim = bounds(
      z$within_reference_covariance
    ),
    ylim = c(
      0.5,
      nrow(
        z
      ) + 1.6
    ),
    yaxt = "n",
    ylab = "",
    xlab = "Weighted log-LC50 covariance (different scale)",
    main = "B  Within publications"
  )

  axis(
    2,
    at = y,
    labels = labels,
    tick = FALSE,
    cex.axis = 0.85
  )

  abline(
    v = 0,
    col = "grey70"
  )

  mtext(
    "Panel B uses its own horizontal scale.",
    side = 3,
    line = 0.3,
    cex = 0.85
  )

  mtext(
    "* Fewer than three publications: limited evidence.",
    side = 1,
    outer = TRUE,
    line = 0.6,
    cex = 0.9
  )

  mtext(
    "Covariance components describe association; they are not causal shares or explained variance.",
    side = 1,
    outer = TRUE,
    line = 1.8,
    cex = 0.85
  )
}

save_covariance_figure(
  "Figure_4_06_covariance_components",
  width_px = 4500,
  height_px = 2700,
  width_in = 15,
  height_in = 9,
  draw = draw_figure_4_06
)


# ----------------------------
# Figure D.7: additional supported associations
# ----------------------------

additional_covariance_pairs <- c(
  "Cu-Hg",
  "Cd-Hg"
)

covariance_additional_support <- covariance_overview %>%
  dplyr::filter(
    Pair %in% additional_covariance_pairs
  ) %>%
  dplyr::mutate(
    Pair = factor(
      Pair,
      levels = additional_covariance_pairs
    )
  ) %>%
  dplyr::arrange(
    Pair
  )

covariance_additional_points <- covariance_pair_values %>%
  dplyr::filter(
    Pair %in% additional_covariance_pairs
  )

write_csv(
  covariance_additional_points,
  file.path(
    plot_data_dir,
    "Figure_D_07_points.csv"
  )
)

write_csv(
  covariance_additional_support,
  file.path(
    plot_data_dir,
    "Figure_D_07_support.csv"
  )
)

draw_figure_D_07 <- function() {

  layout(
    matrix(
      c(
        1,
        2,
        3,
        3
      ),
      nrow = 2,
      byrow = TRUE
    ),
    heights = c(
      4,
      1
    )
  )

  par(
    mar = c(
      4.8,
      4.8,
      3.7,
      1.2
    ),
    oma = c(
      3.4,
      0,
      2.8,
      0
    ),
    cex = 1,
    las = 1
  )

  for (
    j in seq_along(
      additional_covariance_pairs
    )
  ) {

    draw_covariance_pair(
      additional_covariance_pairs[j],
      LETTERS[j]
    )
  }

  par(
    mar = c(
      0,
      0,
      0,
      0
    )
  )

  plot.new()

  shown <- covariance_references[
    covariance_references %in%
      covariance_pair_values$Reference_ID[
        covariance_pair_values$Pair %in%
          additional_covariance_pairs
      ]
  ]

  legend(
    "center",
    legend = shown,
    pt.bg = covariance_colours[
      shown
    ],
    pch = 21,
    pt.cex = 1.4,
    cex = 0.95,
    bty = "n",
    title = "Publications",
    ncol = min(
      4L,
      length(
        shown
      )
    )
  )

  mtext(
    "Additional metal comparisons",
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.2
  )

  mtext(
    "Natural-log LC50; concentrations expressed in micromol/L before transformation.",
    side = 1,
    outer = TRUE,
    line = 0.7,
    cex = 0.85
  )

  mtext(
    "Each point is a matched comparison within one publication. Colours identify publications.",
    side = 1,
    outer = TRUE,
    line = 1.8,
    cex = 0.85
  )
}

save_covariance_figure(
  "Figure_D_07_additional_metal_associations",
  width_px = 3600,
  height_px = 2100,
  width_in = 12,
  height_in = 7,
  draw = draw_figure_D_07
)



# ============================================================
# 21. SAVE FIGURE OBJECTS + FIGURE INDEX
# ============================================================

saveRDS(
  list(
    Figure_4_01 = F00,
    Figure_4_02 = F01,
    Figure_4_03 = F02,
    Figure_4_04 = F04,
    Figure_4_07 = F05,
    Figure_4_08 = F06A,
    Figure_4_09 = F09,
    Figure_4_10 = F07,
    Figure_4_11 = F08A,
    Figure_4_12 = F12,
    Figure_4_13 = F10,
    Figure_4_14 = F11,
    Figure_4_15 = F03,
    Figure_4_16 = F06,
    Figure_4_17 = F08,
    Figure_C_01 = APP_A1
  ),
  file.path(
    object_dir,
    "main_thesis_figure_objects.rds"
  )
)

figure_index <- tibble::tribble(
  ~Thesis_figure, ~Filename_stem, ~Primary_analysis_source,
  "Figure 4.1",  "Figure_4_01_evidence_landscape", "03_audit_metal_support + 04_fit_pooled_metal_model",
  "Figure 4.2",  "Figure_4_02_pooled_results_and_model_estimates", "04_fit_pooled_metal_model",
  "Figure 4.3",  "Figure_4_03_pooled_pairwise_contrasts", "04_fit_pooled_metal_model",
  "Figure 4.4",  "Figure_4_04_within_reference_metal_comparisons", "05_compare_metals_within_references",
  "Figure 4.5",  "Figure_4_05_matched_metal_associations", "14_analyze_metal_covariance",
  "Figure 4.6",  "Figure_4_06_covariance_components", "14_analyze_metal_covariance",
  "Figure 4.7",  "Figure_4_07_stage_results_and_model_estimates", "11_fit_stage_model",
  "Figure 4.8",  "Figure_4_08_stage_contrasts_and_interaction", "11_fit_stage_model",
  "Figure 4.9",  "Figure_4_09_within_reference_stage_comparisons", "09_compare_stages_within_references",
  "Figure 4.10", "Figure_4_10_order_results_and_model_estimates", "13_fit_order_model",
  "Figure 4.11", "Figure_4_11_order_contrasts_and_interaction", "13_fit_order_model",
  "Figure 4.12", "Figure_4_12_exposure_duration_comparisons", "06_compare_exposure_duration",
  "Figure 4.13", "Figure_4_13_salinity_comparisons", "07_assess_salinity_temperature",
  "Figure 4.14", "Figure_4_14_temperature_comparisons", "07_assess_salinity_temperature",
  "Figure 4.15", "Figure_4_15_pooled_model_robustness", "04_fit_pooled_metal_model",
  "Figure 4.16", "Figure_4_16_stage_model_robustness", "11_fit_stage_model",
  "Figure 4.17", "Figure_4_17_order_model_robustness", "13_fit_order_model",
  "Figure C.1",  "Figure_C_01_pH_salinity_confounding", "08_assess_ph_evidence",
  "Figure D.7",  "Figure_D_07_additional_metal_associations", "14_analyze_metal_covariance"
)

write_csv(
  figure_index,
  file.path(output_root, "figure_index.csv")
)

# ============================================================
# 22. FINAL QA + PROVENANCE
# ============================================================

final_stems <- figure_index$Filename_stem

final_pngs <- file.path(
  figure_dir,
  paste0(final_stems, ".png")
)

final_pdfs <- file.path(
  figure_dir,
  paste0(final_stems, ".pdf")
)

stopifnot(
  all(file.exists(final_pngs)),
  all(file.exists(final_pdfs))
)

input_manifest <- tibble::tibble(
  Role = c(
    "Pooled model object",
    "Stage model object",
    "Order model object",
    "Within-Reference metal table directory",
    "Environmental table directory",
    "pH table directory",
    "Covariance overview",
    "Covariance pair values",
    "Covariance components"
  ),
  Path = c(
    pooled_object_file,
    stage_object_file,
    order_object_file,
    within_ref_dir,
    environmental_table_dir,
    ph_table_dir,
    covariance_overview_file,
    covariance_pair_values_file,
    covariance_components_file
  )
)

write_csv(
  input_manifest,
  file.path(output_root, "input_output_manifest.csv")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_root, "sessionInfo.txt")
)

writeLines(
  c(
    "STATUS: MAIN THESIS FIGURE SUITE = PASS",
    paste0("Reader-facing figures: ", nrow(figure_index)),
    "Scientific models and source-verification decisions were not refitted or changed by this script."
  ),
  file.path(output_root, "RUN_COMPLETE.txt")
)

cat("\n\n")
cat("============================================================\n")
cat("MAIN THESIS FIGURE SUITE - COMPLETE\n")
cat("============================================================\n")
cat("\nFigures written to:\n")
cat(
  normalizePath(
    figure_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)
cat("\nSTATUS: MAIN THESIS FIGURE SUITE = PASS\n")

# ============================================================
# END
# ============================================================
