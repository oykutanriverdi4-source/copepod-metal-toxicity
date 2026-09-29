# ============================================================
# 14_analyze_metal_covariance.R
# ============================================================
# Purpose
# Quantify co-variation in matched metal LC50 values across the
# source-verified within-Reference comparison contexts and decompose
# weighted log-LC50 covariance into between-Reference and
# within-Reference components.
#
# Scientific question
# When two metals are measured under comparable matched contexts,
# do their LC50 values tend to vary together across the represented
# evidence?
#
# Analytical unit
# Reference_ID + Source_Context_ID, with a geometric mean per metal
# when more than one retained LC50 Result represents the same metal
# within a matched context.
#
# Weighting
# Each Reference receives equal total weight; that weight is divided
# equally among the matched contexts contributed by that Reference.
# These are descriptive balancing weights, not precision weights.
#
# Uncertainty / sensitivity
# - 20,000 Reference-level bootstrap resamples for supported pairs.
# - LORO (leave-one-Reference-out) influence analysis.
# - Four overlapping evidence scopes retained from the final thesis
#   workflow: approved ledger, strict verification, 96 h, and
#   exclusion of order-only labels.
#
# Important interpretation
# This analysis describes association, not relative toxicity,
# mixture toxicity, shared mechanism or causation. No hypothesis-test
# p-values are used for the covariance analysis.
#
# Inputs
# data/processed/combined_ecotox_wos_harmonized.csv
# data/curated/source_verification/combined_source_verified_rows_used.csv
# data/curated/source_verification/combined_source_adjudication_rows.csv
#
# Main outputs
# outputs/14_analyze_metal_covariance/
#
# Figure note
# Figures produced here are analysis/QA outputs. The thesis-facing
# covariance figures will be regenerated later in the unified
# polished figure-suite script.
#
# Reproducibility rule
# Run from the repository root with copepod-metal-toxicity.Rproj open.
# The scientific calculations, context construction, weighting,
# bootstrap seed/RNG configuration, scope definitions and LORO logic
# are preserved from the final covariance workflow.
# ============================================================


# ============================================================
# 0. SETTINGS
# ============================================================
PROJECT_DIR <- ""
BOOTSTRAP_DRAWS <- 20000L
RANDOM_SEED <- 20260908L
# Frozen for exact reproducibility of the final thesis covariance workflow.
RNG_KIND <- "L'Ecuyer-CMRG"
RNG_NORMAL_KIND <- "Inversion"
RNG_SAMPLE_KIND <- "Rejection"


# ============================================================
# WORD TABLES: FORMAT THE RESULTS COMPUTED BELOW
# ============================================================
write_covariance_word <- function(summary, loro, folder) {
  required <- c("officer", "flextable")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    warning(
      "Optional Word tables skipped because reporting package(s) are not installed: ",
      paste(missing, collapse = ", ")
    )
    return(invisible(NA_character_))
  }
  a <- summary[summary$Subset == "approved_ledger", , drop = FALSE]
  b <- a[a$n_contexts >= 2, , drop = FALSE]
  if (anyDuplicated(summary[c("Subset", "Pair")]) ||
      anyDuplicated(loro[c("Pair", "Omitted_Reference")])) stop("Duplicate result rows.")
  fmt <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "-")
  interpretation <- function(x) ifelse(x$n_contexts == 0, "No verified match",
                                       ifelse(x$n_contexts == 1, "One context: correlation unavailable",
                                              ifelse(x$n_contexts == 2, "Two contexts: correlation uninformative",
                                                     ifelse(x$n_references == 2, "Only two references", "Exploratory association"))))
  make <- function(...) data.frame(..., check.names = FALSE, stringsAsFactors = FALSE)
  main <- make(b$Pair, b$n_contexts, b$n_references, fmt(b$r_equal_reference),
               ifelse(is.finite(b$bootstrap_low), paste(fmt(b$bootstrap_low), "to", fmt(b$bootstrap_high)),
                      "Not estimated"), ifelse(is.finite(b$bootstrap_low), as.character(b$bootstrap_valid), "-"),
               interpretation(b))
  components <- make(b$Pair, b$n_contexts, b$n_references, fmt(b$total_log_covariance),
                     fmt(b$between_reference_covariance), fmt(b$within_reference_covariance), b$n_varying_references)
  scopes <- c("approved_ledger", "strict_verification", "duration_96h", "exclude_order_only_labels")
  scope_table <- make(Pair = a$Pair)
  for (scope in scopes) {
    x <- summary[summary$Subset == scope, , drop = FALSE]
    if (!setequal(x$Pair, a$Pair)) stop("Scope pair coverage differs: ", scope)
    x <- x[match(a$Pair, x$Pair), , drop = FALSE]
    scope_table[[scope]] <- paste0(x$n_contexts, "/", x$n_references, "; ",
                                   ifelse(x$n_contexts > 2, fmt(x$r_equal_reference),
                                          ifelse(x$n_contexts == 2, "not informative", "unavailable")))
  }
  weighting <- make(b$Pair, b$n_contexts, b$n_references, fmt(b$r_equal_reference),
                    fmt(b$r_equal_context), fmt(b$spearman_context), fmt(b$r_within_centered))
  omission <- make(loro$Pair, loro$Omitted_Reference, loro$n_contexts, loro$n_references,
                   fmt(loro$r_equal_reference), fmt(loro$total_log_covariance),
                   fmt(loro$between_reference_covariance), fmt(loro$within_reference_covariance))
  support <- make(a$Pair, a$n_contexts, a$n_references, interpretation(a))
  titles <- c(    "Table C1. Metal LC50 associations in verified contexts.",
                  "Table C2. Decomposition of weighted log-LC50 covariance.",
                  "Table C3. Scope sensitivity of metal associations.",
                  "Table C4. Weighting and rank-based descriptive comparisons.",
                  "Table C5. Leave-one-Reference-out results.",
                  "Table C5 (continued). Leave-one-Reference-out results.",
                  "Table C6. Verified pair coverage and reporting limits.")
  notes <- c(    "Bounds are the 2.5th and 97.5th percentiles of valid equal-reference correlation estimates from 20,000 source-level bootstrap draws (seed 20260908). They are not intervals for covariance components. Very small source counts can yield unstable or degenerate distributions. The computational minimum of three references does not establish adequate support. Cd-Pb has one reference; its two-context r is not informative. C numbering is provisional.",
                 "Covariance uses natural-log molar LC50 values and normalized weights: each reference has equal total weight, distributed equally among its contexts. Total = between-reference + within-reference covariance (rounding may obscure equality). A single-context reference contributes no within-reference variation. A zero between-reference component for Cd-Pb follows from its single-reference design. Components do not identify causes or explained variance.",
                 "Each cell gives contexts/references; equal-reference r. Two-context correlations are suppressed as uninformative. One- or zero-context correlations are unavailable. These are overlapping sensitivity subsets, not independent replications. The order-label scope also excludes missing Species/Order labels in the implemented script. Bootstrap and source omission apply only to the approved-ledger scope.",
                 "All values are descriptive outputs of the covariance analysis. Context weighting gives each context equal weight; reference weighting gives each source equal total weight. Spearman values use context ranks and are not a source-adjusted nonparametric hypothesis test. Two-context correlations are mathematically degenerate; see C1. Within-centred r can depend on very few varying references (C2) and does not adjust for all biological or experimental differences.",
                 "Each row is a source-omission result for the approved-ledger scope. The full table is split over two pages. Ranges across these rows describe sensitivity, not confidence intervals. Correlations based on only two remaining contexts are degenerate; one remaining reference cannot support between-reference inference. A dash denotes an undefined correlation, not a null association.",
                 "Each row is a source-omission result for the approved-ledger scope. The full table is split over two pages. Ranges across these rows describe sensitivity, not confidence intervals. Correlations based on only two remaining contexts are degenerate; one remaining reference cannot support between-reference inference. A dash denotes an undefined correlation, not a null association.",
                 "All 21 candidate pairs and their support are retained. No verified match does not imply absence of biological association or exclusion of the source from the thesis dataset. Single-context computational covariance zeros are not reported as evidence of no covariance. Source decisions were made before this output package; this document does not change the verification ledger.")
  headers <- list(    c("Pair", "Contexts", "References", "Equal-reference r", "Bootstrap bounds", "Valid draws", "Interpretation"),
                      c("Pair", "Contexts", "References", "Total", "Between references", "Within references", "References with within-reference variation"),
                      c("Pair", "Approved ledger", "Strict verification", "96 h", "Exclude order-only labels"),
                      c("Pair", "Contexts", "References", "Equal-reference Pearson r", "Equal-context Pearson r", "Context Spearman rho", "Within-centred r"),
                      c("Pair", "Omitted reference", "Remaining contexts", "Remaining references", "Equal-reference r", "Total covariance", "Between references", "Within references"),
                      c("Pair", "Omitted reference", "Remaining contexts", "Remaining references", "Equal-reference r", "Total covariance", "Between references", "Within references"),
                      c("Pair", "Contexts", "References", "Reporting status"))
  # Source-omission rows are paginated; no row is omitted if a future run grows.
  chunks <- split(seq_len(nrow(omission)), ceiling(seq_len(nrow(omission)) / 18))
  tables <- list(main, components, scope_table, weighting)
  title_ids <- 1:4
  for (k in seq_along(chunks)) {
    tables[[length(tables) + 1L]] <- omission[chunks[[k]], , drop = FALSE]
    title_ids <- c(title_ids, if (k == 1L) 5L else 6L)
  }
  tables[[length(tables) + 1L]] <- support; title_ids <- c(title_ids, 7L)
  doc <- officer::read_docx()
  doc <- officer::body_set_default_section(doc, officer::prop_section(
    page_size = officer::page_size(width = 11.69, height = 8.27, orient = "landscape"),
    page_margins = officer::page_mar(top = .5, bottom = .5, left = .45, right = .45)))
  paragraph <- function(doc, text, bold = FALSE, italic = FALSE) {
    officer::body_add_fpar(doc, officer::fpar(officer::ftext(text,
                                                             officer::fp_text(font.family = "Times New Roman", font.size = 9,
                                                                              bold = bold, italic = italic))))
  }
  for (j in seq_along(tables)) {
    id <- title_ids[j]; x <- tables[[j]]; names(x) <- headers[[id]]
    if (j > 1L) doc <- officer::body_add_break(doc)
    doc <- paragraph(doc, titles[id], bold = TRUE)
    ft <- flextable::flextable(x)
    ft <- flextable::theme_box(ft)
    ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::bg(ft, bg = "#D9D9D9", part = "header")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::align(ft, align = "left", part = "all")
    ft <- flextable::padding(ft, padding = 2, part = "all")
    ft <- flextable::width(ft, width = 10.79 / ncol(x))
    ft <- flextable::set_table_properties(ft, layout = "fixed",
                                          opts_word = list(split = FALSE, repeat_headers = TRUE))
    doc <- flextable::body_add_flextable(doc, ft)
    doc <- paragraph(doc, paste("Note.", notes[id]))
  }
  word_folder <- file.path(dirname(folder), "02_tables")
  if (!dir.exists(word_folder) && !dir.create(word_folder, recursive = TRUE))
    stop("Cannot create Word output folder.")
  target <- file.path(word_folder, "metal_covariance_all_tables.docx")
  print(doc, target = target)
  message("Word tables: ", target)
  invisible(target)
}

# ============================================================
# FIGURE FUNCTIONS: PLOT THE RESULTS COMPUTED BELOW
# ============================================================
render_covariance_figures <- function(summary, pair_values, run_dir) {
  need_summary <- c("Subset", "Pair", "n_contexts", "n_references",
                    "r_equal_reference", "total_log_covariance", "between_reference_covariance",
                    "within_reference_covariance")
  need_pairs <- c("Pair", "Context_Key", "Reference_ID", "A_value", "B_value")
  stopifnot(all(need_summary %in% names(summary)),
            all(need_pairs %in% names(pair_values)))
  s <- summary[summary$Subset == "approved_ledger", , drop = FALSE]
  if (anyDuplicated(s$Pair) || anyDuplicated(pair_values[c("Pair", "Context_Key")]))
    stop("Duplicate pair or context identifiers in saved results.")
  if (!all(pair_values$Pair %in% s$Pair) ||
      anyNA(pair_values$Reference_ID) ||
      any(!is.finite(pair_values$A_value) | pair_values$A_value <= 0) ||
      any(!is.finite(pair_values$B_value) | pair_values$B_value <= 0))
    stop("Invalid saved pair values.")
  for (i in seq_len(nrow(s))) {
    x <- pair_values[pair_values$Pair == s$Pair[i], , drop = FALSE]
    if (nrow(x) != s$n_contexts[i] ||
        length(unique(x$Reference_ID)) != s$n_references[i])
      stop("Saved tables disagree for pair: ", s$Pair[i])
  }
  
  # Groups reflect the thesis metal panel and data support, never p-values.
  core <- c("Cu-Cd", "Cu-Zn", "Cd-Zn")
  eligible <- s$Pair[s$n_contexts > 2 & s$n_references >= 3]
  groups <- list(
    core_metals = intersect(core, eligible),
    additional_metals = setdiff(eligible, core),
    limited_support = s$Pair[s$n_contexts >= 2 & !s$Pair %in% eligible])
  titles <- c(core_metals = "Core metal comparisons",
              additional_metals = "Additional metal comparisons",
              limited_support = "Comparisons supported by few publications")
  s$Figure_group <- ifelse(s$n_contexts == 0, "No matched pair",
                           ifelse(s$n_contexts == 1, "Single context: table only",
                                  ifelse(s$Pair %in% groups$core_metals, "Core metals",
                                         ifelse(s$Pair %in% groups$additional_metals, "Additional metals", "Limited support"))))
  s$Display_note <- ifelse(s$n_contexts == 0, "No correlation available",
                           ifelse(s$n_contexts == 1, "One context: no correlation",
                                  ifelse(s$n_references == 1, "One publication: descriptive only",
                                         ifelse(s$n_contexts == 2, "Two contexts: correlation not informative",
                                                ifelse(s$n_references == 2, "Only two publications", "Exploratory association")))))
  
  folder <- file.path(run_dir, "03_figures")
  if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE))
    stop("Cannot create figure folder.")
  write.csv(s, file.path(folder, "metal_covariance_pair_support.csv"), row.names = FALSE, na = "")
  refs <- sort(unique(pair_values$Reference_ID))
  colors <- setNames(grDevices::hcl.colors(length(refs), "Dark 3"), refs)
  save_png <- function(name, width, height, draw) {
    png(file.path(folder, name), width = width, height = height, res = 300, pointsize = 12, family = "sans")
    on.exit(dev.off())
    draw()
  }
  draw_pair <- function(pair, tag) {
    x <- pair_values[pair_values$Pair == pair, , drop = FALSE]
    z <- s[s$Pair == pair, , drop = FALSE]
    metals <- strsplit(pair, "-", fixed = TRUE)[[1]]
    plot(log(x$A_value), log(x$B_value), pch = 21,
         bg = colors[x$Reference_ID], col = "grey25", cex = 1.4,
         xlab = paste(metals[1], "log-LC50"), ylab = paste(metals[2], "log-LC50"),
         main = paste0(tag, "  ", paste(metals, collapse = " - ")),
         sub = sprintf("%d matched comparisons | %d publications", z$n_contexts, z$n_references))
    limited <- z$n_references < 3 || z$n_contexts <= 2
    label <- if (limited) z$Display_note else if (is.finite(z$r_equal_reference))
      sprintf("Equal-publication r = %.3f", z$r_equal_reference) else "No correlation available"
    mtext(label, side = 3, line = .3, cex = .9,
          col = if (limited) "#8C3B17" else "grey20")
  }
  files <- character()
  for (g in names(groups)) {
    pairs <- groups[[g]]
    if (!length(pairs)) next
    # At most three plots plus a shared reference legend per page.
    pages <- split(pairs, ceiling(seq_along(pairs) / 3))
    for (page in seq_along(pages)) {
      selected <- pages[[page]]
      figure_names <- c(
        core_metals = "metal_covariance_01_Cu_Cd_Zn_associations",
        additional_metals = "metal_covariance_02_other_metals_associations",
        limited_support = "metal_covariance_03_limited_support_associations")
      filename <- paste0(figure_names[[g]],
                         if (length(pages) > 1L) sprintf("_page_%02d", page) else "", ".png")
      save_png(filename, 3600, if (length(selected) == 2L) 2100 else 3000, function() {
        if (length(selected) == 2L) {
          layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4, 1))
          par(mar = c(4.8, 4.8, 3.7, 1.2),
              oma = c(3.4, 0, 2.8, 0), cex = 1, las = 1)
        } else {
          par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3.7, 1.2),
              oma = c(3.4, 0, 2.8, 0), cex = 1, las = 1)
        }
        for (j in seq_along(selected)) draw_pair(selected[j], LETTERS[j])
        if (length(selected) < 3 && length(selected) != 2L)
          for (j in seq_len(3 - length(selected))) plot.new()
        if (length(selected) == 2L) par(mar = c(0, 0, 0, 0))
        plot.new()
        shown <- refs[refs %in% pair_values$Reference_ID[pair_values$Pair %in% selected]]
        legend("center", legend = shown, pt.bg = colors[shown],
               pch = 21, pt.cex = 1.4, cex = .95, bty = "n", title = "Publications",
               ncol = if (length(selected) == 2L) min(4L, length(shown)) else if (length(shown) > 8) 2 else 1)
        mtext(titles[g], side = 3, outer = TRUE, line = 1, font = 2, cex = 1.2)
        mtext("Natural-log LC50; concentrations expressed in micromol/L before transformation.",
              side = 1, outer = TRUE, line = .7, cex = .85)
        mtext("Each point is a matched comparison within one publication. Colours identify publications.",
              side = 1, outer = TRUE, line = 1.8, cex = .85)
      })
      files <- c(files, filename)
    }
  }
  
  # Same stored covariance values, displayed on two explicitly different scales.
  z <- s[s$n_contexts >= 2, , drop = FALSE]
  if (nrow(z)) {
    filename <- "metal_covariance_04_between_within_reference_components.png"
    save_png(filename, 4500, 2700, function() {
      par(mfrow = c(1, 2), oma = c(3.5, 0, 2, 0), las = 1)
      y <- rev(seq_len(nrow(z)))
      labels <- paste0(z$Pair, ifelse(z$n_references < 3, " *", ""),
                       "\n", z$n_contexts, " comparisons / ", z$n_references, " publications")
      label_lines <- unlist(strsplit(labels, "\n", fixed = TRUE))
      label_margin <- max(strwidth(label_lines, units = "inches", cex = .85)) + .35
      par(mai = c(.85, label_margin, .75, .25))
      bounds <- function(x) {
        r <- range(c(0, x), finite = TRUE)
        pad <- if (diff(r) > 0) diff(r) * .12 else .05
        r + c(-pad, pad)
      }
      plot(z$total_log_covariance, y, type = "n",
           xlim = bounds(c(z$total_log_covariance, z$between_reference_covariance)),
           ylim = c(.5, nrow(z) + 1.6), yaxt = "n", ylab = "",
           xlab = "Weighted log-LC50 covariance", main = "A  Total and between publications")
      axis(2, at = y, labels = labels, tick = FALSE, cex.axis = .85)
      abline(v = 0, col = "grey70")
      abline(h = y, col = "grey93")
      points(z$total_log_covariance, y + .12, pch = 16, col = "#234E70", cex = 1.2)
      points(z$between_reference_covariance, y - .12, pch = 15, col = "#269DA3", cex = 1.1)
      legend("top", legend = c("Total", "Between publications"), pch = c(16, 15),
             col = c("#234E70", "#269DA3"), bty = "n", horiz = TRUE, cex = .9)
      plot(z$within_reference_covariance, y, pch = 16, col = "#CE881D", cex = 1.2,
           xlim = bounds(z$within_reference_covariance), ylim = c(.5, nrow(z) + 1.6),
           yaxt = "n", ylab = "", xlab = "Weighted log-LC50 covariance (different scale)",
           main = "B  Within publications")
      axis(2, at = y, labels = labels, tick = FALSE, cex.axis = .85)
      abline(v = 0, col = "grey70")
      mtext("Panel B uses its own horizontal scale.", side = 3, line = .3, cex = .85)
      mtext("* Fewer than three publications: limited evidence.",
            side = 1, outer = TRUE, line = .6, cex = .9)
      mtext("Covariance components describe association; they are not causal shares or explained variance.",
            side = 1, outer = TRUE, line = 1.8, cex = .85)
    })
    files <- c(files, filename)
  }
  writeLines(c(
    "ANALYSIS / QA FIGURE NOTES; POLISHED THESIS FIGURES ARE GENERATED LATER",
    "Core pairs: Cu-Cd, Cu-Zn, Cd-Zn. Groups are defined without p-value selection.",
    "Fewer than three references is a display flag, not a scientific adequacy threshold.",
    "Zero/one-context pairs remain in metal_covariance_pair_support.csv; no scatter panel is drawn for them.",
    "Single-reference and two-context correlations are not highlighted as strong associations.",
    "Reference names/titles are in 01_audit/metal_covariance_reference_index.csv.",
    "The two covariance panels have different horizontal scales; compare numeric axes.",
    "A reference with one context provides no within-reference variation.",
    "Figures use the covariance components and correlations computed in this run.",
    "Colours denote references. No regression lines or confidence ellipses are added.",
    "Main-text versus supplementary placement remains undecided."
  ), file.path(folder, "metal_covariance_figure_notes.txt"))
  cat("\nFIGURES COMPLETE\n", folder, "\n", paste(files, collapse = "\n"), "\n", sep = "")
  invisible(folder)
}

local({
  message("14 | METAL COVARIANCE | EXPLICIT RNG | REPRODUCIBLE OUTPUTS")
  root <- normalizePath(if (nzchar(PROJECT_DIR)) PROJECT_DIR else getwd(),
                        winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "copepod-metal-toxicity.Rproj")))
    stop("Open copepod-metal-toxicity.Rproj or set PROJECT_DIR to the repository root.")
  B <- BOOTSTRAP_DRAWS; seed <- RANDOM_SEED
  stopifnot(length(B) == 1L, is.finite(B), B >= 1000L, B == as.integer(B),
            length(seed) == 1L, is.finite(seed))
  
  # ============================================================
  # 1. EXACT INPUT PATHS AND OUTPUT FOLDERS
  # ============================================================
  main_file <- file.path(root, "data", "processed",
                         "combined_ecotox_wos_harmonized.csv")
  verification_dir <- file.path(root, "data", "curated", "source_verification")
  verified_file <- file.path(verification_dir, "combined_source_verified_rows_used.csv")
  adjudication_file <- file.path(verification_dir, "combined_source_adjudication_rows.csv")
  inputs <- c(Combined = main_file, Verified = verified_file,
              Adjudication = adjudication_file)
  if (!all(file.exists(inputs)))
    stop("Missing current input: ", paste(inputs[!file.exists(inputs)], collapse = " | "),
         ". Locate the approved files; do not blindly regenerate the source ledgers.")
  hash <- function(x) unname(tools::md5sum(x))
  input_hashes <- hash(inputs)
  read_table <- function(path) read.csv(path, stringsAsFactors = FALSE,
                                        check.names = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
  combined <- read_table(main_file)
  v <- read_table(verified_file)
  adjudication <- read_table(adjudication_file)
  required <- c("Result_ID", "Reference_ID", "Species", "Order", "Metal",
                "Duration_days", "LC50_umol_L", "Source_Origin")
  require_columns <- function(x, fields, label) {
    absent <- setdiff(fields, names(x))
    if (length(absent)) stop(label, " lacks: ", paste(absent, collapse = ", "))
  }
  require_columns(combined, c(required, "Harmonization_Status"), "Combined")
  require_columns(v, c(required, "Source_Context_ID", "Final_Pairwise_Use",
                       "Verification_Status"), "Verified ledger")
  require_columns(adjudication, names(v), "Adjudication ledger")
  if (anyNA(combined$Result_ID) || anyDuplicated(combined$Result_ID) ||
      anyNA(adjudication$Result_ID) || anyDuplicated(adjudication$Result_ID))
    stop("Missing or duplicate Result_ID in an input.")
  keep <- !is.na(combined$Harmonization_Status) &
    combined$Harmonization_Status == "HARMONIZED" &
    is.finite(combined$LC50_umol_L) & combined$LC50_umol_L > 0
  d <- combined[keep, , drop = FALSE]
  if (nrow(d) != 304L || nrow(v) != 70L)
    stop("Expected current scope: 304 quantitative records and 70 approved records.")
  accepted <- adjudication[which(adjudication$Final_Pairwise_Use == "YES"), , drop = FALSE]
  if (!setequal(accepted$Result_ID, v$Result_ID)) stop("Ledgers disagree on approved IDs.")
  accepted <- accepted[match(v$Result_ID, accepted$Result_ID), names(v), drop = FALSE]
  for (f in names(v)) {
    a <- as.character(accepted[[f]]); b <- as.character(v[[f]])
    if (!identical(a, b)) stop("Ledgers disagree on field: ", f)
  }
  
  # ============================================================
  # 2. RECORD IDENTITIES, VALUES AND METADATA
  # ============================================================
  stopifnot(!anyNA(d$Result_ID), !anyDuplicated(d$Result_ID),
            !anyNA(v$Result_ID), !anyDuplicated(v$Result_ID),
            all(v$Result_ID %in% d$Result_ID), all(v$Final_Pairwise_Use == "YES"),
            all(is.finite(d$LC50_umol_L) & d$LC50_umol_L > 0),
            all(is.finite(v$LC50_umol_L) & v$LC50_umol_L > 0),
            !anyNA(v$Reference_ID), !anyNA(v$Source_Context_ID))
  m <- match(v$Result_ID, d$Result_ID)
  stopifnot(all(v$Reference_ID == d$Reference_ID[m]), all(v$Metal == d$Metal[m]),
            all(v$Species == d$Species[m]), all(v$Duration_days == d$Duration_days[m]),
            max(abs(log(v$LC50_umol_L) - log(d$LC50_umol_L[m]))) < 1e-10)
  
  # A changed metadata field requires source-decision review, not a silent join.
  metadata <- intersect(c("Order", "Lifestage", "Temperature", "Salinity", "pH",
                          "Exposure_Type", "Conc_Type", "Basis_Interpretation", "External_Stage_Sex"),
                        intersect(names(d), names(v)))
  for (field in metadata) {
    a <- as.character(v[[field]]); b <- as.character(d[[field]][m])
    a[!is.na(a) & a == ""] <- NA_character_; b[!is.na(b) & b == ""] <- NA_character_
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    if (!all(same)) stop("Combined data and source-verification ledger differ in ", field,
                         ". Review the current inputs and source decisions before proceeding.", call. = FALSE)
  }
  
  for (f in c("Result_ID", "Reference_ID", "Source_Context_ID"))
    if (any(!nzchar(trimws(as.character(v[[f]]))))) stop("Blank verified ID: ", f)
  out <- file.path(root, "outputs/14_analyze_metal_covariance")
  dirs <- setNames(file.path(out, c("01_audit", "02_tables", "03_figures", "04_objects")),
                   c("audit", "tables", "figures", "objects"))
  for (folder in dirs)
    if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE)) stop("Cannot create: ", folder)
  write_table <- function(x, name, folder = dirs[["tables"]])
    write.csv(x, file.path(folder, paste0("metal_covariance_", name, ".csv")), row.names = FALSE,
              na = "", fileEncoding = "UTF-8")
  
  # ============================================================
  # 3. VERIFIED CONTEXTS AND METAL VALUES
  # ============================================================
  collapse <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = " | ")
  key <- function(x, cols) do.call(paste, c(x[cols], sep = "::"))
  preferred <- c("Cu", "Cd", "Zn", "Hg", "Ni", "Ag", "Pb")
  metals <- c(intersect(preferred, unique(d$Metal)), setdiff(sort(unique(d$Metal)), preferred))
  pairs <- combn(metals, 2, simplify = FALSE)
  v$Context_Key <- key(v, c("Reference_ID", "Source_Context_ID"))
  strict_status <- c("VERIFIED_FULLTEXT", "VERIFIED_CORRECTED",
                     "VERIFIED_CORRECTED_SPLIT", "VERIFIED_FULLTEXT_WOS")
  context_list <- lapply(split(v, v$Context_Key), function(h) {
    stopifnot(length(unique(h$Species)) == 1L, length(unique(h$Duration_days)) == 1L)
    row <- data.frame(Context_Key = h$Context_Key[1], Reference_ID = h$Reference_ID[1],
                      Source_Context_ID = h$Source_Context_ID[1], Species = h$Species[1], Order = h$Order[1],
                      Duration_days = h$Duration_days[1], n_results = nrow(h),
                      Verification_Status = collapse(h$Verification_Status),
                      Strict_Context = all(h$Verification_Status %in% strict_status), stringsAsFactors = FALSE)
    for (metal in metals) {
      hh <- h[h$Metal == metal, , drop = FALSE]
      row[[metal]] <- if (nrow(hh)) exp(mean(log(hh$LC50_umol_L))) else NA_real_
      row[[paste0("IDs_", metal)]] <- collapse(hh$Result_ID)
    }
    row
  })
  ctx <- do.call(rbind, context_list); rownames(ctx) <- NULL
  write_table(ctx, "verified_contexts", dirs[["audit"]])
  # Presence in this ledger is not the same as inclusion in a primary LMM.
  coverage <- d[required]
  coverage$In_verified_metal_ledger <- coverage$Result_ID %in% v$Result_ID
  coverage$Metal_layer_status <- ifelse(coverage$In_verified_metal_ledger,
                                        "INCLUDED_IN_CURRENT_VERIFIED_LEDGER", "NOT_IN_INCLUDED_LEDGER_REVIEW_EXISTING_DECISIONS")
  write_table(coverage, "record_coverage", dirs[["audit"]])
  
  
  # ============================================================
  # 4. COVARIANCE AND REFERENCE BOOTSTRAP FUNCTIONS
  # ============================================================
  # Equal total weight per Reference, equal context weights within each Reference.
  # These weights define a descriptive target, not inverse-variance precision weights.
  calc <- function(x) {
    n <- nrow(x); refs <- unique(x$Reference_ID); G <- length(refs)
    ans <- list(n_contexts = n, n_references = G, n_varying_references = 0L,
                r_equal_reference = NA_real_, r_equal_context = NA_real_, spearman_context = NA_real_,
                total_log_covariance = NA_real_, between_reference_covariance = NA_real_,
                within_reference_covariance = NA_real_, between_over_total_covariance = NA_real_,
                r_within_centered = NA_real_)
    if (!n) return(ans)
    a <- log(x$A_value); b <- log(x$B_value)
    w <- 1 / as.numeric(table(x$Reference_ID)[x$Reference_ID]); w <- w / sum(w)
    ma <- ave(a, x$Reference_ID, FUN = mean); mb <- ave(b, x$Reference_ID, FUN = mean)
    da <- a - sum(w * a); db <- b - sum(w * b)
    wa <- a - ma; wb <- b - mb
    total <- sum(w * da * db)
    within <- sum(w * wa * wb)
    between <- sum(w * (ma - sum(w*a)) * (mb - sum(w*b)))
    stopifnot(abs(total - within - between) < 1e-9 * max(1, abs(total)))
    ans$total_log_covariance <- total; ans$between_reference_covariance <- between
    ans$within_reference_covariance <- within
    if (abs(total) > 1e-12) ans$between_over_total_covariance <- between / total
    ans$n_varying_references <- sum(vapply(split(seq_len(n), x$Reference_ID), function(i)
      length(i) > 1L && sd(a[i]) > 1e-10 && sd(b[i]) > 1e-10, logical(1)))
    den <- sqrt(sum(w * da^2) * sum(w * db^2))
    if (den > 1e-12) ans$r_equal_reference <- max(-1, min(1, total / den))
    if (n > 1L && sd(a) > 1e-12 && sd(b) > 1e-12) {
      ans$r_equal_context <- cor(a, b)
      ans$spearman_context <- cor(a, b, method = "spearman")
    }
    denw <- sqrt(sum(w * wa^2) * sum(w * wb^2))
    if (denw > 1e-12) ans$r_within_centered <- max(-1, min(1, within / denw))
    ans
  }
  bootstrap_reference <- function(x) {
    G <- length(unique(x$Reference_ID))
    empty <- c(bootstrap_low = NA_real_, bootstrap_high = NA_real_, bootstrap_valid = 0)
    if (G < 3L) return(empty) # operational reporting floor; NOT an adequacy theorem
    moments <- t(vapply(split(x, x$Reference_ID), function(h) {
      a <- log(h$A_value); b <- log(h$B_value)
      c(mean(a), mean(b), mean(a*a), mean(b*b), mean(a*b))
    }, numeric(5)))
    ix <- matrix(sample.int(G, B * G, replace = TRUE), nrow = B)
    mm <- sapply(seq_len(5L), function(j) rowMeans(matrix(moments[as.vector(ix), j], nrow = B)))
    den <- sqrt(pmax(0, (mm[,3]-mm[,1]^2) * (mm[,4]-mm[,2]^2)))
    ok <- is.finite(den) & den > 1e-10
    br <- pmax(-1, pmin(1, (mm[ok,5]-mm[ok,1]*mm[ok,2]) / den[ok]))
    if (!length(br)) return(empty)
    c(bootstrap_low = unname(quantile(br, .025)), bootstrap_high = unname(quantile(br, .975)),
      bootstrap_valid = length(br))
  }
  
  # ============================================================
  # 5. FOUR SCOPES AND LEAVE-ONE-REFERENCE-OUT
  # ============================================================
  set.seed(seed, kind = RNG_KIND, normal.kind = RNG_NORMAL_KIND,
           sample.kind = RNG_SAMPLE_KIND)
  run_rng <- RNGkind()
  stopifnot(identical(unname(run_rng),
                      c(RNG_KIND, RNG_NORMAL_KIND, RNG_SAMPLE_KIND)))
  subsets <- list(approved_ledger = ctx,
                  strict_verification = ctx[ctx$Strict_Context, , drop = FALSE],
                  duration_96h = ctx[ctx$Duration_days == 4, , drop = FALSE],
                  exclude_order_only_labels = ctx[!is.na(ctx$Species) & !is.na(ctx$Order) &
                                                    ctx$Species != ctx$Order, , drop = FALSE])
  summaries <- list(); loros <- list(); pair_data <- list()
  for (sn in names(subsets)) {
    pool <- subsets[[sn]]
    for (ab in pairs) {
      x <- pool[is.finite(pool[[ab[1]]]) & is.finite(pool[[ab[2]]]), , drop = FALSE]
      x$A_value <- x[[ab[1]]]; x$B_value <- x[[ab[2]]]
      pair <- paste(ab, collapse = "-"); st <- calc(x)
      status <- if (st$n_contexts == 0) "NO_VERIFIED_PAIR" else if (st$n_contexts == 1)
        "ONE_CONTEXT_NO_CORRELATION" else if (!is.finite(st$r_equal_reference))
          "ZERO_VARIATION_NO_CORRELATION" else if (st$n_references == 1)
            "ONE_REFERENCE_DESCRIPTIVE" else if (st$n_contexts == 2)
              "TWO_CONTEXTS_CORRELATION_DEGENERATE" else "EXPLORATORY_DESCRIPTIVE"
      boot <- c(bootstrap_low = NA_real_, bootstrap_high = NA_real_, bootstrap_valid = 0)
      if (sn == "approved_ledger" && is.finite(st$r_equal_reference)) boot <- bootstrap_reference(x)
      summaries[[length(summaries)+1L]] <- data.frame(Subset = sn, Pair = pair,
                                                      as.data.frame(st), as.list(boot), Status = status, stringsAsFactors = FALSE)
      if (sn != "approved_ledger" || !nrow(x)) next
      x$Pair <- pair; pair_data[[pair]] <- x
      for (ref in unique(x$Reference_ID)) {
        ll <- calc(x[x$Reference_ID != ref, , drop = FALSE])
        loros[[length(loros)+1L]] <- data.frame(Pair = pair, Omitted_Reference = ref,
                                              as.data.frame(ll), stringsAsFactors = FALSE)
      }
    }
  }
  
  summary <- do.call(rbind, summaries); rownames(summary) <- NULL
  loro <- do.call(rbind, loros)
  pair_values <- do.call(rbind, pair_data)

  # Frozen support checks for the final approved-ledger scope.
  approved_core <- summary[
    summary$Subset == "approved_ledger" &
      summary$Pair %in% c("Cu-Cd", "Cu-Zn", "Cu-Hg", "Cd-Zn", "Cd-Hg"),
    ,
    drop = FALSE
  ]
  approved_core <- approved_core[
    match(c("Cu-Cd", "Cu-Zn", "Cu-Hg", "Cd-Zn", "Cd-Hg"), approved_core$Pair),
    ,
    drop = FALSE
  ]
  stopifnot(
    identical(as.integer(approved_core$n_contexts), c(21L, 6L, 6L, 7L, 5L)),
    identical(as.integer(approved_core$n_references), c(10L, 4L, 4L, 4L, 3L))
  )
  
  # ============================================================
  # 6. RESULT TABLES AND DATA CHECKS
  # ============================================================
  write_table(summary, "covariance_all_pairs_and_scopes")
  write_table(loro, "LORO")
  write_table(pair_values, "pair_context_values")
  status_labels <- c(
    NO_VERIFIED_PAIR = "No verified pair",
    ONE_CONTEXT_NO_CORRELATION = "One context: correlation unavailable",
    ZERO_VARIATION_NO_CORRELATION = "No variation: correlation unavailable",
    ONE_REFERENCE_DESCRIPTIVE = "One reference: descriptive only",
    TWO_CONTEXTS_CORRELATION_DEGENERATE = "Two contexts: correlation is degenerate",
    EXPLORATORY_DESCRIPTIVE = "Exploratory descriptive result")
  overview <- summary[summary$Subset == "approved_ledger",
                      c("Pair", "n_contexts", "n_references", "r_equal_reference",
                        "bootstrap_low", "bootstrap_high", "bootstrap_valid", "Status")]
  overview$Interpretation <- unname(status_labels[overview$Status])
  write_table(overview, "main_scope_overview")
  qa <- data.frame(Check = c("Unique_Result_IDs", "Verified_IDs_in_combined",
                             "Matched_values_and_metadata", "Covariance_decomposition_identity",
                             "All_pairs_in_all_four_scopes", "Approved_ledgers_agree"),
                   Pass = c(TRUE, TRUE, TRUE, TRUE,
                            nrow(summary) == length(pairs) * length(subsets), TRUE))
  stopifnot(all(qa$Pass))
  write_table(qa, "internal_validation", dirs[["audit"]])
  write_table(data.frame(Quantity = c("Combined_rows", "Harmonized_positive_rows",
                                      "Other_combined_rows", "Verified_rows", "Verified_contexts", "Verified_references"),
                         Count = c(nrow(combined), nrow(d), sum(!keep), nrow(v), nrow(ctx),
                                   length(unique(ctx$Reference_ID)))), "input_counts", dirs[["audit"]])
  refs <- sort(unique(ctx$Reference_ID))
  index <- data.frame(Reference_ID = refs, stringsAsFactors = FALSE)
  for (f in intersect(c("Title", "External_Reference_Label", "External_DOI"), names(v)))
    index[[f]] <- vapply(refs, function(ref) collapse(v[[f]][v$Reference_ID == ref]), character(1))
  write_table(index, "reference_index", dirs[["audit"]])
  
  # ============================================================
  # 7. PNG FIGURES
  # ============================================================
  report_folder <- render_covariance_figures(summary, pair_values, out)
  write_covariance_word(summary, loro, report_folder)
  
  # ============================================================
  # 8. PROVENANCE, SAVED RESULTS AND COMPLETION
  # ============================================================
  if (!identical(hash(inputs), input_hashes)) stop("An input changed during this run.")
  write_table(data.frame(Role = names(inputs), File = unname(inputs), MD5 = input_hashes),
              "input_checksums", dirs[["audit"]])
  write_table(data.frame(Bootstrap_draws = B, Seed = seed,
                         RNG_kind = run_rng[1L], RNG_normal_kind = run_rng[2L], RNG_sample_kind = run_rng[3L],
                         Bootstrap_scope = "approved_ledger", Bootstrap_min_references = 3L,
                         LORO_scope = "approved_ledger", Figure_layout = "grouped_pairs_and_separate_covariance_scales"),
              "run_settings", dirs[["audit"]])
  capture.output(sessionInfo(), file = file.path(dirs[["audit"]], "sessionInfo.txt"))
  saveRDS(list(summary = summary, contexts = ctx, LORO = loro,
               pair_context_values = pair_values, seed = seed, B = B, rng = run_rng,
               inputs = inputs, input_hashes = input_hashes),
          file.path(dirs[["objects"]], "metal_covariance_analysis_objects.rds"))
  writeLines(c(
    "Start with 02_tables/metal_covariance_main_scope_overview.csv.",
    "This run writes to outputs/14_analyze_metal_covariance.",
    "RNG settings are explicitly fixed and recorded in 01_audit/metal_covariance_run_settings.csv.",
    "For identical inputs and computation environment, bootstrap sampling no longer inherits the session RNG kind.",
    "Complementary, exploratory covariance analysis; no p-values calculated.",
    "A point is one Reference_ID + Source_Context_ID; within-metal values use geometric means.",
    "Each reference has equal total weight; contexts within a reference share that weight.",
    "Covariance uses normalized weighted moments, not an unbiased sample-covariance denominator.",
    "Bootstrap bounds are 2.5% and 97.5% percentiles for equal-reference correlation r.",
    "They are not intervals for covariance components. Valid resample counts are reported.",
    "Bootstrap and leave-one-reference-out apply to approved_ledger only.",
    "The three-reference bootstrap floor is a reporting rule, not proof of sufficient support.",
    "Very few references can yield unstable or degenerate bootstrap distributions.",
    "A reference with one context contributes no within-reference variation.",
    "Two contexts give a degenerate correlation; one reference supports no between-reference inference.",
    "Between/total covariance may be negative or above one; it is not variance explained.",
    "Strict verification is applied to the entire context, as in the original script.",
    "The order-label sensitivity excludes missing Species/Order as well as equal labels, as before.",
    "Source weighting does not establish independence. Associations do not demonstrate mixture effects.",
    "Pairs with at least two contexts are plotted; all pairs remain in the Word/support tables.",
    "Source names/titles are in 01_audit/metal_covariance_reference_index.csv; author-year labels are not inferred.",
    "Source-verification decisions are frozen in the curated ledgers and are not changed here.",
    "This plan followed earlier exploratory results; it is not a preregistered analysis.",
    "This public-repository script preserves the final thesis covariance calculations and reproducibility settings."
  ), file.path(out, "README_analysis_outputs.txt"))

  input_output_manifest <- data.frame(
    Role = c(
      "Input: combined harmonized dataset",
      "Input: source-verified rows",
      "Input: source-adjudication ledger",
      "Output: main-scope overview",
      "Output: LORO results",
      "Output: reusable covariance objects"
    ),
    Path = c(
      main_file,
      verified_file,
      adjudication_file,
      file.path(dirs[["tables"]], "metal_covariance_main_scope_overview.csv"),
      file.path(dirs[["tables"]], "metal_covariance_LORO.csv"),
      file.path(dirs[["objects"]], "metal_covariance_analysis_objects.rds")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    input_output_manifest,
    file.path(out, "input_output_manifest.csv"),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  writeLines(
    c(
      "STATUS: 14 ANALYZE METAL COVARIANCE = PASS",
      paste0("Approved rows: ", nrow(v)),
      paste0("Verified contexts: ", nrow(ctx)),
      paste0("References: ", length(unique(ctx$Reference_ID))),
      paste0("Bootstrap draws: ", B),
      paste0("Seed: ", seed)
    ),
    file.path(out, "RUN_COMPLETE.txt")
  )

  cat("\nSTATUS: 14 ANALYZE METAL COVARIANCE = PASS\n")
cat("Output folder: ", out, "\n", sep = "")
  print(overview[c("Pair", "n_contexts", "n_references", "r_equal_reference", "Interpretation")],
        row.names = FALSE)
})
