# Called by RUN_REVISION.R in a fresh Rscript --vanilla process.
# Existing generated output folders are moved, not destroyed, before rebuilding.
# Only the explicit generated-path allowlist below is moved. Raw data and .git
# are never moved. Optional database audits are not rerun or removed.
run_revision_checked <- function() {
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  if (!file.exists("copepod-metal-toxicity.Rproj")) stop("Repository root not found.")
  report <- file.path(root, "REVISION_CHECK.txt")
  log_file <- file.path(root, "REVISION_CONSOLE.log")
  writeLines(c("INTEGRATED SOURCE AND REPORTING CORRECTION", paste("Started:", Sys.time()),
               paste("Project:", root), "STATUS: PREFLIGHT", "RUNTIME PASS IS NOT A SCIENTIFIC/RELEASE APPROVAL."), report, useBytes = TRUE)
  append_report <- function(x) cat(paste0(paste(x, collapse = "\n"), "\n"), file = report, append = TRUE)
  console <- file(log_file, open = "wt", encoding = "UTF-8")
  sink_level <- sink.number()
  sink(console, split = TRUE)
  on.exit({ while (sink.number() > sink_level) sink(); close(console) }, add = TRUE)
  stage <- "PREFLIGHT"
  backup <- NA_character_
  result <- tryCatch({
    required_packages <- c("tidyverse", "readxl", "readr", "dplyr", "tibble", "tidyr", "stringr",
                           "lme4", "lmerTest", "emmeans", "pbkrtest", "gt")
    missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop(paste("Missing package(s):", paste(missing, collapse = ", "), "No packages were installed or updated."))
    manifest_path <- "documentation/revision_20261007_input_manifest.csv"
    manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("path", "mode", "md5") %in% names(manifest))) stop("Input manifest schema is invalid.")
    if (anyDuplicated(manifest$path)) stop("Duplicate paths in the input manifest.")
    if (any(!file.exists(manifest$path))) stop(paste("Missing installed input/code files:", paste(manifest$path[!file.exists(manifest$path)], collapse = "; ")))
    file_hash <- function(path, mode) {
      if (mode == "binary") return(unname(tools::md5sum(path)))
      con <- file(path, "rb"); on.exit(close(con), add = TRUE)
      bytes <- readBin(con, "raw", n = file.info(path)$size)
      if (length(bytes) >= 3L && identical(bytes[1:3], as.raw(c(239, 187, 191)))) bytes <- bytes[-(1:3)]
      bytes <- bytes[bytes != as.raw(13L)] # normalize CRLF/LF, not scientific content
      tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
      writeBin(bytes, tmp)
      unname(tools::md5sum(tmp))
    }
    got <- vapply(seq_len(nrow(manifest)), function(i) file_hash(manifest$path[i], manifest$mode[i]), character(1))
    mismatch <- which(tolower(got) != tolower(manifest$md5))
    if (length(mismatch)) stop(paste("Files differ from the assembled correction package:", paste(manifest$path[mismatch], collapse = "; "), "No generated files have been moved yet."))
    append_report(paste("Installed input/code fingerprints matched:", nrow(manifest)))
    r_paths <- manifest$path[grepl("[.]R$", manifest$path)]
    for (p in r_paths) parse(p, encoding = "UTF-8")
    append_report(paste("R syntax parsed:", length(r_paths), "files"))
    test_env <- new.env(parent = globalenv())
    sys.source("tests/test_reporting_helpers.R", envir = test_env)
    unit_results <- test_env$run_reporting_unit_tests(root)
    append_report(c(paste("Reporting unit tests:", nrow(unit_results), "PASS"), capture.output(print(unit_results, row.names = FALSE))))
    # Read the pipeline definition without executing model/data scripts.
    exprs <- parse("scripts/00_run_core_pipeline.R", encoding = "UTF-8")
    assignment <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
      is.symbol(x[[2L]]) && identical(as.character(x[[2L]]), "pipeline"), as.list(exprs))
    if (length(assignment) != 1L) stop("Unique pipeline definition not found.")
    pipeline <- eval(assignment[[1L]][[3L]], envir = baseenv())
    if (length(pipeline) != 16L) stop("Expected 15 core analyses and one figure step.")
    output_dirs <- paste0("outputs/", sub("[.]R$", "", basename(vapply(pipeline[1:15], function(x) x$script, character(1)))))
    generated_paths <- c("data/processed/ecotox_harmonized.csv", "data/processed/combined_ecotox_wos_harmonized.csv",
                         "outputs/00_run_core_pipeline", output_dirs, "outputs/16_make_main_thesis_figures", "results/figures")
    present <- generated_paths[file.exists(generated_paths) | dir.exists(generated_paths)]
    stage <- "BACKUP_EXISTING_GENERATED_FILES"
    backup <- file.path(dirname(root), paste0(basename(root), "_generated_backup_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    if (file.exists(backup) || dir.exists(backup)) stop("Backup path already exists; nothing will be overwritten.")
    if (!dir.create(backup, recursive = FALSE)) stop("Cannot create generated-file backup beside the project.")
    append_report(paste("Generated-file backup:", backup))
    cat("Generated-file backup: ", backup, "\n", sep = "")
    write.csv(data.frame(path = character(), status = character()), file.path(backup, "MOVED_FILES.csv"), row.names = FALSE)
    moved <- character()
    for (p in present) {
      dest <- file.path(backup, p)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      if (!file.rename(file.path(root, p), dest)) stop(paste("Could not move generated path:", p, "Close open result files and send the report. Already moved paths are safely in the backup."))
      moved <- c(moved, p)
      write.csv(data.frame(path = moved, status = "MOVED_TO_BACKUP"), file.path(backup, "MOVED_FILES.csv"), row.names = FALSE)
    }
    append_report(paste("Previous generated paths backed up:", length(moved)))
    # Script 01's optional comparison to an old processed CSV would reject the
    # intentional C6/description correction. The old CSV is backed up above;
    # below we explicitly compare old/new records and permit only those changes.
    stage <- "CORE_PIPELINE"
    append_report("STATUS: RUNNING_CORE_PIPELINE")
    core_env <- new.env(parent = globalenv())
    sys.source("scripts/00_run_core_pipeline.R", envir = core_env, chdir = FALSE)
    stage <- "POST_RUN_CHECKS"
    qa_dir <- "outputs/00_run_core_pipeline"
    write.csv(unit_results, file.path(qa_dir, "reporting_unit_tests.csv"), row.names = FALSE)
    dump_csv <- function(p) {
      append_report(c("", paste("FILE:", p), readLines(p, warn = FALSE, encoding = "UTF-8")))
    }
    read_table <- function(p) read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
    must <- function(condition, msg) if (!isTRUE(condition)) stop(msg, call. = FALSE)
    compare_dataset <- function(rel) {
      old_path <- file.path(backup, rel)
      if (!file.exists(old_path)) return(data.frame(file = rel, changed_cells = NA_integer_, status = "NO_PRIOR_CSV_AVAILABLE"))
      old <- read.csv(old_path, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"))
      new <- read.csv(rel, colClasses = "character", check.names = FALSE, na.strings = c("", "NA"))
      must(identical(names(old), names(new)) && nrow(old) == nrow(new), paste("Unexpected schema/row-count change:", rel))
      key <- if ("Result_ID" %in% names(new)) "Result_ID" else "Result_Number"
      must(!anyNA(new[[key]]) && !anyDuplicated(new[[key]]) && setequal(new[[key]], old[[key]]), paste("Record identity mismatch:", rel))
      old <- old[match(new[[key]], old[[key]]), , drop = FALSE]
      changes <- list()
      for (nm in names(new)) {
        a <- old[[nm]]; b <- new[[nm]]
        same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
        same[is.na(same)] <- FALSE
        an <- suppressWarnings(as.numeric(a)); bn <- suppressWarnings(as.numeric(b))
        num_ok <- is.finite(an) & is.finite(bn) & abs(an - bn) <= 1e-9 * pmax(1, abs(an), abs(bn))
        num_ok[is.na(num_ok)] <- FALSE
        for (i in which(!(same | num_ok))) {
          ref <- if ("Reference_ID" %in% names(new)) new$Reference_ID[i] else new$Reference_Number[i]
          is_obrien <- ref %in% c("E_2977", "2977")
          allowed <- is_obrien && ((nm == "Lifestage" && new[[key]][i] %in% c("E_R115201", "115201") &&
                                    b[i] == "Adult" && a[i] == "Copepodite") || nm == "Harmonization_Rule")
          changes[[length(changes) + 1L]] <- data.frame(file = rel, id = new[[key]][i], column = nm,
            old = a[i], new = b[i], allowed_source_correction = allowed, stringsAsFactors = FALSE)
        }
      }
      if (length(changes)) {
        detail <- do.call(rbind, changes)
        write.csv(detail, file.path(qa_dir, paste0(basename(rel), "_source_correction_changes.csv")), row.names = FALSE, na = "")
        if (!all(detail$allowed_source_correction)) {
          append_report(c("UNEXPECTED DATA CHANGE:", capture.output(print(detail[!detail$allowed_source_correction, ], row.names = FALSE))))
          stop(paste("Unplanned dataset change found:", rel))
        }
      }
      data.frame(file = rel, changed_cells = length(changes), status = "ONLY_ALLOWED_SOURCE_CHANGES")
    }
    dataset_check <- do.call(rbind, lapply(c("data/processed/ecotox_harmonized.csv", "data/processed/combined_ecotox_wos_harmonized.csv"), compare_dataset))
    write.csv(dataset_check, file.path(qa_dir, "source_correction_data_check.csv"), row.names = FALSE)
    append_report(c("", "DATA REBUILD CHECK:", capture.output(print(dataset_check, row.names = FALSE))))
    combined <- read_table("data/processed/combined_ecotox_wos_harmonized.csv")
    c6 <- combined[combined$Result_ID == "E_R115201", , drop = FALSE]
    must(nrow(c6) == 1L && identical(c6$Lifestage, "Adult") && isTRUE(all.equal(as.numeric(c6$LC50_umol_L), 12)), "C6 identity/stage/value check failed.")
    obrien <- combined[combined$Reference_ID == "E_2977", , drop = FALSE]
    must(nrow(obrien) == 13L && all(is.na(obrien$Temperature)), "O'Brien temperature/source scope changed unexpectedly.")
    ledger <- read_table("data/curated/source_verification/OBrien_1988_source_review_20261007.csv")
    must(setequal(ledger$Result_ID, obrien$Result_ID), "O'Brien source ledger membership mismatch.")
    must(all(abs(as.numeric(obrien$LC50_umol_L[match(ledger$Result_ID, obrien$Result_ID)]) - ledger$Source_LC50_umol_L) < 1e-10), "O'Brien LC50 values differ from the source ledger.")
    stage_support <- read_table("outputs/11_fit_stage_model/01_model/stage_model_dataset_overall_support.csv")
    must(stage_support$Results == 70L && stage_support$References == 31L && stage_support$Species == 19L, "Stage support check failed.")
    stage_obj <- readRDS("outputs/11_fit_stage_model/06_objects/stage_model_analysis_objects.rds")
    must(sum(stage_obj$stage_model_data$Result_ID == "E_R115201") == 1L, "C6 is missing/duplicated in the Stage model.")
    for (method in c("LORO", "LOSO")) {
      full <- read_table(paste0("outputs/13_fit_order_model/03_robustness/order_", method, "_all_refits.csv"))
      extremes <- read_table(paste0("outputs/13_fit_order_model/03_robustness/order_", method, "_interaction_extremes.csv"))
      id <- if (method == "LORO") "removed_reference" else "removed_species"
      v <- full$ratio_of_ratios
      valid <- which(!is.na(v))
      chosen <- if (length(valid)) unique(c(valid[which.min(v[valid])], valid[which.max(v[valid])])) else integer()
      must(identical(as.character(extremes[[id]]), as.character(full[[id]][chosen])), paste("Order extremes incorrect:", method))
    }
    index <- read_table("outputs/16_make_main_thesis_figures/figure_index.csv")
    must(nrow(index) == 19L && !anyDuplicated(index$Filename_stem), "Figure index incomplete/duplicated.")
    figure_paths <- as.vector(outer(paste0("results/figures/", index$Filename_stem), c(".pdf", ".png"), paste0))
    must(all(file.exists(figure_paths)) && all(file.info(figure_paths)$size > 0), "One or more regenerated figures is missing/empty.")
    append_report(c("C6: Adult; LC50 = 12 micromol/L; missing numerical assay temperature retained.",
                    "Stage model: 70 records / 31 references / 19 species groups.",
                    "Order min/max row selection: PASS for LORO and LOSO.",
                    "Figure suite regenerated: 19 figures / 38 PDF-PNG files.",
                    "The figure script itself is unchanged in this patch."))
    key_outputs <- c(
      "outputs/00_run_core_pipeline/pipeline_run_log.csv",
      "outputs/04_fit_pooled_metal_model/01_model/metal_pairwise_lc50_ratios.csv",
      "outputs/04_fit_pooled_metal_model/03_robustness/reference_fixed_effects_fit_checks.csv",
      "outputs/04_fit_pooled_metal_model/03_robustness/deletion_direction_and_interval_summary.csv",
      "outputs/11_fit_stage_model/01_model/stage_baseline_metrics.csv",
      "outputs/11_fit_stage_model/01_model/stage_primary_numerical_diagnostics.csv",
      "outputs/11_fit_stage_model/03_robustness/stage_LOSO_summary.csv",
      "outputs/13_fit_order_model/01_model/order_baseline_metrics.csv",
      "outputs/13_fit_order_model/01_model/order_primary_numerical_diagnostics.csv",
      "outputs/13_fit_order_model/03_robustness/order_LORO_interaction_extremes.csv",
      "outputs/13_fit_order_model/03_robustness/order_LOSO_interaction_extremes.csv")
    must(all(file.exists(key_outputs)), "A required correction report was not generated.")
    for (p in key_outputs) dump_csv(p)
    append_report(c("", "R SESSION:", capture.output(sessionInfo()), "", "STATUS: PASS_LOCAL_CORRECTION_RUN",
                    "This confirms execution and the listed regression checks, not all scientific assumptions or source quality.",
                    "No GitHub/Zenodo upload was performed. Word and release metadata have not been updated.",
                    paste("Finished:", Sys.time())))
    cat("\nPASS_LOCAL_CORRECTION_RUN. Send REVISION_CHECK.txt before publication.\n")
    TRUE
  }, error = function(e) {
    append_report(c("", "STATUS: FAILED", paste("Stage:", stage), paste("Error:", conditionMessage(e)),
                    paste("Generated backup (if created):", backup), "Do not push/publish partial outputs. Send this report."))
    cat("\nFAILED at ", stage, ": ", conditionMessage(e), "\nSend REVISION_CHECK.txt.\n", sep = "")
    FALSE
  })
  result
}
ok <- run_revision_checked()
quit(save = "no", status = if (isTRUE(ok)) 0L else 1L, runLast = FALSE)
