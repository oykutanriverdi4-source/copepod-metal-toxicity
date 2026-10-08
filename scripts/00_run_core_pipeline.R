# ============================================================
# 00_run_core_pipeline.R
# ============================================================
# Purpose
# Run the cleaned core thesis workflow in its frozen order:
#   01-15 analytical scripts
#   + the unified thesis figure suite
#
# This master runner does NOT contain scientific calculations.
# It only sources the already-reviewed scripts in sequence and checks
# that each step produces its expected completion marker/output.
#
# Database-audit scripts are intentionally NOT run here because the
# WoS and EnviroTox raw exports are local/external inputs that are not
# redistributed in the public repository. Their workflows remain
# available separately under scripts/database_audits/.
#
# Reproducibility rule
# Open copepod-metal-toxicity.Rproj first and run this script from the
# repository root.
# ============================================================


# ------------------------------------------------------------
# 0. PROJECT ROOT CHECK
# ------------------------------------------------------------

project_file <- "copepod-metal-toxicity.Rproj"

if (!file.exists(project_file)) {
  stop(
    "Repository root not detected. Open ",
    project_file,
    " in RStudio and run this script again."
  )
}

repo_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

cat(
  "Repository root:\n",
  repo_root,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 1. CORE PIPELINE DEFINITION
# ------------------------------------------------------------

pipeline <- list(
  list(
    id = "01",
    name = "Prepare ECOTOX data",
    script = "scripts/analysis/01_prepare_ecotox_data.R",
    expected = "data/processed/ecotox_harmonized.csv"
  ),
  list(
    id = "02",
    name = "Integrate WoS verified records",
    script = "scripts/analysis/02_integrate_wos_verified_records.R",
    expected = "data/processed/combined_ecotox_wos_harmonized.csv"
  ),
  list(
    id = "03",
    name = "Audit metal support",
    script = "scripts/analysis/03_audit_metal_support.R",
    expected = "outputs/03_audit_metal_support/03_objects/metal_support_audit_objects.rds"
  ),
  list(
    id = "04",
    name = "Fit pooled metal model",
    script = "scripts/analysis/04_fit_pooled_metal_model.R",
    expected = "outputs/04_fit_pooled_metal_model/RUN_COMPLETE.txt"
  ),
  list(
    id = "05",
    name = "Compare metals within References",
    script = "scripts/analysis/05_compare_metals_within_references.R",
    expected = "outputs/05_compare_metals_within_references/RUN_COMPLETE.txt"
  ),
  list(
    id = "06",
    name = "Compare exposure duration",
    script = "scripts/analysis/06_compare_exposure_duration.R",
    expected = "outputs/06_compare_exposure_duration/RUN_COMPLETE.txt"
  ),
  list(
    id = "07",
    name = "Assess salinity and temperature",
    script = "scripts/analysis/07_assess_salinity_temperature.R",
    expected = "outputs/07_assess_salinity_temperature/RUN_COMPLETE.txt"
  ),
  list(
    id = "08",
    name = "Assess pH evidence",
    script = "scripts/analysis/08_assess_ph_evidence.R",
    expected = "outputs/08_assess_ph_evidence/RUN_COMPLETE.txt"
  ),
  list(
    id = "09",
    name = "Compare stages within References",
    script = "scripts/analysis/09_compare_stages_within_references.R",
    expected = "outputs/09_compare_stages_within_references/RUN_COMPLETE.txt"
  ),
  list(
    id = "10",
    name = "Audit stage support",
    script = "scripts/analysis/10_audit_stage_support.R",
    expected = "outputs/10_audit_stage_support/RUN_COMPLETE.txt"
  ),
  list(
    id = "11",
    name = "Fit developmental-stage model",
    script = "scripts/analysis/11_fit_stage_model.R",
    expected = "outputs/11_fit_stage_model/RUN_COMPLETE.txt"
  ),
  list(
    id = "12",
    name = "Audit taxonomic-order support",
    script = "scripts/analysis/12_audit_order_support.R",
    expected = "outputs/12_audit_order_support/RUN_COMPLETE.txt"
  ),
  list(
    id = "13",
    name = "Fit taxonomic-order model",
    script = "scripts/analysis/13_fit_order_model.R",
    expected = "outputs/13_fit_order_model/RUN_COMPLETE.txt"
  ),
  list(
    id = "14",
    name = "Analyze matched-metal covariance",
    script = "scripts/analysis/14_analyze_metal_covariance.R",
    expected = "outputs/14_analyze_metal_covariance/RUN_COMPLETE.txt"
  ),
  list(
    id = "15",
    name = "Audit sex-specific evidence",
    script = "scripts/analysis/15_audit_sex_specific_evidence.R",
    expected = "outputs/15_audit_sex_specific_evidence/RUN_COMPLETE.txt"
  ),
  list(
    id = "FIG",
    name = "Generate main thesis figures",
    script = "scripts/figures/01_make_main_thesis_figures.R",
    expected = "outputs/16_make_main_thesis_figures/RUN_COMPLETE.txt"
  )
)


# ------------------------------------------------------------
# 2. PRE-RUN FILE CHECKS
# ------------------------------------------------------------

pipeline_scripts <- vapply(
  pipeline,
  function(x) x$script,
  character(1)
)

missing_scripts <- pipeline_scripts[
  !file.exists(pipeline_scripts)
]

if (length(missing_scripts)) {
  stop(
    "Missing pipeline script(s):\n",
    paste(
      paste0(" - ", missing_scripts),
      collapse = "\n"
    )
  )
}

required_starting_inputs <- c(
  "data/raw/ecotox/ecotox_raw_export.xlsx",
  "data/curated/ecotox/ecotox_working_master_harmonized.xlsx",
  "data/curated/wos/wos_source_verified_records.xlsx",
  "data/curated/source_verification/combined_source_adjudication_rows.csv",
  "data/curated/source_verification/combined_source_verified_rows_used.csv"
)

missing_starting_inputs <- required_starting_inputs[
  !file.exists(required_starting_inputs)
]

if (length(missing_starting_inputs)) {
  stop(
    "Missing required starting input file(s):\n",
    paste(
      paste0(" - ", missing_starting_inputs),
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 3. RUN HELPER
# ------------------------------------------------------------

run_step <- function(step) {

  cat(
    "\n============================================================\n",
    "STEP ",
    step$id,
    " | ",
    step$name,
    "\n",
    "Script: ",
    step$script,
    "\n",
    "============================================================\n",
    sep = ""
  )

  start_time <- Sys.time()

  # Each script gets a clean child environment so temporary objects,
  # helper functions and internal variable names do not leak into the
  # following script. File-based outputs remain the communication
  # layer between pipeline steps.
  step_env <- new.env(
    parent = globalenv()
  )

  tryCatch(
    {
      sys.source(
        step$script,
        envir = step_env,
        chdir = FALSE
      )
    },
    error = function(e) {

      cat(
        "\nPIPELINE FAILED AT STEP ",
        step$id,
        " | ",
        step$name,
        "\n",
        sep = ""
      )

      stop(
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!file.exists(step$expected)) {
    stop(
      "Step ",
      step$id,
      " completed without an R error but the expected output was not found:\n",
      step$expected,
      call. = FALSE
    )
  }

  elapsed <- difftime(
    Sys.time(),
    start_time,
    units = "secs"
  )

  cat(
    "\nSTEP ",
    step$id,
    " PASS | ",
    round(
      as.numeric(elapsed),
      1
    ),
    " seconds\n",
    sep = ""
  )

  data.frame(
    Step = step$id,
    Name = step$name,
    Script = step$script,
    Expected_output = step$expected,
    Status = "PASS",
    Seconds = round(
      as.numeric(elapsed),
      3
    ),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 4. RUN CORE PIPELINE
# ------------------------------------------------------------

pipeline_start <- Sys.time()

run_log <- lapply(
  pipeline,
  run_step
)

run_log <- do.call(
  rbind,
  run_log
)

pipeline_elapsed <- difftime(
  Sys.time(),
  pipeline_start,
  units = "mins"
)


# ------------------------------------------------------------
# 5. FROZEN FINAL BENCHMARK CHECKS
# ------------------------------------------------------------

combined <- read.csv(
  "data/processed/combined_ecotox_wos_harmonized.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  nrow(combined) == 353L,
  sum(
    is.finite(
      suppressWarnings(
        as.numeric(
          combined$LC50_umol_L
        )
      )
    ) &
      suppressWarnings(
        as.numeric(
          combined$LC50_umol_L
        )
      ) > 0
  ) == 304L
)

stage_objects <- readRDS(
  "outputs/10_audit_stage_support/03_objects/stage_support_audit_objects.rds"
)

order_objects <- readRDS(
  "outputs/12_audit_order_support/03_objects/order_support_audit_objects.rds"
)

stopifnot(
  nrow(
    stage_objects$stage_model_data
  ) == 70L,
  nrow(
    order_objects$order_model_data
  ) == 101L
)

pooled_objects <- readRDS(
  "outputs/04_fit_pooled_metal_model/01_model/primary_pooled_analysis_objects.rds"
)

stopifnot(
  nrow(
    pooled_objects$primary_model_data
  ) == 131L
)


# ------------------------------------------------------------
# 6. PIPELINE LOG + COMPLETION MARKER
# ------------------------------------------------------------

pipeline_output_dir <- file.path(
  "outputs",
  "00_run_core_pipeline"
)

dir.create(
  pipeline_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  run_log,
  file.path(
    pipeline_output_dir,
    "pipeline_run_log.csv"
  ),
  row.names = FALSE,
  na = ""
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    pipeline_output_dir,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "STATUS: CORE PIPELINE = PASS",
    "Combined Results: 353",
    "Common-molar quantitative Results: 304",
    "Pooled 96-h Cu-Cd-Zn model Results: 131",
    "Developmental-stage model Results: 70",
    "Taxonomic-order model Results: 101",
    paste0(
      "Pipeline steps completed: ",
      nrow(run_log)
    ),
    paste0(
      "Elapsed minutes: ",
      round(
        as.numeric(
          pipeline_elapsed
        ),
        2
      )
    )
  ),
  file.path(
    pipeline_output_dir,
    "RUN_COMPLETE.txt"
  )
)


# ------------------------------------------------------------
# 7. FINAL CONSOLE SUMMARY
# ------------------------------------------------------------

cat(
  "\n\n============================================================\n",
  "CORE PIPELINE COMPLETE\n",
  "============================================================\n",
  "353 combined Results\n",
  "304 common-molar quantitative Results\n",
  "131 pooled-model Results\n",
  "70 developmental-stage model Results\n",
  "101 taxonomic-order model Results\n",
  nrow(run_log),
  " pipeline steps completed successfully\n",
  "Elapsed time: ",
  round(
    as.numeric(
      pipeline_elapsed
    ),
    2
  ),
  " minutes\n",
  "============================================================\n",
  sep = ""
)
