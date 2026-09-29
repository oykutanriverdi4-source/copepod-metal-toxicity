# ============================================================
# check_repository_before_release.R
# ============================================================
# Purpose
# Inspect the local repository before the first public release.
#
# This script does NOT modify files. It checks:
# - required repository files are present;
# - files ignored by .gitignore are not candidates for commit;
# - local WoS / EnviroTox raw exports are not tracked;
# - obsolete analysis scripts are absent from the commit candidate set;
# - obvious absolute Windows paths and old workflow names are absent
#   from text files that Git would include;
# - core numerical outputs still exist.
#
# Run from the repository root with copepod-metal-toxicity.Rproj open.
# ============================================================


# ------------------------------------------------------------
# 0. PROJECT CHECK
# ------------------------------------------------------------

if (!file.exists("copepod-metal-toxicity.Rproj")) {
  stop(
    "Repository root not detected. Open copepod-metal-toxicity.Rproj ",
    "and run this script from the repository root."
  )
}


# ------------------------------------------------------------
# 1. FILES THAT GIT WOULD CONSIDER
# ------------------------------------------------------------

find_git_executable <- function() {

  git_on_path <- Sys.which("git")

  if (nzchar(git_on_path)) {
    return(unname(git_on_path))
  }

  # GitHub Desktop ships its own Git on Windows. It is often not added
  # to the system PATH, so RStudio may not find `git` even though
  # GitHub Desktop works normally.
  local_appdata <- Sys.getenv(
    "LOCALAPPDATA",
    unset = ""
  )

  if (nzchar(local_appdata)) {

    github_desktop_root <- file.path(
      local_appdata,
      "GitHubDesktop"
    )

    if (dir.exists(github_desktop_root)) {

      app_dirs <- list.dirs(
        github_desktop_root,
        full.names = TRUE,
        recursive = FALSE
      )

      bundled_git <- file.path(
        app_dirs,
        "resources",
        "app",
        "git",
        "cmd",
        "git.exe"
      )

      bundled_git <- bundled_git[
        file.exists(bundled_git)
      ]

      if (length(bundled_git)) {

        # Prefer the most recently modified GitHub Desktop app folder.
        info <- file.info(
          dirname(
            dirname(
              dirname(
                dirname(
                  dirname(
                    bundled_git
                  )
                )
              )
            )
          )
        )

        return(
          bundled_git[
            order(
              info$mtime,
              decreasing = TRUE,
              na.last = TRUE
            )[1]
          ]
        )
      }
    }
  }

  ""
}

git_executable <- find_git_executable()

if (!nzchar(git_executable)) {
  stop(
    "Git could not be located from R. GitHub Desktop may still work, ",
    "but this QA script needs a Git executable to inspect the files ",
    "that would be committed. Send me this message if it persists."
  )
}

cat(
  "Git executable detected:\n",
  git_executable,
  "\n\n",
  sep = ""
)

git_files <- system2(
  git_executable,
  c(
    "ls-files",
    "--cached",
    "--others",
    "--exclude-standard"
  ),
  stdout = TRUE,
  stderr = TRUE
)

git_files <- unique(
  git_files[
    nzchar(git_files)
  ]
)

git_files <- gsub(
  "\\\\",
  "/",
  git_files
)


# ------------------------------------------------------------
# 2. REQUIRED FILES
# ------------------------------------------------------------

required_files <- c(
  "README.md",
  ".gitignore",
  "LICENSE",
  "CITATION.cff",
  "copepod-metal-toxicity.Rproj",
  "scripts/00_run_core_pipeline.R",
  paste0(
    "scripts/analysis/",
    sprintf("%02d", 1:15),
    c(
      "_prepare_ecotox_data.R",
      "_integrate_wos_verified_records.R",
      "_audit_metal_support.R",
      "_fit_pooled_metal_model.R",
      "_compare_metals_within_references.R",
      "_compare_exposure_duration.R",
      "_assess_salinity_temperature.R",
      "_assess_ph_evidence.R",
      "_compare_stages_within_references.R",
      "_audit_stage_support.R",
      "_fit_stage_model.R",
      "_audit_order_support.R",
      "_fit_order_model.R",
      "_analyze_metal_covariance.R",
      "_audit_sex_specific_evidence.R"
    )
  ),
  "scripts/figures/01_make_main_thesis_figures.R",
  "scripts/database_audits/wos/01_screen_wos_records.R",
  "scripts/database_audits/wos/02_check_wos_ecotox_overlap.R",
  "scripts/database_audits/wos/03_prepare_wos_unique_reference_review.R",
  "scripts/database_audits/envirotox/01_screen_envirotox_records.R",
  "scripts/database_audits/envirotox/02_check_envirotox_ecotox_overlap.R",
  "scripts/database_audits/envirotox/03_prepare_envirotox_manual_review.R",
  "documentation/wos_search.md",
  "documentation/envirotox_search.md",
  "documentation/data_provenance.md",
  "data/external/wos/README.md",
  "data/external/envirotox/README.md",
  "data/raw/ecotox/ecotox_raw_export.xlsx",
  "data/curated/ecotox/ecotox_working_master_harmonized.xlsx",
  "data/curated/wos/wos_source_verified_records.xlsx",
  "data/curated/source_verification/combined_source_adjudication_rows.csv",
  "data/curated/source_verification/combined_source_verified_rows_used.csv",
  "data/processed/ecotox_harmonized.csv",
  "data/processed/combined_ecotox_wos_harmonized.csv"
)

missing_required <- required_files[
  !file.exists(required_files)
]


# ------------------------------------------------------------
# 3. FILES THAT MUST NOT BE COMMITTED
# ------------------------------------------------------------

forbidden_exact <- c(
  "data/external/wos/wos_raw_export.xls",
  "data/external/envirotox/envirotox_raw_export.xlsx",
  "outputs/database_audits/wos/01_screen_wos_records/wos_preliminary_screening.xlsx",
  "outputs/database_audits/wos/02_check_wos_ecotox_overlap/wos_ecotox_overlap.xlsx",
  "outputs/database_audits/wos/03_prepare_wos_unique_reference_review/wos_unique_reference_review.xlsx",
  "outputs/database_audits/envirotox/01_screen_envirotox_records/envirotox_preliminary_screening.xlsx",
  "outputs/database_audits/envirotox/02_check_envirotox_ecotox_overlap/envirotox_ecotox_overlap.xlsx",
  "outputs/database_audits/envirotox/03_prepare_envirotox_manual_review/envirotox_manual_review.xlsx"
)

forbidden_tracked <- intersect(
  forbidden_exact,
  git_files
)

obsolete_name_patterns <- c(
  "multivariate",
  "simulation",
  "pilot",
  "common_landscape_plot",
  "final_methods_qa",
  "_updated\\.r$",
  "10_final_figure_suite_polished",
  "11_metal_covariance_full\\.r$",
  "11_metal_covariance\\.r$",
  "/sex\\.r$",
  "00_data_preparation_and_harmonization\\.r$",
  "00_wos_harmonize_and_merge\\.r$"
)

obsolete_files <- git_files[
  vapply(
    git_files,
    function(x) {
      any(
        vapply(
          obsolete_name_patterns,
          function(p) {
            grepl(
              p,
              tolower(x),
              perl = TRUE
            )
          },
          logical(1)
        )
      )
    },
    logical(1)
  )
]


# ------------------------------------------------------------
# 4. TEXT-CONTENT CHECKS
# ------------------------------------------------------------

text_extensions <- c(
  "R",
  "r",
  "md",
  "cff",
  "txt",
  "gitignore",
  "yml",
  "yaml"
)

is_text_candidate <- function(path) {

  if (
    basename(path) == ".gitignore"
  ) {
    return(TRUE)
  }

  ext <- tools::file_ext(path)

  ext %in% text_extensions
}

text_files <- git_files[
  file.exists(git_files) &
    vapply(
      git_files,
      is_text_candidate,
      logical(1)
    )
]

# Do not scan this QA script itself. It contains the literal search
# patterns and would otherwise generate false-positive hits.
text_files <- setdiff(
  text_files,
  "scripts/check_repository_before_release.R"
)

read_text_safely <- function(path) {

  x <- tryCatch(
    readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8"
    ),
    error = function(e) character()
  )

  paste(
    x,
    collapse = "\n"
  )
}

text_content <- setNames(
  lapply(
    text_files,
    read_text_safely
  ),
  text_files
)

find_pattern <- function(pattern) {

  names(
    Filter(
      length,
      lapply(
        text_content,
        function(x) {
          grep(
            pattern,
            x,
            perl = TRUE,
            ignore.case = TRUE
          )
        }
      )
    )
  )
}

absolute_path_hits <- unique(
  c(
    find_pattern("C:/Users/"),
    find_pattern("C:\\\\Users\\\\"),
    find_pattern("OneDrive/")
  )
)

old_workflow_hits <- unique(
  c(
    find_pattern("ECOTOX_Thesis\\.Rproj"),
    find_pattern("00b_WoS_harmonize_and_merge"),
    find_pattern("combined_302"),
    find_pattern("submitted_vs_combined"),
    find_pattern("NEXT STEP:")
  )
)


# ------------------------------------------------------------
# 5. FILE-TYPE WARNINGS
# ------------------------------------------------------------

warn_docx <- git_files[
  grepl(
    "\\.docx$",
    git_files,
    ignore.case = TRUE
  )
]

warn_pdf_outside_figures <- git_files[
  grepl(
    "\\.pdf$",
    git_files,
    ignore.case = TRUE
  ) &
    !grepl(
      "^results/figures/",
      git_files
    )
]

warn_archives <- git_files[
  grepl(
    "\\.(zip|7z|rar)$",
    git_files,
    ignore.case = TRUE
  )
]


# ------------------------------------------------------------
# 6. CORE NUMERICAL OUTPUT CHECKS
# ------------------------------------------------------------

benchmark_messages <- character()

combined_file <- "data/processed/combined_ecotox_wos_harmonized.csv"

if (file.exists(combined_file)) {

  combined <- read.csv(
    combined_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  lc50 <- suppressWarnings(
    as.numeric(
      combined$LC50_umol_L
    )
  )

  quantitative_n <- sum(
    is.finite(lc50) &
      lc50 > 0
  )

  if (
    nrow(combined) != 353L ||
      quantitative_n != 304L
  ) {
    benchmark_messages <- c(
      benchmark_messages,
      paste0(
        "Combined dataset benchmark failed: ",
        nrow(combined),
        " rows / ",
        quantitative_n,
        " quantitative."
      )
    )
  }
}

benchmark_files <- c(
  "outputs/04_fit_pooled_metal_model/RUN_COMPLETE.txt",
  "outputs/11_fit_stage_model/RUN_COMPLETE.txt",
  "outputs/13_fit_order_model/RUN_COMPLETE.txt",
  "outputs/14_analyze_metal_covariance/RUN_COMPLETE.txt",
  "outputs/15_audit_sex_specific_evidence/RUN_COMPLETE.txt",
  "outputs/16_make_main_thesis_figures/RUN_COMPLETE.txt"
)

missing_benchmark_outputs <- benchmark_files[
  !file.exists(
    benchmark_files
  )
]


# ------------------------------------------------------------
# 7. REPORT
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "REPOSITORY PRE-RELEASE QA\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nGit candidate files: ",
  length(git_files),
  "\n",
  sep = ""
)

show_issue <- function(title, x) {

  cat(
    "\n",
    title,
    ": ",
    length(x),
    "\n",
    sep = ""
  )

  if (length(x)) {
    cat(
      paste0(
        " - ",
        x,
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )
  }
}

show_issue(
  "Missing required files",
  missing_required
)

show_issue(
  "Forbidden local/raw files visible to Git",
  forbidden_tracked
)

show_issue(
  "Obsolete analysis files visible to Git",
  obsolete_files
)

show_issue(
  "Absolute/local-path text hits",
  absolute_path_hits
)

show_issue(
  "Old workflow / NEXT STEP text hits",
  old_workflow_hits
)

show_issue(
  "DOCX files visible to Git (review)",
  warn_docx
)

show_issue(
  "PDF files outside results/figures (review)",
  warn_pdf_outside_figures
)

show_issue(
  "Archive files visible to Git (review)",
  warn_archives
)

show_issue(
  "Missing key completed-analysis outputs",
  missing_benchmark_outputs
)

show_issue(
  "Numerical benchmark issues",
  benchmark_messages
)


# ------------------------------------------------------------
# 8. PASS / REVIEW STATUS
# ------------------------------------------------------------

hard_fail <- unique(
  c(
    missing_required,
    forbidden_tracked,
    obsolete_files,
    absolute_path_hits,
    old_workflow_hits,
    missing_benchmark_outputs,
    benchmark_messages
  )
)

if (length(hard_fail)) {

  cat(
    "\nSTATUS: REPOSITORY QA = REVIEW REQUIRED\n"
  )

  stop(
    "Resolve the hard-fail items above before committing/pushing.",
    call. = FALSE
  )
}

cat(
  "\nSTATUS: REPOSITORY QA = PASS\n"
)

if (
  length(warn_docx) ||
    length(warn_pdf_outside_figures) ||
    length(warn_archives)
) {

  cat(
    "Non-blocking file-type warnings remain; review them before release.\n"
  )
} else {
  cat(
    "No non-blocking file-type warnings detected.\n"
  )
}
