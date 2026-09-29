# ============================================================
# 15_audit_sex_specific_evidence.R
# ============================================================
# Purpose
# Audit how explicitly sex is represented in the final combined
# acute-LC50 dataset and identify the limited records/References that
# can support descriptive male-female comparisons.
#
# This is an evidence-availability audit, not a sex-effects model.
# No LC50 pooling, hypothesis tests, or inferential sex-effect
# estimates are performed here.
#
# Scientific interpretation
# - "Not explicit" means that no explicit male/female wording was
#   found in the inspected fields. It does NOT mean the source paper
#   lacked sex information.
# - "Both sexes mentioned" does not establish a pooled or separately
#   analysed male/female design.
# - Adult, sexually mature, gestation, paper titles, and single-letter
#   codes are not used to infer sex.
# - Note-text mentions are review evidence only and never
#   automatically change the record classification.
#
# Inputs
# data/processed/combined_ecotox_wos_harmonized.csv
# data/curated/source_verification/combined_source_adjudication_rows.csv
# data/curated/source_verification/combined_source_verified_rows_used.csv
#
# Main outputs
# outputs/15_audit_sex_specific_evidence/
#
# Reproducibility rule
# Run from the repository root with copepod-metal-toxicity.Rproj open.
# The classification logic is preserved from the final thesis sex
# availability audit; paths, naming and frozen 304-record QA have
# been updated for the public repository.
# ============================================================


# ============================================================
# 0. INPUTS + OUTPUT DIRECTORY
# ============================================================

input_file <- file.path(
  "data", "processed",
  "combined_ecotox_wos_harmonized.csv"
)

ledger_files <- c(
  file.path(
    "data", "curated", "source_verification",
    "combined_source_adjudication_rows.csv"
  ),
  file.path(
    "data", "curated", "source_verification",
    "combined_source_verified_rows_used.csv"
  )
)

expected_quantitative_n <- 304L

required_inputs <- c(input_file, ledger_files)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s): ",
    paste(missing_inputs, collapse = " | ")
  )
}

output_root <- file.path(
  "outputs",
  "15_audit_sex_specific_evidence"
)

if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
}

completion_path <- file.path(output_root, "RUN_COMPLETE.txt")
if (file.exists(completion_path) && !file.remove(completion_path)) {
  stop(
    "Could not remove the previous completion marker. ",
    "Close open output files and retry."
  )
}


# ============================================================
# 1. HELPERS
# ============================================================

read_input <- function(path) {

  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    encoding = "UTF-8",
    na.strings = c("", "NA")
  )

  names(x) <- sub(
    intToUtf8(65279L),
    "",
    names(x),
    fixed = TRUE
  )

  for (name in names(x)) {

    if (is.character(x[[name]])) {

      converted <- iconv(
        x[[name]],
        from = "UTF-8",
        to = "UTF-8",
        sub = NA
      )

      if (any(!is.na(x[[name]]) & is.na(converted))) {
        stop(
          "Invalid UTF-8 in ",
          path,
          ", column ",
          name,
          ". Check source encoding."
        )
      }

      x[[name]] <- converted
    }
  }

  x
}


write_out <- function(x, filename) {

  write.csv(
    x,
    file.path(output_root, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}


clean_text <- function(x) {

  x <- tolower(
    trimws(
      enc2utf8(
        as.character(x)
      )
    )
  )

  x[is.na(x)] <- ""

  x
}


classify_sex_text <- function(x) {

  s <- clean_text(x)

  female <- grepl(
    "\\bfemales?\\b",
    s,
    perl = TRUE
  )

  male <- grepl(
    "\\bmales?\\b",
    s,
    perl = TRUE
  )

  mixed <- grepl(
    "mixed[- ]sex|both sexes|male[s]?\\s*[+/&]\\s*female",
    s,
    perl = TRUE
  )

  negative <- grepl(
    "unknown|unspecified|not (reported|specified|recorded)|undetermined",
    s,
    perl = TRUE
  )

  result <- rep(
    "Not explicit",
    length(s)
  )

  result[female & !male] <- "Female"
  result[male & !female] <- "Male"

  result[
    mixed |
      (female & male)
  ] <- "Both sexes mentioned"

  result[
    negative &
      (female | male | mixed)
  ] <- "Needs review"

  result
}


summarise_records <- function(x, keys) {

  if (!nrow(x)) {
    return(data.frame())
  }

  groups <- x[keys]

  groups[] <- lapply(
    groups,
    function(z) {
      z <- as.character(z)
      z[is.na(z) | z == ""] <- "Not recorded"
      z
    }
  )

  token <- do.call(
    paste,
    c(groups, sep = "\034")
  )

  pieces <- lapply(
    split(seq_len(nrow(x)), token),
    function(index) {

      out <- groups[index[1], , drop = FALSE]

      out$Records <- length(index)

      out$References <- length(
        unique(
          x$Reference_ID[index]
        )
      )

      out$Species_labels <- length(
        unique(
          na.omit(
            x$Species[index]
          )
        )
      )

      out
    }
  )

  result <- do.call(
    rbind,
    pieces
  )

  rownames(result) <- NULL

  result
}


# ============================================================
# 2. READ + VALIDATE FINAL COMBINED DATASET
# ============================================================

d <- read_input(input_file)

required_columns <- c(
  "Reference_ID",
  "Result_ID",
  "Species",
  "Metal",
  "LC50_umol_L"
)

if (!all(required_columns %in% names(d))) {
  stop(
    "Missing columns: ",
    paste(
      setdiff(
        required_columns,
        names(d)
      ),
      collapse = ", "
    )
  )
}

if (anyNA(d$Result_ID) || anyDuplicated(d$Result_ID)) {
  stop(
    "Missing or duplicated Result_ID. ",
    "Check the combined input before counting."
  )
}

if (anyNA(d$Reference_ID)) {
  stop(
    "Missing Reference_ID. ",
    "Check the combined input before counting."
  )
}

lc50_numeric <- suppressWarnings(
  as.numeric(
    d$LC50_umol_L
  )
)

d$Quantitative <- is.finite(lc50_numeric) &
  lc50_numeric > 0

q <- d[
  d$Quantitative,
  ,
  drop = FALSE
]

stopifnot(
  nrow(q) == expected_quantitative_n
)


# ============================================================
# 3. CONSERVATIVE SEX-FIELD CLASSIFICATION
# ============================================================

sex_fields <- intersect(
  c(
    "Sex",
    "sex",
    "Sex_Code",
    "Gender",
    "External_Stage_Sex",
    "Lifestage"
  ),
  names(d)
)

if (!length(sex_fields)) {
  stop(
    "No recognised sex/stage fields were found. ",
    "Inspect the combined-data columns."
  )
}

classes <- lapply(
  d[sex_fields],
  classify_sex_text
)

d$Sex_field_status <- vapply(
  seq_len(nrow(d)),
  function(i) {

    statuses <- unique(
      vapply(
        classes,
        function(x) x[i],
        character(1)
      )
    )

    statuses <- setdiff(
      statuses,
      "Not explicit"
    )

    if (!length(statuses)) {
      "Not explicit"
    } else if (length(statuses) == 1) {
      statuses
    } else {
      "Conflicting fields"
    }
  },
  character(1)
)

d$Sex_field_evidence <- vapply(
  seq_len(nrow(d)),
  function(i) {

    values <- vapply(
      sex_fields,
      function(field) {

        value <- d[[field]][i]

        if (
          is.na(value) ||
          !nzchar(
            trimws(
              as.character(value)
            )
          )
        ) {
          ""
        } else {
          paste0(
            field,
            "=",
            value
          )
        }
      },
      character(1)
    )

    paste(
      values[nzchar(values)],
      collapse = " | "
    )
  },
  character(1)
)

q <- d[
  d$Quantitative,
  ,
  drop = FALSE
]


# ============================================================
# 4. FROZEN FINAL-DATA QA
# ============================================================

status_counts <- table(
  factor(
    q$Sex_field_status,
    levels = c(
      "Not explicit",
      "Female",
      "Male",
      "Both sexes mentioned",
      "Needs review",
      "Conflicting fields"
    )
  )
)

stopifnot(
  unname(status_counts["Not explicit"]) == 279L,
  unname(status_counts["Female"]) == 14L,
  unname(status_counts["Male"]) == 10L,
  unname(status_counts["Both sexes mentioned"]) == 1L,
  unname(status_counts["Needs review"]) == 0L,
  unname(status_counts["Conflicting fields"]) == 0L
)


# ============================================================
# 5. CORE AUDIT OUTPUTS
# ============================================================

write_out(
  d,
  "01_record_audit.csv"
)

write_out(
  summarise_records(
    q,
    "Sex_field_status"
  ),
  "02_quantitative_sex_summary.csv"
)

if ("Source_Origin" %in% names(q)) {

  write_out(
    summarise_records(
      q,
      c(
        "Source_Origin",
        "Sex_field_status"
      )
    ),
    "03_by_source_route.csv"
  )
}

write_out(
  summarise_records(
    q,
    c(
      "Reference_ID",
      "Species",
      "Metal",
      "Sex_field_status"
    )
  ),
  "04_by_reference_species_metal.csv"
)


# ============================================================
# 6. EXPLICIT-SEX RECORDS + MALE/FEMALE SUPPORT
# ============================================================

explicit_sex <- q[
  q$Sex_field_status != "Not explicit",
  ,
  drop = FALSE
]

write_out(
  explicit_sex,
  "05_explicit_sex_records.csv"
)

male_female_rows <- q[
  q$Sex_field_status %in% c(
    "Male",
    "Female"
  ),
  ,
  drop = FALSE
]

support_keys <- c(
  "Reference_ID",
  "Species",
  "Metal"
)

support_groups <- split(
  male_female_rows,
  interaction(
    male_female_rows[support_keys],
    drop = TRUE,
    lex.order = TRUE
  )
)

male_female_support <- do.call(
  rbind,
  lapply(
    support_groups,
    function(x) {

      statuses <- unique(
        x$Sex_field_status
      )

      data.frame(
        Reference_ID = x$Reference_ID[1],
        Species = x$Species[1],
        Metal = x$Metal[1],
        Records = nrow(x),
        Duration_levels = paste(
          sort(
            unique(
              x$Duration_days[
                is.finite(
                  suppressWarnings(
                    as.numeric(
                      x$Duration_days
                    )
                  )
                )
              ]
            )
          ),
          collapse = "; "
        ),
        Has_Male = "Male" %in% statuses,
        Has_Female = "Female" %in% statuses,
        Supports_descriptive_male_female_comparison =
          all(
            c(
              "Male",
              "Female"
            ) %in% statuses
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(male_female_support) <- NULL

write_out(
  male_female_support,
  "06_male_female_support_by_reference_species_metal.csv"
)


# ============================================================
# 7. NOTE-TEXT REVIEW AUDIT
# ============================================================

review <- list()

scan_notes <- function(x, origin) {

  note_columns <- intersect(
    c(
      "External_Source_Note",
      "Decision_Note",
      "Expected_Stage_Sex",
      "Verification_Status",
      "Eligibility_Reason"
    ),
    names(x)
  )

  for (field in note_columns) {

    hit <- grepl(
      "\\b(sex|sexes|male|males|female|females)\\b",
      clean_text(
        x[[field]]
      ),
      perl = TRUE
    )

    if (any(hit)) {

      id_columns <- intersect(
        c(
          "Result_ID",
          "Reference_ID",
          "Species",
          "Metal",
          "Source_Context_ID"
        ),
        names(x)
      )

      y <- x[
        hit,
        id_columns,
        drop = FALSE
      ]

      y$Evidence_file <- origin
      y$Evidence_field <- field
      y$Evidence_text <- x[[field]][hit]

      y$In_current_input <- if (
        "Result_ID" %in% names(y)
      ) {
        y$Result_ID %in% d$Result_ID
      } else {
        NA
      }

      y$In_quantitative_input <- if (
        "Result_ID" %in% names(y)
      ) {
        y$Result_ID %in% q$Result_ID
      } else {
        NA
      }

      review[[length(review) + 1L]] <<- y
    }
  }
}


scan_notes(
  d,
  input_file
)

for (path in ledger_files) {

  scan_notes(
    read_input(path),
    path
  )
}

if (length(review)) {

  all_columns <- unique(
    unlist(
      lapply(
        review,
        names
      )
    )
  )

  review <- lapply(
    review,
    function(x) {

      for (
        field in setdiff(
          all_columns,
          names(x)
        )
      ) {
        x[[field]] <- NA
      }

      x[all_columns]
    }
  )

  write_out(
    unique(
      do.call(
        rbind,
        review
      )
    ),
    "07_note_mentions_for_review.csv"
  )

} else {

  write_out(
    data.frame(
      Message =
        "No sex-related note mentions found in available inputs."
    ),
    "07_note_mentions_for_review.csv"
  )
}


# ============================================================
# 8. INPUT MANIFEST + SUMMARY
# ============================================================

manifest <- data.frame(
  Role = c(
    "Combined harmonized dataset",
    "Source-adjudication ledger",
    "Source-verified ledger"
  ),
  File = required_inputs,
  Exists = file.exists(required_inputs),
  stringsAsFactors = FALSE
)

manifest$MD5 <- vapply(
  manifest$File,
  function(path) {
    unname(
      tools::md5sum(path)
    )
  },
  character(1)
)

write_out(
  manifest,
  "00_input_manifest.csv"
)

explicit_refs <- unique(
  explicit_sex$Reference_ID
)

comparison_rows <- male_female_support[
  male_female_support$Supports_descriptive_male_female_comparison,
  ,
  drop = FALSE
]

report <- c(
  "STATUS: 15 AUDIT SEX-SPECIFIC EVIDENCE = PASS",
  paste(
    "Combined input rows:",
    nrow(d)
  ),
  paste(
    "Positive molar LC50 rows:",
    nrow(q)
  ),
  paste(
    "Quantitative References:",
    length(
      unique(
        q$Reference_ID
      )
    )
  ),
  paste(
    "Fields inspected:",
    paste(
      sex_fields,
      collapse = ", "
    )
  ),
  paste(
    "Explicit-sex quantitative rows:",
    nrow(explicit_sex)
  ),
  paste(
    "References with any explicit-sex quantitative row:",
    length(explicit_refs)
  ),
  paste(
    "Reference × Species × Metal groups containing both explicit male and female rows:",
    sum(
      comparison_rows$Supports_descriptive_male_female_comparison
    )
  ),
  "",
  capture.output(
    print(
      summarise_records(
        q,
        "Sex_field_status"
      ),
      row.names = FALSE
    )
  ),
  "",
  "Interpretation:",
  "Not explicit means no explicit sex in inspected fields; it does not mean the paper lacks sex information.",
  "Both sexes mentioned does not establish pooled or separate male/female experiments.",
  "Explicit male/female labels establish record availability, not a general sex effect.",
  "Note mentions are review candidates and are not automatically verified sex labels.",
  "Reference/species counts overlap across categories and should not be added.",
  "No duration restriction, model, significance test, or pooled sex-effect estimate was applied.",
  "The limited comparative evidence is retained as complementary biological evidence rather than a separate sex-effects model."
)

writeLines(
  report,
  file.path(
    output_root,
    "08_summary.txt"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_root,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "STATUS: 15 AUDIT SEX-SPECIFIC EVIDENCE = PASS",
    paste0(
      "Quantitative Results: ",
      nrow(q)
    ),
    paste0(
      "Explicit-sex quantitative Results: ",
      nrow(explicit_sex)
    ),
    paste0(
      "References with explicit-sex records: ",
      length(explicit_refs)
    ),
    "No sex-effects model was fitted."
  ),
  file.path(
    output_root,
    "RUN_COMPLETE.txt"
  )
)


# ============================================================
# 9. FINAL CONSOLE SUMMARY
# ============================================================

cat(
  paste(
    report,
    collapse = "\n"
  ),
  "\n\n",
  sep = ""
)

cat(
  "Output folder: ",
  output_root,
  "\n",
  sep = ""
)

