# ============================================================
# 02_check_envirotox_ecotox_overlap.R
# ============================================================
# Purpose
# Compare the 54 Saltwater EnviroTox copepod LC50 candidates with
# the ECOTOX copepod screening frame at two levels:
#   1) source/reference overlap;
#   2) Result overlap using species + metal + duration + LC50.
#
# Final thesis benchmark
# - 54 Saltwater EnviroTox candidates
# - 38 exact Result-level ECOTOX overlaps
# - 16 unresolved candidates retained for source review:
#     12 source + species + metal matches without matching duration
#      2 same source/species/metal/duration but different LC50
#      1 source-only match
#      1 source not matched in the ECOTOX screening frame
#
# Inputs
# outputs/database_audits/envirotox/01_screen_envirotox_records/
#   envirotox_preliminary_screening.xlsx
# data/curated/ecotox/ecotox_working_master_harmonized.xlsx
# data/raw/ecotox/ecotox_raw_export.xlsx
#
# Output
# outputs/database_audits/envirotox/02_check_envirotox_ecotox_overlap/
#
# Important interpretation
# A shared source is not automatically a duplicate Result. Non-exact
# matches are retained for original-source review rather than treated
# as additional evidence.
#
# Scientific rule
# Source normalization, bibliographic matching, Result-level matching,
# unit conversion and LC50 numerical-equality logic are preserved from
# the final thesis workflow. Repository cleanup changes paths, naming,
# documentation and QA only.
# ============================================================


# ---- 1. Packages -------------------------------------------------------------
required_packages <- c("readxl", "dplyr", "stringr", "tidyr", "purrr", "tibble", "openxlsx", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(openxlsx)
  library(readr)
})

# ---- 2. Repository file paths -----------------------------------------------

envirotox_file <- file.path(
  "outputs", "database_audits", "envirotox",
  "01_screen_envirotox_records",
  "envirotox_preliminary_screening.xlsx"
)

ecotox_master_file <- file.path(
  "data", "curated", "ecotox",
  "ecotox_working_master_harmonized.xlsx"
)

ecotox_raw_file <- file.path(
  "data", "raw", "ecotox",
  "ecotox_raw_export.xlsx"
)

output_dir <- file.path(
  "outputs", "database_audits", "envirotox",
  "02_check_envirotox_ecotox_overlap"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  output_dir,
  "envirotox_ecotox_overlap.xlsx"
)

completion_path <- file.path(
  output_dir,
  "RUN_COMPLETE.txt"
)

if (
  file.exists(completion_path) &&
  !file.remove(completion_path)
) {
  stop(
    "Could not remove the previous completion marker. ",
    "Close open output files and retry."
  )
}

required_inputs <- c(
  envirotox_file,
  ecotox_master_file,
  ecotox_raw_file
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s): ",
    paste(
      missing_inputs,
      collapse = " | "
    )
  )
}


# ---- 3. Read and validate inputs --------------------------------------------
required_sheets <- list(
  envirotox = "Saltwater_Candidates",
  master = c("Screening_410", "Harmonized_295"),
  raw = "References"
)

check_sheets <- function(path, expected, label) {
  available <- excel_sheets(path)
  absent <- setdiff(expected, available)
  if (length(absent) > 0) {
    stop(
      label, " is missing required sheet(s): ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
}

check_sheets(envirotox_file, required_sheets$envirotox, "EnviroTox workbook")
check_sheets(ecotox_master_file, required_sheets$master, "ECOTOX master workbook")
check_sheets(ecotox_raw_file, required_sheets$raw, "ECOTOX raw workbook")

env <- read_excel(envirotox_file, sheet = "Saltwater_Candidates")
eco_screen <- read_excel(ecotox_master_file, sheet = "Screening_410")
eco_harmonized <- read_excel(ecotox_master_file, sheet = "Harmonized_295")
eco_refs <- read_excel(ecotox_raw_file, sheet = "References")

required_env_cols <- c(
  "EnviroTox_Row", "Metal", "Latin name", "Effect_Value_mg_L",
  "Duration_Days_Numeric", "Source"
)
required_screen_cols <- c(
  "Reference_Number", "Result_Number", "Test_Number", "Metal", "Species",
  "Duration_days", "Analysis_Point", "Analysis_Units", "Eligibility_Status",
  "Eligibility_Reason", "Title"
)
required_ref_cols <- c("Ref. Number", "Author", "Title", "Pub. Year", "Citation")

assert_columns <- function(data, required, label) {
  absent <- setdiff(required, names(data))
  if (length(absent) > 0) {
    stop(label, " is missing column(s): ", paste(absent, collapse = ", "), call. = FALSE)
  }
}

assert_columns(env, required_env_cols, "Saltwater_Candidates")
assert_columns(eco_screen, required_screen_cols, "Screening_410")
assert_columns(eco_refs, required_ref_cols, "References")

if (nrow(env) == 0) stop("Saltwater_Candidates contains no records.", call. = FALSE)

# ---- 4. Normalization helpers ------------------------------------------------
normalize_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- ""
  x <- str_to_lower(x)
  x <- str_replace_all(x, "doi\\s*:?\\s*[^ ]+", " ")
  x <- str_replace_all(x, "ecoref\\s*#?\\s*:?\\s*[0-9]+", " ")
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_squish(x)
}

edit_similarity_one <- function(a, b) {
  if (is.na(a) || is.na(b) || a == "" || b == "") return(0)
  denominator <- max(nchar(a), nchar(b))
  if (denominator == 0) return(0)
  max(0, 1 - as.numeric(adist(a, b)) / denominator)
}

contains_long_title <- function(source_norm, title_norm) {
  if (is.na(source_norm) || is.na(title_norm) ||
      source_norm == "" || title_norm == "" || nchar(title_norm) < 20) {
    return(FALSE)
  }
  str_detect(source_norm, fixed(title_norm)) || str_detect(title_norm, fixed(source_norm))
}

to_mg_L <- function(value, unit) {
  value <- suppressWarnings(as.numeric(value))
  unit_key <- normalize_text(unit)
  unit_key <- str_replace_all(unit_key, " ", "")
  
  case_when(
    is.na(value) ~ NA_real_,
    unit_key %in% c("mg/l", "mgl", "ppm") ~ value,
    unit_key %in% c("ug/l", "ugl", "microg/l", "microgl", "ppb") ~ value / 1000,
    unit_key %in% c("ng/l", "ngl") ~ value / 1000000,
    unit_key %in% c("g/l", "gl") ~ value * 1000,
    TRUE ~ NA_real_
  )
}

relative_difference <- function(a, b) {
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))
  if (is.na(a) || is.na(b)) return(NA_real_)
  denominator <- max(abs(a), abs(b))
  if (denominator == 0) return(0)
  abs(a - b) / denominator
}

# ---- 5. Prepare ECOTOX reference frame --------------------------------------
screen_reference_numbers <- eco_screen %>%
  distinct(Reference_Number) %>%
  pull(Reference_Number)

refs <- eco_refs %>%
  transmute(
    Reference_Number = as.numeric(`Ref. Number`),
    ECOTOX_Author = as.character(Author),
    ECOTOX_Title = as.character(Title),
    ECOTOX_Year = suppressWarnings(as.numeric(`Pub. Year`)),
    ECOTOX_Citation = as.character(Citation)
  ) %>%
  filter(Reference_Number %in% screen_reference_numbers) %>%
  distinct(Reference_Number, .keep_all = TRUE) %>%
  mutate(
    title_norm = normalize_text(ECOTOX_Title),
    citation_norm = normalize_text(ECOTOX_Citation)
  )

env_sources <- env %>%
  count(Source, name = "EnviroTox_Candidate_Records") %>%
  mutate(source_norm = normalize_text(Source))

# ---- 6. Source-level matching ------------------------------------------------
# Every EnviroTox source is compared with every ECOTOX reference represented in
# Screening_410. The full source-citation similarity and title containment are
# retained separately so the match remains auditable.
source_pairs <- tidyr::crossing(
  env_sources %>% select(Source, EnviroTox_Candidate_Records, source_norm),
  refs
) %>%
  rowwise() %>%
  mutate(
    Title_Similarity = edit_similarity_one(source_norm, title_norm),
    Citation_Similarity = edit_similarity_one(source_norm, citation_norm),
    Title_Contained = contains_long_title(source_norm, title_norm),
    Source_Match_Score = max(
      Title_Similarity,
      Citation_Similarity,
      if_else(Title_Contained, 0.97, 0),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

source_ranked <- source_pairs %>%
  group_by(Source) %>%
  arrange(desc(Source_Match_Score), Reference_Number, .by_group = TRUE) %>%
  mutate(Source_Match_Rank = row_number()) %>%
  ungroup()

source_best <- source_ranked %>%
  filter(Source_Match_Rank <= 2) %>%
  select(
    Source, Source_Match_Rank, Reference_Number, ECOTOX_Author, ECOTOX_Title,
    ECOTOX_Year, Title_Similarity, Citation_Similarity, Title_Contained,
    Source_Match_Score
  ) %>%
  pivot_wider(
    names_from = Source_Match_Rank,
    values_from = c(
      Reference_Number, ECOTOX_Author, ECOTOX_Title, ECOTOX_Year,
      Title_Similarity, Citation_Similarity, Title_Contained, Source_Match_Score
    ),
    names_glue = "{.value}_{Source_Match_Rank}"
  ) %>%
  left_join(env_sources %>% select(Source, EnviroTox_Candidate_Records), by = "Source") %>%
  mutate(
    Match_Margin = Source_Match_Score_1 - coalesce(Source_Match_Score_2, 0),
    Source_Overlap_Status = case_when(
      Source_Match_Score_1 >= 0.85 &
        (Title_Contained_1 | Match_Margin >= 0.02 | Citation_Similarity_1 >= 0.95) ~ "MATCHED_AUTO",
      Source_Match_Score_1 >= 0.72 ~ "POSSIBLE_MATCH_REVIEW",
      TRUE ~ "NOT_FOUND_IN_ECOTOX_FRAME"
    ),
    Matched_ECOTOX_Reference = if_else(
      Source_Overlap_Status == "MATCHED_AUTO",
      Reference_Number_1,
      NA_real_
    ),
    Source_Verification_Note = case_when(
      Source_Overlap_Status == "MATCHED_AUTO" ~
        "High-confidence bibliographic match; result-level comparison still required.",
      Source_Overlap_Status == "POSSIBLE_MATCH_REVIEW" ~
        "Possible bibliographic match; verify source manually before using the suggested reference.",
      TRUE ~
        "No reliable source match was found within the ECOTOX copepod screening frame."
    )
  ) %>%
  arrange(Source_Overlap_Status, Source)

# Top three candidates are included as an audit sheet.
source_top_candidates <- source_ranked %>%
  filter(Source_Match_Rank <= 3) %>%
  select(
    Source, EnviroTox_Candidate_Records, Source_Match_Rank,
    Reference_Number, ECOTOX_Author, ECOTOX_Title, ECOTOX_Year,
    Title_Similarity, Citation_Similarity, Title_Contained, Source_Match_Score
  ) %>%
  arrange(Source, Source_Match_Rank)

# ---- 7. Prepare ECOTOX result frame -----------------------------------------
harmonized_status <- eco_harmonized %>%
  transmute(
    Result_Number = as.numeric(Result_Number),
    In_Harmonized_295 = TRUE,
    Harmonization_Status = as.character(Harmonization_Status),
    LC50_target_mg_L_harmonized = suppressWarnings(as.numeric(LC50_target_mg_L))
  ) %>%
  distinct(Result_Number, .keep_all = TRUE)

eco_results <- eco_screen %>%
  transmute(
    ECOTOX_Reference_Number = as.numeric(Reference_Number),
    ECOTOX_Test_Number = as.numeric(Test_Number),
    ECOTOX_Result_Number = as.numeric(Result_Number),
    ECOTOX_Metal = as.character(Metal),
    ECOTOX_Species = as.character(Species),
    ECOTOX_Duration_Days = suppressWarnings(as.numeric(Duration_days)),
    ECOTOX_Analysis_Point = suppressWarnings(as.numeric(Analysis_Point)),
    ECOTOX_Analysis_Units = as.character(Analysis_Units),
    ECOTOX_LC50_mg_L = map2_dbl(Analysis_Point, Analysis_Units, to_mg_L),
    ECOTOX_Eligibility_Status = as.character(Eligibility_Status),
    ECOTOX_Eligibility_Reason = as.character(Eligibility_Reason),
    ECOTOX_Title = as.character(Title),
    species_norm = normalize_text(Species)
  ) %>%
  left_join(harmonized_status, by = c("ECOTOX_Result_Number" = "Result_Number")) %>%
  mutate(
    In_Harmonized_295 = coalesce(In_Harmonized_295, FALSE),
    ECOTOX_LC50_mg_L = coalesce(LC50_target_mg_L_harmonized, ECOTOX_LC50_mg_L)
  )

source_lookup <- source_best %>%
  select(
    Source, Source_Overlap_Status, Matched_ECOTOX_Reference,
    Best_Source_Match_Score = Source_Match_Score_1,
    Source_Match_Margin = Match_Margin,
    Best_ECOTOX_Title = ECOTOX_Title_1
  )

env_prepared <- env %>%
  left_join(source_lookup, by = "Source") %>%
  mutate(
    species_norm = normalize_text(`Latin name`),
    EnviroTox_LC50_mg_L = suppressWarnings(as.numeric(Effect_Value_mg_L)),
    EnviroTox_Duration_Days = suppressWarnings(as.numeric(Duration_Days_Numeric))
  )

# ---- 8. Result-level matching ------------------------------------------------
# Numerical equality uses the documented base-R all.equal default tolerance.
# Explicit scale keeps this scalar comparison relative at every LC50 magnitude.
# This allows computational representation differences, not reporting-rounding
# differences. Unequal or missing values require source review.
# Reference: https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html
lc50_equal <- function(a, b) {
  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(FALSE)
  isTRUE(all.equal(a, b, scale = max(abs(a), abs(b)),
                   check.attributes = FALSE))
}

match_one_record <- function(env_row) {
  source_status <- env_row$Source_Overlap_Status[[1]]
  ref_number <- env_row$Matched_ECOTOX_Reference[[1]]
  
  empty_answer <- tibble(
    ECOTOX_Record_Overlap_Status = "SOURCE_NOT_MATCHED",
    ECOTOX_Reference_Number = NA_real_,
    ECOTOX_Test_Number = NA_real_,
    ECOTOX_Result_Number = NA_real_,
    ECOTOX_Species = NA_character_,
    ECOTOX_Metal = NA_character_,
    ECOTOX_Duration_Days = NA_real_,
    ECOTOX_LC50_mg_L = NA_real_,
    Relative_Value_Difference = NA_real_,
    ECOTOX_Eligibility_Status = NA_character_,
    ECOTOX_Eligibility_Reason = NA_character_,
    In_Harmonized_295 = FALSE,
    Harmonization_Status = NA_character_,
    Candidate_Context_Count = 0L,
    Record_Interpretation = "No automatic ECOTOX source match; inspect as a potential additional source."
  )
  
  if (is.na(ref_number) || source_status != "MATCHED_AUTO") return(empty_answer)
  
  same_source <- eco_results %>%
    filter(ECOTOX_Reference_Number == ref_number)
  
  same_species_metal <- same_source %>%
    filter(
      species_norm == env_row$species_norm[[1]],
      ECOTOX_Metal == env_row$Metal[[1]]
    )
  
  same_context <- same_species_metal %>%
    filter(
      !is.na(ECOTOX_Duration_Days),
      !is.na(env_row$EnviroTox_Duration_Days[[1]]),
      abs(ECOTOX_Duration_Days - env_row$EnviroTox_Duration_Days[[1]]) <= 0.001
    ) %>%
    mutate(
      Relative_Value_Difference = map2_dbl(
        ECOTOX_LC50_mg_L,
        env_row$EnviroTox_LC50_mg_L[[1]],
        relative_difference
      ),
      sort_difference = if_else(is.na(Relative_Value_Difference), Inf, Relative_Value_Difference)
    ) %>%
    arrange(sort_difference, ECOTOX_Result_Number)
  
  if (nrow(same_context) > 0) {
    best <- same_context %>% slice(1)
    difference <- best$Relative_Value_Difference[[1]]
    
    status <- case_when(
      lc50_equal(best$ECOTOX_LC50_mg_L[[1]], env_row$EnviroTox_LC50_mg_L[[1]]) ~ "EXACT_RECORD_OVERLAP",
      TRUE ~ "SAME_CONTEXT_VALUE_DIFF_REVIEW"
    )
    
    interpretation <- case_when(
      status == "EXACT_RECORD_OVERLAP" ~
        "Agreement in matched source, species, metal, duration and numerically equal LC50. This label describes field agreement; it does not resolve multiple experimental contexts.",
      TRUE ~
        "Same source/species/metal/duration but LC50 values are numerically unequal or unavailable; verify in the original source."
    )
    
    return(best %>%
             transmute(
               ECOTOX_Record_Overlap_Status = status,
               ECOTOX_Reference_Number,
               ECOTOX_Test_Number,
               ECOTOX_Result_Number,
               ECOTOX_Species,
               ECOTOX_Metal,
               ECOTOX_Duration_Days,
               ECOTOX_LC50_mg_L,
               Relative_Value_Difference,
               ECOTOX_Eligibility_Status,
               ECOTOX_Eligibility_Reason,
               In_Harmonized_295,
               Harmonization_Status,
               Candidate_Context_Count = nrow(same_context),
               Record_Interpretation = interpretation
             ))
  }
  
  if (nrow(same_species_metal) > 0) {
    best <- same_species_metal %>%
      arrange(ECOTOX_Duration_Days, ECOTOX_Result_Number) %>%
      slice(1)
    
    return(best %>%
             transmute(
               ECOTOX_Record_Overlap_Status = "SOURCE_SPECIES_METAL_ONLY_REVIEW",
               ECOTOX_Reference_Number,
               ECOTOX_Test_Number,
               ECOTOX_Result_Number,
               ECOTOX_Species,
               ECOTOX_Metal,
               ECOTOX_Duration_Days,
               ECOTOX_LC50_mg_L,
               Relative_Value_Difference = NA_real_,
               ECOTOX_Eligibility_Status,
               ECOTOX_Eligibility_Reason,
               In_Harmonized_295,
               Harmonization_Status,
               Candidate_Context_Count = nrow(same_species_metal),
               Record_Interpretation =
                 "Source, species and metal occur in ECOTOX, but the duration does not align; inspect the source manually."
             ))
  }
  
  tibble(
    ECOTOX_Record_Overlap_Status = "SOURCE_ONLY_REVIEW",
    ECOTOX_Reference_Number = ref_number,
    ECOTOX_Test_Number = NA_real_,
    ECOTOX_Result_Number = NA_real_,
    ECOTOX_Species = NA_character_,
    ECOTOX_Metal = NA_character_,
    ECOTOX_Duration_Days = NA_real_,
    ECOTOX_LC50_mg_L = NA_real_,
    Relative_Value_Difference = NA_real_,
    ECOTOX_Eligibility_Status = NA_character_,
    ECOTOX_Eligibility_Reason = NA_character_,
    In_Harmonized_295 = FALSE,
    Harmonization_Status = NA_character_,
    Candidate_Context_Count = nrow(same_source),
    Record_Interpretation =
      "The source exists in ECOTOX, but this species-metal result was not located in the copepod screening frame."
  )
}

record_matches <- map_dfr(
  seq_len(nrow(env_prepared)),
  ~match_one_record(env_prepared[.x, , drop = FALSE])
)

record_overlap <- bind_cols(
  env_prepared %>%
    select(-species_norm) %>%
    mutate(.row_order = row_number()),
  record_matches
) %>%
  mutate(
    ECOTOX_Overlap_Status = ECOTOX_Record_Overlap_Status,
    Thesis_Dataset_Status = case_when(
      ECOTOX_Record_Overlap_Status == "EXACT_RECORD_OVERLAP" &
        In_Harmonized_295 ~ "ALREADY_IN_HARMONIZED_295",
      ECOTOX_Record_Overlap_Status == "EXACT_RECORD_OVERLAP" &
        !In_Harmonized_295 ~ "IN_ECOTOX_BUT_NOT_IN_HARMONIZED_295",
      TRUE ~ "NOT_CONFIRMED_AT_RESULT_LEVEL"
    ),
    Manual_Review_Priority = case_when(
      ECOTOX_Record_Overlap_Status == "SOURCE_NOT_MATCHED" ~ 1L,
      ECOTOX_Record_Overlap_Status == "SOURCE_ONLY_REVIEW" ~ 2L,
      ECOTOX_Record_Overlap_Status == "SAME_CONTEXT_VALUE_DIFF_REVIEW" ~ 3L,
      ECOTOX_Record_Overlap_Status == "SOURCE_SPECIES_METAL_ONLY_REVIEW" ~ 4L,
      TRUE ~ 9L
    )
  ) %>%
  arrange(.row_order) %>%
  select(-.row_order)

# ---- 9. Audit summaries ------------------------------------------------------
record_status_summary <- record_overlap %>%
  count(ECOTOX_Record_Overlap_Status, name = "EnviroTox_Records") %>%
  mutate(
    Percent_of_54 = round(100 * EnviroTox_Records / nrow(record_overlap), 1),
    Interpretation = case_when(
      ECOTOX_Record_Overlap_Status == "EXACT_RECORD_OVERLAP" ~
        "Matching fields and numerically equal LC50; retain the existing result once provenance is established.",
      ECOTOX_Record_Overlap_Status == "SAME_CONTEXT_VALUE_DIFF_REVIEW" ~
        "Same study context but discordant value; inspect original paper/report.",
      ECOTOX_Record_Overlap_Status == "SOURCE_SPECIES_METAL_ONLY_REVIEW" ~
        "Partial overlap only; duration/result provenance requires checking.",
      ECOTOX_Record_Overlap_Status == "SOURCE_ONLY_REVIEW" ~
        "Study is known to ECOTOX, but this result may be additional.",
      TRUE ~ "Potentially additional source; verify bibliographically."
    )
  )

source_status_summary <- source_best %>%
  count(Source_Overlap_Status, name = "EnviroTox_Sources") %>%
  mutate(Percent_of_23 = round(100 * EnviroTox_Sources / nrow(source_best), 1))

thesis_status_summary <- record_overlap %>%
  count(Thesis_Dataset_Status, name = "EnviroTox_Records") %>%
  mutate(Percent_of_54 = round(100 * EnviroTox_Records / nrow(record_overlap), 1))

summary_table <- bind_rows(
  source_status_summary %>%
    transmute(Level = "SOURCE", Status = Source_Overlap_Status,
              Count = EnviroTox_Sources, Percent = Percent_of_23),
  record_status_summary %>%
    transmute(Level = "RESULT", Status = ECOTOX_Record_Overlap_Status,
              Count = EnviroTox_Records, Percent = Percent_of_54),
  thesis_status_summary %>%
    transmute(Level = "THESIS_DATASET", Status = Thesis_Dataset_Status,
              Count = EnviroTox_Records, Percent = Percent_of_54)
)

exact_overlaps <- record_overlap %>%
  filter(ECOTOX_Record_Overlap_Status == "EXACT_RECORD_OVERLAP")

manual_review <- record_overlap %>%
  filter(ECOTOX_Record_Overlap_Status != "EXACT_RECORD_OVERLAP") %>%
  arrange(Manual_Review_Priority, Source, Metal, `Latin name`)

potential_new <- record_overlap %>%
  filter(ECOTOX_Record_Overlap_Status %in% c("SOURCE_NOT_MATCHED", "SOURCE_ONLY_REVIEW")) %>%
  arrange(Manual_Review_Priority, Source, Metal, `Latin name`)

# Guardrails: these checks stop accidental row loss or duplication.
stopifnot(nrow(record_overlap) == nrow(env))
stopifnot(sum(record_status_summary$EnviroTox_Records) == nrow(env))
stopifnot(n_distinct(record_overlap$EnviroTox_Row) == nrow(record_overlap))


status_lookup <- setNames(
  record_status_summary$EnviroTox_Records,
  record_status_summary$ECOTOX_Record_Overlap_Status
)

count_status <- function(name) {
  if (name %in% names(status_lookup)) {
    as.integer(status_lookup[[name]])
  } else {
    0L
  }
}

# Frozen final thesis benchmark.
stopifnot(
  nrow(record_overlap) == 54L,
  n_distinct(record_overlap$Source) == 23L,
  count_status("EXACT_RECORD_OVERLAP") == 38L,
  count_status("SOURCE_SPECIES_METAL_ONLY_REVIEW") == 12L,
  count_status("SAME_CONTEXT_VALUE_DIFF_REVIEW") == 2L,
  count_status("SOURCE_ONLY_REVIEW") == 1L,
  count_status("SOURCE_NOT_MATCHED") == 1L,
  nrow(manual_review) == 16L
)

# ---- 10. README --------------------------------------------------------------
readme <- tribble(
  ~Item, ~Value,
  "Purpose", "Audit overlap between EnviroTox saltwater copepod LC50 candidates and the ECOTOX copepod screening frame.",
  "EnviroTox input", basename(envirotox_file),
  "ECOTOX master input", basename(ecotox_master_file),
  "ECOTOX raw input", basename(ecotox_raw_file),
  "EnviroTox candidate records", as.character(nrow(env)),
  "Unique EnviroTox sources", as.character(n_distinct(env$Source)),
  "Source rule", "High-confidence bibliographic matching uses normalized citation/title similarity plus long-title containment.",
  "Exact result rule", "Same ECOTOX reference, normalized species, metal and duration, with numerically equal positive LC50 values (base-R all.equal default tolerance; scale = larger absolute value).",
  "Unequal value rule", "Any LC50 difference beyond numerical equality requires source review; no percentage-based probable-duplicate class is used.",
  "Important caution", "SOURCE_ONLY and partial-context categories are not automatically new records; verify the original source first.",
  "Harmonized meaning", "In_Harmonized_295 indicates whether the matched ECOTOX Result Number survives in the harmonized master sheet.",
  "Interpretive boundary", "Manual_Review records require original-source verification; overlap status alone does not establish final eligibility.",
  "Created", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

# ---- 11. Write formatted workbook -------------------------------------------
wb <- createWorkbook(creator = "ESS Thesis reproducibility workflow")

header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE,
  border = "Bottom", borderColour = "#FFFFFF"
)
subheader_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#5B9BD5", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
exact_style <- createStyle(fgFill = "#E2F0D9")
review_style <- createStyle(fgFill = "#FFF2CC")
potential_style <- createStyle(fgFill = "#FCE4D6")

write_sheet <- function(wb, sheet_name, data, filter = TRUE) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = header_style, withFilter = filter)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(data)), widths = "auto")
  wide_cols <- which(vapply(data, function(x) {
    is.character(x) && any(nchar(replace(x, is.na(x), "")) > 60)
  }, logical(1)))
  if (length(wide_cols) > 0) {
    setColWidths(wb, sheet_name, cols = wide_cols, widths = 45)
    addStyle(
      wb, sheet_name, createStyle(wrapText = TRUE, valign = "top"),
      rows = 2:(nrow(data) + 1), cols = wide_cols, gridExpand = TRUE, stack = TRUE
    )
  }
  if (nrow(data) > 0) {
    addStyle(
      wb, sheet_name, createStyle(valign = "top"),
      rows = 2:(nrow(data) + 1), cols = seq_len(ncol(data)),
      gridExpand = TRUE, stack = TRUE
    )
  }
}

write_sheet(wb, "README", readme, filter = FALSE)
write_sheet(wb, "Overlap_Summary", summary_table)
write_sheet(wb, "Source_Overlap", source_best)
write_sheet(wb, "Record_Overlap", record_overlap)
write_sheet(wb, "Exact_Overlaps", exact_overlaps)
write_sheet(wb, "Manual_Review", manual_review)
write_sheet(wb, "Potential_New", potential_new)
write_sheet(wb, "Source_Top3_Audit", source_top_candidates)

# Colour-code the status column in the full record audit.
status_col <- which(names(record_overlap) == "ECOTOX_Record_Overlap_Status")
if (length(status_col) == 1 && nrow(record_overlap) > 0) {
  exact_rows <- which(record_overlap$ECOTOX_Record_Overlap_Status == "EXACT_RECORD_OVERLAP") + 1
  review_rows <- which(record_overlap$ECOTOX_Record_Overlap_Status %in% c(
    "SAME_CONTEXT_VALUE_DIFF_REVIEW",
    "SOURCE_SPECIES_METAL_ONLY_REVIEW"
  )) + 1
  potential_rows <- which(record_overlap$ECOTOX_Record_Overlap_Status %in% c(
    "SOURCE_ONLY_REVIEW", "SOURCE_NOT_MATCHED"
  )) + 1
  
  if (length(exact_rows) > 0) addStyle(wb, "Record_Overlap", exact_style, exact_rows, status_col, stack = TRUE)
  if (length(review_rows) > 0) addStyle(wb, "Record_Overlap", review_style, review_rows, status_col, stack = TRUE)
  if (length(potential_rows) > 0) addStyle(wb, "Record_Overlap", potential_style, potential_rows, status_col, stack = TRUE)
}

saveWorkbook(wb, output_file, overwrite = TRUE)


overlap_count_summary <- tibble::tibble(
  Category = c(
    "EnviroTox Saltwater candidates",
    "Exact ECOTOX Result overlaps",
    "Source/species/metal only; duration unresolved",
    "Same context; LC50 value differs",
    "Source-only review",
    "Source not matched",
    "Total unresolved candidates for source review"
  ),
  N = c(
    nrow(record_overlap),
    count_status("EXACT_RECORD_OVERLAP"),
    count_status("SOURCE_SPECIES_METAL_ONLY_REVIEW"),
    count_status("SAME_CONTEXT_VALUE_DIFF_REVIEW"),
    count_status("SOURCE_ONLY_REVIEW"),
    count_status("SOURCE_NOT_MATCHED"),
    nrow(manual_review)
  )
)

readr::write_csv(
  overlap_count_summary,
  file.path(
    output_dir,
    "overlap_count_summary.csv"
  )
)

matching_method <- tibble::tibble(
  Step = c(
    "Source/reference matching",
    "Result context matching",
    "LC50 equality",
    "Interpretation of non-exact matches"
  ),
  Rule = c(
    "Normalized citation/title similarity with long-title containment; high-confidence source matches proceed to Result-level checking.",
    "Within a matched ECOTOX Reference, require normalized species + metal + exposure duration.",
    "Positive mg/L LC50 values compared with base R all.equal using the larger absolute value as scale and default tolerance.",
    "Retain for original-source review; do not automatically classify as additional evidence."
  )
)

readr::write_csv(
  matching_method,
  file.path(
    output_dir,
    "matching_method.csv"
  )
)

input_output_manifest <- tibble::tibble(
  Role = c(
    "Input: EnviroTox preliminary-screening workbook",
    "Input: ECOTOX curated working master",
    "Input: ECOTOX raw workbook",
    "Local output: detailed EnviroTox-ECOTOX overlap workbook",
    "Output: overlap count summary",
    "Output: matching-method documentation"
  ),
  Path = c(
    envirotox_file,
    ecotox_master_file,
    ecotox_raw_file,
    output_file,
    file.path(
      output_dir,
      "overlap_count_summary.csv"
    ),
    file.path(
      output_dir,
      "matching_method.csv"
    )
  )
)

readr::write_csv(
  input_output_manifest,
  file.path(
    output_dir,
    "input_output_manifest.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "STATUS: ENVIROTOX-ECOTOX OVERLAP AUDIT = PASS",
    paste0(
      "Saltwater candidates checked: ",
      nrow(record_overlap)
    ),
    paste0(
      "Exact ECOTOX Result overlaps: ",
      count_status("EXACT_RECORD_OVERLAP")
    ),
    paste0(
      "Unresolved candidates for source review: ",
      nrow(manual_review)
    )
  ),
  file.path(
    output_dir,
    "RUN_COMPLETE.txt"
  )
)

# ---- 12. Console report ------------------------------------------------------
cat("\nEnviroTox-ECOTOX overlap audit completed successfully.\n")
cat("EnviroTox records checked:", nrow(record_overlap), "\n")
cat("EnviroTox sources checked:", nrow(source_best), "\n\n")
print(record_status_summary)
cat(
  "\nDetailed local audit workbook written to:\n",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("\nSTATUS: ENVIROTOX-ECOTOX OVERLAP AUDIT = PASS\n")
