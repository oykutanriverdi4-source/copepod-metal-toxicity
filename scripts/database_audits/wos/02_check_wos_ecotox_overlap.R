# ============================================================
# 02_check_wos_ecotox_overlap.R
# ============================================================
# Purpose
# Compare the 75 candidate Web of Science publications retained
# after preliminary screening with the ECOTOX reference sampling
# frame used in the thesis.
#
# Matching hierarchy
# 1. DOI exact match, when available
# 2. Normalized-title exact match
# 3. Fuzzy-title suggestion for manual verification only
#
# Final thesis benchmark
# - 75 WoS candidate publications
# - 28 publications represented in the ECOTOX retrieval
# - 47 publications not matched to the ECOTOX retrieval
#
# Important interpretation
# - Exact DOI/title matches establish sampling-frame overlap.
# - Fuzzy matches are suggestions only and are NEVER confirmed
#   automatically.
# - "Provisional unique WoS" means that no ECOTOX match was found
#   by this audit; it is not a final full-text eligibility decision.
# - This script does not extract LC50 values and does not make the
#   final 17-publication / 57-result source-verification decision.
#
# Inputs
# outputs/database_audits/wos/01_screen_wos_records/
#   wos_preliminary_screening.xlsx
# data/raw/ecotox/ecotox_raw_export.xlsx
# data/curated/ecotox/ecotox_working_master_harmonized.xlsx
#
# Local detailed output
# outputs/database_audits/wos/02_check_wos_ecotox_overlap/
#   wos_ecotox_overlap.xlsx
#
# The detailed workbook contains bibliographic metadata derived from
# the local WoS export and should remain local/gitignored.
# Compact count summaries and reproducibility metadata are also
# produced for repository documentation.
#
# Scientific rule
# Matching normalization, exact-match rules, fuzzy-similarity rules,
# and overlap classifications are preserved from the final thesis
# workflow. Repository cleanup changes paths, naming, duplicated
# setup code, documentation and frozen QA only.
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(writexl)
library(readr)


# ------------------------------------------------------------
# 2. Define repository file paths
# ------------------------------------------------------------

wos_screening_file <- file.path(
  "outputs",
  "database_audits",
  "wos",
  "01_screen_wos_records",
  "wos_preliminary_screening.xlsx"
)

ecotox_raw_file <- file.path(
  "data",
  "raw",
  "ecotox",
  "ecotox_raw_export.xlsx"
)

ecotox_master_file <- file.path(
  "data",
  "curated",
  "ecotox",
  "ecotox_working_master_harmonized.xlsx"
)

output_root <- file.path(
  "outputs",
  "database_audits",
  "wos",
  "02_check_wos_ecotox_overlap"
)

if (!dir.exists(output_root)) {
  dir.create(
    output_root,
    recursive = TRUE
  )
}

output_file <- file.path(
  output_root,
  "wos_ecotox_overlap.xlsx"
)

completion_path <- file.path(
  output_root,
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


# ------------------------------------------------------------
# 3. Check required inputs
# ------------------------------------------------------------

required_inputs <- c(
  wos_screening_file,
  ecotox_raw_file,
  ecotox_master_file
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


# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

get_optional_column <- function(
    data,
    candidates,
    default = ""
) {
  
  hit <- candidates[
    candidates %in% names(data)
  ]
  
  if (length(hit) == 0) {
    
    return(
      rep(
        default,
        nrow(data)
      )
    )
  }
  
  out <- as.character(
    data[[hit[[1]]]]
  )
  
  out[is.na(out)] <- default
  
  out
}


normalize_title <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  
  x <- iconv(
    x,
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  
  
  x[is.na(x)] <- ""
  
  
  x <- str_to_lower(x)
  
  
  x <- str_replace_all(
    x,
    "&",
    " and "
  )
  
  
  x <- str_replace_all(
    x,
    "[^a-z0-9]+",
    " "
  )
  
  
  x <- str_squish(x)
  
  
  x[
    x == ""
  ] <- NA_character_
  
  
  x
}


normalize_simple_text <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  
  x <- iconv(
    x,
    from = "",
    to = "ASCII//TRANSLIT",
    sub = ""
  )
  
  
  x[is.na(x)] <- ""
  
  
  x <- str_to_lower(x)
  
  
  x <- str_replace_all(
    x,
    "[^a-z0-9]+",
    ""
  )
  
  
  x[
    x == ""
  ] <- NA_character_
  
  
  x
}


extract_first_author_surname <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  
  first_author <- str_replace(
    x,
    ";.*$",
    ""
  )
  
  
  has_comma <- str_detect(
    first_author,
    ","
  )
  
  
  surname <- ifelse(
    has_comma,
    str_extract(
      first_author,
      "^[^,]+"
    ),
    str_extract(
      first_author,
      "^[^[:space:]]+"
    )
  )
  
  
  normalize_simple_text(
    surname
  )
}


extract_year <- function(x) {
  
  x <- as.character(x)
  
  
  suppressWarnings(
    as.integer(
      str_extract(
        x,
        "(?:18|19|20)\\d{2}"
      )
    )
  )
}


extract_doi <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  
  x <- str_to_lower(x)
  
  
  doi <- str_extract(
    x,
    "10\\.\\d{4,9}/[-._;()/:a-z0-9]+"
  )
  
  
  doi <- str_replace(
    doi,
    "[\\.,;]+$",
    ""
  )
  
  
  doi[
    doi == ""
  ] <- NA_character_
  
  
  doi
}


collapse_unique <- function(x) {
  
  x <- as.character(x)
  
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  
  if (length(x) == 0) {
    
    return(
      NA_character_
    )
  }
  
  
  paste(
    sort(
      unique(x)
    ),
    collapse = " | "
  )
}


reference_string_any_in <- function(
    reference_string,
    reference_ids
) {
  
  if (
    is.na(reference_string) ||
    reference_string == "" ||
    length(reference_ids) == 0
  ) {
    
    return(
      FALSE
    )
  }
  
  
  refs <- str_split(
    reference_string,
    "\\s*\\|\\s*"
  )[[1]]
  
  
  any(
    refs %in% reference_ids
  )
}


# ------------------------------------------------------------
# 8. Import WoS screening output
# ------------------------------------------------------------

wos_all <- read_excel(
  wos_screening_file,
  sheet = "All_Records",
  .name_repair = "unique"
)


cat(
  "\nWoS screening records imported:",
  nrow(wos_all),
  "\n"
)


# ------------------------------------------------------------
# 9. Select only candidate references for overlap screening
# ------------------------------------------------------------

candidate_categories <- c(
  "Potentially eligible",
  "Manual review required"
)


wos_candidates <- wos_all %>%
  filter(
    screening_category %in%
      candidate_categories
  )


cat(
  "WoS candidate references selected:",
  nrow(wos_candidates),
  "\n"
)


if (nrow(wos_candidates) == 0) {
  
  stop(
    "No WoS candidate references were found."
  )
}


if (nrow(wos_candidates) != 75) {
  
  warning(
    paste(
      "Expected 75 candidate records from the final preliminary screening,",
      "but found",
      nrow(wos_candidates),
      ". Check whether the screening workbook has changed."
    )
  )
}


# ------------------------------------------------------------
# 10. Reconstruct matching fields for WoS candidates
# ------------------------------------------------------------

wos_candidates <- wos_candidates %>%
  mutate(
    
    wos_title_normalized =
      normalize_title(
        article_title
      ),
    
    wos_first_author =
      extract_first_author_surname(
        authors
      ),
    
    wos_year =
      suppressWarnings(
        as.integer(
          publication_year
        )
      ),
    
    wos_doi =
      extract_doi(
        doi
      )
  )


# ------------------------------------------------------------
# 11. Detect duplicate candidate titles inside the WoS export
# ------------------------------------------------------------

wos_title_counts <- wos_candidates %>%
  count(
    wos_title_normalized,
    name = "wos_same_title_n"
  )


wos_candidates <- wos_candidates %>%
  left_join(
    wos_title_counts,
    by = "wos_title_normalized"
  ) %>%
  mutate(
    
    wos_duplicate_title_flag =
      !is.na(wos_same_title_n) &
      wos_same_title_n > 1
  )


# ------------------------------------------------------------
# 12. Import ECOTOX raw References sheet
# ------------------------------------------------------------

ecotox_reference_raw <- read_excel(
  ecotox_raw_file,
  sheet = "References",
  .name_repair = "unique"
)


cat(
  "\nRaw ECOTOX sampling-frame references imported:",
  nrow(ecotox_reference_raw),
  "\n"
)


# ------------------------------------------------------------
# 13. Extract ECOTOX reference fields
# ------------------------------------------------------------

ecotox_reference_number_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Ref. Number",
      "Reference Number"
    )
  )


ecotox_author_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Author"
    )
  )


ecotox_title_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Title"
    )
  )


ecotox_source_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Source"
    )
  )


ecotox_year_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Pub. Year",
      "Publication Year"
    )
  )


ecotox_citation_vector <-
  get_optional_column(
    ecotox_reference_raw,
    c(
      "Citation"
    )
  )


# ------------------------------------------------------------
# 14. Prepare ECOTOX reference index
# ------------------------------------------------------------

ecotox_refs <- tibble(
  
  ecotox_reference_number =
    as.character(
      ecotox_reference_number_vector
    ),
  
  ecotox_author =
    ecotox_author_vector,
  
  ecotox_title =
    ecotox_title_vector,
  
  ecotox_source =
    ecotox_source_vector,
  
  ecotox_year =
    extract_year(
      ecotox_year_vector
    ),
  
  ecotox_citation =
    ecotox_citation_vector
  
) %>%
  mutate(
    
    ecotox_title_normalized =
      normalize_title(
        ecotox_title
      ),
    
    ecotox_first_author =
      extract_first_author_surname(
        ecotox_author
      ),
    
    ecotox_doi =
      extract_doi(
        ecotox_citation
      )
  )


# ------------------------------------------------------------
# 15. Import Screening_410 and Primary_295 reference IDs
#
# These are optional for the basic raw-reference match,
# but allow us to identify the level at which a matched
# ECOTOX reference entered the thesis workflow.
# ------------------------------------------------------------

screening_410_reference_ids <- character(0)

primary_295_reference_ids <- character(0)


{
  
  screening_410 <- read_excel(
    ecotox_master_file,
    sheet = "Screening_410",
    .name_repair = "unique"
  )
  
  
  primary_295 <- read_excel(
    ecotox_master_file,
    sheet = "Primary_295",
    .name_repair = "unique"
  )
  
  
  screening_410_reference_ids <-
    unique(
      as.character(
        screening_410$Reference_Number
      )
    )
  
  
  primary_295_reference_ids <-
    unique(
      as.character(
        primary_295$Reference_Number
      )
    )
  
  
  screening_410_reference_ids <-
    screening_410_reference_ids[
      !is.na(
        screening_410_reference_ids
      ) &
        screening_410_reference_ids != ""
    ]
  
  
  primary_295_reference_ids <-
    primary_295_reference_ids[
      !is.na(
        primary_295_reference_ids
      ) &
        primary_295_reference_ids != ""
    ]
  
  
  cat(
    "Unique References in Screening_410:",
    length(
      screening_410_reference_ids
    ),
    "\n"
  )
  
  
  cat(
    "Unique References in Primary_295:",
    length(
      primary_295_reference_ids
    ),
    "\n"
  )
}


# ------------------------------------------------------------
# 16. Create DOI exact-match lookup
# ------------------------------------------------------------

doi_lookup <- ecotox_refs %>%
  filter(
    !is.na(ecotox_doi),
    ecotox_doi != ""
  ) %>%
  group_by(
    ecotox_doi
  ) %>%
  summarise(
    
    ecotox_refs_doi =
      collapse_unique(
        ecotox_reference_number
      ),
    
    ecotox_title_doi =
      first(
        ecotox_title
      ),
    
    ecotox_year_doi =
      first(
        ecotox_year
      ),
    
    ecotox_author_doi =
      first(
        ecotox_author
      ),
    
    ecotox_citation_doi =
      first(
        ecotox_citation
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 17. Create normalized-title exact-match lookup
# ------------------------------------------------------------

title_lookup <- ecotox_refs %>%
  filter(
    !is.na(
      ecotox_title_normalized
    ),
    ecotox_title_normalized != ""
  ) %>%
  group_by(
    ecotox_title_normalized
  ) %>%
  summarise(
    
    ecotox_refs_title =
      collapse_unique(
        ecotox_reference_number
      ),
    
    ecotox_title_exact =
      first(
        ecotox_title
      ),
    
    ecotox_year_exact =
      first(
        ecotox_year
      ),
    
    ecotox_author_exact =
      first(
        ecotox_author
      ),
    
    ecotox_citation_exact =
      first(
        ecotox_citation
      ),
    
    ecotox_doi_exact =
      first(
        ecotox_doi
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 18. Apply DOI and exact-title matching
# ------------------------------------------------------------

wos_overlap <- wos_candidates %>%
  left_join(
    doi_lookup,
    by = c(
      "wos_doi" =
        "ecotox_doi"
    )
  ) %>%
  left_join(
    title_lookup,
    by = c(
      "wos_title_normalized" =
        "ecotox_title_normalized"
    )
  )


# ------------------------------------------------------------
# 19. Determine confirmed exact-match method
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    doi_exact_match =
      !is.na(
        ecotox_refs_doi
      ),
    
    title_exact_match =
      !is.na(
        ecotox_refs_title
      ),
    
    
    doi_title_conflict_flag =
      doi_exact_match &
      title_exact_match &
      ecotox_refs_doi !=
      ecotox_refs_title,
    
    
    exact_match_method =
      case_when(
        
        doi_exact_match &
          title_exact_match ~
          "DOI exact + normalized title exact",
        
        doi_exact_match ~
          "DOI exact",
        
        title_exact_match ~
          "Normalized title exact",
        
        TRUE ~
          "No confirmed exact match"
      ),
    
    
    confirmed_exact_match =
      doi_exact_match |
      title_exact_match,
    
    
    matched_ecotox_reference_numbers =
      coalesce(
        ecotox_refs_doi,
        ecotox_refs_title
      ),
    
    
    matched_ecotox_title =
      coalesce(
        ecotox_title_doi,
        ecotox_title_exact
      ),
    
    
    matched_ecotox_year =
      coalesce(
        ecotox_year_doi,
        ecotox_year_exact
      ),
    
    
    matched_ecotox_author =
      coalesce(
        ecotox_author_doi,
        ecotox_author_exact
      ),
    
    
    matched_ecotox_citation =
      coalesce(
        ecotox_citation_doi,
        ecotox_citation_exact
      ),
    
    
    matched_ecotox_doi =
      coalesce(
        wos_doi,
        ecotox_doi_exact
      )
  )


# ------------------------------------------------------------
# 20. Identify ECOTOX workflow level for confirmed matches
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    matched_in_screening_410 =
      map_lgl(
        matched_ecotox_reference_numbers,
        reference_string_any_in,
        reference_ids =
          screening_410_reference_ids
      ),
    
    
    matched_in_primary_295 =
      map_lgl(
        matched_ecotox_reference_numbers,
        reference_string_any_in,
        reference_ids =
          primary_295_reference_ids
      ),
    
    
    confirmed_match_level =
      case_when(
        
        matched_in_primary_295 ~
          "Matched - Primary_295 Reference",
        
        matched_in_screening_410 ~
          "Matched - Screening_410 Reference, not Primary_295",
        
        confirmed_exact_match ~
          "Matched - Raw ECOTOX sampling-frame Reference only",
        
        TRUE ~
          "No confirmed match"
      )
  )


# ============================================================
# 21. FUZZY TITLE MATCHING
# ============================================================


# ------------------------------------------------------------
# 21.1 Find the best fuzzy ECOTOX title for one WoS title
# ------------------------------------------------------------

best_fuzzy_match <- function(
    query_title,
    query_year,
    query_author,
    reference_table
) {
  
  if (
    is.na(query_title) ||
    query_title == "" ||
    nchar(query_title) < 15
  ) {
    
    return(
      tibble(
        
        fuzzy_ecotox_reference_number =
          NA_character_,
        
        fuzzy_ecotox_title =
          NA_character_,
        
        fuzzy_ecotox_year =
          NA_integer_,
        
        fuzzy_ecotox_author =
          NA_character_,
        
        fuzzy_ecotox_doi =
          NA_character_,
        
        fuzzy_ecotox_citation =
          NA_character_,
        
        fuzzy_title_similarity =
          NA_real_,
        
        fuzzy_year_difference =
          NA_integer_,
        
        fuzzy_year_compatible =
          FALSE,
        
        fuzzy_author_match =
          FALSE,
        
        fuzzy_candidate_basis =
          NA_character_
      )
    )
  }
  
  
  candidates <- reference_table %>%
    filter(
      !is.na(
        ecotox_title_normalized
      ),
      ecotox_title_normalized != ""
    )
  
  
  if (nrow(candidates) == 0) {
    
    return(
      tibble(
        
        fuzzy_ecotox_reference_number =
          NA_character_,
        
        fuzzy_ecotox_title =
          NA_character_,
        
        fuzzy_ecotox_year =
          NA_integer_,
        
        fuzzy_ecotox_author =
          NA_character_,
        
        fuzzy_ecotox_doi =
          NA_character_,
        
        fuzzy_ecotox_citation =
          NA_character_,
        
        fuzzy_title_similarity =
          NA_real_,
        
        fuzzy_year_difference =
          NA_integer_,
        
        fuzzy_year_compatible =
          FALSE,
        
        fuzzy_author_match =
          FALSE,
        
        fuzzy_candidate_basis =
          NA_character_
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Prefer candidate sets with compatible author/year metadata
  # when possible.
  # ----------------------------------------------------------
  
  author_available <-
    !is.na(query_author) &&
    query_author != ""
  
  
  year_available <-
    !is.na(query_year)
  
  
  candidate_author_match <-
    rep(
      FALSE,
      nrow(candidates)
    )
  
  
  candidate_year_match <-
    rep(
      FALSE,
      nrow(candidates)
    )
  
  
  if (author_available) {
    
    candidate_author_match <-
      !is.na(
        candidates$ecotox_first_author
      ) &
      candidates$ecotox_first_author ==
      query_author
  }
  
  
  if (year_available) {
    
    candidate_year_match <-
      !is.na(
        candidates$ecotox_year
      ) &
      abs(
        candidates$ecotox_year -
          query_year
      ) <= 1
  }
  
  
  both_match <-
    candidate_author_match &
    candidate_year_match
  
  
  if (any(both_match)) {
    
    candidates <-
      candidates[
        both_match,
      ]
    
    candidate_basis <-
      "First author + publication year"
    
  } else if (any(candidate_year_match)) {
    
    candidates <-
      candidates[
        candidate_year_match,
      ]
    
    candidate_basis <-
      "Publication year"
    
  } else if (any(candidate_author_match)) {
    
    candidates <-
      candidates[
        candidate_author_match,
      ]
    
    candidate_basis <-
      "First author"
    
  } else {
    
    candidate_basis <-
      "All ECOTOX sampling-frame References"
  }
  
  
  # ----------------------------------------------------------
  # Compute normalized edit-distance similarity
  # ----------------------------------------------------------
  
  distances <- as.numeric(
    adist(
      query_title,
      candidates$ecotox_title_normalized
    )
  )
  
  
  denominators <- pmax(
    nchar(
      query_title
    ),
    nchar(
      candidates$ecotox_title_normalized
    )
  )
  
  
  similarities <-
    1 -
    (
      distances /
        denominators
    )
  
  
  best_index <-
    which.max(
      similarities
    )
  
  
  best_row <-
    candidates[
      best_index,
    ]
  
  
  best_year_difference <-
    if (
      !is.na(query_year) &&
      !is.na(
        best_row$ecotox_year
      )
    ) {
      
      abs(
        best_row$ecotox_year -
          query_year
      )
      
    } else {
      
      NA_integer_
    }
  
  
  best_year_compatible <-
    !is.na(
      best_year_difference
    ) &&
    best_year_difference <= 1
  
  
  best_author_match <-
    !is.na(query_author) &&
    query_author != "" &&
    !is.na(
      best_row$ecotox_first_author
    ) &&
    best_row$ecotox_first_author ==
    query_author
  
  
  tibble(
    
    fuzzy_ecotox_reference_number =
      as.character(
        best_row$ecotox_reference_number
      ),
    
    fuzzy_ecotox_title =
      best_row$ecotox_title,
    
    fuzzy_ecotox_year =
      best_row$ecotox_year,
    
    fuzzy_ecotox_author =
      best_row$ecotox_author,
    
    fuzzy_ecotox_doi =
      best_row$ecotox_doi,
    
    fuzzy_ecotox_citation =
      best_row$ecotox_citation,
    
    fuzzy_title_similarity =
      similarities[
        best_index
      ],
    
    fuzzy_year_difference =
      best_year_difference,
    
    fuzzy_year_compatible =
      best_year_compatible,
    
    fuzzy_author_match =
      best_author_match,
    
    fuzzy_candidate_basis =
      candidate_basis
  )
}


# ------------------------------------------------------------
# 21.2 Run fuzzy matching only for records without an exact
#      confirmed match
# ------------------------------------------------------------

unmatched_rows <- which(
  !wos_overlap$confirmed_exact_match
)


if (length(unmatched_rows) > 0) {
  
  fuzzy_results <- map_dfr(
    unmatched_rows,
    function(i) {
      
      fuzzy_result <- best_fuzzy_match(
        
        query_title =
          wos_overlap$wos_title_normalized[
            i
          ],
        
        query_year =
          wos_overlap$wos_year[
            i
          ],
        
        query_author =
          wos_overlap$wos_first_author[
            i
          ],
        
        reference_table =
          ecotox_refs
      )
      
      
      fuzzy_result %>%
        mutate(
          
          record_id =
            wos_overlap$record_id[
              i
            ],
          
          .before = 1
        )
    }
  )
  
  
  wos_overlap <- wos_overlap %>%
    left_join(
      fuzzy_results,
      by = "record_id"
    )
  
} else {
  
  wos_overlap <- wos_overlap %>%
    mutate(
      
      fuzzy_ecotox_reference_number =
        NA_character_,
      
      fuzzy_ecotox_title =
        NA_character_,
      
      fuzzy_ecotox_year =
        NA_integer_,
      
      fuzzy_ecotox_author =
        NA_character_,
      
      fuzzy_ecotox_doi =
        NA_character_,
      
      fuzzy_ecotox_citation =
        NA_character_,
      
      fuzzy_title_similarity =
        NA_real_,
      
      fuzzy_year_difference =
        NA_integer_,
      
      fuzzy_year_compatible =
        FALSE,
      
      fuzzy_author_match =
        FALSE,
      
      fuzzy_candidate_basis =
        NA_character_
    )
}


# ------------------------------------------------------------
# 22. Define possible fuzzy matches
#
# Fuzzy matches are NEVER automatically confirmed.
#
# A record is flagged as a possible fuzzy match when:
# - title similarity >= 0.90 and author/year metadata support
#   the proposed match
#
# OR
#
# - title similarity >= 0.97 even without metadata support
#
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    fuzzy_metadata_support =
      coalesce(
        fuzzy_year_compatible,
        FALSE
      ) |
      coalesce(
        fuzzy_author_match,
        FALSE
      ),
    
    
    possible_fuzzy_match =
      !confirmed_exact_match &
      !is.na(
        fuzzy_title_similarity
      ) &
      (
        (
          fuzzy_title_similarity >= 0.90 &
            fuzzy_metadata_support
        ) |
          fuzzy_title_similarity >= 0.97
      ),
    
    
    fuzzy_review_band =
      case_when(
        
        confirmed_exact_match ~
          "Not applicable - exact match confirmed",
        
        possible_fuzzy_match ~
          "Possible fuzzy match - manual verification required",
        
        !is.na(
          fuzzy_title_similarity
        ) &
          fuzzy_title_similarity >= 0.80 ~
          "Low-confidence fuzzy suggestion",
        
        TRUE ~
          "No useful fuzzy suggestion"
      )
  )


# ------------------------------------------------------------
# 23. Determine workflow level of the best fuzzy suggestion
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    fuzzy_in_screening_410 =
      !is.na(
        fuzzy_ecotox_reference_number
      ) &
      fuzzy_ecotox_reference_number %in%
      screening_410_reference_ids,
    
    
    fuzzy_in_primary_295 =
      !is.na(
        fuzzy_ecotox_reference_number
      ) &
      fuzzy_ecotox_reference_number %in%
      primary_295_reference_ids,
    
    
    fuzzy_match_level =
      case_when(
        
        fuzzy_in_primary_295 ~
          "Suggested match - Primary_295 Reference",
        
        fuzzy_in_screening_410 ~
          "Suggested match - Screening_410 Reference, not Primary_295",
        
        !is.na(
          fuzzy_ecotox_reference_number
        ) ~
          "Suggested match - Raw ECOTOX sampling-frame Reference",
        
        TRUE ~
          NA_character_
      )
  )


# ------------------------------------------------------------
# 24. Assign final overlap status
#
# "Provisional unique" means no exact match and no sufficiently
# supported fuzzy match was identified in the local ECOTOX
# sampling frame.
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    ecotox_sampling_frame_status =
      case_when(
        
        confirmed_exact_match ~
          "Confirmed match in ECOTOX sampling frame",
        
        possible_fuzzy_match ~
          "Possible ECOTOX sampling-frame match - manual check",
        
        TRUE ~
          "No match found in ECOTOX sampling frame"
      ),
    
    
    provisional_unique_wos =
      !confirmed_exact_match &
      !possible_fuzzy_match
  )


# ------------------------------------------------------------
# 25. Add audit notes
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  mutate(
    
    overlap_audit_note =
      case_when(
        
        confirmed_exact_match ~
          paste0(
            "Reference confirmed by ",
            exact_match_method,
            "."
          ),
        
        
        possible_fuzzy_match ~
          paste0(
            "No exact match. Best fuzzy title similarity = ",
            round(
              fuzzy_title_similarity,
              3
            ),
            "; manual verification required."
          ),
        
        
        !is.na(
          fuzzy_title_similarity
        ) ~
          paste0(
            "No confirmed or high-confidence fuzzy match. ",
            "Best title similarity = ",
            round(
              fuzzy_title_similarity,
              3
            ),
            "."
          ),
        
        
        TRUE ~
          "No ECOTOX sampling-frame match identified."
      )
  )


# ------------------------------------------------------------
# 26. Arrange output
# ------------------------------------------------------------

wos_overlap <- wos_overlap %>%
  arrange(
    
    factor(
      ecotox_sampling_frame_status,
      levels = c(
        "Confirmed match in ECOTOX sampling frame",
        "Possible ECOTOX sampling-frame match - manual check",
        "No match found in ECOTOX sampling frame"
      )
    ),
    
    screening_priority,
    
    desc(
      publication_year
    ),
    
    article_title
  )


# ============================================================
# 27. CREATE OUTPUT SUBSETS
# ============================================================


confirmed_ecotox_matches <- wos_overlap %>%
  filter(
    confirmed_exact_match
  )


possible_fuzzy_matches <- wos_overlap %>%
  filter(
    possible_fuzzy_match
  ) %>%
  arrange(
    desc(
      fuzzy_title_similarity
    )
  )


provisional_unique_wos <- wos_overlap %>%
  filter(
    provisional_unique_wos
  ) %>%
  arrange(
    screening_priority,
    desc(
      publication_year
    ),
    article_title
  )


fuzzy_review_all <- wos_overlap %>%
  filter(
    !confirmed_exact_match
  ) %>%
  select(
    
    record_id,
    wos_record_id,
    screening_category,
    
    article_title,
    publication_year,
    authors,
    metals_detected,
    
    fuzzy_ecotox_reference_number,
    fuzzy_ecotox_title,
    fuzzy_ecotox_year,
    fuzzy_ecotox_author,
    fuzzy_ecotox_doi,
    
    fuzzy_title_similarity,
    fuzzy_year_difference,
    fuzzy_year_compatible,
    fuzzy_author_match,
    fuzzy_candidate_basis,
    
    possible_fuzzy_match,
    fuzzy_review_band,
    fuzzy_match_level
  ) %>%
  arrange(
    desc(
      fuzzy_title_similarity
    )
  )


# ============================================================
# 28. CREATE SUMMARY TABLES
# ============================================================


# ------------------------------------------------------------
# 28.1 Main overlap summary
# ------------------------------------------------------------

overlap_summary <- wos_overlap %>%
  count(
    ecotox_sampling_frame_status,
    name = "references"
  )


# ------------------------------------------------------------
# 28.2 Screening-category by overlap status
# ------------------------------------------------------------

screening_category_overlap <- wos_overlap %>%
  count(
    screening_category,
    ecotox_sampling_frame_status,
    name = "references"
  ) %>%
  pivot_wider(
    names_from =
      ecotox_sampling_frame_status,
    values_from =
      references,
    values_fill = 0
  )


# ------------------------------------------------------------
# 28.3 Confirmed ECOTOX workflow-level summary
# ------------------------------------------------------------

confirmed_level_summary <- wos_overlap %>%
  filter(
    confirmed_exact_match
  ) %>%
  count(
    confirmed_match_level,
    name = "references"
  )


# ------------------------------------------------------------
# 28.4 Provisional unique references by metal
#
# Multi-metal papers contribute once to each detected metal.
# ------------------------------------------------------------

unique_metal_summary <- provisional_unique_wos %>%
  select(
    
    record_id,
    
    copper_flag,
    cadmium_flag,
    zinc_flag,
    nickel_flag,
    lead_flag,
    silver_flag,
    chromium_flag,
    mercury_flag
    
  ) %>%
  pivot_longer(
    
    cols = c(
      copper_flag,
      cadmium_flag,
      zinc_flag,
      nickel_flag,
      lead_flag,
      silver_flag,
      chromium_flag,
      mercury_flag
    ),
    
    names_to =
      "metal",
    
    values_to =
      "detected"
  ) %>%
  filter(
    detected
  ) %>%
  mutate(
    
    metal =
      recode(
        metal,
        
        copper_flag = "Cu",
        cadmium_flag = "Cd",
        zinc_flag = "Zn",
        nickel_flag = "Ni",
        lead_flag = "Pb",
        silver_flag = "Ag",
        chromium_flag = "Cr",
        mercury_flag = "Hg"
      )
  ) %>%
  count(
    metal,
    name = "provisional_unique_references"
  ) %>%
  arrange(
    factor(
      metal,
      levels = c(
        "Cu",
        "Cd",
        "Zn",
        "Ni",
        "Pb",
        "Ag",
        "Cr",
        "Hg"
      )
    )
  )


# ------------------------------------------------------------
# 28.5 Provisional unique references by publication year
# ------------------------------------------------------------

unique_year_summary <- provisional_unique_wos %>%
  count(
    publication_year,
    name =
      "provisional_unique_references"
  ) %>%
  arrange(
    desc(
      publication_year
    )
  )


# ------------------------------------------------------------
# 29. Document match logic
# ------------------------------------------------------------

match_methods <- tibble(
  
  step = c(
    1,
    2,
    3,
    4
  ),
  
  
  method = c(
    "DOI exact match",
    "Normalized title exact match",
    "Fuzzy title suggestion",
    "Provisional unique classification"
  ),
  
  
  rule = c(
    paste0(
      "Exact DOI equality after DOI normalization. ",
      "Used when DOI metadata are available."
    ),
    
    paste0(
      "Exact equality after lowercasing, ASCII transliteration, ",
      "punctuation removal, whitespace normalization, ",
      "and replacement of '&' with 'and'."
    ),
    
    paste0(
      "Best normalized-title edit-distance match. ",
      "A fuzzy suggestion is flagged when similarity >= 0.90 ",
      "with compatible first-author or publication-year metadata, ",
      "or similarity >= 0.97 without metadata support."
    ),
    
    paste0(
      "No DOI/title exact match and no sufficiently supported ",
      "fuzzy match in the ECOTOX sampling frame. ",
      "This remains provisional until manual/source verification."
    )
  ),
  
  
  interpretation = c(
    "Confirmed sampling-frame overlap",
    "Confirmed sampling-frame overlap",
    "Manual verification only; never automatically confirmed",
    "Candidate for later full-text eligibility assessment"
  )
)


# ------------------------------------------------------------
# 30. Create an ECOTOX reference index for audit
# ------------------------------------------------------------

ecotox_reference_index <- ecotox_refs %>%
  mutate(
    
    in_screening_410 =
      ecotox_reference_number %in%
      screening_410_reference_ids,
    
    in_primary_295 =
      ecotox_reference_number %in%
      primary_295_reference_ids,
    
    
    thesis_workflow_level =
      case_when(
        
        in_primary_295 ~
          "Primary_295",
        
        in_screening_410 ~
          "Screening_410, not Primary_295",
        
        TRUE ~
          "Raw ECOTOX sampling frame only"
      )
  ) %>%
  select(
    
    ecotox_reference_number,
    ecotox_author,
    ecotox_title,
    ecotox_source,
    ecotox_year,
    ecotox_doi,
    ecotox_citation,
    
    in_screening_410,
    in_primary_295,
    thesis_workflow_level,
    
    ecotox_title_normalized,
    ecotox_first_author
  )


# ------------------------------------------------------------
# 31. Integrity checks
# ------------------------------------------------------------

overlap_count_check <-
  nrow(
    confirmed_ecotox_matches
  ) +
  nrow(
    possible_fuzzy_matches
  ) +
  nrow(
    provisional_unique_wos
  )


if (
  overlap_count_check !=
  nrow(
    wos_overlap
  )
) {
  
  stop(
    paste(
      "Overlap classification check failed:",
      overlap_count_check,
      "classified out of",
      nrow(
        wos_overlap
      )
    )
  )
}


# ------------------------------------------------------------
# 32. Export workbook
# ------------------------------------------------------------

write_xlsx(
  list(
    
    Overlap_Summary =
      overlap_summary,
    
    Screen_Category_Overlap =
      screening_category_overlap,
    
    Confirmed_Level_Summary =
      confirmed_level_summary,
    
    Unique_Metal_Summary =
      unique_metal_summary,
    
    Unique_Year_Summary =
      unique_year_summary,
    
    Match_Methods =
      match_methods,
    
    Confirmed_ECOTOX =
      confirmed_ecotox_matches,
    
    Possible_Fuzzy_Matches =
      possible_fuzzy_matches,
    
    Provisional_Unique_WoS =
      provisional_unique_wos,
    
    Fuzzy_Review_All =
      fuzzy_review_all,
    
    All_Candidates_Overlap =
      wos_overlap,
    
    ECOTOX_Reference_Index =
      ecotox_reference_index
  ),
  
  output_file
)


# ------------------------------------------------------------
# 33. Final console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "WOS - ECOTOX OVERLAP AUDIT\n"
)

cat(
  "============================================\n"
)


cat(
  "\nWoS candidate references:",
  nrow(
    wos_overlap
  ),
  "\n"
)


cat(
  "Confirmed matches in ECOTOX sampling frame:",
  nrow(
    confirmed_ecotox_matches
  ),
  "\n"
)


cat(
  "Possible fuzzy matches requiring manual check:",
  nrow(
    possible_fuzzy_matches
  ),
  "\n"
)


cat(
  "Provisional unique WoS references:",
  nrow(
    provisional_unique_wos
  ),
  "\n"
)


cat(
  "\nOverlap classification check:",
  overlap_count_check,
  "records classified out of",
  nrow(
    wos_overlap
  ),
  "\n"
)


{
  
  cat(
    "\nConfirmed matches represented in Primary_295:",
    sum(
      confirmed_ecotox_matches$matched_in_primary_295,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "Confirmed matches represented in Screening_410 but not Primary_295:",
    sum(
      confirmed_ecotox_matches$matched_in_screening_410 &
        !confirmed_ecotox_matches$matched_in_primary_295,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "Confirmed matches present only in the raw ECOTOX sampling frame:",
    sum(
      confirmed_ecotox_matches$confirmed_exact_match &
        !confirmed_ecotox_matches$matched_in_screening_410,
      na.rm = TRUE
    ),
    "\n"
  )
}


cat(
  "\nOutput saved as:",
  output_file,
  "\n"
)


cat(
  "\nIMPORTANT:\n"
)


cat(
  paste0(
    "The 'Provisional_Unique_WoS' sheet does not yet represent ",
    "final eligible external evidence. It contains candidate WoS ",
    "references for which no confirmed or sufficiently supported ",
    "match was identified in the ECOTOX sampling frame. ",
    "Possible fuzzy matches must be manually resolved before ",
    "full-text eligibility assessment.\n"
  )
)