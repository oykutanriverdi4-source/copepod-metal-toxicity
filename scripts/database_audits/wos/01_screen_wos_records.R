# ============================================================
# 01_screen_wos_records.R
# ============================================================
# Purpose
# Reproduce the preliminary title/abstract screening of the
# Web of Science Core Collection search used to identify
# supplementary copepod acute-metal toxicity literature.
#
# Search represented by the local raw export
# Database: Web of Science Core Collection
# Search date: 5 September 2026
# Field: Topic
# Query:
# copepod* AND ("LC50" OR "LC 50" OR "median lethal concentration") AND
# (copper OR cadmium OR zinc OR nickel OR lead OR silver OR chromium OR mercury)
#
# Final thesis benchmark
# - 117 bibliographic records in the raw export
# - 75 records retained as preliminary candidates for ECOTOX-overlap
#   assessment (Potentially eligible + Manual review required)
#
# Important scope notes
# - This is title/abstract/author-keyword screening only.
# - Keywords Plus is retained for audit purposes because it is part
#   of the WoS Topic search, but it is NOT used for the subsequent
#   screening decision.
# - This script does NOT compare records with ECOTOX.
# - This script does NOT make final full-text eligibility decisions.
# - Missing information is retained for manual review rather than
#   treated automatically as exclusion evidence.
#
# Local-only input
# data/external/wos/wos_raw_export.xls
#
# The WoS raw export is not redistributed in the public repository.
# Users with WoS access can reproduce this step by running the exact
# search above and placing their export at the expected local path.
#
# Local audit output
# outputs/database_audits/wos/01_screen_wos_records/
#
# NOTE:
# The detailed screening workbook contains bibliographic metadata
# from the local WoS export and should remain local/gitignored.
# A compact count summary and run metadata are also written.
#
# Scientific rule
# Screening logic, text parsing, classification categories and the
# treatment of uncertain records are preserved from the final thesis
# workflow. Repository cleanup changes only paths, naming,
# documentation and reproducibility checks.
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
# 2. Define file paths
# ------------------------------------------------------------

wos_file <- file.path(
  "data",
  "external",
  "wos",
  "wos_raw_export.xls"
)

output_root <- file.path(
  "outputs",
  "database_audits",
  "wos",
  "01_screen_wos_records"
)

if (!dir.exists(output_root)) {
  dir.create(output_root, recursive = TRUE)
}

output_file <- file.path(
  output_root,
  "wos_preliminary_screening.xlsx"
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
# 3. Check input file
# ------------------------------------------------------------

if (!file.exists(wos_file)) {
  
  stop(
    paste0(
      "Local WoS export not found: ",
      wos_file,
      ". Run the documented Web of Science search and place the ",
      "export at this path. The raw export is intentionally not ",
      "redistributed in the public repository."
    )
  )
}


# ------------------------------------------------------------
# 4. Define reusable regex patterns
# ------------------------------------------------------------

# Recognises:
# LC50
# LC 50
# LC-50
# LC(50)
# LC (50)

lc50_pattern <- paste0(
  "\\bLC\\s*[-_ ]?\\(?\\s*50\\s*\\)?",
  "|median lethal concentration",
  "|lethal concentration\\s*\\(?\\s*50\\s*\\)?"
)


# ASCII-safe representations of optional dash characters

dash_pattern <- "(?:-|\\x{2013}|\\x{2014})"


# ------------------------------------------------------------
# 5. Helper functions
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


extract_duration_hours <- function(x) {
  
  x <- as.character(x)
  
  if (
    is.na(x) ||
    x == ""
  ) {
    
    return(
      numeric(0)
    )
  }
  
  
  duration_regex <- paste0(
    "\\b",
    "(\\d+(?:\\.\\d+)?)",
    "\\s*",
    dash_pattern,
    "?",
    "\\s*",
    "(h|hr|hrs|hour|hours|d|day|days)",
    "\\b"
  )
  
  
  matches <- str_match_all(
    x,
    regex(
      duration_regex,
      ignore_case = TRUE
    )
  )[[1]]
  
  
  if (nrow(matches) == 0) {
    
    return(
      numeric(0)
    )
  }
  
  
  values <- suppressWarnings(
    as.numeric(
      matches[, 2]
    )
  )
  
  
  units <- str_to_lower(
    matches[, 3]
  )
  
  
  day_units <- units %in%
    c(
      "d",
      "day",
      "days"
    )
  
  
  values[
    day_units
  ] <- values[
    day_units
  ] * 24
  
  
  values <- values[
    !is.na(values) &
      is.finite(values)
  ]
  
  
  sort(
    unique(values)
  )
}


extract_lc50_duration_hours <- function(x) {
  
  x <- as.character(x)
  
  if (
    is.na(x) ||
    x == ""
  ) {
    
    return(
      numeric(0)
    )
  }
  
  
  # Split the text into sentence-like units.
  # Only units explicitly containing LC50 are used for
  # LC50-specific duration screening.
  
  text_units <- unlist(
    str_split(
      x,
      "(?<=[.!?;])\\s+"
    )
  )
  
  
  lc50_units <- text_units[
    str_detect(
      text_units,
      regex(
        lc50_pattern,
        ignore_case = TRUE
      )
    )
  ]
  
  
  if (length(lc50_units) == 0) {
    
    return(
      numeric(0)
    )
  }
  
  
  values <- unlist(
    map(
      lc50_units,
      extract_duration_hours
    ),
    use.names = FALSE
  )
  
  
  values <- values[
    !is.na(values) &
      is.finite(values)
  ]
  
  
  sort(
    unique(values)
  )
}


format_duration_values <- function(x) {
  
  if (length(x) == 0) {
    
    return("")
  }
  
  paste(
    x,
    collapse = ", "
  )
}


collapse_flags <- function(...) {
  
  x <- unlist(
    list(...),
    use.names = FALSE
  )
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  if (length(x) == 0) {
    
    return("")
  }
  
  paste(
    unique(x),
    collapse = "; "
  )
}


# ------------------------------------------------------------
# 6. Import Web of Science export
# ------------------------------------------------------------

wos_raw <- suppressWarnings(
  read_excel(
    wos_file,
    col_types = "text",
    .name_repair = "unique"
  )
)


cat(
  "\nWeb of Science records imported:",
  nrow(wos_raw),
  "\n"
)


cat(
  "Number of variables imported:",
  ncol(wos_raw),
  "\n"
)


# Frozen search-export benchmark from the final thesis workflow.
stopifnot(
  nrow(wos_raw) == 117L
)


# ------------------------------------------------------------
# 7. Extract relevant WoS fields
# ------------------------------------------------------------

title_vector <- get_optional_column(
  wos_raw,
  c(
    "Article Title",
    "Title"
  )
)


abstract_vector <- get_optional_column(
  wos_raw,
  c(
    "Abstract"
  )
)


author_keywords_vector <- get_optional_column(
  wos_raw,
  c(
    "Author Keywords"
  )
)


keywords_plus_vector <- get_optional_column(
  wos_raw,
  c(
    "Keywords Plus"
  )
)


authors_vector <- get_optional_column(
  wos_raw,
  c(
    "Authors",
    "Author Full Names"
  )
)


source_vector <- get_optional_column(
  wos_raw,
  c(
    "Source Title"
  )
)


year_vector <- get_optional_column(
  wos_raw,
  c(
    "Publication Year"
  )
)


document_type_vector <- get_optional_column(
  wos_raw,
  c(
    "Document Type"
  )
)


doi_vector <- get_optional_column(
  wos_raw,
  c(
    "DOI Link",
    "DOI"
  )
)


wos_id_vector <- get_optional_column(
  wos_raw,
  c(
    "UT (Unique WOS ID)",
    "UT"
  )
)


# ------------------------------------------------------------
# 8. Prepare screening dataset
#
# Keywords Plus is intentionally excluded from screening_text.
# ------------------------------------------------------------

wos_screen <- wos_raw %>%
  transmute(
    
    record_id = row_number(),
    
    wos_record_id = if_else(
      wos_id_vector == "",
      paste0(
        "WOS_",
        record_id
      ),
      wos_id_vector
    ),
    
    authors =
      authors_vector,
    
    article_title =
      title_vector,
    
    publication_year =
      extract_year(
        year_vector
      ),
    
    source_title =
      source_vector,
    
    document_type =
      document_type_vector,
    
    author_keywords =
      author_keywords_vector,
    
    keywords_plus =
      keywords_plus_vector,
    
    abstract =
      abstract_vector,
    
    doi =
      extract_doi(
        doi_vector
      ),
    
    normalized_title =
      normalize_title(
        title_vector
      )
    
  ) %>%
  mutate(
    
    title_lower =
      str_to_lower(
        article_title
      ),
    
    author_keywords_lower =
      str_to_lower(
        author_keywords
      ),
    
    title_keyword_text =
      str_squish(
        paste(
          article_title,
          author_keywords,
          sep = " | "
        )
      ),
    
    screening_text =
      str_squish(
        paste(
          article_title,
          author_keywords,
          abstract,
          sep = " | "
        )
      ),
    
    screening_text_lower =
      str_to_lower(
        screening_text
      )
  )


# ------------------------------------------------------------
# 9. Define copepod terms
# ------------------------------------------------------------

copepod_pattern <- paste(
  c(
    "copepod",
    "copepoda",
    "calanoid",
    "calanoida",
    "harpacticoid",
    "harpacticoida",
    "cyclopoid",
    "cyclopoida",
    "acartia",
    "eurytemora",
    "tigriopus",
    "tisbe",
    "nitocra",
    "amphiascus",
    "pseudodiaptomus",
    "paracyclopina",
    "paracalanus",
    "calanus",
    "temora",
    "scutellidium",
    "robertsonia",
    "acanthocyclops"
  ),
  collapse = "|"
)


# ------------------------------------------------------------
# 10. Identify copepod relevance
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    copepod_core =
      str_detect(
        title_keyword_text,
        regex(
          copepod_pattern,
          ignore_case = TRUE
        )
      ),
    
    copepod_anywhere =
      str_detect(
        screening_text,
        regex(
          copepod_pattern,
          ignore_case = TRUE
        )
      ),
    
    copepod_relevance =
      case_when(
        
        copepod_core ~
          "Core title/keyword relevance",
        
        copepod_anywhere ~
          "Copepod identified in abstract",
        
        TRUE ~
          "No clear copepod relevance"
      )
  )


# ------------------------------------------------------------
# 11. Identify the eight target metals
#
# Lead is treated conservatively because "lead" can also be
# an ordinary English verb.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    copper_flag =
      str_detect(
        screening_text_lower,
        "\\bcopper\\b|\\bcupric\\b"
      ) |
      str_detect(
        screening_text,
        "\\bCu\\b"
      ),
    
    
    cadmium_flag =
      str_detect(
        screening_text_lower,
        "\\bcadmium\\b"
      ) |
      str_detect(
        screening_text,
        "\\bCd\\b"
      ),
    
    
    zinc_flag =
      str_detect(
        screening_text_lower,
        "\\bzinc\\b"
      ) |
      str_detect(
        screening_text,
        "\\bZn\\b"
      ),
    
    
    nickel_flag =
      str_detect(
        screening_text_lower,
        "\\bnickel\\b"
      ) |
      str_detect(
        screening_text,
        "\\bNi\\b"
      ),
    
    
    lead_symbol_flag =
      str_detect(
        screening_text,
        "\\bPb\\b"
      ),
    
    
    lead_keyword_flag =
      str_detect(
        author_keywords_lower,
        "(^|[;,|])\\s*lead\\s*($|[;,|])"
      ),
    
    
    lead_chemical_context_flag =
      str_detect(
        screening_text_lower,
        paste0(
          
          "\\blead\\s+",
          "(?:toxicity|toxicities|exposure|",
          "concentration|concentrations|",
          "contamination|levels?|ions?|",
          "nitrate|chloride|acetate|sulfate|metal)",
          "\\b",
          
          "|",
          
          "\\b(?:toxicity|toxicities|exposure|",
          "concentration|concentrations|",
          "contamination)\\s+",
          "(?:of|to|with)\\s+lead\\b",
          
          "|",
          
          "\\blead\\s*(?:,|;|/|&|and)\\s*",
          "(?:copper|cadmium|zinc|nickel|silver|chromium|mercury)\\b",
          
          "|",
          
          "\\b(?:copper|cadmium|zinc|nickel|silver|chromium|mercury)",
          "\\s*(?:,|;|/|&|and)\\s*lead\\b"
        )
      ),
    
    
    lead_flag =
      lead_symbol_flag |
      lead_keyword_flag |
      lead_chemical_context_flag,
    
    
    silver_flag =
      str_detect(
        screening_text_lower,
        "\\bsilver\\b"
      ) |
      str_detect(
        screening_text,
        "\\bAg\\b"
      ),
    
    
    chromium_flag =
      str_detect(
        screening_text_lower,
        "\\bchromium\\b"
      ) |
      str_detect(
        screening_text,
        "\\bCr\\b"
      ),
    
    
    mercury_flag =
      str_detect(
        screening_text_lower,
        "\\bmercury\\b|\\bmercuric\\b"
      ) |
      str_detect(
        screening_text,
        "\\bHg\\b"
      ),
    
    
    # Generic metal wording is retained as a manual-review
    # signal even when the individual target metal is not
    # explicitly identified in the WoS metadata.
    
    generic_metal_scope_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\btrace metals?\\b|",
          "\\bheavy metals?\\b|",
          "\\bmetal toxicity\\b|",
          "\\bmetal toxicities\\b|",
          "\\btoxicity of metals?\\b"
        )
      )
  )


# ------------------------------------------------------------
# 12. Summarise metal detection
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  rowwise() %>%
  mutate(
    
    metals_detected =
      paste(
        c(
          if (copper_flag) "Cu" else NULL,
          if (cadmium_flag) "Cd" else NULL,
          if (zinc_flag) "Zn" else NULL,
          if (nickel_flag) "Ni" else NULL,
          if (lead_flag) "Pb" else NULL,
          if (silver_flag) "Ag" else NULL,
          if (chromium_flag) "Cr" else NULL,
          if (mercury_flag) "Hg" else NULL
        ),
        collapse = ", "
      ),
    
    metal_count =
      sum(
        c(
          copper_flag,
          cadmium_flag,
          zinc_flag,
          nickel_flag,
          lead_flag,
          silver_flag,
          chromium_flag,
          mercury_flag
        )
      ),
    
    target_metal =
      metal_count > 0,
    
    metal_scope_relevant =
      target_metal |
      generic_metal_scope_flag
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 13. Identify LC50 and acute lethal relevance
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    lc50_flag =
      str_detect(
        screening_text,
        regex(
          lc50_pattern,
          ignore_case = TRUE
        )
      ),
    
    
    mortality_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bmortalit|",
          "\\blethal|",
          "\\bsurvival|",
          "\\bsurvivability"
        )
      ),
    
    
    acute_flag =
      str_detect(
        screening_text_lower,
        "\\bacute\\b"
      ),
    
    
    acute_lethal_relevance =
      lc50_flag |
      (
        acute_flag &
          mortality_flag
      ),
    
    
    numeric_lc50_in_abstract =
      str_detect(
        abstract,
        regex(
          paste0(
            "(?:",
            lc50_pattern,
            ")",
            ".{0,120}",
            "\\d+(?:\\.\\d+)?"
          ),
          ignore_case = TRUE
        )
      ),
    
    
    # Recognises:
    # LC50 >
    # LC50 <
    # LC50 >=
    # LC50 <=
    # and Unicode >= / <= if present in the source text.
    
    censored_lc50_flag =
      str_detect(
        screening_text,
        regex(
          paste0(
            "(?:",
            lc50_pattern,
            ")",
            "\\s*[:=]?\\s*",
            "(?:[<>]=?|\\x{2265}|\\x{2264})"
          ),
          ignore_case = TRUE
        )
      ),
    
    
    endpoint_status =
      case_when(
        
        lc50_flag ~
          "LC50 identified",
        
        acute_flag &
          mortality_flag ~
          "Possible acute lethal endpoint",
        
        TRUE ~
          "No clear acute LC50 endpoint"
      )
  )


# ------------------------------------------------------------
# 14. Identify marine or estuarine relevance
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    marine_estuarine_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bmarine\\b|",
          "\\bestuar|",
          "\\bsalt[- ]?water\\b|",
          "\\bsea[- ]?water\\b|",
          "\\bbrackish\\b|",
          "\\bcoastal\\b|",
          "\\bintertidal\\b|",
          "\\boceanic\\b|",
          "\\bharbou?r\\b"
        )
      ),
    
    
    freshwater_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bfresh[- ]?water\\b|",
          "\\blimnetic\\b"
        )
      ),
    
    
    freshwater_only_flag =
      freshwater_flag &
      !marine_estuarine_flag,
    
    
    habitat_status =
      case_when(
        
        marine_estuarine_flag &
          freshwater_flag ~
          "Mixed marine/freshwater context",
        
        marine_estuarine_flag ~
          "Marine/estuarine evidence",
        
        freshwater_only_flag ~
          "Freshwater-only evidence",
        
        TRUE ~
          "Habitat requires verification"
      )
  )


# ------------------------------------------------------------
# 15. Identify laboratory-test evidence
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    laboratory_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\blaborator|",
          "\\bbioassay|",
          "\\btoxicity test|",
          "\\bexperimental|",
          "\\bexpos(?:ed|ure)"
        )
      ),
    
    
    laboratory_status =
      if_else(
        laboratory_flag,
        "Laboratory/experimental evidence",
        "Laboratory setting requires verification"
      )
  )


# ------------------------------------------------------------
# 16. Extract all reported durations
#
# These values are descriptive only.
# A chronic duration elsewhere in the abstract must NOT
# automatically exclude an acute LC50 record.
# ------------------------------------------------------------

all_duration_values <- map(
  wos_screen$screening_text,
  extract_duration_hours
)


wos_screen$all_durations_hours <-
  map_chr(
    all_duration_values,
    format_duration_values
  )


# ------------------------------------------------------------
# 17. Extract LC50-specific durations
#
# Only durations from text units explicitly containing LC50
# are allowed to generate duration-based scope decisions.
# ------------------------------------------------------------

lc50_duration_values <- map(
  wos_screen$screening_text,
  extract_lc50_duration_hours
)


wos_screen$lc50_durations_hours <-
  map_chr(
    lc50_duration_values,
    format_duration_values
  )


wos_screen$lc50_duration_detected <-
  map_lgl(
    lc50_duration_values,
    ~ length(.x) > 0
  )


wos_screen$lc50_duration_24_96h_flag <-
  map_lgl(
    lc50_duration_values,
    ~ any(
      .x >= 24 &
        .x <= 96
    )
  )


wos_screen <- wos_screen %>%
  mutate(
    
    lc50_duration_status =
      case_when(
        
        lc50_duration_24_96h_flag ~
          "24-96 h LC50 duration identified",
        
        lc50_duration_detected &
          !lc50_duration_24_96h_flag ~
          "Only LC50 duration(s) outside 24-96 h identified",
        
        lc50_flag &
          !lc50_duration_detected ~
          "LC50 duration requires verification",
        
        TRUE ~
          "LC50 duration not applicable"
      )
  )


# ------------------------------------------------------------
# 18. Identify waterborne-exposure evidence
#
# Waterborne evidence can be explicit wording or an aqueous
# concentration unit such as mg/L or ug/L.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    explicit_waterborne_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bwaterborne\\b|",
          "\\baqueous exposure\\b|",
          "\\bwater exposure\\b|",
          "\\btest solution\\b|",
          "\\bexposed in seawater\\b|",
          "\\bexposed in sea water\\b|",
          "\\bdissolved metal\\b|",
          "\\bdissolved metals\\b"
        )
      ),
    
    
    aqueous_concentration_unit_flag =
      str_detect(
        screening_text,
        regex(
          paste0(
            
            "(?:",
            
            # ug, microgram, mg, ng, etc.
            "(?:u|\\x{00B5}|\\x{03BC})?g",
            "|mg",
            "|ng",
            "|(?:u|\\x{00B5}|\\x{03BC})mol",
            "|mmol",
            "|mol",
            
            ")",
            
            "\\s*(?:[A-Za-z]+\\s*)?",
            
            "(?:",
            "/\\s*L",
            "|L\\s*(?:-|\\x{2212})?1",
            ")",
            
            "|",
            
            "\\b(?:uM|mM|nM)\\b"
          ),
          ignore_case = TRUE
        )
      ),
    
    
    waterborne_evidence_flag =
      explicit_waterborne_flag |
      aqueous_concentration_unit_flag,
    
    
    dietary_exposure_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bdietary exposure\\b|",
          "\\bdietary uptake\\b|",
          "\\bfoodborne\\b|",
          "\\bvia (?:the )?diet\\b|",
          "\\bcontaminated food\\b|",
          "\\bcontaminated algae\\b"
        )
      )
  )


# ------------------------------------------------------------
# 19. Identify sediment-related exposure
#
# Simple mention of sediment is separated from stronger
# evidence that sediment itself was the test matrix.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    sediment_mention_flag =
      str_detect(
        screening_text_lower,
        "\\bsediment"
      ),
    
    
    sediment_matrix_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bsediment toxicity\\b|",
          "\\bsediment bioassay\\b|",
          "\\bsediment exposure\\b|",
          "\\bspiked sediment\\b|",
          "\\bsediment-associated\\b|",
          "\\bpore[- ]?water\\b|",
          "\\bsediment-dwelling\\b|",
          "\\bresting eggs? in (?:the )?sediment\\b|",
          "\\bkg\\s+(?:dry\\s+)?sediment\\b|",
          
          "\\bsediment.{0,40}",
          "(?:mg|ug|ng|mmol|mol)\\s*/?\\s*kg\\b",
          
          "|",
          
          "(?:mg|ug|ng|mmol|mol).{0,30}",
          "kg.{0,20}sediment"
        )
      ),
    
    
    sediment_matrix_only_flag =
      sediment_matrix_flag &
      !waterborne_evidence_flag
  )


# ------------------------------------------------------------
# 20. Summarise exposure route
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    exposure_route_status =
      case_when(
        
        waterborne_evidence_flag &
          sediment_matrix_flag ~
          "Waterborne and sediment-related evidence both present",
        
        waterborne_evidence_flag ~
          "Waterborne exposure evidence identified",
        
        dietary_exposure_flag ~
          "Dietary exposure mentioned",
        
        sediment_matrix_flag ~
          "Sediment exposure matrix identified",
        
        sediment_mention_flag ~
          "Sediment mentioned; exposure route requires verification",
        
        TRUE ~
          "Exposure route requires verification"
      )
  )


# ------------------------------------------------------------
# 21. Identify mixture and multi-metal study structure
#
# Multi-metal papers are NOT automatically excluded because
# metals may have been tested separately.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    mixture_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bmixture|",
          "\\bcombined toxic|",
          "\\bcombined effect|",
          "\\bjoint action|",
          "\\bjoint toxic|",
          "\\bbinary mixture|",
          "\\bternary mixture|",
          "\\bco[- ]?exposure"
        )
      ),
    
    
    isolated_single_arm_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bisolated\\b|",
          "\\bindividual(?:ly)?\\b|",
          "\\bsingle[- ]metal\\b|",
          "\\bseparately\\b|",
          "\\bmetal alone\\b|",
          "\\bindividual metals\\b"
        )
      ),
    
    
    exposure_structure =
      case_when(
        
        mixture_flag &
          isolated_single_arm_flag ~
          "Single-metal and mixture arms mentioned",
        
        mixture_flag ~
          "Mixture study - verify separate single-metal arms",
        
        metal_count > 1 ~
          "Multiple target metals - verify separate exposures",
        
        metal_count == 1 ~
          "Single target metal identified",
        
        generic_metal_scope_flag ~
          "Metal identity requires verification",
        
        TRUE ~
          "Exposure structure unclear"
      )
  )


# ------------------------------------------------------------
# 22. Identify deliberate experimental modifiers
#
# Ordinary salinity and temperature variation is not itself
# treated as an exclusion criterion.
#
# Explicit co-stressors or modifier treatments are sent to
# manual review.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    deliberate_modifier_flag =
      str_detect(
        screening_text_lower,
        paste0(
          
          "dissolved organic matter|",
          "dissolved organic carbon|",
          "\\bhumic\\b|",
          "\\bfulvic\\b|",
          
          "\\bultraviolet\\b|",
          "\\buv[- ]?[ab]?\\b|",
          
          "\\bocean acidification\\b|",
          "\\bseawater acidification\\b|",
          "\\bacidified seawater\\b|",
          "\\belevated pco2\\b|",
          "\\bhigh pco2\\b|",
          "\\bco2-driven\\b|",
          
          "\\boceanic warming\\b|",
          "\\bnear-future warming\\b|",
          "\\bheat exposure\\b|",
          "\\bheat stress\\b|",
          "\\bthermal stress\\b|",
          
          "\\bquicklime\\b|",
          "\\bcalcium oxide\\b|",
          
          "\\bfood factor\\b|",
          "\\bfood availability\\b|",
          "\\bnutrition\\b|",
          "\\bnutritional\\b|",
          "\\bfeeding regime\\b|",
          "\\bfeeding condition\\b|",
          "\\bwith food\\b|",
          "\\bwithout food\\b|",
          "\\bstarvation\\b|",
          "\\bdiet manipulation\\b|",
          
          "\\bpredator cue|",
          
          "\\bpre[- ]?exposure\\b|",
          "\\bpretreated\\b|",
          "\\bpre-treatment\\b"
        )
      )
  )


# ------------------------------------------------------------
# 23. Identify chemical-form issues
#
# These generate manual-review flags rather than automatic
# exclusions because an article may also contain an eligible
# plain-metal arm.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    nanoparticle_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "\\bnanopart|",
          "\\bnano[- ]?(?:zinc|copper|silver)|",
          "\\bzinc oxide nanoparticle|",
          "\\bzno nanoparticle"
        )
      ),
    
    
    zinc_pyrithione_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "zinc pyrithione|",
          "\\bznpt\\b"
        )
      ),
    
    
    chromium_iii_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "trivalent chromium|",
          "chromium\\s*\\(iii\\)|",
          "\\bcr\\s*\\(iii\\)"
        )
      ),
    
    
    chromium_vi_flag =
      str_detect(
        screening_text_lower,
        paste0(
          "hexavalent chromium|",
          "chromium\\s*\\(vi\\)|",
          "\\bcr\\s*\\(vi\\)|",
          "\\bdichromate\\b"
        )
      ),
    
    
    leaching_flag =
      str_detect(
        screening_text_lower,
        "\\bleach(?:ing|ate)"
      ),
    
    
    chemistry_review_flag =
      nanoparticle_flag |
      zinc_pyrithione_flag |
      chromium_iii_flag |
      chromium_vi_flag |
      leaching_flag
  )


# ------------------------------------------------------------
# 24. Identify secondary literature
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    secondary_source_flag =
      str_detect(
        str_to_lower(
          paste(
            document_type,
            article_title
          )
        ),
        paste0(
          "\\breview\\b|",
          "\\bmeta-analysis\\b|",
          "\\bmeta analysis\\b|",
          "\\bsystematic review\\b"
        )
      )
  )


# ------------------------------------------------------------
# 25. Define explicit scope failures
#
# IMPORTANT:
# Only LC50-specific durations can generate a duration-based
# exclusion signal.
#
# Other durations elsewhere in the abstract are descriptive
# and cannot exclude a record.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    explicit_lc50_duration_failure =
      lc50_duration_detected &
      !lc50_duration_24_96h_flag,
    
    
    explicit_scope_failure =
      freshwater_only_flag |
      explicit_lc50_duration_failure |
      sediment_matrix_only_flag,
    
    
    complex_design_flag =
      mixture_flag |
      dietary_exposure_flag |
      sediment_matrix_flag |
      deliberate_modifier_flag |
      chemistry_review_flag |
      censored_lc50_flag
  )


# ------------------------------------------------------------
# 26. Construct review flags
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  rowwise() %>%
  mutate(
    
    review_flags =
      collapse_flags(
        
        if (
          !target_metal &
          generic_metal_scope_flag
        )
          "Generic metal terminology; target-metal identity requires verification"
        else "",
        
        
        if (
          !copepod_core &
          copepod_anywhere
        )
          "Copepod relevance mainly identified from abstract"
        else "",
        
        
        if (
          !lc50_flag &
          acute_lethal_relevance
        )
          "Acute lethal response identified but LC50 is not explicit"
        else "",
        
        
        if (
          !marine_estuarine_flag &
          !freshwater_only_flag
        )
          "Marine/estuarine setting not explicit"
        else "",
        
        
        if (!laboratory_flag)
          "Laboratory setting not explicit"
        else "",
        
        
        if (
          lc50_flag &
          !lc50_duration_detected
        )
          "LC50 duration not explicit"
        else "",
        
        
        if (!waterborne_evidence_flag)
          "Waterborne exposure not explicit"
        else "",
        
        
        if (mixture_flag)
          "Mixture or combined exposure mentioned"
        else "",
        
        
        if (
          metal_count > 1 &
          !mixture_flag
        )
          "Multiple target metals reported; verify that exposures were separate"
        else "",
        
        
        if (dietary_exposure_flag)
          "Dietary exposure mentioned"
        else "",
        
        
        if (sediment_matrix_flag)
          "Sediment exposure matrix mentioned"
        else if (sediment_mention_flag)
          "Sediment mentioned in study context"
        else "",
        
        
        if (deliberate_modifier_flag)
          "Deliberate experimental modifier or co-stressor mentioned"
        else "",
        
        
        if (nanoparticle_flag)
          "Nanoparticle exposure mentioned"
        else "",
        
        
        if (zinc_pyrithione_flag)
          "Zinc pyrithione mentioned"
        else "",
        
        
        if (chromium_iii_flag)
          "Cr(III) chemistry mentioned"
        else "",
        
        
        if (chromium_vi_flag)
          "Cr(VI) chemistry requires concentration-basis verification"
        else "",
        
        
        if (leaching_flag)
          "Leaching exposure mentioned"
        else "",
        
        
        if (censored_lc50_flag)
          "Censored LC50 value may be present"
        else ""
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 27. Construct fields requiring source verification
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  rowwise() %>%
  mutate(
    
    criteria_to_verify =
      collapse_flags(
        
        if (!copepod_core)
          "Tested organism"
        else "",
        
        
        if (
          !target_metal &
          generic_metal_scope_flag
        )
          "Target-metal identity"
        else "",
        
        
        if (!marine_estuarine_flag)
          "Marine/estuarine habitat"
        else "",
        
        
        if (!laboratory_flag)
          "Laboratory test"
        else "",
        
        
        if (!waterborne_evidence_flag)
          "Waterborne exposure"
        else "",
        
        
        if (
          mixture_flag |
          metal_count > 1
        )
          "Eligible single-metal exposure arm"
        else "",
        
        
        if (
          lc50_flag &
          !lc50_duration_detected
        )
          "24-96 h LC50 duration"
        else "",
        
        
        if (!numeric_lc50_in_abstract)
          "Usable point LC50"
        else "",
        
        
        if (chemistry_review_flag)
          "Chemical identity/form and concentration basis"
        else "",
        
        
        if (deliberate_modifier_flag)
          "Plain-metal control or unmodified exposure arm"
        else ""
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 28. Define high-confidence preliminary candidates
#
# "Potentially eligible" is deliberately conservative.
# Missing information is sent to Manual Review.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    high_confidence_candidate =
      copepod_core &
      target_metal &
      lc50_flag &
      marine_estuarine_flag &
      laboratory_flag &
      waterborne_evidence_flag &
      lc50_duration_24_96h_flag &
      !explicit_scope_failure &
      !complex_design_flag
  )


# ------------------------------------------------------------
# 29. Assign preliminary screening categories
#
# Order matters.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    screening_category =
      case_when(
        
        # Secondary literature
        
        secondary_source_flag &
          copepod_anywhere ~
          "Secondary source - citation audit",
        
        
        # Explicit freshwater-only scope failure
        
        copepod_anywhere &
          metal_scope_relevant &
          freshwater_only_flag ~
          "Likely outside thesis scope",
        
        
        # Explicit LC50 duration outside 24-96 h
        
        copepod_anywhere &
          metal_scope_relevant &
          explicit_lc50_duration_failure ~
          "Likely outside thesis scope",
        
        
        # Explicit sediment-only exposure matrix
        
        copepod_anywhere &
          metal_scope_relevant &
          sediment_matrix_only_flag ~
          "Likely outside thesis scope",
        
        
        # High-confidence title/abstract candidate
        
        high_confidence_candidate ~
          "Potentially eligible",
        
        
        # Acute lethal evidence, but one or more criteria
        # require verification
        
        copepod_anywhere &
          metal_scope_relevant &
          acute_lethal_relevance &
          !explicit_scope_failure ~
          "Manual review required",
        
        
        # Copepod + target/generic metal relevance but no
        # explicit acute lethal LC50 evidence
        
        copepod_anywhere &
          metal_scope_relevant &
          !acute_lethal_relevance ~
          "Contextual relevance only",
        
        
        # No copepod relevance
        
        !copepod_anywhere ~
          "Not relevant",
        
        
        # No target or generic metal relevance
        
        !metal_scope_relevant ~
          "Not relevant",
        
        
        TRUE ~
          "Manual review required"
      ),
    
    
    screening_priority =
      case_when(
        
        screening_category ==
          "Potentially eligible" ~
          1L,
        
        screening_category ==
          "Manual review required" ~
          2L,
        
        screening_category ==
          "Secondary source - citation audit" ~
          3L,
        
        screening_category ==
          "Contextual relevance only" ~
          4L,
        
        screening_category ==
          "Likely outside thesis scope" ~
          5L,
        
        TRUE ~
          6L
      )
  )


# ------------------------------------------------------------
# 30. Construct preliminary screening reasons
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  rowwise() %>%
  mutate(
    
    preliminary_reason =
      case_when(
        
        screening_category ==
          "Potentially eligible" ~
          paste0(
            "Main thesis-scope criteria are supported by the available ",
            "title/abstract metadata; final eligibility still requires ",
            "source verification."
          ),
        
        
        screening_category ==
          "Manual review required" ~
          if_else(
            review_flags == "",
            paste0(
              "Potentially relevant acute metal-toxicity evidence identified, ",
              "but at least one eligibility criterion requires verification."
            ),
            review_flags
          ),
        
        
        screening_category ==
          "Secondary source - citation audit" ~
          paste0(
            "Secondary synthesis rather than a primary toxicity study; ",
            "retain for citation and coverage audit."
          ),
        
        
        screening_category ==
          "Contextual relevance only" ~
          paste0(
            "Copepod and metal relevance identified, but no clear acute ",
            "mortality LC50 evidence was detected in the available metadata."
          ),
        
        
        screening_category ==
          "Likely outside thesis scope" &
          freshwater_only_flag ~
          "Freshwater-only evidence identified.",
        
        
        screening_category ==
          "Likely outside thesis scope" &
          explicit_lc50_duration_failure ~
          paste0(
            "Only LC50 duration(s) outside the 24-96 h thesis window ",
            "were identified."
          ),
        
        
        screening_category ==
          "Likely outside thesis scope" &
          sediment_matrix_only_flag ~
          paste0(
            "The available metadata indicate a sediment exposure matrix ",
            "rather than a waterborne exposure."
          ),
        
        
        !copepod_anywhere ~
          paste0(
            "No copepod relevance identified in title, ",
            "author keywords, or abstract."
          ),
        
        
        !metal_scope_relevant ~
          paste0(
            "No target-metal or generic metal-toxicity relevance ",
            "identified in title, author keywords, or abstract."
          ),
        
        
        TRUE ~
          "Manual verification required."
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 31. Add blank fields for later full-text screening
#
# These are intentionally left blank at this stage.
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  mutate(
    
    full_text_decision = "",
    
    full_text_exclusion_reason = "",
    
    full_text_notes = ""
  )


# ------------------------------------------------------------
# 32. Arrange records by screening priority
# ------------------------------------------------------------

wos_screen <- wos_screen %>%
  arrange(
    screening_priority,
    desc(publication_year),
    article_title
  )


# ------------------------------------------------------------
# 33. Create main screening export
# ------------------------------------------------------------

screening_export <- wos_screen %>%
  select(
    
    record_id,
    wos_record_id,
    
    screening_priority,
    screening_category,
    preliminary_reason,
    review_flags,
    criteria_to_verify,
    
    full_text_decision,
    full_text_exclusion_reason,
    full_text_notes,
    
    authors,
    article_title,
    publication_year,
    source_title,
    document_type,
    doi,
    
    author_keywords,
    keywords_plus,
    abstract,
    
    copepod_relevance,
    copepod_core,
    copepod_anywhere,
    
    metals_detected,
    metal_count,
    target_metal,
    generic_metal_scope_flag,
    metal_scope_relevant,
    
    copper_flag,
    cadmium_flag,
    zinc_flag,
    nickel_flag,
    lead_flag,
    silver_flag,
    chromium_flag,
    mercury_flag,
    
    endpoint_status,
    lc50_flag,
    mortality_flag,
    acute_flag,
    acute_lethal_relevance,
    numeric_lc50_in_abstract,
    censored_lc50_flag,
    
    habitat_status,
    marine_estuarine_flag,
    freshwater_flag,
    freshwater_only_flag,
    
    laboratory_status,
    laboratory_flag,
    
    all_durations_hours,
    lc50_durations_hours,
    lc50_duration_detected,
    lc50_duration_24_96h_flag,
    lc50_duration_status,
    
    exposure_route_status,
    explicit_waterborne_flag,
    aqueous_concentration_unit_flag,
    waterborne_evidence_flag,
    dietary_exposure_flag,
    
    sediment_mention_flag,
    sediment_matrix_flag,
    sediment_matrix_only_flag,
    
    exposure_structure,
    mixture_flag,
    isolated_single_arm_flag,
    
    deliberate_modifier_flag,
    
    nanoparticle_flag,
    zinc_pyrithione_flag,
    chromium_iii_flag,
    chromium_vi_flag,
    leaching_flag,
    chemistry_review_flag,
    
    secondary_source_flag,
    
    explicit_lc50_duration_failure,
    explicit_scope_failure,
    complex_design_flag,
    high_confidence_candidate,
    
    normalized_title
  )


# ------------------------------------------------------------
# 34. Create category subsets
# ------------------------------------------------------------

potentially_eligible <- screening_export %>%
  filter(
    screening_category ==
      "Potentially eligible"
  )


manual_review <- screening_export %>%
  filter(
    screening_category ==
      "Manual review required"
  )


secondary_sources <- screening_export %>%
  filter(
    screening_category ==
      "Secondary source - citation audit"
  )


contextual_relevance <- screening_export %>%
  filter(
    screening_category ==
      "Contextual relevance only"
  )


likely_outside_scope <- screening_export %>%
  filter(
    screening_category ==
      "Likely outside thesis scope"
  )


not_relevant <- screening_export %>%
  filter(
    screening_category ==
      "Not relevant"
  )


# ------------------------------------------------------------
# 35. Create overall screening summary
# ------------------------------------------------------------

screening_summary <- screening_export %>%
  count(
    screening_priority,
    screening_category,
    name = "records"
  ) %>%
  arrange(
    screening_priority
  )


# ------------------------------------------------------------
# 36. Create metal-level screening summary
# ------------------------------------------------------------

metal_summary <- screening_export %>%
  select(
    record_id,
    screening_category,
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
    names_to = "metal",
    values_to = "detected"
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
    screening_category,
    name = "records"
  ) %>%
  pivot_wider(
    names_from = screening_category,
    values_from = records,
    values_fill = 0
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
# 37. Create publication-year summary
# ------------------------------------------------------------

year_summary <- screening_export %>%
  filter(
    screening_category %in%
      c(
        "Potentially eligible",
        "Manual review required"
      )
  ) %>%
  count(
    publication_year,
    screening_category,
    name = "records"
  ) %>%
  arrange(
    desc(publication_year),
    screening_category
  )


# ------------------------------------------------------------
# 38. Create generic-metal review subset
#
# These records mention generic trace/heavy-metal toxicity but
# do not explicitly identify one of the eight target metals in
# title, author keywords, or abstract.
# ------------------------------------------------------------

generic_metal_review <- screening_export %>%
  filter(
    generic_metal_scope_flag &
      !target_metal
  ) %>%
  select(
    record_id,
    screening_category,
    article_title,
    publication_year,
    source_title,
    author_keywords,
    abstract,
    review_flags,
    criteria_to_verify
  )


# ------------------------------------------------------------
# 39. Create screening-criteria documentation
# ------------------------------------------------------------

criteria_table <- tibble(
  
  criterion = c(
    "Taxonomy",
    "Target metals",
    "Generic metal terminology",
    "Test setting",
    "Exposure environment",
    "Exposure route",
    "Exposure structure",
    "Biological response",
    "Endpoint",
    "Exposure duration",
    "LC50 value",
    "Sediment exposure",
    "Chemical interpretation",
    "Deliberate modifiers",
    "Environmental variation",
    "Keywords Plus"
  ),
  
  
  thesis_scope = c(
    "Copepods, including Calanoida, Harpacticoida, and Cyclopoida",
    "Cu, Cd, Zn, Ni, Pb, Ag, Cr, and Hg",
    "Generic trace-metal/heavy-metal records require metal-identity verification",
    "Laboratory toxicity tests",
    "Marine or estuarine salt/brackish-water systems",
    "Waterborne exposure",
    "Single-metal exposure",
    "Mortality",
    "LC50",
    "24-96 h LC50",
    "Usable point LC50 estimate",
    "Sediment exposure is outside the waterborne Primary scope",
    "Chemical identity and concentration basis must be interpretable",
    "Deliberate co-stressor or modifier treatments require separate assessment",
    "Ordinary variation in temperature and salinity is not an automatic exclusion criterion",
    "Not used for automatic screening decisions"
  ),
  
  
  preliminary_screening_treatment = c(
    "Title, author keywords, and abstract are screened; uncertain cases are retained for review",
    "All eight target metals are screened independently",
    "Generic metal wording can trigger manual review but cannot establish final metal eligibility",
    "Explicit laboratory evidence is recorded; missing information is sent to manual review",
    "Marine/estuarine evidence is recorded; ambiguous habitat is sent to manual review",
    "Explicit waterborne wording and aqueous concentration units are treated as waterborne evidence",
    "Multi-metal and mixture papers are retained for review because separate single-metal arms may exist",
    "Mortality, lethal, and survival terminology is screened",
    "LC50, LC 50, LC-50, and LC(50) variants are recognised",
    "Only durations linked to LC50-containing text units are used for duration-based scope decisions",
    "Numeric LC50 values are flagged where visible; source verification remains required",
    "Strong sediment-matrix evidence is separated from simple contextual mentions of sediment",
    "Nanoparticles, zinc pyrithione, chromium forms, leaching, and related issues are flagged for review",
    "Explicit co-stressors or experimental modifiers are sent to manual review rather than automatically discarded",
    "Salinity and temperature remain part of experimental context unless combined with a separately identified co-stressor treatment",
    "Keywords Plus is preserved in the workbook but excluded from screening_text and automated eligibility decisions"
  )
)


# ------------------------------------------------------------
# 40. Classification integrity check
# ------------------------------------------------------------

classification_check <-
  nrow(potentially_eligible) +
  nrow(manual_review) +
  nrow(secondary_sources) +
  nrow(contextual_relevance) +
  nrow(likely_outside_scope) +
  nrow(not_relevant)


if (
  classification_check !=
  nrow(screening_export)
) {
  
  stop(
    paste(
      "Classification check failed:",
      classification_check,
      "records classified out of",
      nrow(screening_export)
    )
  )
}


# Frozen candidate benchmark used in the thesis search-flow summary.
candidate_n <-
  nrow(potentially_eligible) +
  nrow(manual_review)

stopifnot(
  nrow(screening_export) == 117L,
  candidate_n == 75L
)


# ------------------------------------------------------------
# 41. Export screening workbook
# ------------------------------------------------------------

write_xlsx(
  list(
    
    Screening_Summary =
      screening_summary,
    
    Metal_Summary =
      metal_summary,
    
    Year_Summary =
      year_summary,
    
    Screening_Criteria =
      criteria_table,
    
    Potentially_Eligible =
      potentially_eligible,
    
    Manual_Review =
      manual_review,
    
    Generic_Metal_Review =
      generic_metal_review,
    
    Secondary_Sources =
      secondary_sources,
    
    Contextual_Relevance =
      contextual_relevance,
    
    Likely_Outside_Scope =
      likely_outside_scope,
    
    Not_Relevant =
      not_relevant,
    
    All_Records =
      screening_export
  ),
  
  output_file
)


screening_count_summary <- tibble::tibble(
  Category = c(
    "Total WoS records",
    "Potentially eligible",
    "Manual review required",
    "Candidate records for ECOTOX overlap audit",
    "Secondary sources",
    "Contextual relevance only",
    "Likely outside thesis scope",
    "Not relevant",
    "Generic-metal identity review"
  ),
  N = c(
    nrow(screening_export),
    nrow(potentially_eligible),
    nrow(manual_review),
    candidate_n,
    nrow(secondary_sources),
    nrow(contextual_relevance),
    nrow(likely_outside_scope),
    nrow(not_relevant),
    nrow(generic_metal_review)
  )
)

readr::write_csv(
  screening_count_summary,
  file.path(
    output_root,
    "screening_count_summary.csv"
  )
)

search_metadata <- tibble::tibble(
  Field = c(
    "Database",
    "Search date",
    "Search field",
    "Query",
    "Keywords Plus screening treatment"
  ),
  Value = c(
    "Web of Science Core Collection",
    "2026-09-05",
    "Topic",
    paste0(
      'copepod* AND ("LC50" OR "LC 50" OR "median lethal concentration") AND ',
      '(copper OR cadmium OR zinc OR nickel OR lead OR silver OR chromium OR mercury)'
    ),
    paste(
      "Included in the Topic search but not used in the subsequent",
      "title/abstract/author-keyword screening"
    )
  )
)

readr::write_csv(
  search_metadata,
  file.path(
    output_root,
    "search_metadata.csv"
  )
)

input_output_manifest <- tibble::tibble(
  Role = c(
    "Local input: raw WoS export",
    "Local output: detailed screening workbook",
    "Output: screening count summary",
    "Output: search metadata"
  ),
  Path = c(
    wos_file,
    output_file,
    file.path(
      output_root,
      "screening_count_summary.csv"
    ),
    file.path(
      output_root,
      "search_metadata.csv"
    )
  )
)

readr::write_csv(
  input_output_manifest,
  file.path(
    output_root,
    "input_output_manifest.csv"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_root,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "STATUS: WOS PRELIMINARY SCREENING = PASS",
    paste0(
      "Raw WoS records: ",
      nrow(screening_export)
    ),
    paste0(
      "Candidate records for ECOTOX overlap audit: ",
      candidate_n
    ),
    "Detailed screening workbook is local-only and should remain gitignored."
  ),
  file.path(
    output_root,
    "RUN_COMPLETE.txt"
  )
)


# ------------------------------------------------------------
# 42. Final console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "WEB OF SCIENCE PRELIMINARY SCREENING\n"
)

cat(
  "============================================\n"
)


cat(
  "\nTotal WoS records:",
  nrow(screening_export),
  "\n"
)


cat(
  "Potentially eligible:",
  nrow(potentially_eligible),
  "\n"
)


cat(
  "Manual review required:",
  nrow(manual_review),
  "\n"
)


cat(
  "Candidate records for ECOTOX overlap audit:",
  candidate_n,
  "\n"
)


cat(
  "Secondary sources:",
  nrow(secondary_sources),
  "\n"
)


cat(
  "Contextual relevance only:",
  nrow(contextual_relevance),
  "\n"
)


cat(
  "Likely outside thesis scope:",
  nrow(likely_outside_scope),
  "\n"
)


cat(
  "Not relevant:",
  nrow(not_relevant),
  "\n"
)


cat(
  "Generic-metal records requiring identity review:",
  nrow(generic_metal_review),
  "\n"
)


cat(
  "\nClassification check:",
  classification_check,
  "records classified out of",
  nrow(screening_export),
  "\n"
)


cat(
  "\nDetailed local screening workbook saved as:",
  output_file,
  "\n"
)

cat(
  "STATUS: WOS PRELIMINARY SCREENING = PASS\n"
)


cat(
  "\nIMPORTANT:\n"
)


cat(
  paste0(
    "This workbook represents preliminary title/abstract screening only. ",
    "Potentially eligible and manual-review records have not yet been ",
    "compared with ECOTOX and have not yet received final full-text ",
    "eligibility decisions.\n"
  )
)

