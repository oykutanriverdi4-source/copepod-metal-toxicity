# ============================================================
# 03_prepare_wos_unique_reference_review.R
# ============================================================
# Purpose
# Prepare the 47 Web of Science publications not matched to the
# ECOTOX retrieval for structured full-text eligibility review.
#
# This is a REFERENCE-level triage only.
# It does NOT:
# - make final eligibility decisions,
# - extract individual LC50 Results,
# - decide relevance from the later Cu-Cd-Zn modelling subset,
# - replace review of the original publication.
#
# Final thesis workflow benchmark
# - 75 candidate WoS publications entered the overlap audit
# - 28 were represented in the ECOTOX retrieval after manual
#   verification of four high-confidence fuzzy matches
# - 47 publications remained for full-text review
# - after original-publication review and harmonization,
#   17 publications contributed 57 quantitative LC50 Results
#   to the final combined dataset
#
# Inputs
# outputs/database_audits/wos/02_check_wos_ecotox_overlap/
#   wos_ecotox_overlap.xlsx
# data/external/wos/wos_raw_export.xls
#
# Local detailed output
# outputs/database_audits/wos/03_prepare_wos_unique_reference_review/
#   wos_unique_reference_review.xlsx
#
# The detailed workbook contains WoS bibliographic metadata and
# blank/manual review fields. It should remain local/gitignored.
# Compact count summaries, review criteria and run metadata can be
# retained in the public reproducibility archive.
#
# Scientific rule
# The original triage criteria and category logic are preserved.
# Repository cleanup changes paths, naming, documentation and
# reproducibility checks only.
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(writexl)
library(readr)


# ------------------------------------------------------------
# 2. Define repository file paths
# ------------------------------------------------------------

overlap_file <- file.path(
  "outputs",
  "database_audits",
  "wos",
  "02_check_wos_ecotox_overlap",
  "wos_ecotox_overlap.xlsx"
)

wos_raw_file <- file.path(
  "data",
  "external",
  "wos",
  "wos_raw_export.xls"
)

output_root <- file.path(
  "outputs",
  "database_audits",
  "wos",
  "03_prepare_wos_unique_reference_review"
)

if (!dir.exists(output_root)) {
  dir.create(
    output_root,
    recursive = TRUE
  )
}

output_file <- file.path(
  output_root,
  "wos_unique_reference_review.xlsx"
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
# 3. Check required files
# ------------------------------------------------------------

if (!file.exists(overlap_file)) {
  
  stop(
    paste0("WoS overlap workbook not found: ", overlap_file)
  )
}


if (!file.exists(wos_raw_file)) {
  
  stop(
    paste0(
      "Local WoS raw export not found: ",
      wos_raw_file,
      ". The raw export is intentionally not redistributed in the public repository."
    )
  )
}


# ------------------------------------------------------------
# 4. Check required input sheet
# ------------------------------------------------------------

overlap_sheets <- excel_sheets(
  overlap_file
)


if (
  !"Provisional_Unique_WoS" %in%
  overlap_sheets
) {
  
  stop(
    paste0(
      "Sheet 'Provisional_Unique_WoS' was not found in ",
      overlap_file,
      "."
    )
  )
}


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


as_flag <- function(x) {
  
  if (is.logical(x)) {
    
    x[is.na(x)] <- FALSE
    
    return(x)
  }
  
  
  x_chr <- str_to_lower(
    str_trim(
      as.character(x)
    )
  )
  
  
  x_chr %in%
    c(
      "true",
      "t",
      "1",
      "yes",
      "y"
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


safe_text <- function(x) {
  
  x <- as.character(x)
  
  x[is.na(x)] <- ""
  
  x
}


# ------------------------------------------------------------
# 6. Import provisional unique WoS references
# ------------------------------------------------------------

unique_wos <- read_excel(
  overlap_file,
  sheet = "Provisional_Unique_WoS",
  .name_repair = "unique"
)


cat(
  "\nProvisional unique WoS references imported:",
  nrow(unique_wos),
  "\n"
)


if (nrow(unique_wos) == 0) {
  
  stop(
    "No provisional unique WoS references were found."
  )
}


stopifnot(
  nrow(unique_wos) == 47L
)


# ------------------------------------------------------------
# 7. Check essential columns
# ------------------------------------------------------------

required_columns <- c(
  "record_id",
  "screening_category",
  "article_title",
  "publication_year",
  "authors",
  "source_title",
  "abstract",
  "copepod_core",
  "copepod_anywhere",
  "metals_detected",
  "metal_count",
  "target_metal",
  "generic_metal_scope_flag",
  "metal_scope_relevant",
  "lc50_flag",
  "acute_lethal_relevance",
  "numeric_lc50_in_abstract",
  "censored_lc50_flag",
  "marine_estuarine_flag",
  "freshwater_only_flag",
  "laboratory_flag",
  "lc50_duration_detected",
  "lc50_duration_24_96h_flag",
  "waterborne_evidence_flag",
  "dietary_exposure_flag",
  "sediment_matrix_flag",
  "sediment_matrix_only_flag",
  "mixture_flag",
  "deliberate_modifier_flag",
  "chemistry_review_flag",
  "explicit_lc50_duration_failure",
  "explicit_scope_failure",
  "complex_design_flag",
  "high_confidence_candidate",
  "ecotox_sampling_frame_status",
  "provisional_unique_wos"
)


missing_columns <- setdiff(
  required_columns,
  names(unique_wos)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Required column(s) missing from overlap workbook:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# 8. Convert imported flag columns to logical values
# ------------------------------------------------------------

flag_columns <- c(
  "copepod_core",
  "copepod_anywhere",
  "target_metal",
  "generic_metal_scope_flag",
  "metal_scope_relevant",
  
  "copper_flag",
  "cadmium_flag",
  "zinc_flag",
  "nickel_flag",
  "lead_flag",
  "silver_flag",
  "chromium_flag",
  "mercury_flag",
  
  "lc50_flag",
  "mortality_flag",
  "acute_flag",
  "acute_lethal_relevance",
  "numeric_lc50_in_abstract",
  "censored_lc50_flag",
  
  "marine_estuarine_flag",
  "freshwater_flag",
  "freshwater_only_flag",
  
  "laboratory_flag",
  
  "lc50_duration_detected",
  "lc50_duration_24_96h_flag",
  
  "explicit_waterborne_flag",
  "aqueous_concentration_unit_flag",
  "waterborne_evidence_flag",
  "dietary_exposure_flag",
  
  "sediment_mention_flag",
  "sediment_matrix_flag",
  "sediment_matrix_only_flag",
  
  "mixture_flag",
  "isolated_single_arm_flag",
  
  "deliberate_modifier_flag",
  
  "nanoparticle_flag",
  "zinc_pyrithione_flag",
  "chromium_iii_flag",
  "chromium_vi_flag",
  "leaching_flag",
  "chemistry_review_flag",
  
  "secondary_source_flag",
  
  "explicit_lc50_duration_failure",
  "explicit_scope_failure",
  "complex_design_flag",
  "high_confidence_candidate",
  
  "confirmed_exact_match",
  "possible_fuzzy_match",
  "provisional_unique_wos"
)


flag_columns <- intersect(
  flag_columns,
  names(unique_wos)
)


unique_wos <- unique_wos %>%
  mutate(
    across(
      all_of(flag_columns),
      as_flag
    )
  )


# ------------------------------------------------------------
# 9. Standardize selected numeric fields
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  mutate(
    
    record_id =
      suppressWarnings(
        as.integer(record_id)
      ),
    
    publication_year =
      suppressWarnings(
        as.integer(publication_year)
      ),
    
    metal_count =
      suppressWarnings(
        as.integer(metal_count)
      )
  )


# ------------------------------------------------------------
# 10. Confirm that the input contains provisional unique records
# ------------------------------------------------------------

if (
  any(
    !unique_wos$provisional_unique_wos,
    na.rm = TRUE
  )
) {
  
  warning(
    paste0(
      "At least one row in 'Provisional_Unique_WoS' is not ",
      "flagged as provisional_unique_wos = TRUE."
    )
  )
}


if (
  "possible_fuzzy_match" %in%
  names(unique_wos) &&
  any(
    unique_wos$possible_fuzzy_match,
    na.rm = TRUE
  )
) {
  
  warning(
    paste0(
      "At least one possible fuzzy match remains in the ",
      "provisional unique sheet. Review the overlap audit."
    )
  )
}


# ------------------------------------------------------------
# 11. Recover DOI values from the raw WoS export
#
# The preliminary script used DOI Link where available.
# The raw WoS export also contains a plain-text DOI field.
# record_id corresponds to the row order in the raw WoS export.
# ------------------------------------------------------------

wos_raw <- suppressWarnings(
  read_excel(
    wos_raw_file,
    col_types = "text",
    .name_repair = "unique"
  )
)


raw_doi_vector <- get_optional_column(
  wos_raw,
  c(
    "DOI"
  )
)


raw_doi_index <- tibble(
  
  record_id =
    seq_len(
      nrow(wos_raw)
    ),
  
  recovered_wos_doi =
    extract_doi(
      raw_doi_vector
    )
)


unique_wos <- unique_wos %>%
  left_join(
    raw_doi_index,
    by = "record_id"
  ) %>%
  mutate(
    
    existing_doi =
      extract_doi(
        safe_text(doi)
      ),
    
    best_doi =
      coalesce(
        recovered_wos_doi,
        existing_doi
      ),
    
    doi_url =
      if_else(
        !is.na(best_doi) &
          best_doi != "",
        paste0(
          "https://doi.org/",
          best_doi
        ),
        ""
      )
  )


# ------------------------------------------------------------
# 12. Calculate transparent core-criteria support
#
# Seven metadata criteria are counted:
# 1. Copepod relevance
# 2. Target-metal identity
# 3. LC50 endpoint
# 4. Marine/estuarine context
# 5. Laboratory evidence
# 6. Waterborne evidence
# 7. 24-96 h LC50 duration
#
# This count is used ONLY to prioritize full-text review.
# It is NOT an eligibility score.
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  rowwise() %>%
  mutate(
    
    core_criteria_supported_n =
      sum(
        c(
          copepod_anywhere,
          target_metal,
          lc50_flag,
          marine_estuarine_flag,
          laboratory_flag,
          waterborne_evidence_flag,
          lc50_duration_24_96h_flag
        ),
        na.rm = TRUE
      ),
    
    core_criteria_missing_n =
      7L -
      core_criteria_supported_n
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 13. Define explicit outside-scope reason
#
# These are deliberately narrow.
# Missing information is NOT treated as an exclusion.
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  rowwise() %>%
  mutate(
    
    clear_outside_scope_reason =
      collapse_flags(
        
        if (freshwater_only_flag)
          "Freshwater-only evidence"
        else "",
        
        
        if (sediment_matrix_only_flag)
          "Sediment-only exposure matrix"
        else "",
        
        
        if (explicit_lc50_duration_failure)
          "Only LC50 duration(s) outside the 24-96 h scope"
        else ""
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 14. Define complex-design reason
#
# These records are NOT excluded.
# They are separated because an eligible plain single-metal
# arm may coexist with an ineligible or secondary treatment.
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  rowwise() %>%
  mutate(
    
    complex_design_reason =
      collapse_flags(
        
        if (mixture_flag)
          "Mixture or combined-exposure design"
        else "",
        
        
        if (dietary_exposure_flag)
          "Dietary exposure component"
        else "",
        
        
        if (
          sediment_matrix_flag &
          !sediment_matrix_only_flag
        )
          "Sediment-related component alongside other exposure evidence"
        else "",
        
        
        if (deliberate_modifier_flag)
          "Deliberate modifier or co-stressor treatment"
        else "",
        
        
        if (nanoparticle_flag)
          "Nanoparticle exposure"
        else "",
        
        
        if (zinc_pyrithione_flag)
          "Zinc pyrithione exposure"
        else "",
        
        
        if (chromium_iii_flag)
          "Cr(III) chemistry"
        else "",
        
        
        if (chromium_vi_flag)
          "Cr(VI) chemistry / concentration-basis issue"
        else "",
        
        
        if (leaching_flag)
          "Leaching exposure"
        else "",
        
        
        if (censored_lc50_flag)
          "Potentially censored LC50"
        else ""
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 15. Construct targeted full-text checklist
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  rowwise() %>%
  mutate(
    
    full_text_checklist =
      collapse_flags(
        
        "Confirm source-level study details and identify every candidate LC50 Result",
        
        
        if (!copepod_core)
          "Confirm that the tested organism is an eligible copepod"
        else "",
        
        
        if (!target_metal)
          "Confirm that at least one of the eight target metals was tested"
        else "",
        
        
        if (!marine_estuarine_flag)
          "Confirm marine/estuarine saltwater or brackish-water exposure"
        else "",
        
        
        if (!laboratory_flag)
          "Confirm laboratory test setting"
        else "",
        
        
        if (!waterborne_evidence_flag)
          "Confirm waterborne exposure route"
        else "",
        
        
        if (
          mixture_flag |
          metal_count > 1
        )
          "Identify separate single-metal exposure arm(s)"
        else "",
        
        
        if (!lc50_flag)
          "Confirm mortality LC50 endpoint"
        else "",
        
        
        if (!lc50_duration_24_96h_flag)
          "Confirm that candidate LC50 duration is within 24-96 h"
        else "",
        
        
        if (!numeric_lc50_in_abstract)
          "Locate usable point LC50 estimate in full text, table, figure, or supplement"
        else "",
        
        
        if (chemistry_review_flag)
          "Confirm chemical form and concentration basis"
        else "",
        
        
        if (deliberate_modifier_flag)
          "Separate plain exposure arm from deliberate modifier treatment(s)"
        else "",
        
        
        if (censored_lc50_flag)
          "Confirm whether an uncensored point LC50 is available"
        else ""
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 16. Assign reference-level triage categories
#
# Priority:
# 1. High-priority full-text candidate
# 2. Manual eligibility review
# 3. Modifier / complex-design candidate
# 4. Clear outside scope from metadata
#
# IMPORTANT:
# No category here is a final eligibility decision.
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  mutate(
    
    triage_category =
      case_when(
        
        # Explicit, narrow scope failures
        
        explicit_scope_failure |
          freshwater_only_flag |
          sediment_matrix_only_flag |
          explicit_lc50_duration_failure ~
          "Clear outside scope from metadata",
        
        
        # Complex studies may still contain eligible plain arms
        
        complex_design_flag |
          mixture_flag |
          dietary_exposure_flag |
          sediment_matrix_flag |
          deliberate_modifier_flag |
          chemistry_review_flag |
          censored_lc50_flag ~
          "Modifier / complex-design candidate",
        
        
        # Strong metadata support and no complex-design flags
        
        target_metal &
          lc50_flag &
          acute_lethal_relevance &
          core_criteria_supported_n >= 6 &
          !explicit_scope_failure ~
          "High-priority full-text candidate",
        
        
        # Everything else remains unresolved rather than excluded
        
        TRUE ~
          "Manual eligibility review"
      ),
    
    
    triage_priority =
      case_when(
        
        triage_category ==
          "High-priority full-text candidate" ~
          1L,
        
        triage_category ==
          "Manual eligibility review" ~
          2L,
        
        triage_category ==
          "Modifier / complex-design candidate" ~
          3L,
        
        triage_category ==
          "Clear outside scope from metadata" ~
          4L,
        
        TRUE ~
          5L
      )
  )


# ------------------------------------------------------------
# 17. Construct triage reasons
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  rowwise() %>%
  mutate(
    
    triage_reason =
      case_when(
        
        triage_category ==
          "High-priority full-text candidate" ~
          paste0(
            "Strong title/abstract support for the thesis scope ",
            "(",
            core_criteria_supported_n,
            "/7 core metadata criteria supported) and no major ",
            "complex-design flag. Full text should be checked first."
          ),
        
        
        triage_category ==
          "Manual eligibility review" ~
          paste0(
            "Potentially relevant evidence remains, but one or more ",
            "core eligibility elements require source verification. ",
            "Supported metadata criteria: ",
            core_criteria_supported_n,
            "/7."
          ),
        
        
        triage_category ==
          "Modifier / complex-design candidate" ~
          if_else(
            complex_design_reason != "",
            paste0(
              "Potentially eligible evidence may coexist with a ",
              "complex or secondary treatment: ",
              complex_design_reason,
              "."
            ),
            paste0(
              "Complex study design requires source-level separation ",
              "of potentially eligible exposure arms."
            )
          ),
        
        
        triage_category ==
          "Clear outside scope from metadata" ~
          if_else(
            clear_outside_scope_reason != "",
            clear_outside_scope_reason,
            "Explicit title/abstract metadata indicate a scope mismatch."
          ),
        
        
        TRUE ~
          "Manual source review required."
      )
    
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 18. Add blank manual full-text decision fields
#
# A reference can contain several Tests/Results.
# Therefore a paper-level eligibility decision only indicates
# whether eligible evidence exists in the source.
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  mutate(
    
    full_text_checked = "",
    
    full_text_source = "",
    
    full_text_decision = "",
    
    full_text_exclusion_reason = "",
    
    eligible_metals = "",
    
    species = "",
    
    taxonomic_order = "",
    
    life_stage = "",
    
    laboratory_confirmed = "",
    
    marine_estuarine_confirmed = "",
    
    waterborne_confirmed = "",
    
    single_metal_arm_confirmed = "",
    
    eligible_duration_hours = "",
    
    usable_point_lc50 = "",
    
    chemical_basis_interpretable = "",
    
    modifier_status = "",
    
    number_of_usable_lc50_results = "",
    
    result_level_extraction_needed = "",
    
    analytical_relevance = "",
    
    full_text_notes = ""
  )


# ------------------------------------------------------------
# 19. Arrange review order
# ------------------------------------------------------------

unique_wos <- unique_wos %>%
  arrange(
    triage_priority,
    desc(core_criteria_supported_n),
    desc(publication_year),
    article_title
  )


# ------------------------------------------------------------
# 20. Create compact review workbook columns
# ------------------------------------------------------------

review_export <- unique_wos %>%
  select(
    
    # Review control
    
    record_id,
    wos_record_id,
    
    triage_priority,
    triage_category,
    triage_reason,
    
    core_criteria_supported_n,
    core_criteria_missing_n,
    
    clear_outside_scope_reason,
    complex_design_reason,
    full_text_checklist,
    
    
    # Manual full-text fields
    
    full_text_checked,
    full_text_source,
    full_text_decision,
    full_text_exclusion_reason,
    
    eligible_metals,
    species,
    taxonomic_order,
    life_stage,
    
    laboratory_confirmed,
    marine_estuarine_confirmed,
    waterborne_confirmed,
    single_metal_arm_confirmed,
    
    eligible_duration_hours,
    usable_point_lc50,
    chemical_basis_interpretable,
    modifier_status,
    
    number_of_usable_lc50_results,
    result_level_extraction_needed,
    analytical_relevance,
    full_text_notes,
    
    
    # Bibliographic information
    
    authors,
    article_title,
    publication_year,
    source_title,
    document_type,
    
    best_doi,
    doi_url,
    
    author_keywords,
    keywords_plus,
    abstract,
    
    
    # Previous screening information
    
    screening_category,
    preliminary_reason,
    review_flags,
    criteria_to_verify,
    
    
    # Biological / chemical scope
    
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
    
    
    # Endpoint / duration
    
    endpoint_status,
    lc50_flag,
    acute_lethal_relevance,
    numeric_lc50_in_abstract,
    censored_lc50_flag,
    
    all_durations_hours,
    lc50_durations_hours,
    lc50_duration_detected,
    lc50_duration_24_96h_flag,
    lc50_duration_status,
    
    
    # Environment / exposure
    
    habitat_status,
    marine_estuarine_flag,
    freshwater_flag,
    freshwater_only_flag,
    
    laboratory_status,
    laboratory_flag,
    
    exposure_route_status,
    explicit_waterborne_flag,
    aqueous_concentration_unit_flag,
    waterborne_evidence_flag,
    dietary_exposure_flag,
    
    sediment_mention_flag,
    sediment_matrix_flag,
    sediment_matrix_only_flag,
    
    
    # Exposure complexity
    
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
    
    
    # Previous scope / overlap audit
    
    explicit_lc50_duration_failure,
    explicit_scope_failure,
    complex_design_flag,
    high_confidence_candidate,
    
    ecotox_sampling_frame_status,
    provisional_unique_wos,
    overlap_audit_note,
    
    fuzzy_ecotox_reference_number,
    fuzzy_ecotox_title,
    fuzzy_ecotox_year,
    fuzzy_title_similarity
  )


# ------------------------------------------------------------
# 21. Create category-specific review sheets
# ------------------------------------------------------------

high_priority_fulltext <- review_export %>%
  filter(
    triage_category ==
      "High-priority full-text candidate"
  )


manual_eligibility <- review_export %>%
  filter(
    triage_category ==
      "Manual eligibility review"
  )


modifier_complex_design <- review_export %>%
  filter(
    triage_category ==
      "Modifier / complex-design candidate"
  )


clear_outside_scope <- review_export %>%
  filter(
    triage_category ==
      "Clear outside scope from metadata"
  )


full_text_queue <- review_export %>%
  filter(
    triage_category !=
      "Clear outside scope from metadata"
  ) %>%
  arrange(
    triage_priority,
    desc(core_criteria_supported_n),
    desc(publication_year),
    article_title
  )


# ------------------------------------------------------------
# 22. Create overall triage summary
# ------------------------------------------------------------

triage_summary <- review_export %>%
  count(
    triage_priority,
    triage_category,
    name = "references"
  ) %>%
  arrange(
    triage_priority
  ) %>%
  mutate(
    
    percent_of_unique_pool =
      round(
        100 *
          references /
          sum(references),
        1
      )
  )


# ------------------------------------------------------------
# 23. Create triage-by-metal summary
#
# Multi-metal references contribute once to each detected metal.
# ------------------------------------------------------------

triage_by_metal <- review_export %>%
  select(
    
    record_id,
    triage_category,
    
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
    triage_category,
    name = "references"
  ) %>%
  pivot_wider(
    names_from = triage_category,
    values_from = references,
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
# 24. Create unresolved generic-metal summary
# ------------------------------------------------------------

generic_metal_summary <- review_export %>%
  filter(
    generic_metal_scope_flag &
      !target_metal
  ) %>%
  count(
    triage_category,
    name = "references"
  )


# ------------------------------------------------------------
# 25. Create triage-by-year summary
# ------------------------------------------------------------

triage_by_year <- review_export %>%
  count(
    publication_year,
    triage_category,
    name = "references"
  ) %>%
  mutate(
    triage_order =
      factor(
        triage_category,
        levels = c(
          "High-priority full-text candidate",
          "Manual eligibility review",
          "Modifier / complex-design candidate",
          "Clear outside scope from metadata"
        )
      )
  ) %>%
  arrange(
    desc(publication_year),
    triage_order
  ) %>%
  select(
    publication_year,
    triage_category,
    references
  )

# ------------------------------------------------------------
# 26. Create eligibility-criteria documentation
# ------------------------------------------------------------

eligibility_criteria <- tibble(
  
  criterion = c(
    "Taxonomy",
    "Target metals",
    "Environment",
    "Test setting",
    "Exposure route",
    "Exposure structure",
    "Biological response",
    "Endpoint",
    "Exposure duration",
    "LC50 value",
    "Chemical interpretation",
    "Modifier treatment",
    "Reference-level decision",
    "Analytical relevance"
  ),
  
  
  primary_scope = c(
    "Copepod",
    "Cu, Cd, Zn, Ni, Pb, Ag, Cr, or Hg",
    "Marine or estuarine salt/brackish water",
    "Laboratory",
    "Waterborne",
    "At least one separate single-metal exposure arm",
    "Mortality",
    "LC50",
    "24-96 h",
    "Usable point estimate rather than only a censored/range value",
    "Chemical identity and concentration basis sufficiently interpretable",
    "Plain Primary evidence separated from deliberate modifier treatments",
    "At least one eligible LC50 Result in the source is sufficient to move the Reference to result-level extraction",
    "Assessed only after eligible Results are extracted across all eight metals"
  ),
  
  
  triage_rule = c(
    "Missing metadata triggers source verification, not automatic exclusion",
    "Generic trace/heavy-metal wording alone requires metal-identity verification",
    "Freshwater-only evidence can be classified outside the thesis scope",
    "Missing laboratory wording requires source verification",
    "Missing waterborne wording requires source verification",
    "Mixture or multi-metal studies are retained if separate single-metal arms may exist",
    "Mortality or lethal response must be confirmed",
    "LC50 must be confirmed from the source",
    "Only explicit LC50 durations outside 24-96 h constitute a metadata-level duration failure",
    "The full paper, table, figure, or supplement may provide the usable point estimate",
    "Nanoparticles, Cr forms, zinc pyrithione, leaching, and ambiguous bases require source review",
    "Modifier studies are not treated as invalid; plain and modifier evidence are separated",
    "Reference-level eligibility does not equal one Result and does not automatically add data to the thesis dataset",
    "Do not use current Cu-Cd-Zn model membership as an eligibility rule"
  )
)


# ------------------------------------------------------------
# 27. Create manual decision guide
# ------------------------------------------------------------

decision_guide <- tibble(
  
  field = c(
    "full_text_checked",
    "full_text_decision",
    "laboratory_confirmed",
    "marine_estuarine_confirmed",
    "waterborne_confirmed",
    "single_metal_arm_confirmed",
    "usable_point_lc50",
    "chemical_basis_interpretable",
    "modifier_status",
    "result_level_extraction_needed",
    "analytical_relevance"
  ),
  
  
  suggested_values = c(
    "Yes | No",
    
    paste0(
      "Eligible - Primary evidence present | ",
      "Eligible - Primary + Modifier evidence present | ",
      "Eligible - Modifier evidence only | ",
      "Exclude | Unresolved"
    ),
    
    "Yes | No | Unresolved",
    
    "Yes | No | Unresolved",
    
    "Yes | No | Unresolved",
    
    "Yes | No | Unresolved",
    
    "Yes | No | Censored only | Range only | Unresolved",
    
    "Yes | No | Unresolved",
    
    "Plain only | Modifier only | Plain + Modifier | Not applicable | Unresolved",
    
    "Yes | No | Unresolved",
    
    paste0(
      "Leave blank until eligible LC50 Results are extracted. ",
      "Later describe whether the source changes metal support, ",
      "cross-metal connectivity, biological support, duration support, ",
      "or only contextual interpretation."
    )
  ),
  
  
  interpretation = c(
    "Whether the original paper or adequate source material was inspected",
    
    paste0(
      "Reference-level conclusion only. ",
      "A positive decision means at least one potentially usable ",
      "Result exists and should be extracted separately."
    ),
    
    "Confirm laboratory test setting",
    
    "Confirm marine/estuarine saltwater or brackish-water context",
    
    "Confirm the candidate LC50 was generated by waterborne exposure",
    
    "Confirm at least one eligible single-metal arm exists",
    
    "Confirm a usable point LC50 can be recovered",
    
    "Confirm metal form and concentration basis are interpretable for the intended evidence layer",
    
    "Distinguish ordinary Primary evidence from deliberate modifier treatments",
    
    "Set to Yes when at least one eligible Result should be entered into a separate result-level extraction table",
    
    "Do not decide this from the paper title alone"
  )
)


# ------------------------------------------------------------
# 28. Create coverage summary
# ------------------------------------------------------------

coverage_summary <- tibble(
  
  metric = c(
    "Provisional unique WoS References entering triage",
    "High-priority full-text candidates",
    "Manual eligibility reviews",
    "Modifier / complex-design candidates",
    "Clear outside-scope References from metadata",
    "References remaining in the active full-text queue"
  ),
  
  
  value = c(
    nrow(review_export),
    nrow(high_priority_fulltext),
    nrow(manual_eligibility),
    nrow(modifier_complex_design),
    nrow(clear_outside_scope),
    nrow(full_text_queue)
  )
)


# ------------------------------------------------------------
# 29. Integrity checks
# ------------------------------------------------------------

triage_count_check <-
  nrow(high_priority_fulltext) +
  nrow(manual_eligibility) +
  nrow(modifier_complex_design) +
  nrow(clear_outside_scope)


if (
  triage_count_check !=
  nrow(review_export)
) {
  
  stop(
    paste(
      "Triage classification check failed:",
      triage_count_check,
      "records classified out of",
      nrow(review_export)
    )
  )
}


if (
  n_distinct(review_export$record_id) !=
  nrow(review_export)
) {
  
  stop(
    paste0(
      "Duplicate record_id values were detected in the ",
      "provisional unique WoS pool."
    )
  )
}


# ------------------------------------------------------------
# 30. Export workbook
# ------------------------------------------------------------

write_xlsx(
  list(
    
    Coverage_Summary =
      coverage_summary,
    
    Triage_Summary =
      triage_summary,
    
    Triage_By_Metal =
      triage_by_metal,
    
    Triage_By_Year =
      triage_by_year,
    
    Generic_Metal_Summary =
      generic_metal_summary,
    
    Eligibility_Criteria =
      eligibility_criteria,
    
    Decision_Guide =
      decision_guide,
    
    FullText_Queue =
      full_text_queue,
    
    High_Priority_FullText =
      high_priority_fulltext,
    
    Manual_Eligibility =
      manual_eligibility,
    
    Modifier_Complex_Design =
      modifier_complex_design,
    
    Clear_Outside_Scope =
      clear_outside_scope,
    
    All_Unique_Triage =
      review_export
  ),
  
  output_file
)


# ------------------------------------------------------------
# 31. REPOSITORY-SAFE SUMMARY + RUN METADATA
# ------------------------------------------------------------

triage_count_summary <- coverage_summary %>%
  dplyr::rename(
    Metric = metric,
    N = value
  )

readr::write_csv(
  triage_count_summary,
  file.path(
    output_root,
    "triage_count_summary.csv"
  )
)

readr::write_csv(
  eligibility_criteria,
  file.path(
    output_root,
    "eligibility_criteria.csv"
  )
)

readr::write_csv(
  decision_guide,
  file.path(
    output_root,
    "manual_decision_guide.csv"
  )
)

input_output_manifest <- tibble::tibble(
  Role = c(
    "Input: WoS-ECOTOX overlap workbook",
    "Local input: WoS raw export",
    "Local output: detailed unique-reference review workbook",
    "Output: triage count summary",
    "Output: eligibility criteria",
    "Output: manual decision guide"
  ),
  Path = c(
    overlap_file,
    wos_raw_file,
    output_file,
    file.path(
      output_root,
      "triage_count_summary.csv"
    ),
    file.path(
      output_root,
      "eligibility_criteria.csv"
    ),
    file.path(
      output_root,
      "manual_decision_guide.csv"
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
    "STATUS: WOS UNIQUE-REFERENCE TRIAGE = PASS",
    paste0(
      "Publications entering triage: ",
      nrow(review_export)
    ),
    paste0(
      "Active full-text review queue: ",
      nrow(full_text_queue)
    ),
    paste0(
      "Clear outside-scope from metadata: ",
      nrow(clear_outside_scope)
    ),
    "This triage does not make final source-eligibility decisions."
  ),
  file.path(
    output_root,
    "RUN_COMPLETE.txt"
  )
)


# ------------------------------------------------------------
# 32. FINAL CONSOLE SUMMARY
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "WOS UNIQUE-REFERENCE ELIGIBILITY TRIAGE\n"
)

cat(
  "============================================\n"
)


cat(
  "\nProvisional unique WoS References:",
  nrow(review_export),
  "\n"
)


cat(
  "High-priority full-text candidates:",
  nrow(high_priority_fulltext),
  "\n"
)


cat(
  "Manual eligibility reviews:",
  nrow(manual_eligibility),
  "\n"
)


cat(
  "Modifier / complex-design candidates:",
  nrow(modifier_complex_design),
  "\n"
)


cat(
  "Clear outside-scope References from metadata:",
  nrow(clear_outside_scope),
  "\n"
)


cat(
  "Active full-text review queue:",
  nrow(full_text_queue),
  "\n"
)


cat(
  "\nTriage classification check:",
  triage_count_check,
  "records classified out of",
  nrow(review_export),
  "\n"
)


cat(
  "\nRecovered DOI values:",
  sum(
    !is.na(review_export$best_doi) &
      review_export$best_doi != ""
  ),
  "out of",
  nrow(review_export),
  "\n"
)


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
    "This workbook is a reference-level full-text triage only. ",
    "No reference has been automatically added to or excluded from ",
    "the thesis dataset. Final source decisions must be recorded ",
    "manually, and eligible LC50 Results must later be extracted at ",
    "the Result level before any eight-metal support/comparability ",
    "assessment is repeated.\n"
  )
)