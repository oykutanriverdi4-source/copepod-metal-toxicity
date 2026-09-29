# ============================================================
# 01_screen_envirotox_records.R
# ============================================================
# Purpose
# Reproduce the preliminary screening of the EnviroTox export for
# acute LC50 records on the eight prespecified metals in
# marine/estuarine copepods.
#
# Database search represented by the local export
# EnviroTox version 2.0.0
# Search date: 5 September 2026
# Advanced Search conditions (AND):
# - Heavy Metals contains "1"
# - Trophic Level contains "INVERT"
# - Test type = "A"
# - Test statistic = "LC50"
#
# Final thesis screening benchmarks
# - 3,375 raw test records
# - 106 copepod-order records
# - 76 records for the eight prespecified metals
# - 22 Freshwater exclusions
# - 54 Saltwater candidate records
# - 17 species in the Saltwater candidate set
# - 23 source citations in the Saltwater candidate set
#
# Important scope notes
# - Copepods are identified by Taxonomic order because EnviroTox v2
#   stores these taxa under the older class label "Maxillopoda".
# - Saltwater is used as the database-level proxy for the thesis
#   marine/estuarine scope; original test conditions are reviewed
#   subsequently.
# - Candidate records are not treated as final eligible records.
# - Possible within-EnviroTox duplicates are flagged, not deleted.
#
# Local-only input
# data/external/envirotox/envirotox_raw_export.xlsx
#
# Local detailed output
# outputs/database_audits/envirotox/01_screen_envirotox_records/
#   envirotox_preliminary_screening.xlsx
#
# The raw EnviroTox export and detailed audit workbook are third-party
# database exports / derived audit materials and are kept local rather
# than redistributed in the public repository. Compact screening
# summaries and search metadata are written separately.
#
# Scientific rule
# Taxonomic filtering, eight-metal filtering, Medium-based screening,
# technical eligibility checks, duplicate flags and frozen numerical
# benchmarks are preserved from the final thesis workflow.
# ============================================================


# ============================================================
# 0. PACKAGES, FILES, AND OUTPUT DIRECTORY
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tibble",
  "writexl",
  "readr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  if (identical(getOption("repos")[["CRAN"]], "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  install.packages(missing_packages)
}

library(readxl)
library(dplyr)
library(tibble)
library(writexl)
library(readr)

options(width = 200)

input_file <- file.path(
  "data",
  "external",
  "envirotox",
  "envirotox_raw_export.xlsx"
)

output_dir <- file.path(
  "outputs",
  "database_audits",
  "envirotox",
  "01_screen_envirotox_records"
)

output_file <- file.path(
  output_dir,
  "envirotox_preliminary_screening.xlsx"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Local EnviroTox export not found: ",
      input_file,
      ". Place the EnviroTox workbook exported from the documented ",
      "Advanced Search at this path."
    )
  )
}

cat("\n============================================\n")
cat("ENVIROTOX PRELIMINARY SCREENING\n")
cat("============================================\n")
cat("Input:", normalizePath(input_file), "\n")


# ============================================================
# 1. READ THE THREE ENVIROTOX SHEETS AS TEXT
# ============================================================

available_sheets <- excel_sheets(input_file)
required_sheets <- c("test", "substance", "taxonomy")

if (!all(required_sheets %in% available_sheets)) {
  stop(
    paste0(
      "The workbook must contain these sheets: ",
      paste(required_sheets, collapse = ", "),
      ".\nAvailable sheets: ",
      paste(available_sheets, collapse = ", ")
    )
  )
}

read_envirotox_sheet <- function(sheet_name) {
  suppressWarnings(
    read_excel(
      input_file,
      sheet = sheet_name,
      col_types = "text",
      .name_repair = "minimal"
    )
  )
}

test_raw <- read_envirotox_sheet("test") %>%
  mutate(
    EnviroTox_Row = row_number() + 1L,
    .before = 1
  )

substance_raw <- read_envirotox_sheet("substance") %>%
  mutate(
    Substance_Row = row_number() + 1L,
    .before = 1
  )

taxonomy_raw <- read_envirotox_sheet("taxonomy") %>%
  mutate(
    Taxonomy_Row = row_number() + 1L,
    .before = 1
  )

required_test_columns <- c(
  "CAS",
  "Chemical name",
  "Latin name",
  "Trophic Level",
  "Effect",
  "Effect value",
  "Unit",
  "Test type",
  "Test statistic",
  "Duration",
  "Duration (days)",
  "Duration (hours)",
  "Effect is 5X above water solubility",
  "Source",
  "version",
  "Reported chemical name",
  "original CAS"
)

required_taxonomy_columns <- c(
  "Latin name",
  "Trophic Level",
  "Medium",
  "Taxonomic kingdom",
  "Taxonomic phylum or division",
  "Taxonomic subphylum",
  "Taxonomic superclass",
  "Taxonomic class",
  "Taxonomic order",
  "Taxonomic family"
)

if (!all(required_test_columns %in% names(test_raw))) {
  stop(
    "The test sheet does not contain all expected EnviroTox columns."
  )
}

if (!all(required_taxonomy_columns %in% names(taxonomy_raw))) {
  stop(
    "The taxonomy sheet does not contain all expected EnviroTox columns."
  )
}

if (anyDuplicated(taxonomy_raw$`Latin name`) > 0) {
  stop(
    "The taxonomy sheet contains repeated Latin names. Review the taxonomy join before continuing."
  )
}


# ============================================================
# 2. DEFINE THE PRESPECIFIED SCOPE
# ============================================================

target_orders <- c(
  "Calanoida",
  "Harpacticoida",
  "Cyclopoida"
)

metal_lookup <- tribble(
  ~CAS,           ~Metal,
  "Metalgrp.Ag", "Ag",
  "Metalgrp.Cd", "Cd",
  "Metalgrp.Cr", "Cr",
  "Metalgrp.Cu", "Cu",
  "Metalgrp.Hg", "Hg",
  "Metalgrp.Ni", "Ni",
  "Metalgrp.Pb", "Pb",
  "Metalgrp.Zn", "Zn"
)

target_metal_groups <- metal_lookup$CAS

normalize_text <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

parse_number_safely <- function(x) {
  x <- normalize_text(x)
  suppressWarnings(as.numeric(x))
}

taxonomy_for_join <- taxonomy_raw %>%
  transmute(
    `Latin name` = normalize_text(`Latin name`),
    Taxonomy_Row,
    Trophic_Level_Taxonomy = normalize_text(`Trophic Level`),
    Medium = normalize_text(Medium),
    Taxonomic_Kingdom = normalize_text(`Taxonomic kingdom`),
    Taxonomic_Phylum = normalize_text(`Taxonomic phylum or division`),
    Taxonomic_Subphylum = normalize_text(`Taxonomic subphylum`),
    Taxonomic_Superclass = normalize_text(`Taxonomic superclass`),
    Taxonomic_Class = normalize_text(`Taxonomic class`),
    Taxonomic_Order = normalize_text(`Taxonomic order`),
    Taxonomic_Family = normalize_text(`Taxonomic family`)
  )


# ============================================================
# 3. JOIN TEST RECORDS TO TAXONOMY AND CREATE QA FLAGS
# ============================================================

screened_all <- test_raw %>%
  mutate(
    across(all_of(required_test_columns), normalize_text),
    Effect_Value_mg_L = parse_number_safely(`Effect value`),
    Duration_Days_Numeric = parse_number_safely(`Duration (days)`),
    Duration_Hours_Numeric = parse_number_safely(`Duration (hours)`),
    Solubility_5X_Flag = parse_number_safely(
      `Effect is 5X above water solubility`
    )
  ) %>%
  left_join(
    taxonomy_for_join,
    by = "Latin name"
  ) %>%
  left_join(
    metal_lookup,
    by = "CAS"
  ) %>%
  mutate(
    Taxonomy_Match_Flag = !is.na(Taxonomy_Row),
    Copepod_Order_Flag = Taxonomic_Order %in% target_orders,
    Target_Metal_Flag = CAS %in% target_metal_groups,
    Saltwater_Flag = Medium == "Saltwater",
    Acute_Test_Flag = toupper(trimws(`Test type`)) == "A",
    LC50_Flag = toupper(trimws(`Test statistic`)) == "LC50",
    Mortality_Flag = grepl(
      "mortality",
      tolower(`Effect`),
      fixed = TRUE
    ),
    Positive_Effect_Value_Flag =
      !is.na(Effect_Value_mg_L) & Effect_Value_mg_L > 0,
    Unit_mg_L_Flag = tolower(gsub("[[:space:]]", "", Unit)) == "mg/l",
    Positive_Duration_Flag =
      !is.na(Duration_Hours_Numeric) & Duration_Hours_Numeric > 0
  )

if (any(!screened_all$Taxonomy_Match_Flag)) {
  warning(
    sum(!screened_all$Taxonomy_Match_Flag),
    " test rows did not match the taxonomy sheet by Latin name."
  )
}


# ============================================================
# 4. APPLY THE DATABASE-LEVEL SCREENING FLOW
# ============================================================

# Step A: all copepod-order records returned by the EnviroTox query.
all_copepod_heavy_metals <- screened_all %>%
  filter(Copepod_Order_Flag)

# Step B: retain the eight metals prespecified before analysis.
target8_all_media <- all_copepod_heavy_metals %>%
  filter(Target_Metal_Flag)

# Step C: document records excluded only because they are freshwater.
freshwater_excluded <- target8_all_media %>%
  filter(Medium == "Freshwater") %>%
  mutate(
    Preliminary_Exclusion_Reason =
      "Freshwater record outside the marine/estuarine thesis scope"
  )

# Step D: preliminary marine/estuarine candidate set.
# These conditions reproduce the intended database search and retain
# only rows that can potentially be harmonized with the thesis data.
saltwater_candidates <- target8_all_media %>%
  filter(
    Saltwater_Flag,
    Acute_Test_Flag,
    LC50_Flag,
    Mortality_Flag,
    Positive_Effect_Value_Flag,
    Unit_mg_L_Flag,
    Positive_Duration_Flag
  )

# Any target-metal copepod rows that fail a technical eligibility check
# are preserved here rather than silently discarded.
technical_exclusions <- target8_all_media %>%
  filter(
    !Saltwater_Flag |
      !Acute_Test_Flag |
      !LC50_Flag |
      !Mortality_Flag |
      !Positive_Effect_Value_Flag |
      !Unit_mg_L_Flag |
      !Positive_Duration_Flag
  ) %>%
  mutate(
    Preliminary_Exclusion_Reason = "",
    Preliminary_Exclusion_Reason = ifelse(
      !Saltwater_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Outside Saltwater medium; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !Acute_Test_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Test type is not A; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !LC50_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Test statistic is not LC50; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !Mortality_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Effect is not mortality; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !Positive_Effect_Value_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Effect value is missing, nonnumeric, or nonpositive; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !Unit_mg_L_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Unit is not mg/L; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = ifelse(
      !Positive_Duration_Flag,
      paste0(
        Preliminary_Exclusion_Reason,
        "Duration is missing, nonnumeric, or nonpositive; "
      ),
      Preliminary_Exclusion_Reason
    ),
    Preliminary_Exclusion_Reason = sub(
      "; $",
      "",
      Preliminary_Exclusion_Reason
    )
  )

# Copepod records for heavy metals outside the prespecified eight-metal
# universe are retained for audit but are not eligible for this thesis.
non_target_metals <- all_copepod_heavy_metals %>%
  filter(!Target_Metal_Flag) %>%
  mutate(
    Preliminary_Exclusion_Reason =
      "Metal outside the eight prespecified thesis metals"
  )


# ============================================================
# 5. FLAG POSSIBLE WITHIN-ENVIROTOX DUPLICATES
# ============================================================

# This is only a flag. No record is deleted at this stage.
saltwater_candidates <- saltwater_candidates %>%
  mutate(
    Duplicate_Key = paste(
      tolower(trimws(`Latin name`)),
      Metal,
      format(Effect_Value_mg_L, scientific = FALSE, trim = TRUE),
      tolower(trimws(Unit)),
      format(Duration_Hours_Numeric, scientific = FALSE, trim = TRUE),
      tolower(trimws(`original CAS`)),
      tolower(trimws(Source)),
      sep = " | "
    )
  ) %>%
  group_by(Duplicate_Key) %>%
  mutate(
    Within_EnviroTox_Group_Size = n(),
    Within_EnviroTox_Duplicate_Flag = n() > 1
  ) %>%
  ungroup() %>%
  mutate(
    Source_Verification_Status = "Not checked",
    ECOTOX_Overlap_Status = "Not checked",
    Final_Eligibility_Decision = "Pending",
    Final_Decision_Reason = "",
    Reviewer_Notes = ""
  ) %>%
  arrange(Metal, `Latin name`, Duration_Hours_Numeric, Effect_Value_mg_L)

within_envirotox_duplicates <- saltwater_candidates %>%
  filter(Within_EnviroTox_Duplicate_Flag) %>%
  arrange(Duplicate_Key, EnviroTox_Row)


# ============================================================
# 6. CREATE AUDIT SUMMARIES
# ============================================================

screening_flow <- tribble(
  ~Step, ~Dataset, ~Records, ~Interpretation,
  1L,
  "Raw EnviroTox test export",
  nrow(test_raw),
  "All records returned by the Advanced Search query",
  2L,
  "Copepod orders",
  nrow(all_copepod_heavy_metals),
  "Calanoida, Harpacticoida, or Cyclopoida",
  3L,
  "Eight target metals",
  nrow(target8_all_media),
  "Ag, Cd, Cr, Cu, Hg, Ni, Pb, or Zn",
  4L,
  "Freshwater exclusions",
  nrow(freshwater_excluded),
  "Outside the marine/estuarine thesis scope",
  5L,
  "Saltwater preliminary candidates",
  nrow(saltwater_candidates),
  "Requires source verification and ECOTOX overlap assessment"
)

metal_summary <- tibble(
  Metal = c("Ag", "Cd", "Cr", "Cu", "Hg", "Ni", "Pb", "Zn")
) %>%
  left_join(
    target8_all_media %>%
      count(Metal, name = "All_Media_Records"),
    by = "Metal"
  ) %>%
  left_join(
    freshwater_excluded %>%
      count(Metal, name = "Freshwater_Excluded"),
    by = "Metal"
  ) %>%
  left_join(
    saltwater_candidates %>%
      count(Metal, name = "Saltwater_Candidates"),
    by = "Metal"
  ) %>%
  mutate(
    across(
      c(
        All_Media_Records,
        Freshwater_Excluded,
        Saltwater_Candidates
      ),
      ~ ifelse(is.na(.x), 0L, as.integer(.x))
    )
  )

species_summary <- saltwater_candidates %>%
  count(
    Taxonomic_Order,
    Taxonomic_Family,
    `Latin name`,
    Metal,
    name = "Candidate_Records"
  ) %>%
  arrange(Taxonomic_Order, `Latin name`, Metal)

source_summary <- saltwater_candidates %>%
  count(Source, name = "Candidate_Records", sort = TRUE) %>%
  mutate(
    Source_Verification_Status = "Not checked",
    ECOTOX_Overlap_Status = "Not checked",
    Reviewer_Notes = ""
  )

copepod_taxonomy <- taxonomy_for_join %>%
  filter(Taxonomic_Order %in% target_orders) %>%
  arrange(Medium, Taxonomic_Order, `Latin name`)

expected_counts <- tribble(
  ~Metric, ~Expected, ~Observed,
  "Raw test rows", 3375L, nrow(test_raw),
  "Copepod-order rows", 106L, nrow(all_copepod_heavy_metals),
  "Eight-metal rows, all media", 76L, nrow(target8_all_media),
  "Freshwater exclusions", 22L, nrow(freshwater_excluded),
  "Saltwater preliminary candidates", 54L, nrow(saltwater_candidates),
  "Saltwater unique species", 17L,
  n_distinct(saltwater_candidates$`Latin name`),
  "Saltwater unique sources", 23L,
  n_distinct(saltwater_candidates$Source)
) %>%
  mutate(
    Check = ifelse(Expected == Observed, "PASS", "REVIEW")
  )

stopifnot(
  all(
    expected_counts$Check == "PASS"
  )
)


# ============================================================
# 7. README / METHODS NOTES FOR THE OUTPUT WORKBOOK
# ============================================================

readme <- tribble(
  ~Item, ~Value,
  "Purpose",
  paste(
    "Preliminary audit of EnviroTox records for acute metal",
    "toxicity in marine/estuarine copepods."
  ),
  "Input file",
  "Local EnviroTox Advanced Search export (not redistributed)",
  "Input sheets",
  "test; substance; taxonomy",
  "Database query",
  paste(
    "Heavy Metals contains 1 AND Trophic Level contains INVERT",
    "AND Test type = A AND Test statistic = LC50"
  ),
  "Copepod definition",
  "Taxonomic order is Calanoida, Harpacticoida, or Cyclopoida",
  "Why order was used",
  paste(
    "EnviroTox v2 uses the older class label Maxillopoda for these",
    "copepods; order-level filtering is more reliable."
  ),
  "Target metals",
  "Ag, Cd, Cr, Cu, Hg, Ni, Pb, and Zn",
  "Environment rule",
  paste(
    "Medium = Saltwater is retained as the database-level proxy",
    "for the marine/estuarine thesis scope."
  ),
  "Freshwater rule",
  paste(
    "Freshwater records are documented but excluded because they",
    "fall outside the predefined environmental scope."
  ),
  "Outcome rule",
  "Acute test type A; LC50; mortality; positive numeric value in mg/L",
  "Duration rule",
  "A positive numeric duration in hours is required for harmonization",
  "Solubility flag",
  paste(
    "Effect is 5X above water solubility is retained for QA;",
    "flagged rows must be reviewed and are not silently removed."
  ),
  "Duplicate rule",
  paste(
    "Possible within-EnviroTox duplicates are flagged, not deleted,",
    "before source-level adjudication."
  ),
  "Interpretation",
  paste(
    "Saltwater_Candidates is a source-verification set, not a final",
    "addition to the ECOTOX analysis dataset."
  ),
  "Interpretive boundary",
  paste(
    "This workbook documents database-level screening only;",
    "source verification and ECOTOX overlap are separate audit steps."
  )
)


# ============================================================
# 8. SELECT AND ORDER USER-FACING COLUMNS
# ============================================================

candidate_columns <- c(
  "EnviroTox_Row",
  "Metal",
  "CAS",
  "original CAS",
  "Chemical name",
  "Reported chemical name",
  "Latin name",
  "Medium",
  "Taxonomic_Order",
  "Taxonomic_Family",
  "Effect",
  "Effect_Value_mg_L",
  "Unit",
  "Test type",
  "Test statistic",
  "Duration",
  "Duration_Days_Numeric",
  "Duration_Hours_Numeric",
  "Solubility_5X_Flag",
  "Source",
  "version",
  "Within_EnviroTox_Duplicate_Flag",
  "Within_EnviroTox_Group_Size",
  "Duplicate_Key",
  "Source_Verification_Status",
  "ECOTOX_Overlap_Status",
  "Final_Eligibility_Decision",
  "Final_Decision_Reason",
  "Reviewer_Notes"
)

saltwater_candidates_export <- saltwater_candidates %>%
  select(all_of(candidate_columns))

within_envirotox_duplicates_export <- within_envirotox_duplicates %>%
  select(all_of(candidate_columns))

freshwater_excluded_export <- freshwater_excluded %>%
  select(
    EnviroTox_Row,
    Metal,
    CAS,
    `original CAS`,
    `Chemical name`,
    `Reported chemical name`,
    `Latin name`,
    Medium,
    Taxonomic_Order,
    Taxonomic_Family,
    Effect,
    Effect_Value_mg_L,
    Unit,
    `Test type`,
    `Test statistic`,
    Duration,
    Duration_Hours_Numeric,
    Source,
    Preliminary_Exclusion_Reason
  ) %>%
  arrange(Metal, `Latin name`, Duration_Hours_Numeric)

technical_exclusions_export <- technical_exclusions %>%
  select(
    EnviroTox_Row,
    Metal,
    CAS,
    `original CAS`,
    `Latin name`,
    Medium,
    Taxonomic_Order,
    Effect,
    `Effect value`,
    Unit,
    `Test type`,
    `Test statistic`,
    Duration,
    `Duration (hours)`,
    Source,
    Preliminary_Exclusion_Reason
  ) %>%
  arrange(Metal, `Latin name`, EnviroTox_Row)

non_target_metals_export <- non_target_metals %>%
  select(
    EnviroTox_Row,
    CAS,
    `original CAS`,
    `Chemical name`,
    `Reported chemical name`,
    `Latin name`,
    Medium,
    Taxonomic_Order,
    Effect,
    `Effect value`,
    Unit,
    `Test type`,
    `Test statistic`,
    Duration,
    `Duration (hours)`,
    Source,
    Preliminary_Exclusion_Reason
  ) %>%
  arrange(CAS, `Latin name`, EnviroTox_Row)

all_copepod_heavy_metals_export <- all_copepod_heavy_metals %>%
  select(
    EnviroTox_Row,
    Metal,
    CAS,
    `original CAS`,
    `Chemical name`,
    `Reported chemical name`,
    `Latin name`,
    Medium,
    Taxonomic_Order,
    Taxonomic_Family,
    Effect,
    `Effect value`,
    Unit,
    `Test type`,
    `Test statistic`,
    Duration,
    `Duration (hours)`,
    `Effect is 5X above water solubility`,
    Source,
    version
  ) %>%
  arrange(CAS, `Latin name`, EnviroTox_Row)


# ============================================================
# 9. WRITE THE AUDIT WORKBOOK
# ============================================================

output_sheets <- list(
  README = readme,
  Screening_Flow = screening_flow,
  Expected_Counts = expected_counts,
  Metal_Summary = metal_summary,
  Species_Summary = species_summary,
  Source_Summary = source_summary,
  Saltwater_Candidates = saltwater_candidates_export,
  Freshwater_Excluded = freshwater_excluded_export,
  Technical_Exclusions = technical_exclusions_export,
  Within_DB_Duplicates = within_envirotox_duplicates_export,
  NonTarget_Metals = non_target_metals_export,
  All_Copepod_HM = all_copepod_heavy_metals_export,
  Copepod_Taxonomy = copepod_taxonomy
)

write_xlsx(
  x = output_sheets,
  path = output_file
)


screening_count_summary <- tibble::tibble(
  Metric = expected_counts$Metric,
  N = expected_counts$Observed
)

readr::write_csv(
  screening_count_summary,
  file.path(
    output_dir,
    "screening_count_summary.csv"
  )
)

readr::write_csv(
  metal_summary,
  file.path(
    output_dir,
    "metal_summary.csv"
  )
)

search_metadata <- tibble::tibble(
  Field = c(
    "Database",
    "Version",
    "Search date",
    "Heavy Metals",
    "Trophic Level",
    "Test type",
    "Test statistic",
    "Copepod orders",
    "Marine/estuarine screening proxy"
  ),
  Value = c(
    "EnviroTox",
    "2.0.0",
    "2026-09-05",
    'contains "1"',
    'contains "INVERT"',
    "A",
    "LC50",
    "Calanoida; Harpacticoida; Cyclopoida",
    "Medium = Saltwater"
  )
)

readr::write_csv(
  search_metadata,
  file.path(
    output_dir,
    "search_metadata.csv"
  )
)

input_output_manifest <- tibble::tibble(
  Role = c(
    "Local input: EnviroTox raw export",
    "Local output: detailed preliminary-screening workbook",
    "Output: screening count summary",
    "Output: metal summary",
    "Output: search metadata"
  ),
  Path = c(
    input_file,
    output_file,
    file.path(
      output_dir,
      "screening_count_summary.csv"
    ),
    file.path(
      output_dir,
      "metal_summary.csv"
    ),
    file.path(
      output_dir,
      "search_metadata.csv"
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
    "STATUS: ENVIROTOX PRELIMINARY SCREENING = PASS",
    paste0(
      "Raw test rows: ",
      nrow(test_raw)
    ),
    paste0(
      "Copepod-order rows: ",
      nrow(all_copepod_heavy_metals)
    ),
    paste0(
      "Eight-metal rows, all media: ",
      nrow(target8_all_media)
    ),
    paste0(
      "Freshwater exclusions: ",
      nrow(freshwater_excluded)
    ),
    paste0(
      "Saltwater preliminary candidates: ",
      nrow(saltwater_candidates)
    )
  ),
  file.path(
    output_dir,
    "RUN_COMPLETE.txt"
  )
)


# ============================================================
# 10. FINAL CONSOLE REPORT
# ============================================================

cat("\nScreening completed.\n")
cat("Raw test rows:", nrow(test_raw), "\n")
cat("Copepod-order rows:", nrow(all_copepod_heavy_metals), "\n")
cat("Eight-metal rows, all media:", nrow(target8_all_media), "\n")
cat("Freshwater exclusions:", nrow(freshwater_excluded), "\n")
cat("Saltwater preliminary candidates:", nrow(saltwater_candidates), "\n")
cat(
  "Saltwater unique species:",
  n_distinct(saltwater_candidates$`Latin name`),
  "\n"
)
cat(
  "Saltwater unique sources:",
  n_distinct(saltwater_candidates$Source),
  "\n"
)
cat(
  "Possible within-EnviroTox duplicate rows:",
  nrow(within_envirotox_duplicates),
  "\n"
)
cat("\nDetailed local audit workbook written to:\n", normalizePath(output_file), "\n")
cat("\nSTATUS: ENVIROTOX PRELIMINARY SCREENING = PASS\n")

