# ============================================================
# 03_prepare_envirotox_manual_review.R
# ============================================================
# Purpose
# Create an editable, source-first manual-review workbook for the
# 16 EnviroTox records that were not resolved by the automatic
# EnviroTox-ECOTOX overlap audit.
#
# This script PREPARES the review workbook; it does not make or
# reproduce the final source-verification decisions by itself.
#
# Final thesis workflow context
# - 54 Saltwater EnviroTox candidates entered the overlap audit.
# - 38 were exact Result-level ECOTOX overlaps.
# - 16 required original-source review across 8 source citations.
# - Subsequent source review confirmed no separate EnviroTox-derived
#   quantitative records for the combined dataset.
# - The review did, however, contribute to checking overlapping and
#   unresolved records, including reinstatement of one ECOTOX record
#   after original-source review.
#
# Input
# outputs/database_audits/envirotox/02_check_envirotox_ecotox_overlap/
#   envirotox_ecotox_overlap.xlsx
#
# Local output
# outputs/database_audits/envirotox/03_prepare_envirotox_manual_review/
#   envirotox_manual_review.xlsx
#
# The workbook is intentionally editable and local-only. It contains
# database-derived row-level audit material and manual-review fields,
# so it should not be committed to the public repository.
#
# Safety
# By default, an existing manual-review workbook is NOT overwritten.
# Set OVERWRITE_EXISTING <- TRUE only when replacement is intentional.
#
# Scientific rule
# The source-first review structure, controlled vocabularies, evidence
# fields and decision categories are preserved from the final thesis
# workflow. Repository cleanup changes paths, naming, documentation
# and frozen QA only.
# ============================================================


# ---- 1. User setting ---------------------------------------------------------
OVERWRITE_EXISTING <- FALSE

# ---- 2. Packages -------------------------------------------------------------
required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "tibble",
  "openxlsx",
  "readr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Install required package(s) before running this script: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(openxlsx)
  library(readr)
})

# ---- 3. Repository file paths -----------------------------------------------

input_file <- file.path(
  "outputs",
  "database_audits",
  "envirotox",
  "02_check_envirotox_ecotox_overlap",
  "envirotox_ecotox_overlap.xlsx"
)

output_dir <- file.path(
  "outputs",
  "database_audits",
  "envirotox",
  "03_prepare_envirotox_manual_review"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

output_file <- file.path(
  output_dir,
  "envirotox_manual_review.xlsx"
)

if (!file.exists(input_file)) {
  stop(
    "EnviroTox overlap workbook not found: ",
    input_file
  )
}

if (
  file.exists(output_file) &&
  !isTRUE(OVERWRITE_EXISTING)
) {
  stop(
    "A manual-review workbook already exists and was not overwritten: ",
    output_file,
    ". This protects any decisions entered manually. ",
    "Set OVERWRITE_EXISTING <- TRUE only if replacement is intended."
  )
}

cat(
  "\nOverlap input: ",
  input_file,
  "\n",
  "Manual-review output: ",
  output_file,
  "\n\n",
  sep = ""
)


# ---- 4. Read and validate the automatic audit -------------------------------
available_sheets <- excel_sheets(input_file)
if (!"Manual_Review" %in% available_sheets) {
  stop("The input workbook does not contain a Manual_Review sheet.", call. = FALSE)
}

manual_input <- read_excel(input_file, sheet = "Manual_Review")

required_columns <- c(
  "EnviroTox_Row", "Metal", "CAS", "Chemical name", "Reported chemical name",
  "Latin name", "Medium", "Taxonomic_Order", "Taxonomic_Family", "Effect",
  "Effect_Value_mg_L", "Unit", "Test type", "Test statistic", "Duration",
  "Duration_Days_Numeric", "Duration_Hours_Numeric", "Source",
  "Source_Overlap_Status", "Matched_ECOTOX_Reference", "Best_ECOTOX_Title",
  "ECOTOX_Record_Overlap_Status", "ECOTOX_Reference_Number",
  "ECOTOX_Test_Number", "ECOTOX_Result_Number", "ECOTOX_Species",
  "ECOTOX_Duration_Days", "ECOTOX_LC50_mg_L", "Relative_Value_Difference",
  "ECOTOX_Eligibility_Status", "ECOTOX_Eligibility_Reason",
  "In_Harmonized_295", "Harmonization_Status", "Record_Interpretation",
  "Thesis_Dataset_Status", "Manual_Review_Priority"
)

missing_columns <- setdiff(required_columns, names(manual_input))
if (length(missing_columns) > 0) {
  stop(
    "Manual_Review is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(manual_input) == 0) {
  stop("Manual_Review contains no records.", call. = FALSE)
}

# Frozen final thesis manual-review queue.
stopifnot(
  nrow(manual_input) == 16L,
  dplyr::n_distinct(manual_input$Source) == 8L
)

cat("Manual records found:", nrow(manual_input), "\n")
cat("Unique sources found:", dplyr::n_distinct(manual_input$Source), "\n\n")

# ---- 5. Build the source-level review table ---------------------------------
source_review <- manual_input %>%
  group_by(Source) %>%
  summarise(
    Review_Order = min(Manual_Review_Priority, na.rm = TRUE),
    Candidate_Record_Count = n(),
    Automated_Source_Status = first(Source_Overlap_Status),
    Matched_ECOTOX_Reference = first(Matched_ECOTOX_Reference),
    Best_ECOTOX_Title = first(Best_ECOTOX_Title),
    Record_Statuses = paste(sort(unique(ECOTOX_Record_Overlap_Status)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Review_Order, Source) %>%
  mutate(
    Source_ID = sprintf("SRC%02d", row_number()),
    Verification_Task = case_when(
      Automated_Source_Status == "NOT_FOUND_IN_ECOTOX_FRAME" ~
        "Verify bibliography, obtain the source, and determine whether the study and LC50 record are genuinely absent from ECOTOX.",
      str_detect(Record_Statuses, fixed("SOURCE_ONLY_REVIEW")) ~
        "The study is known to ECOTOX; locate the specific species-metal endpoint and determine why the result is absent.",
      str_detect(Record_Statuses, fixed("SAME_CONTEXT_VALUE_DIFF_REVIEW")) ~
        "Check the original table for the exact LC50 value, unit, concentration basis, and whether multiple tests are reported.",
      TRUE ~
        "Check the original table for exposure duration and determine whether EnviroTox and ECOTOX describe the same experiment."
    ),
    Full_Text_Status = "NOT_SEARCHED",
    Full_Text_or_URL = "",
    Bibliographic_Identity = "UNCERTAIN",
    Source_Evidence_Location = "",
    Source_Level_Notes = "",
    Source_Review_Complete = "NO"
  ) %>%
  select(
    Source_ID, Review_Order, Source, Candidate_Record_Count,
    Automated_Source_Status, Matched_ECOTOX_Reference, Best_ECOTOX_Title,
    Record_Statuses, Verification_Task, Full_Text_Status, Full_Text_or_URL,
    Bibliographic_Identity, Source_Evidence_Location, Source_Level_Notes,
    Source_Review_Complete
  )

source_lookup <- source_review %>% select(Source_ID, Source)

# ---- 6. Build the record-level review table ---------------------------------
record_review <- manual_input %>%
  left_join(source_lookup, by = "Source") %>%
  arrange(Manual_Review_Priority, Source_ID, Metal, `Latin name`, EnviroTox_Row) %>%
  mutate(
    Review_ID = sprintf("REC%02d", row_number()),
    Review_Status = "NOT_STARTED",
    Verified_Species = "UNCERTAIN",
    Verified_Medium = "UNCLEAR",
    Verified_Endpoint = "UNCLEAR",
    Verified_Exposure = "UNCLEAR",
    Verified_Duration_Hours = NA_real_,
    Verified_LC50_Value = NA_real_,
    Verified_Original_Unit = "",
    Concentration_Basis = "UNCLEAR",
    Same_Experiment_as_ECOTOX = "UNCERTAIN",
    Original_Source_Page_Table = "",
    Evidence_Note = "",
    Final_Decision = "",
    Final_Decision_Reason = "",
    Reviewer_Notes = "",
    Reviewed_By = "",
    Review_Date = as.Date(NA)
  ) %>%
  transmute(
    Review_ID,
    Source_ID,
    Priority = Manual_Review_Priority,
    EnviroTox_Row,
    Metal,
    Latin_Name = `Latin name`,
    Taxonomic_Order,
    Taxonomic_Family,
    EnviroTox_Medium = Medium,
    Chemical_Name = `Chemical name`,
    Reported_Chemical_Name = `Reported chemical name`,
    EnviroTox_Effect = Effect,
    EnviroTox_Test_Type = `Test type`,
    EnviroTox_Test_Statistic = `Test statistic`,
    EnviroTox_Duration_Text = Duration,
    EnviroTox_Duration_Hours = Duration_Hours_Numeric,
    EnviroTox_LC50_mg_L = Effect_Value_mg_L,
    EnviroTox_Unit = Unit,
    Source,
    Automatic_Overlap_Status = ECOTOX_Record_Overlap_Status,
    Matched_ECOTOX_Reference,
    ECOTOX_Test_Number,
    ECOTOX_Result_Number,
    ECOTOX_Species,
    ECOTOX_Duration_Hours = ECOTOX_Duration_Days * 24,
    ECOTOX_LC50_mg_L,
    Relative_Value_Difference,
    ECOTOX_Eligibility_Status,
    ECOTOX_Eligibility_Reason,
    In_Harmonized_295,
    Harmonization_Status,
    Thesis_Dataset_Status,
    Automatic_Interpretation = Record_Interpretation,
    Review_Status,
    Verified_Species,
    Verified_Medium,
    Verified_Endpoint,
    Verified_Exposure,
    Verified_Duration_Hours,
    Verified_LC50_Value,
    Verified_Original_Unit,
    Concentration_Basis,
    Same_Experiment_as_ECOTOX,
    Original_Source_Page_Table,
    Evidence_Note,
    Final_Decision,
    Final_Decision_Reason,
    Reviewer_Notes,
    Reviewed_By,
    Review_Date
  )

# ---- 7. Dropdown lists -------------------------------------------------------
decision_lists <- list(
  Full_Text_Status = c(
    "NOT_SEARCHED", "FULL_TEXT_FOUND", "ABSTRACT_ONLY",
    "METADATA_ONLY", "NOT_FOUND"
  ),
  Bibliographic_Identity = c(
    "SAME_SOURCE_AS_ECOTOX", "SOURCE_ABSENT_FROM_ECOTOX", "UNCERTAIN"
  ),
  Yes_No = c("YES", "NO"),
  Yes_No_Uncertain = c("YES", "NO", "UNCERTAIN"),
  Verified_Medium = c(
    "MARINE", "ESTUARINE", "BRACKISH", "FRESHWATER", "UNCLEAR"
  ),
  Verified_Endpoint = c(
    "LC50_MORTALITY", "OTHER_LETHAL_ENDPOINT", "SUBLETHAL_ENDPOINT", "UNCLEAR"
  ),
  Verified_Exposure = c(
    "SINGLE_TARGET_METAL", "METAL_MIXTURE", "OTHER_CHEMICAL", "UNCLEAR"
  ),
  Concentration_Basis = c(
    "TARGET_METAL", "TEST_COMPOUND", "DISSOLVED_METAL",
    "TOTAL_METAL", "UNCLEAR"
  ),
  Review_Status = c("NOT_STARTED", "IN_PROGRESS", "COMPLETE", "BLOCKED"),
  Final_Decision = c(
    "DUPLICATE_ECOTOX",
    "NEW_ELIGIBLE_RECORD",
    "PRESENT_IN_ECOTOX_BUT_EXCLUDED",
    "EXCLUDE_OUT_OF_SCOPE",
    "EXCLUDE_METHOD",
    "UNRESOLVED"
  )
)

max_list_length <- max(vapply(decision_lists, length, integer(1)))
lists_sheet <- tibble(Row = seq_len(max_list_length))
for (list_name in names(decision_lists)) {
  values <- decision_lists[[list_name]]
  lists_sheet[[list_name]] <- c(values, rep(NA_character_, max_list_length - length(values)))
}
lists_sheet <- lists_sheet %>% select(-Row)

# ---- 8. README and live progress sheet --------------------------------------
readme <- tribble(
  ~Item, ~Instruction,
  "Purpose", "Manually verify the 16 records not resolved by the automatic EnviroTox-ECOTOX overlap audit.",
  "Automatic audit", "Keep the automatic overlap workbook unchanged; enter manual decisions only in this review workbook.",
  "Start here", "Open Source_Review and work from the smallest Review_Order upward.",
  "Source first", "Review each source once; multiple linked records share the same Source_ID.",
  "Evidence", "Record a page, table, figure, appendix, report section, or stable URL for every completed decision.",
  "Record review", "Enter source-verified species, medium, endpoint, exposure, duration, value, unit, and concentration basis.",
  "Duplicate", "Use DUPLICATE_ECOTOX only when the original source shows that both databases represent the same experimental result.",
  "New eligible", "Use NEW_ELIGIBLE_RECORD only when the result is within scope, independently identifiable, and absent from ECOTOX.",
  "Excluded", "Use PRESENT_IN_ECOTOX_BUT_EXCLUDED when the result exists in ECOTOX but was already removed by an explicit thesis rule.",
  "Unresolved", "Keep UNRESOLVED when evidence is insufficient; do not force a decision.",
  "Completion", "Set Review_Status to COMPLETE only after Final_Decision, reason, and evidence location have been entered.",
  "Output protection", "This script refuses to overwrite an existing manual-review workbook unless OVERWRITE_EXISTING is deliberately set to TRUE.",
  "Historical field name", "In_Harmonized_295 is retained from the original audit workbook as a provenance field; it does not redefine the final 304-result combined dataset.",
  "Input workbook", basename(input_file),
  "Created", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

progress_labels <- c(
  "Total sources",
  "Sources complete",
  "Total records",
  "Records complete",
  "Duplicate ECOTOX",
  "New eligible record",
  "Present in ECOTOX but excluded",
  "Exclude out of scope",
  "Exclude method",
  "Unresolved"
)

progress <- tibble(Metric = progress_labels, Value = NA_character_)

# ---- 9. Workbook styles ------------------------------------------------------
wb <- createWorkbook(creator = "ESS Thesis reproducibility workflow")

header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE,
  border = "Bottom", borderColour = "#FFFFFF"
)

instruction_style <- createStyle(
  fontColour = "#1F1F1F", fgFill = "#D9EAF7", wrapText = TRUE,
  valign = "top"
)

manual_style <- createStyle(
  fgFill = "#FFF2CC", wrapText = TRUE, valign = "top",
  border = "LeftRightTopBottom", borderColour = "#E6B800"
)

automatic_style <- createStyle(
  fgFill = "#F2F2F2", wrapText = TRUE, valign = "top"
)

complete_style <- createStyle(fgFill = "#E2F0D9")
blocked_style <- createStyle(fgFill = "#F4CCCC")
date_style <- createStyle(numFmt = "yyyy-mm-dd")
decimal_style <- createStyle(numFmt = "0.000000")
percent_style <- createStyle(numFmt = "0.0%")

write_table_sheet <- function(wb, sheet_name, data, tab_colour = NULL) {
  addWorksheet(wb, sheet_name, tabColour = tab_colour)
  writeData(wb, sheet_name, data, headerStyle = header_style, withFilter = TRUE)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(data)), widths = "auto")
  if (nrow(data) > 0) {
    addStyle(
      wb, sheet_name, createStyle(valign = "top"),
      rows = 2:(nrow(data) + 1), cols = seq_len(ncol(data)),
      gridExpand = TRUE, stack = TRUE
    )
  }
}

# ---- 10. Write README --------------------------------------------------------
addWorksheet(wb, "README", tabColour = "#1F4E78")
writeData(wb, "README", readme, headerStyle = header_style)
freezePane(wb, "README", firstRow = TRUE)
setColWidths(wb, "README", cols = 1, widths = 24)
setColWidths(wb, "README", cols = 2, widths = 95)
addStyle(
  wb, "README", instruction_style,
  rows = 2:(nrow(readme) + 1), cols = 1:2,
  gridExpand = TRUE, stack = TRUE
)
setRowHeights(wb, "README", rows = 2:(nrow(readme) + 1), heights = 34)

# Create the Progress sheet here so the workbook opens in the intended order.
# Its formulas are written after Source_Review and Record_Review exist.
addWorksheet(wb, "Progress", tabColour = "#4472C4")
writeData(wb, "Progress", progress, headerStyle = header_style)
setColWidths(wb, "Progress", cols = 1, widths = 36)
setColWidths(wb, "Progress", cols = 2, widths = 18)
freezePane(wb, "Progress", firstRow = TRUE)

# ---- 11. Write source and record review sheets ------------------------------
write_table_sheet(wb, "Source_Review", source_review, tab_colour = "#70AD47")
write_table_sheet(wb, "Record_Review", record_review, tab_colour = "#FFC000")

source_manual_columns <- match(
  c(
    "Full_Text_Status", "Full_Text_or_URL", "Bibliographic_Identity",
    "Source_Evidence_Location", "Source_Level_Notes", "Source_Review_Complete"
  ),
  names(source_review)
)

record_manual_columns <- match(
  c(
    "Review_Status", "Verified_Species", "Verified_Medium",
    "Verified_Endpoint", "Verified_Exposure", "Verified_Duration_Hours",
    "Verified_LC50_Value", "Verified_Original_Unit", "Concentration_Basis",
    "Same_Experiment_as_ECOTOX", "Original_Source_Page_Table", "Evidence_Note",
    "Final_Decision", "Final_Decision_Reason", "Reviewer_Notes",
    "Reviewed_By", "Review_Date"
  ),
  names(record_review)
)

source_auto_columns <- setdiff(seq_len(ncol(source_review)), source_manual_columns)
record_auto_columns <- setdiff(seq_len(ncol(record_review)), record_manual_columns)

addStyle(
  wb, "Source_Review", automatic_style,
  rows = 2:(nrow(source_review) + 1), cols = source_auto_columns,
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Source_Review", manual_style,
  rows = 2:(nrow(source_review) + 1), cols = source_manual_columns,
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Record_Review", automatic_style,
  rows = 2:(nrow(record_review) + 1), cols = record_auto_columns,
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "Record_Review", manual_style,
  rows = 2:(nrow(record_review) + 1), cols = record_manual_columns,
  gridExpand = TRUE, stack = TRUE
)

# Practical widths: automatic evidence remains visible without producing
# extremely wide sheets; yellow manual-entry fields receive usable space.
setColWidths(wb, "Source_Review", cols = match("Source", names(source_review)), widths = 65)
setColWidths(wb, "Source_Review", cols = match("Best_ECOTOX_Title", names(source_review)), widths = 45)
setColWidths(wb, "Source_Review", cols = match("Verification_Task", names(source_review)), widths = 55)
setColWidths(wb, "Source_Review", cols = match(c("Full_Text_or_URL", "Source_Level_Notes"), names(source_review)), widths = 35)

setColWidths(wb, "Record_Review", cols = match("Source", names(record_review)), widths = 60)
setColWidths(wb, "Record_Review", cols = match("Automatic_Interpretation", names(record_review)), widths = 48)
setColWidths(
  wb, "Record_Review",
  cols = match(c("Evidence_Note", "Final_Decision_Reason", "Reviewer_Notes"), names(record_review)),
  widths = 35
)
setColWidths(wb, "Record_Review", cols = record_manual_columns, widths = 22)
setRowHeights(wb, "Source_Review", rows = 2:(nrow(source_review) + 1), heights = 68)
setRowHeights(wb, "Record_Review", rows = 2:(nrow(record_review) + 1), heights = 54)

# Numeric formatting.
source_ref_col <- match("Matched_ECOTOX_Reference", names(source_review))
addStyle(
  wb, "Source_Review", createStyle(numFmt = "0"),
  rows = 2:(nrow(source_review) + 1), cols = source_ref_col,
  gridExpand = TRUE, stack = TRUE
)

decimal_columns <- match(
  c("EnviroTox_LC50_mg_L", "ECOTOX_LC50_mg_L", "Verified_LC50_Value"),
  names(record_review)
)
decimal_columns <- decimal_columns[!is.na(decimal_columns)]
addStyle(
  wb, "Record_Review", decimal_style,
  rows = 2:(nrow(record_review) + 1), cols = decimal_columns,
  gridExpand = TRUE, stack = TRUE
)

relative_col <- match("Relative_Value_Difference", names(record_review))
addStyle(
  wb, "Record_Review", percent_style,
  rows = 2:(nrow(record_review) + 1), cols = relative_col,
  gridExpand = TRUE, stack = TRUE
)

date_col <- match("Review_Date", names(record_review))
addStyle(
  wb, "Record_Review", date_style,
  rows = 2:(nrow(record_review) + 1), cols = date_col,
  gridExpand = TRUE, stack = TRUE
)

# ---- 12. Write dropdown-list sheet and apply validations --------------------
write_table_sheet(wb, "Decision_Lists", lists_sheet, tab_colour = "#A5A5A5")
writeData(
  wb, "Decision_Lists",
  "These controlled vocabularies feed the dropdowns. Do not rename their headers.",
  startRow = nrow(lists_sheet) + 4, startCol = 1
)

list_range <- function(column_name) {
  column_number <- match(column_name, names(lists_sheet))
  final_row <- length(decision_lists[[column_name]]) + 1
  paste0(
    "'Decision_Lists'!$", int2col(column_number), "$2:$",
    int2col(column_number), "$", final_row
  )
}

apply_dropdown <- function(sheet, data, column_name, list_name) {
  column_number <- match(column_name, names(data))
  if (is.na(column_number)) stop("Dropdown target column not found: ", column_name)
  dataValidation(
    wb, sheet,
    cols = column_number,
    rows = 2:(nrow(data) + 1),
    type = "list",
    value = list_range(list_name),
    allowBlank = TRUE
  )
}

apply_dropdown("Source_Review", source_review, "Full_Text_Status", "Full_Text_Status")
apply_dropdown("Source_Review", source_review, "Bibliographic_Identity", "Bibliographic_Identity")
apply_dropdown("Source_Review", source_review, "Source_Review_Complete", "Yes_No")

apply_dropdown("Record_Review", record_review, "Review_Status", "Review_Status")
apply_dropdown("Record_Review", record_review, "Verified_Species", "Yes_No_Uncertain")
apply_dropdown("Record_Review", record_review, "Verified_Medium", "Verified_Medium")
apply_dropdown("Record_Review", record_review, "Verified_Endpoint", "Verified_Endpoint")
apply_dropdown("Record_Review", record_review, "Verified_Exposure", "Verified_Exposure")
apply_dropdown("Record_Review", record_review, "Concentration_Basis", "Concentration_Basis")
apply_dropdown("Record_Review", record_review, "Same_Experiment_as_ECOTOX", "Yes_No_Uncertain")
apply_dropdown("Record_Review", record_review, "Final_Decision", "Final_Decision")

# Highlight completed and blocked record statuses.
review_status_col <- match("Review_Status", names(record_review))
review_status_letter <- int2col(review_status_col)
conditionalFormatting(
  wb, "Record_Review",
  cols = seq_len(ncol(record_review)),
  rows = 2:(nrow(record_review) + 1),
  type = "expression",
  rule = paste0("$", review_status_letter, "2=\"COMPLETE\""),
  style = complete_style
)
conditionalFormatting(
  wb, "Record_Review",
  cols = seq_len(ncol(record_review)),
  rows = 2:(nrow(record_review) + 1),
  type = "expression",
  rule = paste0("$", review_status_letter, "2=\"BLOCKED\""),
  style = blocked_style
)

# ---- 13. Live progress summary ----------------------------------------------
source_complete_col <- int2col(match("Source_Review_Complete", names(source_review)))
record_complete_col <- int2col(match("Review_Status", names(record_review)))
decision_col <- int2col(match("Final_Decision", names(record_review)))

source_first_row <- 2
source_last_row <- nrow(source_review) + 1
record_first_row <- 2
record_last_row <- nrow(record_review) + 1

progress_formulas <- c(
  paste0("=COUNTA('Source_Review'!$A$", source_first_row, ":$A$", source_last_row, ")"),
  paste0("=COUNTIF('Source_Review'!$", source_complete_col, "$", source_first_row,
         ":$", source_complete_col, "$", source_last_row, ",\"YES\")"),
  paste0("=COUNTA('Record_Review'!$A$", record_first_row, ":$A$", record_last_row, ")"),
  paste0("=COUNTIF('Record_Review'!$", record_complete_col, "$", record_first_row,
         ":$", record_complete_col, "$", record_last_row, ",\"COMPLETE\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"DUPLICATE_ECOTOX\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"NEW_ELIGIBLE_RECORD\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"PRESENT_IN_ECOTOX_BUT_EXCLUDED\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"EXCLUDE_OUT_OF_SCOPE\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"EXCLUDE_METHOD\")"),
  paste0("=COUNTIF('Record_Review'!$", decision_col, "$", record_first_row,
         ":$", decision_col, "$", record_last_row, ",\"UNRESOLVED\")")
)

writeFormula(wb, "Progress", x = progress_formulas, startCol = 2, startRow = 2)
addStyle(
  wb, "Progress", instruction_style,
  rows = 2:(nrow(progress) + 1), cols = 1:2,
  gridExpand = TRUE, stack = TRUE
)

# ---- 14. Save ---------------------------------------------------------------
saveWorkbook(wb, output_file, overwrite = isTRUE(OVERWRITE_EXISTING))


preparation_summary <- tibble::tibble(
  Metric = c(
    "Manual-review records",
    "Unique source citations",
    "Priority 1 records",
    "Priority 2 records",
    "Priority 3 records",
    "Priority 4 records"
  ),
  N = c(
    nrow(record_review),
    nrow(source_review),
    sum(record_review$Priority == 1, na.rm = TRUE),
    sum(record_review$Priority == 2, na.rm = TRUE),
    sum(record_review$Priority == 3, na.rm = TRUE),
    sum(record_review$Priority == 4, na.rm = TRUE)
  )
)

readr::write_csv(
  preparation_summary,
  file.path(
    output_dir,
    "manual_review_queue_summary.csv"
  )
)

input_output_manifest <- tibble::tibble(
  Role = c(
    "Input: automatic EnviroTox-ECOTOX overlap workbook",
    "Local editable output: manual-review workbook",
    "Output: manual-review queue summary"
  ),
  Path = c(
    input_file,
    output_file,
    file.path(
      output_dir,
      "manual_review_queue_summary.csv"
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
    "STATUS: ENVIROTOX MANUAL-REVIEW WORKBOOK PREPARATION = PASS",
    paste0("Records prepared for review: ", nrow(record_review)),
    paste0("Sources prepared for review: ", nrow(source_review)),
    "This status confirms workbook preparation only; it does not indicate that manual source review is complete."
  ),
  file.path(
    output_dir,
    "WORKBOOK_PREPARED.txt"
  )
)

cat("\nEnviroTox manual-review workbook prepared successfully.\n")
cat("Sources prepared for review:", nrow(source_review), "\n")
cat("Records prepared for review:", nrow(record_review), "\n")
cat(
  "Local editable workbook written to:\n",
  normalizePath(
    output_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat(
  "STATUS: ENVIROTOX MANUAL-REVIEW WORKBOOK PREPARATION = PASS\n"
)
