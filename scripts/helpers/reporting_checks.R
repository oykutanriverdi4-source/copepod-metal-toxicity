# Reporting helpers only: no model fitting, data exclusion or inference changes.
# Raw lme4 messages are retained. A singularity message is not re-labelled as
# an independent optimizer failure. "No alert recorded" is not model approval.
# See https://lme4.github.io/lme4/reference/isSingular.html and convergence.html.

report_text_values <- function(x) {
  out <- as.character(unlist(x, recursive = TRUE, use.names = FALSE))
  out <- out[!is.na(out)]
  out[nzchar(trimws(out))]
}

report_collapse_text <- function(x) {
  out <- report_text_values(x)
  if (!length(out)) NA_character_ else paste(out, collapse = " | ")
}

report_classify_optinfo <- function(optinfo, singular) {
  stopifnot(is.list(optinfo), is.logical(singular), length(singular) == 1L,
            !is.na(singular))
  messages <- report_text_values(optinfo$conv$lme4$messages)
  # Exact known messages only. Unknown wording stays in "other" rather than
  # being silently dismissed as a singularity message.
  known_singular_messages <- c(
    "boundary (singular) fit: see help('isSingular')",
    'boundary (singular) fit: see help("isSingular")',
    "boundary (singular) fit: see ?isSingular"
  )
  singular_message <- trimws(messages) %in% known_singular_messages
  other_messages <- messages[!singular_message]
  codes_raw <- unlist(optinfo$conv$opt, recursive = TRUE, use.names = FALSE)
  code_available <- length(codes_raw) > 0L
  codes_text <- as.character(codes_raw)
  codes <- suppressWarnings(as.numeric(codes_text))
  # Missing/empty entries in a supplied optimizer code are not discarded.
  codes_text[is.na(codes_text)] <- "NA"
  codes_text[!nzchar(trimws(codes_text))] <- "EMPTY"
  code_uninterpretable <- code_available && any(!is.finite(codes))
  nonzero_code <- any(is.finite(codes) & codes != 0)
  optimizer_warnings <- report_text_values(optinfo$warnings)
  other_alert <- length(other_messages) > 0L || nonzero_code ||
    code_uninterpretable || length(optimizer_warnings) > 0L
  singular_evidence <- singular || any(singular_message)
  status <- if (singular_evidence && other_alert) {
    "SINGULARITY_AND_OTHER_ALERT_RECORDED"
  } else if (singular_evidence) {
    "SINGULARITY_ONLY_RECORDED"
  } else if (other_alert) {
    "OTHER_NUMERICAL_ALERT_RECORDED"
  } else {
    "NO_NUMERICAL_ALERT_RECORDED"
  }
  data.frame(
    singular = singular,
    lme4_messages = report_collapse_text(messages),
    singularity_message = report_collapse_text(messages[singular_message]),
    other_lme4_message = report_collapse_text(other_messages),
    optimizer_code = report_collapse_text(codes_text),
    optimizer_code_available = code_available,
    optimizer_nonzero_code = nonzero_code,
    optimizer_code_uninterpretable = code_uninterpretable,
    optimizer_warning = report_collapse_text(optimizer_warnings),
    other_numerical_alert = other_alert,
    numerical_report_status = status,
    stringsAsFactors = FALSE
  )
}

report_fit_diagnostics <- function(model) {
  report_classify_optinfo(
    model@optinfo,
    singular = lme4::isSingular(model, tol = 1e-4)
  )
}

# Counts refer to reported ratios and interval endpoints, not to new tests.
# Missing, invalid or absent refits are counted explicitly, not as agreement.
report_ratio_support <- function(dat, ratio, lower, upper, id,
                                 expected_refits, interval_adjustment) {
  required <- c(ratio, lower, upper, id)
  if (!is.data.frame(dat) || !all(required %in% names(dat))) {
    stop("Missing required columns for the ratio-support report.")
  }
  if (length(expected_refits) != 1L || !is.finite(expected_refits) ||
      expected_refits < 0 || expected_refits != as.integer(expected_refits)) {
    stop("expected_refits must be one non-negative integer.")
  }
  ids <- as.character(dat[[id]])
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("Each deletion refit must have one non-missing, unique identifier.")
  }
  if (nrow(dat) > expected_refits) {
    stop("More refits supplied than the recorded number of groups.")
  }
  values <- dat[[ratio]]
  lo <- dat[[lower]]
  hi <- dat[[upper]]
  if (!is.numeric(values) || !is.numeric(lo) || !is.numeric(hi)) {
    stop("Ratio and interval columns must be numeric.")
  }
  valid_point <- is.finite(values) & values > 0
  valid_interval <- valid_point & is.finite(lo) & is.finite(hi) &
    lo > 0 & lo <= hi & lo <= values & values <= hi
  n_valid <- sum(valid_interval)
  below <- sum(hi[valid_interval] < 1)
  above <- sum(lo[valid_interval] > 1)
  includes <- sum(lo[valid_interval] <= 1 & hi[valid_interval] >= 1)
  all_intervals_available <- expected_refits > 0L &&
    nrow(dat) == expected_refits && n_valid == expected_refits
  data.frame(
    expected_refits = as.integer(expected_refits),
    returned_refits = nrow(dat),
    absent_refits = as.integer(expected_refits - nrow(dat)),
    valid_point_estimates = sum(valid_point),
    point_estimates_below_1 = sum(values[valid_point] < 1),
    point_estimates_above_1 = sum(values[valid_point] > 1),
    point_estimates_equal_1 = sum(values[valid_point] == 1),
    valid_intervals = n_valid,
    intervals_below_1 = below,
    intervals_above_1 = above,
    intervals_including_1 = includes,
    unavailable_or_invalid_intervals = nrow(dat) - n_valid,
    all_expected_intervals_below_1 = if (all_intervals_available) {
      below == expected_refits
    } else NA,
    all_expected_intervals_above_1 = if (all_intervals_available) {
      above == expected_refits
    } else NA,
    interval_level = 0.95,
    interval_adjustment = interval_adjustment,
    interpretation = paste(
      "Overlapping deletion checks, not independent confirmations;",
      "interval inclusion is reported separately from point direction."
    ),
    stringsAsFactors = FALSE
  )
}
