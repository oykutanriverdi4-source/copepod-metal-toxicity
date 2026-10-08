# Isolated reporting tests. This file never sources an analysis script,
# fits a model, changes data, installs packages, or contacts GitHub/Zenodo.
# Execute only from a deliberately selected review/test project.
run_reporting_unit_tests <- function(root = ".") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  helper_file <- file.path(root, "scripts", "helpers", "reporting_checks.R")
  if (!file.exists(helper_file)) stop("Reporting test package root not found.")
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr is required for replaying the two existing row-selection expressions.")
  }
  helper_env <- new.env(parent = baseenv())
  sys.source(helper_file, envir = helper_env)
  classify <- helper_env$report_classify_optinfo
  counts <- helper_env$report_ratio_support
  passed <- character()
  expect <- function(name, condition) {
    if (!isTRUE(condition)) stop(paste("FAILED:", name))
    passed <<- c(passed, name)
  }
  expect_error <- function(name, expression) {
    got_error <- tryCatch({ force(expression); FALSE }, error = function(e) TRUE)
    expect(name, got_error)
  }
  opt <- function(messages = NULL, code = 0, warnings = NULL) {
    list(conv = list(lme4 = list(messages = messages), opt = code), warnings = warnings)
  }
  sing <- "boundary (singular) fit: see help('isSingular')"
  x <- classify(opt(), FALSE)
  expect("No alert with zero optimizer code", !x$other_numerical_alert && !x$singular)
  x <- classify(opt(sing), TRUE)
  expect("Singularity is preserved separately", x$singular && !x$other_numerical_alert && !is.na(x$singularity_message))
  expect("Raw singularity message remains available", identical(x$lme4_messages, sing))
  x <- classify(opt(c(sing, "Model failed to converge with max|grad|")), TRUE)
  expect("Singularity plus gradient alert retains both", x$other_numerical_alert && !is.na(x$other_lme4_message) && !is.na(x$singularity_message))
  x <- classify(opt(code = 1), FALSE)
  expect("Nonzero optimizer code is retained without lme4 message", x$other_numerical_alert && x$optimizer_nonzero_code)
  x <- classify(opt(warnings = "optimizer warning"), FALSE)
  expect("Optimizer warning is retained", x$other_numerical_alert && !is.na(x$optimizer_warning))
  x <- classify(opt(code = NULL), FALSE)
  expect("Absent optimizer code is explicitly marked", !x$optimizer_code_available && !x$optimizer_nonzero_code)
  x <- classify(opt(code = c(0, NA_real_)), FALSE)
  expect("Missing supplied optimizer code is not silently discarded", x$optimizer_code_uninterpretable && x$other_numerical_alert)
  x <- classify(opt(code = "not a code"), FALSE)
  expect("Unexpected optimizer code is not called success", x$optimizer_code_uninterpretable && x$other_numerical_alert)
  x <- classify(opt("boundary (singular) fit plus other unknown text"), TRUE)
  expect("Unknown message is not suppressed by broad matching", x$other_numerical_alert)
  x <- classify(opt(c(NA_character_, "", sing)), TRUE)
  expect("Empty messages do not create a separate alert", !x$other_numerical_alert && !is.na(x$singularity_message))

  z <- data.frame(id = letters[1:3], r = c(0.5, 2, 1), lo = c(0.2, 1.1, 0.8), hi = c(0.8, 3, 1.2))
  y <- counts(z, "r", "lo", "hi", "id", 3, "test only")
  expect("Direction and interval categories are counted separately", y$point_estimates_below_1 == 1 && y$point_estimates_above_1 == 1 && y$point_estimates_equal_1 == 1 && y$intervals_below_1 == 1 && y$intervals_above_1 == 1 && y$intervals_including_1 == 1)
  z$hi[1] <- 1
  y <- counts(z, "r", "lo", "hi", "id", 3, "test only")
  expect("An interval ending at 1 includes 1", y$intervals_below_1 == 0 && y$intervals_including_1 == 2)
  z$lo[2] <- NA_real_
  y <- counts(z, "r", "lo", "hi", "id", 3, "test only")
  expect("Missing intervals are not treated as agreement", y$valid_intervals == 2 && y$unavailable_or_invalid_intervals == 1 && is.na(y$all_expected_intervals_above_1))
  z$lo[2] <- 4
  y <- counts(z, "r", "lo", "hi", "id", 3, "test only")
  expect("Reversed bounds are invalid", y$unavailable_or_invalid_intervals == 1)
  y <- counts(z[0, ], "r", "lo", "hi", "id", 3, "test only")
  expect("Empty refits never become universal agreement", y$absent_refits == 3 && is.na(y$all_expected_intervals_below_1))
  y <- counts(z[1, ], "r", "lo", "hi", "id", 3, "test only")
  expect("Absent refits are explicit", y$absent_refits == 2 && is.na(y$all_expected_intervals_below_1))
  expect_error("Duplicate refit identifiers are rejected", counts(z[c(1, 1), ], "r", "lo", "hi", "id", 3, "test only"))
  expect_error("More refits than expected are rejected", counts(z, "r", "lo", "hi", "id", 2, "test only"))
  expect_error("Noninteger expected counts are rejected", counts(z, "r", "lo", "hi", "id", 3.5, "test only"))

  get_assignment <- function(exprs, name) {
    selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]) && as.character(x[[2L]]) == name, as.list(exprs))
    if (length(selected) != 1L) stop(paste("Unique assignment not found:", name))
    selected[[1L]]
  }
  expect("Helper parses in R", length(parse(helper_file, encoding = "UTF-8")) > 0L)
  # Evaluate ONLY the two patched row-selection assignments, never the model script.
  order_exprs <- parse(file.path(root, "scripts", "analysis", "13_fit_order_model.R"), encoding = "UTF-8")
  for (method in c("loro", "loso")) {
    env <- new.env(parent = asNamespace("dplyr"))
    name <- paste0(method, "_interaction_extremes")
    id <- if (method == "loro") "removed_reference" else "removed_species"
    for (vals in list(c(3, 1, 2), c(1, 1, 2), c(2, 2, 2), 1, numeric(), c(NA_real_, 3, 1))) {
      synthetic <- data.frame(ratio_of_ratios = vals)
      synthetic[[id]] <- if (length(vals)) paste0("id", seq_along(vals)) else character()
      assign(paste0(method, "_results"), synthetic, envir = env)
      eval(get_assignment(order_exprs, name), envir = env)
      got <- get(name, envir = env)
      valid <- which(!is.na(vals))
      expected_ids <- if (length(valid)) unique(c(valid[which.min(vals[valid])], valid[which.max(vals[valid])])) else integer()
      expect(paste("Order", method, "synthetic selection", paste(vals, collapse = ",")), identical(as.character(got[[id]]), synthetic[[id]][expected_ids]))
    }
  }
  expect("All reporting unit tests reached completion", TRUE)
  data.frame(test = passed, status = "PASS", stringsAsFactors = FALSE)
}
