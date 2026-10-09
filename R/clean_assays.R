# Turning the workbooks into the two tables the fits consume.
#
# Each output row is one experiment: a control arm and a treated arm from the
# same infectious feed, under one set of conditions. Every row carries the
# workbook, worksheet and columns it came from, so a number in a table can
# always be traced back to a cell in a spreadsheet.

#' Read the dissection day from a worksheet name
#'
#' Sporozoite worksheets end in `spz10` or `spz13`, naming the day after the
#' infectious blood meal on which the salivary glands were dissected. This is
#' taken from the worksheet name rather than a typed column, since the typed
#' column carried a transcription slip.
#'
#' @param sheet_name A worksheet name.
#'
#' @return The dissection day, or `NA` for a worksheet that names none.
#'
#' @export
dissection_day <- function(sheet_name) {
  m <- regmatches(sheet_name, regexpr("spz[0-9]+", sheet_name))
  if (!length(m)) {
    return(NA_integer_)
  }
  as.integer(sub("spz", "", m))
}

#' Whether exposure came before or after the infectious blood meal
#'
#' Timing is counted forward from the infectious blood meal throughout, so
#' exposure after it is positive and exposure before it is negative. This
#' follows the supplementary information and the fitting code, where a positive
#' delay means the mosquitoes were already infected when they met the net.
#'
#' Beware that the oocyst file the blocking model was previously fitted from
#' used the opposite sign, counting backwards from the blood meal.
#'
#' @param sheet_name A worksheet name.
#'
#' @return `1` where the worksheet says exposure followed the blood meal, `-1`
#'   where it preceded it, and `NA` where it says neither.
#'
#' @keywords internal
exposure_sign <- function(sheet_name) {
  if (grepl("post", sheet_name, ignore.case = TRUE)) {
    1
  } else if (grepl("pre", sheet_name, ignore.case = TRUE)) {
    -1
  } else {
    NA_real_
  }
}

#' Pair the control and treated arms of a tabulated worksheet
#'
#' Arms are paired within a worksheet by exposure time, resting duration, number
#' of blood meals and position down the column, which together identify one
#' experiment.
#'
#' @param tab A table as returned by [tabulate_assay_sheet()].
#' @param concentration Fit only treated arms at this concentration, in mg per
#'   square metre.
#'
#' @return A data frame with one row per experiment.
#'
#' @keywords internal
pair_arms <- function(tab, concentration = 200) {
  if (!nrow(tab)) {
    return(empty_paired_table())
  }
  tab$duration_key <- ifelse(is.na(tab$duration_min), -1, tab$duration_min)
  tab$exposure_key <- ifelse(is.na(tab$exposure_h), -1, tab$exposure_h)
  key <- paste(tab$sheet, tab$exposure_key, tab$duration_key,
               tab$n_bloodmeals, tab$subgroup, sep = "\r")

  out <- lapply(split(tab, key), function(g) {
    ctl <- g[g$arm == "control", ]
    trt <- g[g$arm == "treated" & !is.na(g$concentration) &
               g$concentration == concentration, ]
    if (nrow(ctl) != 1L || nrow(trt) != 1L) {
      return(NULL)
    }
    data.frame(
      worksheet = ctl$sheet,
      control_column = ctl$column,
      treated_column = trt$column,
      exposure_h = ctl$exposure_h,
      duration_min = ctl$duration_min,
      n_bloodmeals = ctl$n_bloodmeals,
      replicate = ctl$subgroup,
      n_control = ctl$n, n_positive_control = ctl$n_positive,
      n_treated = trt$n, n_positive_treated = trt$n_positive,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) empty_paired_table() else out[order(out$worksheet, out$replicate), ]
}

empty_paired_table <- function() {
  data.frame(
    worksheet = character(0), control_column = character(0),
    treated_column = character(0), exposure_h = numeric(0), duration_min = numeric(0),
    n_bloodmeals = integer(0), replicate = integer(0), n_control = integer(0),
    n_positive_control = integer(0), n_treated = integer(0),
    n_positive_treated = integer(0), stringsAsFactors = FALSE
  )
}

#' Find columns holding the same data twice
#'
#' A column whose sub-groups repeat the same values has been entered twice.
#' Columns that take only one value, such as an arm in which nothing was
#' infected, are ignored, since those match one another trivially.
#'
#' @param values A list of character vectors, one per sub-group of a column.
#'
#' @return The indices of sub-groups that repeat an earlier one.
#'
#' @export
repeated_subgroups <- function(values) {
  if (length(values) < 2L) {
    return(integer(0))
  }
  informative <- vapply(values, function(v) length(unique(v)) > 1L, logical(1))
  keys <- vapply(values, function(v) paste(sort(v), collapse = "\r"), character(1))
  which(duplicated(keys) & informative)
}
