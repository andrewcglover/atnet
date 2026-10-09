# Building the two tables the fits consume, from the workbooks in
# data_private/fitting/.
#
# Each output row is one experiment: a control arm and a 200 mg/m2 arm from the
# same infectious feed under one set of conditions, with the workbook, worksheet
# and columns it came from. The corrections listed by assay_corrections() are
# applied here and nowhere else.

#' Collapse a worksheet of the per-replicate workbook to one group per column
#'
#' In that workbook a worksheet is one replicate, so the blank cells within a
#' column do not separate experiments and the sub-groups belong together. Any
#' sub-group repeating an earlier one is dropped first, since a repeat is the
#' same mosquitoes entered twice rather than more of them.
#'
#' @param path Path to the workbook.
#' @param sheet A worksheet name.
#'
#' @return A table in the form of [tabulate_assay_sheet()], with one row per
#'   column and `subgroup` set to one.
#'
#' @keywords internal
collapse_replicate_sheet <- function(path, sheet) {
  raw <- suppressMessages(readxl::read_excel(path, sheet = sheet,
                                             col_names = FALSE,
                                             col_types = "text",
                                             .name_repair = "minimal"))
  tab <- tabulate_assay_sheet(raw, sheet)
  if (!nrow(tab)) {
    return(tab)
  }

  # Matched against tab$column, which is normalised, so normalise here too.
  headers <- gsub("\\s+", " ", trimws(as.character(unlist(raw[1L, ],
                                                          use.names = FALSE))))
  body <- raw[-1L, , drop = FALSE]
  drop <- lapply(seq_along(headers), function(j) {
    groups <- split_at_blanks(as.character(body[[j]]))
    if (length(groups) < 2L) {
      return(NULL)
    }
    rep_idx <- repeated_subgroups(groups)
    if (!length(rep_idx)) {
      return(NULL)
    }
    data.frame(column = headers[[j]], subgroup = rep_idx, stringsAsFactors = FALSE)
  })
  drop <- do.call(rbind, drop)
  if (!is.null(drop)) {
    keep <- !paste(tab$column, tab$subgroup) %in% paste(drop$column, drop$subgroup)
    tab <- tab[keep, , drop = FALSE]
  }

  key <- paste(tab$column, tab$arm, tab$concentration, tab$exposure_h,
               tab$duration_min, tab$n_bloodmeals, sep = "\r")
  out <- lapply(split(tab, key), function(g) {
    first <- g[1L, ]
    first$subgroup <- 1L
    first$n <- sum(g$n)
    first$n_positive <- sum(g$n_positive)
    first
  })
  do.call(rbind, out)
}

#' Build the sporozoite table for the extrinsic incubation period fit
#'
#' @param dir Folder holding the workbooks, normally `data_private/fitting`.
#' @param concentration Fit only treated arms at this concentration.
#' @param checks The counts [missing_day10_replicate()] expects, kept with the
#'   data since they are data.
#'
#' @return A data frame with one row per experiment.
#'
#' @export
build_sporozoite_table <- function(dir = "data_private/fitting",
                                   concentration = 200,
                                   checks = file.path(dirname(dir),
                                                      "expected_day10_replicate.csv")) {
  per_replicate <- find_workbook(dir, "per_replicate")
  pooled <- find_workbook(dir, "pooled")

  sheets <- readxl::excel_sheets(per_replicate)
  sheets <- sheets[grepl("spz", sheets)]
  # Holds a copy of another worksheet; the genuine replicate is taken from the
  # pooled workbook below.
  sheets <- setdiff(sheets, "72h 2BF post 09.07.26 spz10")

  rows <- lapply(sheets, function(s) {
    paired <- pair_arms(collapse_replicate_sheet(per_replicate, s), concentration)
    if (!nrow(paired)) {
      return(NULL)
    }
    paired$workbook <- basename(per_replicate)
    paired$dissection_day <- dissection_day(s)
    paired
  })
  out <- do.call(rbind, rows)
  out <- rbind(out, missing_day10_replicate(pooled, checks))

  # Counted forward from the infectious blood meal, as the supplementary
  # information and the fitting code do.
  out$hours_after_infection <- out$exposure_h *
    vapply(out$worksheet, exposure_sign, numeric(1))
  out$exposure_h <- NULL
  out$replicate <- NULL

  out <- out[order(out$dissection_day, out$hours_after_infection,
                   out$n_bloodmeals, out$worksheet), ]
  finalise_experiment_table(out, one_row_each = TRUE)
}

#' Number the experiments, label them, and write the timing both ways round
#'
#' The fits index experiments by a number from one upwards. The label repeats
#' the fields that define an experiment, so a row can be read without going back
#' to the other columns, and is named for the fields it holds.
#'
#' The timing is also written counting backwards from the infectious blood meal,
#' as `hours_before_infection`, because the file the blocking model was
#' previously fitted from carried both. It is derived here rather than recorded,
#' so the two can never disagree.
#'
#' @param x A table from [build_sporozoite_table()] or [build_oocyst_table()].
#' @param one_row_each Whether each row is already one experiment.
#'
#' @return `x` with `experiment_id`, `hours_before_infection` and a labelled
#'   column added.
#'
#' @keywords internal
finalise_experiment_table <- function(x, one_row_each) {
  x$hours_before_infection <- -x$hours_after_infection
  stopifnot(all(x$hours_after_infection + x$hours_before_infection == 0))
  label <- paste(x$workbook, x$worksheet, x$hours_after_infection,
                 x$n_bloodmeals, if (one_row_each) 1L else x$subgroup,
                 sep = " | ")
  x$experiment_id <- if (one_row_each) {
    seq_len(nrow(x))
  } else {
    match(label, unique(label))
  }
  x[[paste0("experiment (workbook | worksheet | hours_after_infection | ",
            "n_bloodmeals | subgroup)")]] <- label
  rownames(x) <- NULL
  x
}

#' The day-10 replicate held only in the pooled workbook
#'
#' One replicate was entered in the pooled workbook but not in the
#' per-replicate one, whose corresponding worksheet holds a copy of an earlier
#' replicate instead. It is the last group down each column of
#' `#7 72h post spz10`. The counts are checked against the values that make the
#' two workbooks reconcile, so a changed workbook raises an error rather than
#' passing silently. Those values are data, so they are read from a file kept
#' with the workbooks rather than written here.
#'
#' @param pooled Path to the pooled workbook.
#' @param checks Path to a CSV with one row per number of blood meals, giving
#'   `n_bloodmeals`, `n_control`, `n_positive_control`, `n_treated` and
#'   `n_positive_treated` as this replicate should have them.
#'
#' @return A one-worksheet data frame in the form of [build_sporozoite_table()].
#'
#' @keywords internal
missing_day10_replicate <- function(pooled, checks) {
  if (!file.exists(checks)) {
    stop("The expected counts for the day-10 replicate are missing: ", checks,
         call. = FALSE)
  }
  expected <- utils::read.csv(checks)
  stopifnot(nrow(expected) > 0L)
  tab <- tabulate_assay_workbook(pooled, "#7 72h post spz10")
  last <- do.call(rbind, lapply(
    split(tab, paste(tab$column, tab$n_bloodmeals)),
    function(g) g[which.max(g$subgroup), ]))

  out <- lapply(seq_len(nrow(expected)), function(i) {
    e <- expected[i, ]
    ctl <- last[last$arm == "control" & last$n_bloodmeals == e$n_bloodmeals, ]
    trt <- last[last$arm == "treated" & last$n_bloodmeals == e$n_bloodmeals, ]
    stopifnot(nrow(ctl) == 1L, nrow(trt) == 1L)
    if (!identical(c(ctl$n, ctl$n_positive, trt$n, trt$n_positive),
                   as.integer(c(e$n_control, e$n_positive_control,
                                e$n_treated, e$n_positive_treated)))) {
      stop("The last group of '#7 72h post spz10' is not the replicate this ",
           "correction was written for. Re-check the workbook against ",
           "data_private/DATA_SOURCES.md before fitting.", call. = FALSE)
    }
    data.frame(
      worksheet = "#7 72h post spz10",
      control_column = ctl$column, treated_column = trt$column,
      exposure_h = 72, duration_min = NA_real_,
      n_bloodmeals = as.integer(e$n_bloodmeals), replicate = 4L,
      n_control = ctl$n, n_positive_control = ctl$n_positive,
      n_treated = trt$n, n_positive_treated = trt$n_positive,
      workbook = basename(pooled), dissection_day = 10L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
