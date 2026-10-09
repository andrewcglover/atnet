# Reading the sporozoite prevalence assays from the laboratory workbooks.
#
# A sheet holds one column per experimental group, with a header naming the
# group and one row per dissected mosquito below it. Within a column, blank
# cells separate sub-groups (typically replicates); the columns are separated
# independently of one another, so a blank in one column does not break another.

#' Interpret an assay column header
#'
#' Headers name the arm, the time of net exposure relative to the infectious
#' blood meal, and the number of blood meals, for example `"CTL 72h 1BF"` or
#' `"200mg 6d"`. Where the number of blood meals is absent it is taken to be
#' one.
#'
#' @param header A single column header, which may contain line breaks.
#'
#' @return A list with elements `arm` (`"control"`, `"treated"` or `NA` for a
#'   column that is not an assay group), `exposure_h` (hours after the
#'   infectious blood meal, negative for exposure before it) and
#'   `n_bloodmeals`.
#'
#' @keywords internal
parse_assay_header <- function(header) {
  tokens <- strsplit(gsub("\\s+", " ", trimws(header)), " ", fixed = TRUE)[[1]]
  blank <- list(arm = NA_character_, exposure_h = NA_real_, n_bloodmeals = NA_integer_)
  if (!length(tokens)) {
    return(blank)
  }

  arm <- if (any(grepl("^CTL$", tokens, ignore.case = TRUE))) {
    "control"
  } else if (any(grepl("^[0-9.]+mg$", tokens, ignore.case = TRUE))) {
    "treated"
  } else {
    NA_character_
  }
  if (is.na(arm)) {
    return(blank)
  }

  time_token <- tokens[grepl("^-?[0-9.]+[hd]$", tokens)]
  exposure_h <- if (!length(time_token)) {
    NA_real_
  } else {
    value <- as.numeric(sub("[hd]$", "", time_token[[1]]))
    if (grepl("d$", time_token[[1]])) value * 24 else value
  }

  bf_token <- tokens[grepl("^[0-9]BF$", tokens, ignore.case = TRUE)]
  n_bloodmeals <- if (length(bf_token)) {
    as.integer(sub("BF$", "", bf_token[[1]], ignore.case = TRUE))
  } else {
    1L
  }

  list(arm = arm, exposure_h = exposure_h, n_bloodmeals = n_bloodmeals)
}

#' Split a column into sub-groups at its blank cells
#'
#' Trailing blanks are discarded; runs of blanks within the column separate
#' consecutive sub-groups.
#'
#' @param x A character vector holding one column of a sheet, blanks as `NA`.
#'
#' @return A list of character vectors, one per sub-group, in order. A column
#'   that is entirely blank gives an empty list.
#'
#' @keywords internal
split_at_blanks <- function(x) {
  filled <- !is.na(x)
  if (!any(filled)) {
    return(list())
  }
  x <- x[seq_len(max(which(filled)))]
  blank <- is.na(x)
  starts_run <- blank & !c(FALSE, utils::head(blank, -1))
  group <- cumsum(starts_run)
  unname(split(x[!blank], group[!blank]))
}

#' Count dissected and infected mosquitoes in a sub-group
#'
#' A cell counts as a dissected mosquito when it holds a number, and as an
#' infected mosquito when that number is above zero. Trailing annotation such
#' as the asterisk in `"0*"` is ignored, so an annotated zero still counts as a
#' dissected, uninfected mosquito.
#'
#' @param x A character vector of cells from one sub-group.
#'
#' @return A named integer vector with elements `n` and `n_positive`.
#'
#' @keywords internal
count_group <- function(x) {
  value <- suppressWarnings(as.numeric(sub("[^0-9.]+$", "", x)))
  c(n = sum(!is.na(value)), n_positive = sum(value > 0, na.rm = TRUE))
}

#' Tabulate the assay groups on one sheet
#'
#' @param sheet A data frame or matrix holding the sheet exactly as stored,
#'   with the column headers in the first row and no row names.
#' @param sheet_name Recorded in the result so that rows can be traced back.
#'
#' @return A data frame with one row per sub-group, with columns `sheet`,
#'   `arm`, `exposure_h`, `n_bloodmeals`, `subgroup` (its position within the
#'   column, counting from one), `n` and `n_positive`.
#'
#' @export
tabulate_assay_sheet <- function(sheet, sheet_name = NA_character_) {
  sheet <- as.data.frame(sheet, stringsAsFactors = FALSE)
  if (nrow(sheet) < 2L) {
    return(empty_assay_table())
  }
  headers <- as.character(unlist(sheet[1L, ], use.names = FALSE))
  body <- sheet[-1L, , drop = FALSE]

  out <- lapply(seq_along(headers), function(j) {
    meta <- parse_assay_header(headers[[j]])
    if (is.na(meta$arm)) {
      return(NULL)
    }
    groups <- split_at_blanks(as.character(body[[j]]))
    if (!length(groups)) {
      return(NULL)
    }
    counts <- vapply(groups, count_group, c(n = 0L, n_positive = 0L))
    data.frame(
      sheet = sheet_name,
      arm = meta$arm,
      exposure_h = meta$exposure_h,
      n_bloodmeals = meta$n_bloodmeals,
      subgroup = seq_along(groups),
      n = as.integer(counts["n", ]),
      n_positive = as.integer(counts["n_positive", ]),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  if (is.null(out)) empty_assay_table() else out
}

empty_assay_table <- function() {
  data.frame(
    sheet = character(0), arm = character(0), exposure_h = numeric(0),
    n_bloodmeals = integer(0), subgroup = integer(0), n = integer(0),
    n_positive = integer(0), stringsAsFactors = FALSE
  )
}

#' Tabulate the assay groups in a workbook
#'
#' @param path Path to an `.xlsx` workbook.
#' @param sheets Sheets to read. The default reads every sheet except `README`.
#'
#' @return A data frame as returned by [tabulate_assay_sheet()], stacked over
#'   the sheets read.
#'
#' @export
tabulate_assay_workbook <- function(path, sheets = NULL) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("readxl is needed to read assay workbooks.", call. = FALSE)
  }
  if (is.null(sheets)) {
    sheets <- setdiff(readxl::excel_sheets(path), "README")
  }
  out <- lapply(sheets, function(s) {
    raw <- readxl::read_excel(path, sheet = s, col_names = FALSE,
                              col_types = "text", .name_repair = "minimal")
    tabulate_assay_sheet(raw, sheet_name = s)
  })
  do.call(rbind, c(out, list(empty_assay_table())))
}
