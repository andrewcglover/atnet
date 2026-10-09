# The oocyst table for the transmission-blocking fit.
#
# That fit works on the oocyst burden of each mosquito rather than on counts of
# infected mosquitoes, so this table has one row per mosquito.

#' Read a worksheet as one row per mosquito
#'
#' @param path Path to the workbook.
#' @param sheet A worksheet name.
#'
#' @return A data frame with one row per mosquito, carrying the worksheet and
#'   column it came from, the conditions, the sub-group it sat in, and its
#'   oocyst burden.
#'
#' @keywords internal
assay_values <- function(path, sheet) {
  raw <- suppressMessages(readxl::read_excel(path, sheet = sheet,
                                             col_names = FALSE,
                                             col_types = "text",
                                             .name_repair = "minimal"))
  if (nrow(raw) < 2L) {
    return(NULL)
  }
  headers <- as.character(unlist(raw[1L, ], use.names = FALSE))
  body <- raw[-1L, , drop = FALSE]

  out <- lapply(seq_along(headers), function(j) {
    meta <- parse_assay_header(headers[[j]])
    if (is.na(meta$arm)) {
      return(NULL)
    }
    groups <- split_at_blanks(as.character(body[[j]]))
    if (!length(groups)) {
      return(NULL)
    }
    do.call(rbind, lapply(seq_along(groups), function(k) {
      v <- suppressWarnings(as.numeric(sub("[^0-9.]+$", "", groups[[k]])))
      v <- v[!is.na(v)]
      if (!length(v)) {
        return(NULL)
      }
      data.frame(
        worksheet = sheet, column = gsub("\\s+", " ", trimws(headers[[j]])),
        arm = meta$arm,
        concentration = meta$concentration, exposure_h = meta$exposure_h,
        duration_min = meta$duration_min, n_bloodmeals = meta$n_bloodmeals,
        subgroup = k, oocysts = v, stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, out)
}

#' Build the oocyst table for the transmission-blocking fit
#'
#' @param dir Folder holding the workbooks, normally `data_private/fitting`.
#' @param concentration Fit only treated arms at this concentration.
#' @param duration Fit only the standard rest on the net, in minutes.
#'
#' @return A data frame with one row per mosquito.
#'
#' @export
build_oocyst_table <- function(dir = "data_private/fitting",
                               concentration = 200, duration = 6) {
  wb <- assay_workbooks()
  pooled <- file.path(dir, wb$file[grepl("pooledreps", wb$file)])
  extra <- file.path(dir, wb$file[grepl("0h_exposure", wb$file)])

  sheets <- readxl::excel_sheets(pooled)
  sheets <- sheets[grepl(" int$| exposure pre$", sheets)]

  rows <- lapply(sheets, function(s) {
    d <- assay_values(pooled, s)
    if (is.null(d)) {
      return(NULL)
    }
    # The timing is in the header where it varies, and in the worksheet name
    # where it does not.
    sign <- exposure_sign(s)
    d$exposure_h[is.na(d$exposure_h)] <- 0
    d$hours_after_infection <- d$exposure_h * ifelse(is.na(sign), 1, sign)
    d$exposure_h <- NULL
    # The first two sub-groups of the 72h one-blood-meal columns of this
    # worksheet are one experiment (Saxena, 27 Apr 2026 16:29).
    if (s == "#5 6h 72h post int") {
      merge_me <- d$hours_after_infection == 72 & d$n_bloodmeals == 1L &
        d$subgroup <= 2L
      d$subgroup[merge_me] <- 1L
    }
    d$workbook <- basename(pooled)
    d
  })
  out <- do.call(rbind, rows)
  out <- rbind(out, zero_hour_experiments(extra))

  keep <- out$arm == "control" |
    (!is.na(out$concentration) & out$concentration == concentration)
  out <- out[keep, , drop = FALSE]
  # Only the standard rest on the net; the worksheet varying it is otherwise
  # at the same timing, so a missing duration means the standard one.
  out <- out[is.na(out$duration_min) | out$duration_min == duration, ,
             drop = FALSE]

  out <- out[order(out$worksheet, out$hours_after_infection, out$subgroup,
                   out$arm), ]
  label_experiments(out, one_row_each = FALSE)
}

#' The concurrent-exposure experiments held in a separate workbook
#'
#' This workbook carries one column pair per experiment, headed by its date.
#' Its first two experiments are already in the pooled workbook and are dropped
#' (Saxena, 27 Apr 2026 15:23).
#'
#' @param path Path to the workbook.
#'
#' @return A data frame with one row per mosquito.
#'
#' @keywords internal
zero_hour_experiments <- function(path) {
  raw <- suppressMessages(readxl::read_excel(path, sheet = "Data By Experiment",
                                             col_names = FALSE,
                                             col_types = "text",
                                             .name_repair = "minimal"))
  dates <- as.character(unlist(raw[1L, ], use.names = FALSE))
  groups <- as.character(unlist(raw[2L, ], use.names = FALSE))
  body <- raw[-(1:2), , drop = FALSE]

  # Each date heads a pair of columns, control then treated.
  starts <- which(!is.na(dates) & toupper(dates) != "DATE")
  starts <- starts[-(1:2)]   # the first two experiments are already pooled

  out <- lapply(starts, function(j) {
    pair <- lapply(c(j, j + 1L), function(k) {
      arm <- if (grepl("^CTL", groups[[k]], ignore.case = TRUE)) "control" else "treated"
      v <- suppressWarnings(as.numeric(sub("[^0-9.]+$", "",
                                           as.character(body[[k]]))))
      v <- v[!is.na(v)]
      if (!length(v)) {
        return(NULL)
      }
      data.frame(
        worksheet = "Data By Experiment", column = paste(dates[[j]], groups[[k]]),
        arm = arm, concentration = if (arm == "treated") 200 else NA_real_,
        hours_after_infection = 0, duration_min = NA_real_, n_bloodmeals = 1L,
        subgroup = match(j, starts), oocysts = v, stringsAsFactors = FALSE
      )
    })
    do.call(rbind, pair)
  })
  out <- do.call(rbind, out)
  out$workbook <- basename(path)
  out
}
