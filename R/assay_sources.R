# Where the fitted data come from, and the corrections applied on the way.
#
# Three workbooks feed the two fits. They live in data_private/fitting/, which
# is not in this repository, each prefixed with the time it was received
# (YYYYMMDDhhmm), so the data behind a fit can be traced to the message that
# carried it. Who sent each workbook, and the correspondence behind each
# correction, are recorded in data_private/DATA_SOURCES.md, kept with the data.

#' The workbooks the fits are built from
#'
#' @return A data frame with one row per workbook, giving its role, its name
#'   without the time it was received, its layout, and which fits it feeds.
#'
#' @export
assay_workbooks <- function() {
  data.frame(
    role = c("per_replicate", "pooled", "zero_hour"),
    name = c("200mgPCLnetdata_individualreps.xlsx",
             "200mgPCLnetdata_pooledreps.xlsx",
             "200mgm2_ELQ-453_0h_exposure.xlsx"),
    layout = c("one worksheet per replicate",
               "one worksheet per experimental design, replicates stacked",
               "one column pair per experiment, dated"),
    feeds = c("sporozoite", "sporozoite and oocyst", "oocyst"),
    stringsAsFactors = FALSE
  )
}

#' Find a workbook the fits are built from
#'
#' @param dir Folder holding the workbooks, normally `data_private/fitting`.
#' @param role One of the roles in [assay_workbooks()].
#'
#' @return The path to the one file in `dir` with that workbook's name,
#'   prefixed with the time it was received.
#'
#' @export
find_workbook <- function(dir, role) {
  wb <- assay_workbooks()
  name <- wb$name[wb$role == role]
  if (length(name) != 1L) {
    stop("Unknown workbook role '", role, "'.", call. = FALSE)
  }
  files <- list.files(dir)
  hits <- files[endsWith(files, paste0("_", name)) &
                  grepl("^[0-9]{12}_", files)]
  if (length(hits) != 1L) {
    stop(sprintf("Expected one '<YYYYMMDDhhmm>_%s' in %s, found %d.",
                 name, dir, length(hits)), call. = FALSE)
  }
  file.path(dir, hits)
}

#' Corrections applied to the assay data
#'
#' Every departure from a plain reading of the workbooks, with the evidence for
#' it. Each is applied by [build_sporozoite_table()] or [build_oocyst_table()],
#' and `analysis/02_build_fitting_data.R` reports them all, so that no
#' correction is silent.
#'
#' @return A data frame with one row per correction.
#'
#' @export
assay_corrections <- function() {
  data.frame(
    id = c("merge_10_9_25", "merge_oocyst_s5", "drop_dup_11_6_25",
           "drop_dup_09_07_26", "day_from_worksheet", "treated_200_only",
           "standard_duration_only"),
    applies_to = c("sporozoite", "oocyst", "sporozoite", "sporozoite",
                   "both", "both", "oocyst"),
    description = c(
      "Worksheet '6h 72h 2BF post 10.9.25 spz13' is one experiment, not three. Its staggered blank cells read as three bands when the sheet is read across.",
      "In worksheet '#5 6h 72h post int', the first two groups of the 72h one-blood-meal columns are one experiment.",
      "In worksheet '6h 72h 2BF post 11.6.25 spz13', the one-blood-meal columns hold the same 24 values twice. Count them once.",
      "Worksheet '72h 2BF post 09.07.26 spz10' holds a copy of '6h 72h 2BF post 11.3.25 spz10'. Ignore it and take the genuine replicate from the pooled workbook.",
      "The dissection day comes from the worksheet name, not from a typed column, which carried a transcription slip.",
      "Only nets at 200 mg per square metre are fitted, so the 100 mg column of '#1 0h int' is excluded.",
      "Only the standard six-minute rest on the net is fitted, so the three and one minute arms of '#2 6,3,1 min exposure pre' are excluded. The model has no dimension for how long a mosquito rests."
    ),
    evidence = c(
      "confirmed by the laboratory in writing",
      "confirmed by the laboratory in writing",
      "identical values, verified",
      "identical values, verified; the pooled workbook reconciles only with this reading",
      "colour coding in the pooled workbook links the rows to a day 13 worksheet",
      "decision, A C Glover, 9 Oct 2026",
      "model scope; six minutes is the laboratory's standard exposure"
    ),
    stringsAsFactors = FALSE
  )
}
