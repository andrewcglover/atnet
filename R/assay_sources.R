# Where the fitted data come from, and the corrections applied on the way.
#
# Three workbooks feed the two fits. Each is named for when the email carrying
# it was received, so that the data behind a fit can always be traced to a
# message. They live in data_private/fitting/, which is not in this repository.

#' The workbooks the fits are built from
#'
#' @return A data frame with one row per workbook, giving the file name, when it
#'   was received, who sent it, and which fits it feeds.
#'
#' @export
assay_workbooks <- function() {
  data.frame(
    file = c("202609181722_200mgPCLnetdata_individualreps.xlsx",
             "202609181722_200mgPCLnetdata_pooledreps.xlsx",
             "202604271323_200mgm2_ELQ-453_0h_exposure.xlsx"),
    received = c("2026-09-18 17:22", "2026-09-18 17:22", "2026-04-27 13:23"),
    sender = c("Aditi Saxena", "Aditi Saxena", "forwarded by A C Glover"),
    layout = c("one worksheet per replicate",
               "one worksheet per experimental design, replicates stacked",
               "one column pair per experiment, dated"),
    feeds = c("sporozoite", "sporozoite and oocyst", "oocyst"),
    stringsAsFactors = FALSE
  )
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
      "Saxena, 27 Apr 2026 15:23",
      "Saxena, 27 Apr 2026 16:29",
      "identical values, verified",
      "identical values, verified; the pooled workbook reconciles only with this reading",
      "colour coding in the pooled workbook links the rows to a day 13 worksheet",
      "decision, A C Glover, 9 Oct 2026",
      "model scope; Adams manuscript draft, 'our standard 6-minute exposure time'"
    ),
    stringsAsFactors = FALSE
  )
}
