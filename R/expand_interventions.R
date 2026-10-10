# Extending a site file's intervention history into the projection window.
#
# Moved unchanged from the malariasimulation fork's dev/InterventionExpansion.R
# with the projection pipeline (analysis/07_run_projections.R).

#' Carry a site's last year of interventions forward
#'
#' Appends `expand_year` copies of the last historical year's intervention
#' rows, one per future year, with only the year changed. The projection
#' pipeline then overwrites the future rows it controls (nets, and any
#' intervention it switches off).
#'
#' @param site_data A single-site object from `site::subset_site()`, with an
#'   `interventions` data frame carrying a `year` column.
#' @param expand_year The number of future years to append.
#' @param delay,counterfactual Not implemented; anything other than the
#'   defaults gives a warning and is ignored.
#'
#' @return `site_data` with the extended `interventions` data frame.
#'
#' @export
expand_interventions <- function(site_data, expand_year, delay = 0, counterfactual = FALSE) {
  if (delay != 0 || isTRUE(counterfactual))
    warning("expand_interventions: delay != 0 / counterfactual != FALSE are not implemented; ignoring.")
  idf      <- site_data$interventions
  last_yr  <- max(idf$year)
  last_rows <- idf[idf$year == last_yr, , drop = FALSE]
  future <- do.call(rbind, lapply(seq_len(expand_year), function(dy) {
    r <- last_rows; r$year <- last_yr + dy; r
  }))
  rownames(future) <- NULL
  site_data$interventions <- rbind(idf, future)
  site_data
}
