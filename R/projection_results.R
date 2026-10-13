# Reading and summarising the projections written by
# analysis/07_run_projections.R, which saves one file per region and
# draws-table row in <run folder>/jobs/. Each file holds the 2022-24 segment
# (`pre_projection`) and the arms from 31 Dec 2024 (`results`), as daily
# values.
#
# Summaries follow one order throughout. Firstly, values are worked out for
# each region, row and arm. Secondly, arms are compared within a row, since
# all arms of a row start from the same end-2024 state and random-number
# state. Finally, medians and intervals are taken across rows. National
# values add up population-scaled counts before any rate is formed, so they
# are never averages of regional rates or percentages.

#' The arms of the paper, in plotting order
#'
#' @return A data frame with `arm` (as in the results) and `label`.
#'
#' @export
projection_arms <- function() {
  data.frame(
    arm   = c("none", "pyr", "pyr_pbo", "pyr_cfp", "atn", "pyr_atn",
              "pyr_cfp_atn"),
    label = c("No future nets", "Pyr", "Pyr-PBO", "Pyr-CFP", "ATN", "Pyr-ATN",
              "Pyr-CFP-ATN")
  )
}

# The region tag used in file names, as in analysis/07_run_projections.R
# ("Ségou" becomes "Segou").
region_tag <- function(region) {
  gsub("[^A-Za-z0-9]", "", iconv(region, to = "ASCII//TRANSLIT"))
}

#' List the job files of a projection run
#'
#' @param dir The run folder, e.g. `outputs/projections/mli_hldrawn_dummy`.
#' @param regions Optional regions to keep, with or without accents.
#' @param rows Optional draws-table rows to keep.
#'
#' @return A data frame with one line per file (`region_tag`, `draw`,
#'   `file`), ordered by region and row. Requested rows with no file are
#'   reported in a warning and left out. Files whose names start with
#'   "FAKE_" (made-up data for developing plots) are read, with a warning.
#'
#' @export
projection_job_files <- function(dir, regions = NULL, rows = NULL) {
  files <- list.files(file.path(dir, "jobs"), pattern = "_row[0-9]+\\.rds$",
                      full.names = TRUE)
  if (!length(files)) {
    stop("No job files in ", file.path(dir, "jobs"), call. = FALSE)
  }
  parts <- regmatches(basename(files),
                      regexec("^(FAKE_)?[a-z]+_(.+)_row([0-9]+)\\.rds$",
                              basename(files)))
  jobs <- data.frame(
    region_tag = vapply(parts, `[`, "", 3),
    draw       = as.integer(vapply(parts, `[`, "", 4)),
    file       = files
  )
  if (any(startsWith(basename(files), "FAKE_"))) {
    warning("Reading FAKE data from ", dir, ": for plot layout only.",
            call. = FALSE)
  }
  if (!is.null(regions)) jobs <- jobs[jobs$region_tag %in% region_tag(regions), ]
  if (!is.null(rows)) {
    tags <- if (is.null(regions)) unique(jobs$region_tag) else region_tag(regions)
    wanted <- expand.grid(region_tag = tags, draw = rows,
                          stringsAsFactors = FALSE)
    absent <- wanted[!paste(wanted$region_tag, wanted$draw) %in%
                       paste(jobs$region_tag, jobs$draw), ]
    if (nrow(absent)) {
      warning(nrow(absent), " requested job file(s) not found, e.g. ",
              absent$region_tag[1], " row ", absent$draw[1], call. = FALSE)
    }
    jobs <- jobs[jobs$draw %in% rows, ]
  }
  if (!nrow(jobs)) {
    stop("No job files match the regions and rows asked for.", call. = FALSE)
  }
  jobs <- jobs[order(jobs$region_tag, jobs$draw), ]
  rownames(jobs) <- NULL
  jobs
}

#' Read daily projection results
#'
#' Reads the job files of a run into one data frame of daily values.
#'
#' @inheritParams projection_job_files
#' @param columns Optional columns to keep besides the identifiers
#'   (`timestep`, `region`, `arm`, `draw`, `seed`, `year_rel`). Keeping only
#'   what a plot needs saves memory: all 33 columns for 9 regions and 100
#'   rows come to about 4.5 GB.
#' @param pre If `TRUE`, the 2022-24 segment is included, with `arm`
#'   "2022-2024".
#'
#' @return A data frame of daily values, one line per region, row, arm and
#'   day.
#'
#' @export
read_projection_jobs <- function(dir, regions = NULL, rows = NULL,
                                 columns = NULL, pre = TRUE) {
  ids <- c("timestep", "region", "arm", "draw", "seed", "year_rel")
  jobs <- projection_job_files(dir, regions, rows)
  out <- lapply(jobs$file, function(f) {
    job <- readRDS(f)
    x <- if (pre) rbind(job$pre_projection, job$results) else job$results
    if (!is.null(columns)) {
      absent <- setdiff(columns, names(x))
      if (length(absent)) {
        stop("Column(s) not in ", basename(f), ": ",
             paste(absent, collapse = ", "), call. = FALSE)
      }
      x <- x[, union(intersect(ids, names(x)), columns), drop = FALSE]
    }
    x
  })
  do.call(rbind, out)
}

#' The days of one campaign cycle
#'
#' From one mass campaign up to the day before the next. The default, the
#' cycle starting at the second campaign (2028 in the Mali runs), is the
#' period quoted in the paper.
#'
#' @param campaign_days The future campaign days, `future_campaign_days` in
#'   the run's `settings.rds`.
#' @param cycle Which cycle; 1 starts at the first campaign.
#'
#' @return A named integer vector: `from` and `to`, timesteps, both included.
#'
#' @export
campaign_cycle <- function(campaign_days, cycle = 2L) {
  campaign_days <- sort(campaign_days)
  if (cycle < 1L || cycle >= length(campaign_days)) {
    stop("Cycle ", cycle, " needs campaigns at both ends; there are ",
         length(campaign_days), " campaigns.", call. = FALSE)
  }
  c(from = as.integer(campaign_days[cycle]),
    to   = as.integer(campaign_days[cycle + 1L] - 1L))
}

# The count columns of period_counts(): totals over the period, which can be
# scaled by population and added across regions.
count_columns <- function() {
  c("person_years", "clin_cases", "severe_cases", "detect_2_10", "age_2_10",
    "net_user_days", "infectious_bites", "mosquito_days",
    "exposed_mosquito_days")
}

#' Totals over a period for each region, row and arm
#'
#' @param x Daily results, as from [read_projection_jobs()] with all columns.
#' @param from,to The first and last timestep of the period, both included.
#'
#' @return One line per region, row (`draw`) and arm, with the number of
#'   `days` and the totals over the period: `person_years`; `clin_cases` and
#'   `severe_cases` (all ages); `detect_2_10` and `age_2_10` (person-days
#'   aged 2-10, detected by microscopy and in all); `net_user_days`;
#'   `infectious_bites`; `mosquito_days` and `exposed_mosquito_days` (adult
#'   females, all and in the ATN-exposed classes). The rates of
#'   [projection_rates()] are added.
#'
#' @export
period_counts <- function(x, from, to) {
  x <- x[x$timestep >= from & x$timestep <= to, , drop = FALSE]
  if (!nrow(x)) stop("No days between ", from, " and ", to, call. = FALSE)
  sum_cols <- function(cols) rowSums(x[, cols, drop = FALSE])
  eir_cols <- setdiff(grep("^EIR_", names(x), value = TRUE), "EIR_total_pp")
  daily <- data.frame(
    region = x$region, draw = x$draw, arm = x$arm,
    days = 1,
    # site_parameters() forces these three age groups, which cover everyone.
    person_years = sum_cols(c("n_age_0_1824", "n_age_1825_5474",
                              "n_age_5475_36499")) / 365,
    clin_cases = sum_cols(grep("^n_inc_clinical_", names(x), value = TRUE)),
    severe_cases = sum_cols(grep("^n_inc_severe_", names(x), value = TRUE)),
    detect_2_10 = x$n_detect_lm_730_3649,
    age_2_10 = x$n_age_730_3649,
    net_user_days = x$n_use_net,
    infectious_bites = sum_cols(eir_cols),
    mosquito_days = sum_cols(grep("^total_M_", names(x), value = TRUE)),
    exposed_mosquito_days = sum_cols(grep("^M_exposed_", names(x), value = TRUE))
  )
  out <- stats::aggregate(daily[c("days", count_columns())],
                          by = daily[c("region", "draw", "arm")], FUN = sum)
  projection_rates(out)
}

#' Rates from period totals
#'
#' @param v Totals over a period, as from [period_counts()].
#'
#' @return `v` with `pfpr2to10_pct` (microscopy, ages 2-10, %),
#'   `clin_inc_per1000` and `severe_inc_per1000` (cases per 1,000 people per
#'   year), `net_use_pct` (% of people using a net), `annual_eir` (infectious
#'   bites per person per year) and `atn_exposed_pct` (% of adult female
#'   mosquitoes in the ATN-exposed classes).
#'
#' @export
projection_rates <- function(v) {
  v$pfpr2to10_pct <- 100 * v$detect_2_10 / v$age_2_10
  v$clin_inc_per1000 <- 1000 * v$clin_cases / v$person_years
  v$severe_inc_per1000 <- 1000 * v$severe_cases / v$person_years
  v$net_use_pct <- 100 * v$net_user_days / (365 * v$person_years)
  v$annual_eir <- v$infectious_bites / v$person_years
  v$atn_exposed_pct <- 100 * v$exposed_mosquito_days / v$mosquito_days
  v
}

#' Per-row totals for a whole run
#'
#' Applies [period_counts()] to each job file in turn, so that only one file
#' is held in memory at a time.
#'
#' @inheritParams projection_job_files
#' @inheritParams period_counts
#'
#' @return As [period_counts()], for every file read.
#'
#' @export
projection_row_values <- function(dir, from, to, regions = NULL, rows = NULL) {
  jobs <- projection_job_files(dir, regions, rows)
  out <- lapply(jobs$file, function(f) {
    job <- readRDS(f)
    period_counts(rbind(job$pre_projection, job$results), from, to)
  })
  do.call(rbind, out)
}

#' Regional populations from a site file
#'
#' @param site A site file, e.g. `readRDS("data_private/site_files/MLI.rds")`.
#' @param years The years to average over; by default 2029 and 2030, the
#'   middle of the 2028-31 cycle.
#'
#' @return A data frame with `region` and `pop`.
#'
#' @export
region_populations <- function(site, years = 2029:2030) {
  p <- site$population$population_total
  p <- p[p$year %in% years, ]
  absent <- setdiff(years, p$year)
  if (length(absent)) {
    stop("No population for year(s) ", paste(absent, collapse = ", "),
         call. = FALSE)
  }
  by_year <- stats::aggregate(pop ~ name_1 + year, data = p, FUN = sum)
  out <- stats::aggregate(pop ~ name_1, data = by_year, FUN = mean)
  names(out) <- c("region", "pop")
  out
}

#' National values for each row and arm
#'
#' Each region's totals are scaled by its real population over its simulated
#' population, then added across regions, and the rates are formed from the
#' national totals. Rows lacking any region are reported in a warning and
#' left out, so a missing region cannot shrink a national total unnoticed.
#'
#' @param v Per-row totals, as from [projection_row_values()].
#' @param pop Real populations, as from [region_populations()]. Every region
#'   in `pop` counts towards the national total.
#' @param name What to call the nation in `region`.
#'
#' @return As [period_counts()], with one line per row and arm.
#'
#' @export
national_values <- function(v, pop, name = "National") {
  tags <- region_tag(pop$region)
  v <- v[region_tag(v$region) %in% tags, ]
  key <- paste(v$draw, v$arm)
  complete <- tapply(region_tag(v$region), key,
                     function(r) all(tags %in% r))
  if (!all(complete)) {
    warning(sum(!complete), " row x arm combination(s) lack a region and ",
            "are left out of the national values.", call. = FALSE)
    v <- v[key %in% names(complete)[complete], ]
  }
  if (!nrow(v)) stop("No row has every region.", call. = FALSE)
  simulated_pop <- 365 * v$person_years / v$days
  scale <- pop$pop[match(region_tag(v$region), tags)] / simulated_pop
  v[count_columns()] <- v[count_columns()] * scale
  out <- stats::aggregate(v[count_columns()], by = v[c("draw", "arm")],
                          FUN = sum)
  out$days <- stats::aggregate(v["days"], by = v[c("draw", "arm")],
                               FUN = max)$days
  out <- cbind(region = name, out[c("draw", "arm", "days", count_columns())])
  projection_rates(out)
}

#' Clinical cases averted, paired within rows
#'
#' Compares each arm with the comparator arm of the same region and row.
#' Applied to [national_values()], the percentage averted is therefore the
#' national ratio, not an average of regional percentages.
#'
#' @param v Per-row totals, as from [projection_row_values()] or
#'   [national_values()].
#' @param comparator The arm compared against.
#'
#' @return One line per region, row and arm (other than the comparator):
#'   `cases_averted` (over the period), `cases_averted_per1000` (per 1,000
#'   people per year) and `pct_cases_averted`.
#'
#' @export
arm_contrasts <- function(v, comparator = "none") {
  ref <- v[v$arm == comparator, c("region", "draw", "clin_cases",
                                  "person_years")]
  if (!nrow(ref)) stop("No lines for arm ", comparator, call. = FALSE)
  names(ref)[3:4] <- c("ref_cases", "ref_person_years")
  out <- merge(v[v$arm != comparator, c("region", "draw", "arm", "clin_cases")],
               ref, by = c("region", "draw"))
  out$comparator <- comparator
  out$cases_averted <- out$ref_cases - out$clin_cases
  out$cases_averted_per1000 <- 1000 * out$cases_averted / out$ref_person_years
  out$pct_cases_averted <- 100 * out$cases_averted / out$ref_cases
  out[c("region", "draw", "arm", "comparator", "cases_averted",
        "cases_averted_per1000", "pct_cases_averted")]
}

#' Median and interval across rows
#'
#' @param x A data frame with one line per row (`draw`) in each group.
#' @param cols The columns to summarise.
#' @param by The grouping columns.
#' @param probs The interval's quantiles; by default the 95% interval.
#' @param digits Decimal places in `text`: one number, or a vector named by
#'   column.
#'
#' @return One line per group and column: `quantity`, `n_rows`, `median`,
#'   `lower`, `upper` and `text`, the last written "XX \[XX, XX\]".
#'
#' @export
summarise_rows <- function(x, cols, by = c("region", "arm"),
                           probs = c(0.025, 0.975), digits = 1) {
  groups <- split(x, x[by], drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) {
    if (anyDuplicated(g$draw)) {
      stop("A group has more than one line for the same row.", call. = FALSE)
    }
    do.call(rbind, lapply(cols, function(col) {
      q <- stats::quantile(g[[col]], c(0.5, probs), names = FALSE)
      data.frame(g[1, by, drop = FALSE], quantity = col, n_rows = nrow(g),
                 median = q[1], lower = q[2], upper = q[3])
    }))
  }))
  d <- if (is.null(names(digits))) rep(digits, nrow(out)) else digits[out$quantity]
  out$text <- format_interval(out$median, out$lower, out$upper, d)
  rownames(out) <- NULL
  out
}

#' Write a median and interval as "XX \[XX, XX\]"
#'
#' @param median,lower,upper Numbers.
#' @param digits Decimal places.
#'
#' @return A character vector.
#'
#' @export
format_interval <- function(median, lower, upper, digits = 1) {
  f <- function(v) mapply(formatC, v, digits = digits,
                          MoreArgs = list(format = "f", big.mark = ","))
  sprintf("%s [%s, %s]", f(median), f(lower), f(upper))
}
