# Refit the site file's yearly ITN distributions for the projections' net retention.
#
# The site file's itn_input_dist (the share of people given a net on 1 January each year,
# 2000-2024) was fitted with netz to the site file's annual net-use estimates (itn_use),
# read at mid-year, under the site file's own mean retention (about 2014 days). This
# script reproduces that fit first, as a check that the method is the same, and then
# repeats it at the retention 07 uses (1.8 years = 657 days), so that simulated use
# follows the same estimates with the shorter retention. Distributions stay on 1 January;
# a year whose estimate cannot be reached even by giving everyone a net is capped at 1.
#
# Output: data_private/itn_distribution/itn_input_dist_ret<days>.csv, one row per region
# and year, read by 07 in place of the site file's itn_input_dist. It stays in
# data_private/ (not on GitHub): its values are derived directly from the site file.
#
#   Rscript analysis/06b_fit_itn_distribution.R            (working directory = atnet root)
# Optional: NET_RETENTION_DAYS (default 657, as in 07).

SITE_FILE <- file.path("data_private", "site_files", "MLI.rds")
retention_days <- as.numeric(Sys.getenv("NET_RETENTION_DAYS", "657"))
stopifnot(is.finite(retention_days), retention_days > 0)
start_year <- 2000L        # as in 07: day t is day t of 365-day years from 1 January 2000
use_day    <- 183L         # the site file's fit reads each year's use at mid-year

site <- readRDS(SITE_FILE)
iv   <- site$interventions

fit_region <- function(x, retention) {
  x  <- x[order(x$year), ]
  ut <- (x$year - start_year) * 365L + use_day
  dt <- (x$year - start_year) * 365L + x$itn_distribution_day
  d  <- netz::usage_to_model_distribution(
    usage = x$itn_use, usage_timesteps = ut, distribution_timesteps = dt,
    distribution_lower = rep(0, nrow(x)), distribution_upper = rep(1, nrow(x)),
    net_loss_function = netz::net_loss_exp, mean_retention = retention)
  u  <- netz::model_distribution_to_usage(ut, d, dt, net_loss_function = netz::net_loss_exp,
                                          mean_retention = retention)
  data.frame(year = x$year, timestep = dt, target_use = x$itn_use, coverage = d, fitted_use = u,
             site_file_coverage = x$itn_input_dist)
}

out <- do.call(rbind, lapply(unique(iv$name_1), function(r) {
  x <- iv[iv$name_1 == r, ]
  # The check: at the site file's own retention, netz must give back its distributions.
  site_ret <- unique(x$mean_retention)
  stopifnot(length(site_ret) == 1L)
  back <- fit_region(x, site_ret)
  err  <- max(abs(back$coverage - back$site_file_coverage))
  if (err > 1e-3) stop(r, ": netz does not reproduce the site file's distributions (max error ", err, ")")
  cbind(region = iconv(r, to = "ASCII//TRANSLIT"), fit_region(x, retention_days),
        site_file_retention_days = site_ret, reproduction_error = err)
}))
out$capped         <- out$coverage > 1 - 1e-9
out$retention_days <- retention_days
out$netz_version   <- as.character(utils::packageVersion("netz"))
out$site_file_md5  <- unname(tools::md5sum(SITE_FILE))

dir <- file.path("data_private", "itn_distribution")
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
f <- file.path(dir, sprintf("itn_input_dist_ret%d.csv", round(retention_days)))
write.csv(out, f, row.names = FALSE)

message(sprintf("Wrote %s: %d regions x %d years; site file reproduced to within %.1e",
                f, length(unique(out$region)), length(unique(out$year)), max(out$reproduction_error)))
k <- out$year >= 2010 & out$year <= 2021
s <- do.call(rbind, lapply(split(out, out$region), function(d) {
  kk <- d$year >= 2010 & d$year <= 2021
  data.frame(region = d$region[1], target_use_2010_21 = mean(d$target_use[kk]),
             fitted_use_2010_21 = mean(d$fitted_use[kk]), years_capped = sum(d$capped),
             first_capped = if (any(d$capped)) min(d$year[d$capped]) else NA_integer_,
             max_shortfall = max(d$target_use - d$fitted_use))
}))
print(s, row.names = FALSE, digits = 3)
