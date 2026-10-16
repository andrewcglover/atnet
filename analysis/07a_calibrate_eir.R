# Recalibrate each region's baseline EIR for this model, with mrc-ide's cali package.
#
# The site file's EIRs were fitted for upstream malariasimulation. The fork differs in
# its EIP (Erlang stages rather than a fixed delay), and the projections also use a
# 1.8-year net retention (the site file's is about 2014 days) and the compartment counts
# set in 07. Run at the site file's EIR, the history would therefore give a different
# prevalence from the one it was fitted to; this finds the EIR that restores it.
#
# What is fitted: the annual mean PfPR 2-10 of the site file, 2010 to branch_year - 1
# (2021), the part of the history every realisation shares. 2022-2024 are not used: their
# nets depend on each realisation's ITN draw, and the site file's PfPR for those years is a
# single value carried forward. cali::calibrate() is used as it stands. It solves for one
# EIR per region by running the history and comparing output with target (it minimises
# sum(output - target)); the fitting window is inside summary_fn below, because cali's own
# summary function averages the whole run and reads a column name (730_3650) that
# site_parameters() does not produce (730_3649).
#
# The model is built by 07's own build_params(region, "history", NULL), so calibration and
# projections cannot differ: same fork, site package, retention, IRS rule, compartment
# counts and solver settings, with ATN off. 07 is sourced with country_test_mode on, which
# stops it after setup.
#
# One region per process, so regions run side by side:
#   ATN_FORK_PATH=<fork> Rscript analysis/07a_calibrate_eir.R Mopti
# Optional environment variables:
#   ATNET_SITE_LIB  library holding site 1.0.10 (default: the one 08's pre-flight uses).
#   CALI_POP        human populations to try, in order (default 10000,100000). cali moves
#                   to the next when a run eliminates malaria, then raises the EIR.
#   CALI_VALIDATE   1 (default) re-runs the history at the calibrated EIR with HUMAN_POP
#                   people and records how closely it matches the target, year by year.
#   CALI_DRY        1 runs one evaluation at the site file's EIR and the first population,
#                   prints it against the target and stops: a timing and a check of the
#                   summary function, with no calibration.
#   CALI_OUT        output folder (default outputs/calibration/<date>_dq<>_sl<>).
# Writes <CALI_OUT>/<region>.rds. Run with the working directory at the package root.
# Afterwards, 07b_collect_calibrated_eir.R <folder> writes the tracked table that 07 reads.
#
# ORDER MATTERS: cali imports malariasimulation, so loading cali first would load whatever
# is installed under that name. 07 loads the fork first (ATN_FORK_PATH or an installed
# fork); cali is only touched afterwards, and the fork is checked again below.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Give exactly one region, e.g. Rscript analysis/07a_calibrate_eir.R Mopti")

# The library holding site 1.0.10: the same setting, and the same default, as 08's pre-flight.
site_lib <- Sys.getenv("ATNET_SITE_LIB", "C:/Users/ag4218/AppData/Local/R/atnet-lib/4.4")
if (nzchar(site_lib)) .libPaths(c(normalizePath(site_lib), .libPaths()))
if (utils::packageVersion("site") != "1.0.10") {
  stop("site ", utils::packageVersion("site"), " is first on the library path; the Mali site file ",
       "needs 1.0.10. Set ATNET_SITE_LIB to a library holding it.", call. = FALSE)
}

# CALIBRATED_EIR_FILE = none: 07 would otherwise replace the site file's EIRs with a table
# this script is about to replace, and stop if that table was fitted for another model.
Sys.setenv(SWEEP_REGIONS = args[[1]], CALIBRATED_EIR_FILE = "none")
options(country_test_mode = TRUE)
source("analysis/07_run_projections.R")
stopifnot(length(sweep_regions) == 1L)
region <- sweep_regions

if (!"omega_atn" %in% names(malariasimulation::get_parameters())) {
  stop("malariasimulation is not the ATN fork.", call. = FALSE)
}
if (!requireNamespace("cali", quietly = TRUE)) stop("cali is not installed.", call. = FALSE)
# cali's check_fit() calls knitr::kable() although knitr is only a Suggests of cali.
if (!requireNamespace("knitr", quietly = TRUE)) stop("knitr is needed by cali's check_fit().", call. = FALSE)
if (!isNamespaceLoaded("malariasimulation") ||
    !"omega_atn" %in% names(getNamespace("malariasimulation")$get_parameters())) {
  stop("cali would run the wrong malariasimulation.", call. = FALSE)
}

# ---- target and window ------------------------------------------------------------------
fit_first_year <- 2010L
fit_last_year  <- as.integer(branch_year) - 1L
fit_years      <- fit_first_year:fit_last_year

site_row <- site_obj$sites[site_obj$sites$name_1 == region, , drop = FALSE]
ms       <- site::subset_site(site_obj, site_row)
pf       <- ms$prevalence
target   <- setNames(pf$pfpr[match(fit_years, pf$year)], fit_years)
stopifnot(!anyNA(target), length(target) == length(fit_years))

# Annual mean of the daily PfPR 2-10 over the fitting window. Timestep t of a run is day
# t of 365-day years counted from start_year, as in 07.
first_step <- (fit_first_year - start_year) * 365L
last_step  <- (fit_last_year + 1L - start_year) * 365L
summary_fn <- function(x) {
  d <- x[x$timestep > first_step & x$timestep <= last_step, , drop = FALSE]
  yr <- start_year + (d$timestep - 1L) %/% 365L
  out <- tapply(d$n_detect_lm_730_3649 / d$n_age_730_3649, yr, mean)
  stopifnot(length(out) == length(fit_years))
  out
}

# ---- parameters: the history, ATN off, ending where 07's history ends ---------------------
p <- build_params(region, "history", NULL)
p$timesteps <- history_end_day

pops <- as.numeric(strsplit(Sys.getenv("CALI_POP", "10000,100000"), ",")[[1]])
tx_2000 <- ms$interventions$tx_cov[ms$interventions$year == start_year]
site_eir <- ms$eir$eir

out_dir <- Sys.getenv("CALI_OUT", file.path("outputs", "calibration",
                      sprintf("%s_dq%d_sl%d", format(Sys.Date(), "%Y%m%d"), deltaq_use, spor_len_use)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("\nCalibrating %s: %d-%d, mean PfPR2-10 %.4f, site-file EIR %.2f, populations %s",
                region, fit_first_year, fit_last_year, mean(target), site_eir, toString(pops)))

# CALI_DRY=1: one evaluation at the site file's own EIR and the first population, printed
# against the target, then stop. It times one run and shows the summary function working,
# which is all a calibration is made of.
if (identical(Sys.getenv("CALI_DRY", "0"), "1")) {
  pd <- p
  pd$human_population <- pops[1]
  pd <- malariasimulation::set_equilibrium(pd, init_EIR = site_eir)
  set.seed(SEED_BASE)
  td <- proc.time()[["elapsed"]]
  sim <- malariasimulation::run_simulation(timesteps = pd$timesteps, parameters = pd)
  mins <- (proc.time()[["elapsed"]] - td) / 60
  fit <- summary_fn(sim)
  message(sprintf("\nDRY RUN %s: one evaluation, %d people, EIR %.2f: %.1f min", region, pops[1],
                  site_eir, mins))
  print(data.frame(year = fit_years, simulated = round(unname(fit), 4),
                   target = round(unname(target), 4)), row.names = FALSE)
  message(sprintf("mean simulated %.4f vs target %.4f", mean(fit), mean(target)))
  quit(save = "no", status = 0)
}

set.seed(SEED_BASE)
t0 <- proc.time()[["elapsed"]]
eir <- cali::calibrate(
  parameters       = p,
  target           = target,
  summary_function = summary_fn,
  eq_prevalence    = mean(target),
  eq_ft            = tx_2000,
  human_population = pops
)
minutes <- (proc.time()[["elapsed"]] - t0) / 60
message(sprintf("%s: calibrated EIR %.3f (site file %.3f) in %.1f min", region, eir, site_eir, minutes))

# ---- validation: how close is the history at this EIR, at the projections' population? ----
validation <- NULL
if (!identical(Sys.getenv("CALI_VALIDATE", "1"), "0")) {
  pv <- p
  pv$human_population <- human_pop
  pv <- malariasimulation::set_equilibrium(pv, init_EIR = eir)
  set.seed(SEED_BASE)
  tv <- proc.time()[["elapsed"]]
  sim <- malariasimulation::run_simulation(timesteps = pv$timesteps, parameters = pv)
  fit <- summary_fn(sim)
  validation <- list(human_pop = human_pop, simulated = fit, target = target,
                     mean_diff = mean(fit - target), minutes = (proc.time()[["elapsed"]] - tv) / 60)
  message(sprintf("%s validation at %d people: simulated mean %.4f vs target %.4f (difference %+.4f)",
                  region, human_pop, mean(fit), mean(target), mean(fit - target)))
}

saveRDS(list(region = region, eir = eir, site_file_eir = site_eir, target = target,
             window = c(fit_first_year, fit_last_year), populations = pops,
             minutes = minutes, validation = validation, seed = SEED_BASE,
             deltaq = deltaq_use, spor_len = spor_len_use, retention_days = retention_time,
             fork_version = fork_version, sim_identity = sim_identity,
             cali = as.character(utils::packageVersion("cali")), code_md5 = code_md5),
        file.path(out_dir, sprintf("%s.rds", region_tag(region))))
message("Saved to ", out_dir)
