# =====================================================================
# 07_run_projections.R
#   Country-level admin-1 (no urban/rural split): 7-year forward projections
#   (2025-2031 for Mali) under SEVEN future net arms, one simulation per row of
#   the draws table (inst/extdata/projection_draws.csv, analysis/06_make_draws.R)
#   for each region and arm. Mass campaigns every three years (2025, 2028,
#   2031), each at the start of the calendar month two months before the month
#   of peak rainfall (1 June for Mali, whose regions peak 8-13 August), with
#   monthly continuous distribution between:
#
#     none        — no future nets
#     pyr         — pyrethroid-only (Pyr)
#     pyr_pbo     — pyrethroid + PBO (Pyr-PBO)
#     pyr_cfp     — pyrethroid + chlorfenapyr (Pyr-CFP)
#     atn         — antimalarial net, no insecticide (ATN)
#     pyr_atn     — pyrethroid + antimalarial (Pyr-ATN)
#     pyr_cfp_atn — Pyr-CFP + antimalarial (Pyr-CFP-ATN)
#
#   An eighth arm, pyr_cfp_mc_atn_cd (Pyr-CFP campaigns, ATN continuous
#   distribution), is not in the paper but can still be run via SWEEP_ARMS.
#
#   Each simulation runs in three segments, resuming from the saved state:
#     1. 2000-2021, ONE run per region shared by every realisation (seed
#        SEED_BASE), with the site file's net efficacy. Saved under HISTORY_DIR
#        and reused by later batches, so all realisations share the same history.
#     2. 2022-2024, one run per region and realisation i (draws-table row i,
#        seed SEED_BASE + i): nets distributed in these years take the efficacy
#        of ITN draw i for the site's net type, at that year's resistance. Nets
#        handed out earlier keep the site values while still in use.
#     3. 2025-2031, every arm resumes from realisation i's end-2024 state
#        (30 December, the day before the first future distribution), with its
#        random-number state restored, so the arms of a realisation are
#        compared on the same random numbers.
#   Every segment and arm uses deltaq = deltaq_use, so the saved states fit.
#   Simulation i of every region and arm uses row i of the draws table:
#     - the EIP and blocking kernel parameters, omega_atn and the antimalarial
#       half-life (ATN arms only);
#     - p_A: p_a_atn for the ATN arm, p_a_aitn for Pyr-ATN and Pyr-CFP-ATN;
#     - itn_draw: the draw of the Churcher et al. (2024) net parameters
#       (read_itn_draws()) used for the nets distributed from 2022, at each
#       year's pyrethroid resistance in the region.
#   Arms without an antimalarial use only itn_draw and the seed.
#
#   Settings (environment variables):
#     DRAW_ROWS         rows of the draws table to run, e.g. 1-100 (required)
#     ANTIMAL_HL_YEARS  fix the antimalarial half-life (years) in every
#                       realisation, e.g. 1, 2.5 or 5; unset = each realisation's
#                       drawn value
#     SWEEP_ARMS        comma-separated subset of arms (e.g. the three ATN
#                       arms for the fixed half-life runs)
#     SWEEP_REGIONS     comma-separated subset of regions (default all)
#     STOP_BY           e.g. "2026-10-12 04:00" (local time): no job is started
#                       that would not finish by then, judged by the slowest job
#                       already saved; the sweep then ends and saves its summary
#     STAGE             unset = everything in one go (a local run); on the
#                       cluster (analysis/08_submit_cluster.R) "histories" writes
#                       settings.rds and builds the histories, then
#                       "realisations" runs the jobs, which need both to exist
#     HUMAN_POP         simulated population (default 100000); smaller for tests
#     END_DAY           stop the arms on this day instead of the last one (tests)
#     HISTORY_DIR       folder of the saved histories (tests use their own)
#     SWEEP_CORES, NET_RETENTION_DAYS, OUT_SUFFIX, ATN_FORK_PATH, ATNET_CODE_DIR,
#     ATNET_VERSION: see below
#   A "realisation" is simulation i: row i of the draws table, with seed
#   SEED_BASE + i. Results go to the folder outputs/projections/<iso>_hl<tag>
#   (hl tag "drawn", or the fixed half-life), one file per region and
#   realisation in jobs/ holding its 2022-2024 segment and every arm, with the
#   columns in keep_columns and each segment's run time. A job whose file
#   exists is skipped, so a stopped sweep is restarted by launching it again;
#   nothing is overwritten. A job in which any segment failed is saved in
#   jobs_failed/ instead, for diagnosis, and is run again by the next launch.
#   Each sweep also writes sweep_rows<a>-<b>.rds, or sweep_<regions>_rows<a>-<b>.rds
#   for a subset of regions (settings, run times, failures); a sweep that left jobs
#   failed or not started is finished by launching it again. The folder's
#   settings.rds stops a later batch of realisations from mixing settings or
#   versions (the fork, the model's packages, R, the site file, this code).
#
#   Change COUNTRY_ISO to run any available country. Default: MLI (Mali).
#
#   Moved from the malariasimulation fork's dev/c24med_projection_run.R
#   (atn-dev at 967d2bd), which used median parameters throughout. Inputs:
#     site file         data_private/site_files/<ISO>.rds
#                       (copied from the fork's dev/site_files/without_split/)
#     draws table       inst/extdata/projection_draws.csv
#     ITN draws         inst/extdata/churcher2024_itn_draws/
#     results           outputs/projections/ (histories in history/)
#   and expand_interventions() is a function of this package.
#
#   Run with the working directory at the root of this package. The
#   malariasimulation fork is loaded from source when ATN_FORK_PATH is set
#   (a local run), and otherwise from the installed library (on the cluster,
#   where the fork is installed at a fixed tag):
#       ATN_FORK_PATH=C:/Users/ag4218/Local/github/malariasimulation \
#         DRAW_ROWS=1-100 Rscript analysis/07_run_projections.R
#   Either way the run stops unless the loaded malariasimulation is the fork.
#   On the cluster the code (this script, R/ and inst/) is run from a frozen
#   copy, ATNET_CODE_DIR, made at submission, so editing the package cannot
#   change a batch already queued; data_private/ and outputs/ stay at the root.
# =====================================================================

# The package root (data_private/, outputs/) and the code to run, which is the
# root itself unless ATNET_CODE_DIR names a frozen copy.
ATNET_PATH <- normalizePath(".")
CODE_PATH  <- normalizePath(Sys.getenv("ATNET_CODE_DIR", ATNET_PATH), mustWork = TRUE)
pkgload::load_all(CODE_PATH, quiet = TRUE)   # atnet: expand_interventions(), read_itn_draws()
FORK_PATH  <- Sys.getenv("ATN_FORK_PATH", "")
if (nzchar(FORK_PATH)) {
  FORK_PATH <- normalizePath(FORK_PATH)
  pkgload::load_all(FORK_PATH, quiet = TRUE)
} else {
  library(malariasimulation)
}
if (!"omega_atn" %in% names(malariasimulation::get_parameters())) {
  stop("The loaded malariasimulation is not the ATN fork. Set ATN_FORK_PATH to the fork's ",
       "folder, or install the fork.", call. = FALSE)
}
git_out <- function(dir, args) {
  tryCatch(system(sprintf('git -C "%s" %s', dir, args), intern = TRUE),
           error = function(e) NA_character_, warning = function(w) NA_character_)
}
# Which fork produced the results. From a source folder: its commit, with
# "-dirty" if its code (R/, src/, DESCRIPTION, NAMESPACE, inst/) has uncommitted
# changes; installed: the package's version and GitHub commit.
if (nzchar(FORK_PATH)) {
  fork_dirty   <- length(git_out(FORK_PATH, paste("status --porcelain --untracked-files=no --",
                                                  "R src DESCRIPTION NAMESPACE inst"))) > 0
  fork_sha     <- paste0(git_out(FORK_PATH, "rev-parse HEAD"), if (fork_dirty) "-dirty")
  fork_version <- paste0(git_out(FORK_PATH, "describe --tags --always"), if (fork_dirty) "-dirty")
} else {
  d <- utils::packageDescription("malariasimulation")
  fork_sha     <- if (is.null(d$RemoteSha)) NA_character_ else d$RemoteSha
  fork_version <- paste(c(d$Version, d$RemoteRef, d$RemoteSha), collapse = " ")
}
# And which version of this package ("-dirty" = uncommitted changes). On the
# cluster, where git may be missing, the submission passes it in ATNET_VERSION.
atnet_version <- Sys.getenv("ATNET_VERSION", "")
if (!nzchar(atnet_version)) atnet_version <- git_out(ATNET_PATH, "describe --always --dirty")
library(site)
library(dplyr); library(tidyr)
library(parallel)

# What the simulations themselves depend on besides their inputs: the fork, the
# packages that run the model, and R. A saved history is reused only if it was
# built with the same; settings.rds adds the site file and this pipeline's code
# (code_md5, below), so that one results folder never mixes versions.
pkg_id <- function(p) {
  d <- utils::packageDescription(p)
  paste(c(d$Version, d$RemoteSha), collapse = " ")
}
sim_identity <- list(
  fork_sha  = fork_sha,
  packages  = vapply(c("malariasimulation", "individual", "dqrng", "malariaEquilibrium", "site"),
                     pkg_id, character(1)),
  r_version = R.version$version.string
)
# hipercow's id for the installation of the cluster's library (none locally).
conan_id <- local({
  d <- list.files(file.path(.libPaths(), ".conan"))
  if (length(d)) max(d) else NA_character_
})

# ---------------------------------------------------------------------
# 0. Paths & top-level settings
# ---------------------------------------------------------------------
COUNTRY_ISO <- "MLI"   # ISO3 of site file under data_private/site_files/
SITE_FILE   <- file.path("data_private", "site_files", sprintf("%s.rds", COUNTRY_ISO))

# Settings that no longer exist. Each is refused if set, so that an old command
# cannot appear to choose something the run then ignores.
removed_settings <- c(
  ATN_OMEGA     = "omega_atn is now drawn per realisation, from the draws table",
  ATN_CHEM_DOSE = "chem_dose_atn was replaced by omega_atn on 2026-10-09",
  FIT_STAMP     = "the kernel parameters now come from the draws table, whose fit is pinned in analysis/06_make_draws.R"
)
for (v in names(removed_settings)) {
  if (nzchar(Sys.getenv(v, ""))) {
    stop(v, " is no longer used: ", removed_settings[[v]], ".", call. = FALSE)
  }
}

# The draws table: row i holds every uncertain input of simulation i (realisation i).
DRAWS_FILE <- file.path(CODE_PATH, "inst", "extdata", "projection_draws.csv")

# STAGE: unset runs the histories and then the realisations in one go. On the
# cluster the two run as separate submissions: "histories" (one job, writing
# settings.rds and every region's history) and then "realisations" (many jobs,
# which require settings.rds and the histories to exist and never write them,
# so jobs running side by side cannot race to write the same file).
STAGE <- Sys.getenv("STAGE", "")
if (!STAGE %in% c("", "histories", "realisations")) {
  stop("STAGE must be unset, 'histories' or 'realisations', not '", STAGE, "'.", call. = FALSE)
}

# DRAW_ROWS: the realisations to run, as rows of the draws table, first-last
# (e.g. 1-100) or a single one. Required, so a large run cannot start by
# accident; a setup-only check (country_test_mode) and the histories stage, which
# runs no realisation, default to realisation 1.
draw_rows_env <- Sys.getenv("DRAW_ROWS", "")
if (!nzchar(draw_rows_env)) {
  if (!isTRUE(getOption("country_test_mode")) && STAGE != "histories") {
    stop("Set DRAW_ROWS to the rows of the draws table to run, e.g. DRAW_ROWS=1-100.",
         call. = FALSE)
  }
  draw_rows_env <- "1"
}
draw_bounds <- suppressWarnings(as.integer(strsplit(draw_rows_env, "-", fixed = TRUE)[[1]]))
if (length(draw_bounds) == 1L) draw_bounds <- rep(draw_bounds, 2L)
if (length(draw_bounds) != 2L || anyNA(draw_bounds) || draw_bounds[1] < 1L ||
    draw_bounds[2] < draw_bounds[1]) {
  stop("DRAW_ROWS must be first-last (e.g. 1-100) or a single row, not '", draw_rows_env, "'.",
       call. = FALSE)
}
draw_rows <- seq.int(draw_bounds[1], draw_bounds[2])

# Random numbers: the shared 2000-2021 history of each region uses seed SEED_BASE;
# from 2022 simulation i uses seed SEED_BASE + i, the same for every region, and
# its arms all continue from the same random-number state.
SEED_BASE <- 20261011L

# First year in which nets take the realisation's ITN draw, and so the first year
# that differs between realisations; the history segment runs up to the end of the
# year before.
branch_year <- 2022L
# Saved 2000-2021 histories, one file per region, reused by every later batch of
# realisations (and by the fixed half-life runs). A history is reused only if its
# parameters are identical to the ones this run would build, and it was built with
# the same fork, model packages and R (sim_identity). HISTORY_DIR
# overrides the folder, so that tests never share the final runs' histories.
HISTORY_DIR <- Sys.getenv("HISTORY_DIR", file.path(ATNET_PATH, "outputs", "projections", "history"))

# ANTIMAL_HL_YEARS: fix the antimalarial half-life (years) in every realisation,
# for the sensitivity runs at 1, 2.5 and 5 years; everything else still comes from
# the realisation's row of the draws table. Unset = each realisation's drawn
# half-life (half_life_atn_years).
#   Launch as: ANTIMAL_HL_YEARS=5 SWEEP_ARMS=atn,pyr_atn,pyr_cfp_atn DRAW_ROWS=1-100 \
#                Rscript analysis/07_run_projections.R
hl_env <- Sys.getenv("ANTIMAL_HL_YEARS", "")
hl_fixed_years <- if (nzchar(hl_env)) as.numeric(hl_env) else NULL
if (!is.null(hl_fixed_years) && !(is.finite(hl_fixed_years) && hl_fixed_years > 0)) {
  stop("ANTIMAL_HL_YEARS must be a positive number of years, not '", hl_env, "'.",
       call. = FALSE)
}
hl_tag <- if (is.null(hl_fixed_years)) "drawn" else
  gsub("\\.", "p", format(hl_fixed_years, trim = TRUE))  # 2.5 -> "2p5", 5 -> "5"

# OUT_SUFFIX: extra tag appended to the results folder name, e.g. for retention
# sweeps or add-on arm runs. Nothing in it is ever overwritten.
out_extra_suffix <- Sys.getenv("OUT_SUFFIX", "")
OUT_DIR    <- file.path(ATNET_PATH, "outputs", "projections",
                        sprintf("%s_hl%s%s", tolower(COUNTRY_ISO), hl_tag, out_extra_suffix))
JOB_DIR    <- file.path(OUT_DIR, "jobs")
# Jobs in which any segment failed, kept for diagnosis; they are not skipped.
FAILED_DIR <- file.path(OUT_DIR, "jobs_failed")

# SWEEP_CORES: number of parallel workers. Override via env var for concurrent
# runs. Required on the cluster, where it must equal the cores the job was given
# (the default below would start up to 16 workers whatever the job's share).
sweep_cores_env <- Sys.getenv("SWEEP_CORES", "")
if (nzchar(STAGE) && !nzchar(sweep_cores_env)) {
  stop("Set SWEEP_CORES to the cores the job was given when STAGE is set.", call. = FALSE)
}
N_CORES <- if (nchar(sweep_cores_env) > 0L) {
  as.integer(sweep_cores_env)
} else {
  min(16L, max(1L, parallel::detectCores() - 1L))
}
if (is.na(N_CORES) || N_CORES < 1L) {
  stop("SWEEP_CORES must be a positive whole number, not '", sweep_cores_env, "'.", call. = FALSE)
}

stop_by_env <- Sys.getenv("STOP_BY", "")
stop_by <- if (nzchar(stop_by_env)) as.POSIXct(stop_by_env, format = "%Y-%m-%d %H:%M") else NULL
if (!is.null(stop_by) && (is.na(stop_by) || stop_by <= Sys.time())) {
  stop("STOP_BY must be a future local time as YYYY-MM-DD HH:MM, not '", stop_by_env, "'.",
       call. = FALSE)
}

# HUMAN_POP: simulated population. The default is the paper's; tests use fewer.
# It is part of each history's file name, so tests never reuse a full history.
human_pop      <- as.integer(Sys.getenv("HUMAN_POP", "100000"))
if (is.na(human_pop) || human_pop < 1L) {
  stop("HUMAN_POP must be a positive whole number, not '", Sys.getenv("HUMAN_POP"), "'.", call. = FALSE)
}
n_future_years <- 7L      # 2025-2031: campaigns in 2025, 2028 and 2031
# Columns kept from 2022 on (regular expressions; daily values). The history files
# keep every column. M_exposed_<species> is the species' ATN-exposed adults
# (S, E and I) and total_M_<species> all its adults; see annotate().
keep_columns <- c("^timestep$", "^n_detect_lm_", "^n_age_", "^n_inc_clinical_",
                  "^n_inc_severe_", "^EIR_", "^total_M_", "^M_exposed_", "^n_use_net$",
                  "^n_use_atn$", "^(region|arm|draw|seed|year_rel|pfpr2to10|clin_inc)$")
# Mass campaigns are this many calendar months before the month of peak rainfall,
# at the start of the month (so 60-90 days before the peak); see section 5.
campaign_months_before_peak <- 2L
# The seven arms in the paper. pyr_cfp_mc_atn_cd can still be run via SWEEP_ARMS.
arms           <- c("none", "pyr", "pyr_pbo", "pyr_cfp", "atn", "pyr_atn", "pyr_cfp_atn")
# SWEEP_ARMS: comma-separated subset to run (default = the seven above). Used for sensitivity
# runs that only need the affected arms (e.g. the three ATN arms for a fixed half-life, since
# the other arms do not use it and would reproduce the main run exactly).
sweep_arms_env <- Sys.getenv("SWEEP_ARMS", "")
if (nzchar(sweep_arms_env)) {
  arms <- trimws(strsplit(sweep_arms_env, ",")[[1]])
}
# Net retention: the mean duration of ITN use in Mali, 1.8 years (national median of Glover et
# al., eLife 2025, 14:RP108745, Figure 5A). That estimate already allows for nets being given at
# random, including to people who still hold one (its Eq. 6), which is how malariasimulation
# distributes them, so it is entered as it stands. One national value, and no urban/rural split,
# by design. It replaces the site file's mean_retention (about 2014 days) for every net, including
# 2000-2024, whose distributions were fitted with the site value, so the history needs
# recalibrating before the large run.
retention_years <- 1.8
# NET_RETENTION_DAYS: override the retention (days). Empty = retention_years.
ret_env            <- Sys.getenv("NET_RETENTION_DAYS", "")
retention_override <- if (nzchar(ret_env)) as.numeric(ret_env) else NULL
deltaq_use     <- 10L

# Future non-net site interventions (same defaults as previous scripts).
future_interventions <- list(
  case_management = TRUE,
  smc             = TRUE,
  irs             = FALSE,
  rtss            = FALSE,
  r21             = FALSE,
  pmc             = FALSE,
  lsm             = FALSE
)

site_obj      <- readRDS(SITE_FILE)
site_file_md5 <- unname(tools::md5sum(SITE_FILE))
country_name  <- unique(site_obj$country)[1]
message(sprintf("Country: %s (%s)", country_name, COUNTRY_ISO))
message(sprintf("Fork %s | site package %s | %s", fork_version, sim_identity$packages[["site"]],
                sim_identity$r_version))

# This script and the package code it calls, ignoring comments and layout: any
# other change starts a new results folder (settings.rds).
code_md5 <- local({
  f <- tempfile()
  on.exit(unlink(f))
  code <- file.path(CODE_PATH, c("analysis/07_run_projections.R", "R/expand_interventions.R",
                                 "R/itn_draws.R"))
  writeLines(unlist(lapply(code, function(x) deparse(parse(x, keep.source = FALSE, encoding = "UTF-8")))), f)
  unname(tools::md5sum(f))
})
message(sprintf("Realisations (draws-table rows): %d-%d (%d) | seeds %d-%d | cores: %d", draw_bounds[1], draw_bounds[2],
                length(draw_rows), SEED_BASE + draw_bounds[1], SEED_BASE + draw_bounds[2], N_CORES))
message(sprintf("Antimalarial half-life: %s | arms: %s",
                if (is.null(hl_fixed_years)) "drawn per realisation" else sprintf("fixed at %g yr", hl_fixed_years),
                paste(arms, collapse = ",")))
message(sprintf("OUT_DIR: %s", OUT_DIR))

retention_time <- if (!is.null(retention_override)) retention_override else retention_years * 365
message(sprintf("Net retention = %.1f days (%s)", retention_time,
                if (is.null(retention_override)) sprintf("%g years", retention_years) else "manual override"))

# IRS efficacy must wane with time since spraying. malariasimulation keeps each
# person's last spray date for life, so a slope of the wrong sign drives an effect
# to 0 or 1 for good, years after the product was last used. Some fitted slopes in
# the site package's IRS table have that sign: deterrence of DDT (+0.004, the Mali
# site file, 2000-2016) and Vectron, and mortality of Vectron. Each such slope is
# replaced by Actellic 300CS's, whose deterrence slope Sherrard-Smith et al. (2018,
# Nat Commun 9:4982, Supplementary Table 3) set equal to its mortality slope.
# Every replacement is reported and recorded in settings.rds; the site file itself
# is never edited. An IRS product used in any year whose curves still give an
# impossible probability (over 20 years after spraying) stops the run.
irs_gamma_actellic <- c(ls_gamma = -0.009, ms_gamma = -0.009)
iv <- site_obj$interventions
irs_corrections <- list()
for (par in names(irs_gamma_actellic)) {
  bad <- which(iv[[par]] > 0)
  for (ins in unique(iv$irs_insecticide[bad])) {
    k <- bad[iv$irs_insecticide[bad] == ins]
    irs_corrections[[length(irs_corrections) + 1L]] <- data.frame(
      insecticide = ins, parameter = par, from = paste(unique(iv[[par]][k]), collapse = ","),
      to = irs_gamma_actellic[[par]], years = paste(range(iv$year[k]), collapse = "-"),
      region_years = length(k), stringsAsFactors = FALSE)
  }
  iv[[par]][bad] <- irs_gamma_actellic[[par]]
}
irs_corrections <- do.call(rbind, irs_corrections)   # NULL if none
for (k in seq_len(NROW(irs_corrections))) {
  x <- irs_corrections[k, ]
  message(sprintf("%sIRS %s: %s %s -> %g (Actellic's) in %d region-years, %s%s",
                  if (x$parameter == "ls_gamma") "WARNING: " else "", x$insecticide, x$parameter,
                  x$from, x$to, x$region_years, x$years,
                  if (x$parameter == "ls_gamma") "; this sets the product's residual life" else ""))
}
# The problems, if any, of one IRS parameter set over 20 years after spraying, with
# the probabilities computed as malariasimulation does (R/vector_control.R).
irs_par <- c("ls_theta", "ls_gamma", "ks_theta", "ks_gamma", "ms_theta", "ms_gamma")
irs_problems <- function(x, k0 = malariasimulation::get_parameters()$k0, t = 0:(20 * 365)) {
  if (anyNA(unlist(x[irs_par]))) return("a parameter is missing")
  ls <- plogis(x$ls_theta + x$ls_gamma * t)
  ks <- k0 * plogis(x$ks_theta + x$ks_gamma * t)
  ms <- plogis(x$ms_theta + x$ms_gamma * t)
  js <- 1 - ls - ks
  ls_ <- ls * (1 - ms); ks_ <- ks * (1 - ms); js_ <- js * (1 - ms) + ms
  rs <- (1 - ks_ / k0) * js_ / (ls_ + js_)
  ss <- ks_ / k0
  problem <- c("feeding success falls with time" = x$ks_gamma < 0,
               "killed plus fed exceeds 1"        = any(js < 0),
               "P(repelled) leaves 0-1"           = any(rs < 0 | rs > 1),
               "P(survives) leaves 0-1"           = any(ss < 0 | ss > 1))
  names(problem)[problem]
}
irs_used <- unique(iv[!is.na(iv$irs_cov) & iv$irs_cov > 0, c("irs_insecticide", irs_par)])
for (k in seq_len(nrow(irs_used))) {
  problem <- irs_problems(irs_used[k, ])
  if (length(problem)) {
    stop(sprintf("IRS %s in the site file: %s. Decide how to correct it in this script ",
                 irs_used$irs_insecticide[k], paste(problem, collapse = "; ")),
         "(next to irs_gamma_actellic); never edit the site file.", call. = FALSE)
  }
}
site_obj$interventions <- iv
rm(iv, irs_used)

# ---------------------------------------------------------------------
# 1. CD coverage math
#    Mass campaigns reach campaign_cov of the population. Continuous
#    distribution (CD) maintains a floor of use equal to cd_floor_frac of the
#    post-campaign use (Mali-specific, Glover et al. 2025). The monthly top-up
#    coverage is the one that holds use at that floor when nets are given at
#    random and lost at the retention rate.
# ---------------------------------------------------------------------
campaign_cov  <- 0.95
cd_floor_frac <- 0.3
cd_interval   <- 365 / 12

d_cd     <- exp(-cd_interval / retention_time)
cd_floor <- cd_floor_frac * campaign_cov
cd_cov   <- cd_floor * (1 - d_cd) / (1 - cd_floor * d_cd)
message(sprintf("CD: campaign_cov=%.2f -> floor=%.4f, monthly_top-up_cov=%.4f",
                campaign_cov, cd_floor, cd_cov))

# ---------------------------------------------------------------------
# 2. Net-efficacy draws of Churcher et al. (2024): 1000 draws of dn0, rn0
#    and gamman (years) per net type, each a whole curve across resistance
#    levels 0, 0.01, ..., 1 (see ?read_itn_draws).
# ---------------------------------------------------------------------
itn_draws <- read_itn_draws()

# The draws for each insecticidal arm's net. Pyr-ATN and Pyr-CFP-ATN carry the
# same insecticide as Pyr and Pyr-CFP, and the same draw, so each pair differs
# only by the antimalarial.
net_of_arm <- c(pyr = "only", pyr_atn = "only", pyr_pbo = "pbo",
                pyr_cfp = "cfp", pyr_cfp_atn = "cfp")
# The draws for the net type the site file records for each distribution year.
net_of_site_type <- c(pyrethroid_only = "only", pyrethroid_pbo = "pbo",
                      pyrethroid_pyrrole = "cfp")

# dn0, rn0 and gamman (days) of draw k of a net type at each resistance in res,
# each taken at the nearest resistance level.
net_draw <- function(net, k, res) {
  x <- itn_draws[[net]]
  x <- x[x$draw == k, ]
  j <- vapply(res, function(r) which.min(abs(x$resistance - r)), integer(1))
  list(dn0 = x$dn0[j], rn0 = x$rn0[j], gamman = x$gamman[j] * 365)
}

# ---------------------------------------------------------------------
# 3. The draws table: one row per simulation
# ---------------------------------------------------------------------
draws_tab <- read.csv(DRAWS_FILE)
if (max(draw_rows) > nrow(draws_tab) || !identical(draws_tab$draw, seq_len(nrow(draws_tab)))) {
  stop(sprintf("DRAW_ROWS asks for rows up to %d, but %s has rows 1 to %d.",
               max(draw_rows), DRAWS_FILE, nrow(draws_tab)), call. = FALSE)
}
if (max(draws_tab$itn_draw) > min(vapply(itn_draws, function(x) max(x$draw), numeric(1)))) {
  stop("The draws table asks for more ITN draws than the ITN files hold.", call. = FALSE)
}

# The value of p_A each antimalarial arm uses: the ATN prior for a net with no
# insecticide, the AITN prior for one co-treated with insecticide. The ATN rows of
# the mixed-delivery arm are ATNs.
p_a_of_arm <- c(atn = "p_a_atn", pyr_atn = "p_a_aitn", pyr_cfp_atn = "p_a_aitn",
                pyr_cfp_mc_atn_cd = "p_a_atn")

# The inputs of simulation i (row i of the draws table).
row_inputs <- function(i) {
  d <- draws_tab[i, ]
  hl_years <- if (is.null(hl_fixed_years)) d$half_life_atn_years else hl_fixed_years
  list(
    kern = list(
      s_half_eip  = d$s_half_eip,
      nH_eip      = d$nH_eip,
      rho_frac    = d$rho_frac,
      s_half_pre  = d$s_half_pre,
      nH_pre      = d$nH_pre,
      B_max_post  = d$B_max_post,
      s_half_post = d$s_half_post,
      nH_post     = d$nH_post
    ),
    omega_atn = d$omega_atn,
    p_a       = c(p_a_atn = d$p_a_atn, p_a_aitn = d$p_a_aitn),
    hl_days   = hl_years * 365,   # antimalarial half-life
    itn_draw  = d$itn_draw
  )
}

render_overrides <- list(
  human_population                      = human_pop,
  prevalence_rendering_min_ages         = 2 * 365,
  prevalence_rendering_max_ages         = 10 * 365,
  age_group_rendering_min_ages          = 2 * 365,
  age_group_rendering_max_ages          = 10 * 365,
  clinical_incidence_rendering_min_ages = 0,
  clinical_incidence_rendering_max_ages = 100 * 365
)
form_overrides <- list(use_eip_hill = TRUE, use_bompard = TRUE)

# ---------------------------------------------------------------------
# 5. Load site; establish the calendar
# ---------------------------------------------------------------------
regions    <- site_obj$sites$name_1
# SWEEP_REGIONS: run only these regions (names with or without accents, e.g.
# Segou); the calendar still uses every region.
sweep_regions_env <- Sys.getenv("SWEEP_REGIONS", "")
sweep_regions <- regions
if (nzchar(sweep_regions_env)) {
  ascii <- function(x) iconv(x, to = "ASCII//TRANSLIT")
  pick  <- trimws(strsplit(sweep_regions_env, ",")[[1]])
  hit   <- match(ascii(pick), ascii(regions))
  if (anyNA(hit)) {
    stop("SWEEP_REGIONS: no region called ", paste(pick[is.na(hit)], collapse = ", "),
         ". Regions: ", paste(regions, collapse = ", "), call. = FALSE)
  }
  sweep_regions <- regions[hit]
}
message(sprintf("Regions: %s", paste(sweep_regions, collapse = ", ")))

# The sweep's summary file, named by its realisations and, for a subset of
# regions, by those regions too, so that the cluster's jobs for the same
# realisations in different regions do not collide.
region_tag <- function(region) gsub("[^A-Za-z0-9]", "", iconv(region, to = "ASCII//TRANSLIT"))
sweep_tag  <- if (setequal(sweep_regions, regions)) "" else
  paste0(paste(region_tag(sweep_regions), collapse = "-"), "_")
SWEEP_FILE <- file.path(OUT_DIR, sprintf("sweep_%srows%d-%d.rds", sweep_tag, draw_bounds[1], draw_bounds[2]))
# A finished sweep is never run again. One that left jobs failed or not started
# (STOP_BY) is finished by launching it again: its saved jobs are skipped, and its
# summary is replaced by the new one, which covers every job.
if (file.exists(SWEEP_FILE) && !isTRUE(getOption("country_test_mode")) && STAGE != "histories") {
  prev <- readRDS(SWEEP_FILE)
  if (!NROW(prev$failures) && all(prev$jobs$started | prev$jobs$skipped)) {
    stop(SWEEP_FILE, " already exists: these realisations are done. Move it, or set OUT_SUFFIX.",
         call. = FALSE)
  }
  message("Finishing the jobs that ", basename(SWEEP_FILE), " records as failed or not started")
}
start_year <- min(site_obj$interventions$year)
hist_last  <- max(site_obj$interventions$year)
future_yr0 <- hist_last + 1L
n_years    <- (hist_last - start_year + 1L) + n_future_years
n_steps    <- n_years * 365L
future_start_day <- (future_yr0 - start_year) * 365L
# Last day of the shared history segment (31 December of branch_year - 1).
history_end_day  <- (branch_year - start_year) * 365L
# Last day of the per-realisation segment: the day before the first future
# distribution, which is on future_start_day (31 December 2024, the January
# slot of the distribution grid), so every arm's own nets fall in its segment.
pre_end_day      <- future_start_day - 1L
stopifnot(branch_year > start_year, branch_year <= future_yr0, history_end_day < pre_end_day)

# END_DAY: stop the arms on this day instead of the last (n_steps), for tests and
# timing; the net schedules are unchanged.
end_day_env <- Sys.getenv("END_DAY", "")
END_DAY <- if (nzchar(end_day_env)) suppressWarnings(as.integer(end_day_env)) else n_steps
if (is.na(END_DAY) || END_DAY <= pre_end_day || END_DAY > n_steps) {
  stop(sprintf("END_DAY must be a day from %d to %d, not '%s'.", pre_end_day + 1L, n_steps, end_day_env),
       call. = FALSE)
}
if (END_DAY < n_steps) message(sprintf("END_DAY: arms stop on day %d of %d", END_DAY, n_steps))

# Mass campaign dates. Peak rainfall is found for each region from the site's
# seasonality with malariasimulation's own peak_season_offset() (day of a
# 365-day year), and the national peak is the median across regions. The
# campaign month is campaign_months_before_peak calendar months before the peak
# month, wrapping into November or December for a January or February peak, so
# the campaigns stay in the first, fourth and seventh future years and each
# precedes the following peak. Each campaign takes that month's slot on the
# monthly continuous-distribution grid (evenly spaced, so within two days of the
# 1st), replacing that month's top-up; the spacing the top-up coverage was
# calculated for is therefore unchanged. The expression below is the grid's
# own (see build_future_schedule()), so the dates match it exactly.
seas <- site_obj$seasonality$seasonality_parameters
rainfall_floor <- malariasimulation::get_parameters()$rainfall_floor
peak_days <- setNames(vapply(seq_len(nrow(seas)), function(i) {
  malariasimulation::peak_season_offset(list(
    model_seasonality = TRUE,
    g0 = seas$g0[i],
    g  = c(seas$g1[i], seas$g2[i], seas$g3[i]),
    h  = c(seas$h1[i], seas$h2[i], seas$h3[i]),
    rainfall_floor = rainfall_floor
  ))
}, numeric(1)), seas$name_1)
peak_day       <- median(peak_days)
month_starts   <- cumsum(c(1L, 31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L))
peak_month     <- findInterval(peak_day, month_starts)
campaign_month <- (peak_month - 1L - campaign_months_before_peak) %% 12L + 1L
future_campaign_days <- round(future_start_day +
  (campaign_month - 1L + 12L * seq(0L, n_future_years - 1L, by = 3L)) * cd_interval)
message(sprintf("Peak rainfall: day %g-%g by region, national median day %g (%s); campaigns in %s",
                min(peak_days), max(peak_days), peak_day, month.name[peak_month],
                month.name[campaign_month]))

# ---------------------------------------------------------------------
# 6. Net schedule builders
# ---------------------------------------------------------------------
# The site file's net distributions up to the projection window. With a
# realisation's inputs (inp), the nets distributed from branch_year take the
# efficacy of its ITN draw for the site's net type, at that year's resistance in the region
# (from the same resistance table as the future nets); earlier nets, and all nets
# when inp is NULL (the shared history), keep the site file's values.
build_past_schedule <- function(idf, region, inp = NULL) {
  pd <- idf[!is.na(idf$itn_input_dist) & idf$itn_input_dist > 0, ]
  if (nrow(pd) == 0L) return(NULL)
  sch <- list(
    timesteps = (pd$year - start_year) * 365L + pd$itn_distribution_day,
    coverages = pd$itn_input_dist,
    dn0 = pd$dn0,
    rn  = pd$rn0,
    rnm = pd$rnm,
    gam = pd$gamman * 365
  )
  drawn <- which(pd$year >= branch_year)
  if (!is.null(inp) && length(drawn)) {
    rtab <- site_obj$vectors$pyrethroid_resistance
    rtab <- rtab[rtab$name_1 == region, ]
    for (k in drawn) {
      net <- net_of_site_type[[pd$net_type[k]]]
      res <- rtab$pyrethroid_resistance[match(pd$year[k], rtab$year)]
      d   <- net_draw(net, inp$itn_draw, res)
      sch$dn0[k] <- d$dn0; sch$rn[k] <- d$rn0; sch$gam[k] <- d$gamman
    }
  }
  sch
}

# inp: the realisation's inputs, from row_inputs(); its ITN draw sets the efficacy of the
# insecticidal nets, and its antimalarial half-life the decay of an ATN's netting.
build_future_schedule <- function(arm, region, inp) {
  if (arm %in% c("none", "history")) return(NULL)

  grid    <- unique(round(seq(future_start_day, n_steps, by = cd_interval)))
  is_camp <- vapply(grid, function(t) any(abs(t - future_campaign_days) < 1L), logical(1))
  if (sum(is_camp) != length(future_campaign_days)) {
    # Cannot happen while the campaign dates are built from the grid (section 5);
    # without this check a campaign off the grid would silently not take place.
    stop(sprintf("%d of %d mass campaigns fall on the distribution grid; the rest would be lost.",
                 sum(is_camp), length(future_campaign_days)), call. = FALSE)
  }
  cov     <- ifelse(is_camp, campaign_cov, cd_cov)
  n       <- length(grid)

  rtab      <- site_obj$vectors$pyrethroid_resistance
  rtab      <- rtab[rtab$name_1 == region, ]
  grid_year <- start_year + (grid %/% 365L)
  res_grid  <- rtab$pyrethroid_resistance[match(grid_year, rtab$year)]

  if (arm == "atn") {
    # Non-insecticidal ATN: no pyrethroid mortality; rn=rnm so repellency is
    # mechanical-only; gamman set to the antimalarial half-life (inconsequential
    # since rn≈rnm, but makes intent explicit).
    dn0 <- rep(0, n)
    rn  <- rep(0.24, n)
    rnm <- rep(0.24 - 1e-9, n)
    gam <- rep(inp$hl_days, n)
  } else if (arm == "pyr_cfp_mc_atn_cd") {
    # Mixed-delivery: Pyr-CFP (insecticidal) for the mass-campaign distributions;
    # non-insecticidal ATN for the CD top-ups.
    p_mc <- net_draw("cfp", inp$itn_draw, res_grid[is_camp])
    dn0 <- rn <- gam <- numeric(n)
    rnm <- numeric(n)
    # Mass campaigns: Pyr-CFP efficacy, resistance-projected
    dn0[is_camp]  <- p_mc$dn0
    rn [is_camp]  <- p_mc$rn0
    rnm[is_camp]  <- 0.24
    gam[is_camp]  <- p_mc$gamman
    # CD top-ups: non-insecticidal ATN
    dn0[!is_camp] <- 0
    rn [!is_camp] <- 0.24
    rnm[!is_camp] <- 0.24 - 1e-9
    gam[!is_camp] <- inp$hl_days
  } else {
    # The realisation's draw of the arm's net, at each year's projected resistance.
    p_net <- net_draw(net_of_arm[[arm]], inp$itn_draw, res_grid)
    dn0   <- p_net$dn0
    rn    <- p_net$rn0
    rnm   <- rep(0.24, n)
    gam   <- p_net$gamman
  }
  list(timesteps = grid, coverages = cov, dn0 = dn0, rn = rn, rnm = rnm, gam = gam,
       is_camp = is_camp)
}

# ---------------------------------------------------------------------
# 7. Build parameters for one region x arm x realisation (draw = its row of the
#    draws table). arm = "history" builds the shared 2000-2021 segment: site net
#    efficacy throughout, no future nets, ATN off, and no realisation (draw = NULL).
# ---------------------------------------------------------------------
build_params <- function(region, arm, draw) {

  inp <- if (arm == "history") NULL else row_inputs(draw)

  site_row <- site_obj$sites[site_obj$sites$name_1 == region, , drop = FALSE]
  ms       <- site::subset_site(site_obj, site_row)

  # (IRS slopes of the wrong sign were corrected in site_obj when it was read.)
  ms_ext <- expand_interventions(ms, expand_year = n_future_years)

  fut    <- ms_ext$interventions$year >= future_yr0

  net_zero_cols <- c("itn_input_dist", "itn_use")
  for (col in net_zero_cols) {
    if (col %in% names(ms_ext$interventions)) {
      ms_ext$interventions[[col]][fut] <- 0
    }
  }

  toggle_cols <- list(
    case_management = "tx_cov",  smc  = "smc_cov", irs = "irs_cov",
    rtss            = "rtss_cov", r21 = "r21_cov", pmc = "pmc_cov",
    lsm             = "lsm_cov"
  )
  for (nm in names(toggle_cols)) {
    if (!isTRUE(future_interventions[[nm]])) {
      col <- toggle_cols[[nm]]
      if (col %in% names(ms_ext$interventions)) {
        ms_ext$interventions[[col]][fut] <- 0
      }
    }
  }

  fsch <- build_future_schedule(arm, region, inp)

  # ATN construction-time overrides for antimalarial arms, all from the realisation:
  # p_A by net type, the potency decay rate from the antimalarial half-life,
  # omega_atn and the kernel parameters. Other arms keep the defaults (ATN off).
  atn_overrides <- list()
  if (arm %in% c("atn", "pyr_atn", "pyr_cfp_atn")) {
    atn_overrides <- c(list(
      p_atn         = inp$p_a[[p_a_of_arm[[arm]]]],
      gamma_atn     = log(2) / inp$hl_days,
      omega_atn     = inp$omega_atn,
      Q0_atn        = fsch$coverages,
      t0_atn        = fsch$timesteps,
      n_atn         = length(fsch$timesteps)
    ), inp$kern)
  } else if (arm == "pyr_cfp_mc_atn_cd") {
    # Mixed-delivery: drug events = CD top-ups only; mass campaigns are non-drug displacement
    # events that overwrite ATN holders in the IBM but deliver no antimalarial.
    # atn_displace_t0/Q0 lets compute_atn_kernels collapse Q_atn_t at each campaign.
    cd_rows   <- !fsch$is_camp
    mc_rows   <-  fsch$is_camp
    cd_t0  <- fsch$timesteps[cd_rows]
    cd_Q0  <- fsch$coverages[cd_rows]
    mc_t0  <- fsch$timesteps[mc_rows]
    mc_Q0  <- fsch$coverages[mc_rows]
    atn_overrides <- c(list(
      p_atn             = inp$p_a[[p_a_of_arm[[arm]]]],
      gamma_atn         = log(2) / inp$hl_days,
      omega_atn         = inp$omega_atn,
      t0_atn            = cd_t0,
      Q0_atn            = cd_Q0,
      n_atn             = length(cd_t0),
      atn_displace_t0   = mc_t0,
      atn_displace_Q0   = mc_Q0
    ), inp$kern)
  }

  p <- site::site_parameters(
    interventions = ms_ext$interventions,
    demography    = ms_ext$demography,
    vectors       = ms_ext$vectors$vector_species,
    seasonality   = ms_ext$seasonality$seasonality_parameters,
    eir           = ms$eir$eir,
    # deltaq is the same in every segment and arm (ATN or not), so that a saved
    # state can be resumed by any arm; without an antimalarial the extra
    # exposure classes stay empty.
    overrides     = c(render_overrides, form_overrides, list(deltaq = deltaq_use),
                      atn_overrides, list(ode_max_steps = 1e8, a_tol = 0.1))
  )

  past <- build_past_schedule(ms$interventions, region, inp)

  ts  <- c(past$timesteps,  fsch$timesteps)
  cov <- c(past$coverages,  fsch$coverages)
  dn0 <- c(past$dn0,        fsch$dn0)
  rn  <- c(past$rn,         fsch$rn)
  rnm <- c(past$rnm,        fsch$rnm)
  gam <- c(past$gam,        fsch$gam)

  ord <- order(ts)
  ts <- ts[ord]; cov <- cov[ord]; dn0 <- dn0[ord]
  rn <- rn[ord]; rnm <- rnm[ord]; gam <- gam[ord]

  n_sp <- length(p$species)
  p <- set_bednets(p,
    timesteps = ts,
    coverages = cov,
    retention = retention_time,
    dn0       = matrix(rep(dn0, n_sp), ncol = n_sp),
    rn        = matrix(rep(rn,  n_sp), ncol = n_sp),
    rnm       = matrix(rep(rnm, n_sp), ncol = n_sp),
    gamman    = gam
  )
  # site_parameters() has already called set_bednets() with the site file's own
  # retention, and the fork sets the ATN coverage loss rate from the first such
  # call only, so set it here from the retention the nets now use.
  p$lambda_atn <- 1 / retention_time

  p
}

# ---------------------------------------------------------------------
# 8. The three segments
#    The history keeps ALL rendered columns; from 2022 only keep_columns.
#    A seed is set immediately before each segment that does not restore the
#    saved random-number state; run_resumable_simulation() seeds the model's
#    own C++ generator from R's, so this fixes every random draw. Resuming was
#    checked to reproduce an uninterrupted run exactly (7a.1, 2026-10-11).
# ---------------------------------------------------------------------
annotate <- function(d, region, arm, draw, seed) {
  d <- as.data.frame(d)
  for (sp in sub("^total_M_", "", grep("^total_M_", names(d), value = TRUE))) {
    d[[paste0("M_exposed_", sp)]] <-
      rowSums(d[grep(sprintf("^(Sv|Ev|Iv)_exposed_%s_count$", sp), names(d))])
  }
  d |>
    mutate(
      region    = region,
      arm       = arm,
      draw      = draw,
      seed      = seed,
      year_rel  = (timestep - future_start_day) / 365,
      pfpr2to10 = n_detect_lm_730_3649 / n_age_730_3649,
      clin_inc  = n_inc_clinical_0_1824 + n_inc_clinical_1825_5474 + n_inc_clinical_5475_36499,
      EIR_total_pp = rowSums(across(starts_with("EIR_"))) / human_pop
    )
}

trim_columns <- function(d) d[grepl(paste(keep_columns, collapse = "|"), names(d))]

history_file <- function(region) {
  file.path(HISTORY_DIR, sprintf("%s_%s_pop%d_ret%d.rds", tolower(COUNTRY_ISO), region_tag(region),
                                 human_pop, round(retention_time)))
}

job_file <- function(job, dir = JOB_DIR) {
  file.path(dir, sprintf("%s_%s_row%04d.rds", tolower(COUNTRY_ISO), region_tag(job$region), job$draw))
}

# Written under a temporary name and then renamed, so a run stopped mid-write
# leaves no half-written file behind.
save_whole <- function(x, f) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(f, ".tmp")
  saveRDS(x, tmp)
  if (!file.rename(tmp, f)) stop("Could not rename ", tmp, " to ", f, call. = FALSE)
}

minutes_since <- function(t0) (proc.time()[["elapsed"]] - t0) / 60

# Segment 1: the region's 2000-2021 history, run once and saved. An existing
# file is reused only if it was built from identical parameters and seed.
run_history <- function(region) {
  f <- history_file(region)
  p <- build_params(region, "history", NULL)
  if (file.exists(f)) {
    h <- readRDS(f)
    if (!identical(h$parameters, p) || !identical(h$seed, SEED_BASE) ||
        !identical(h$end_day, history_end_day) || !identical(h$sim_identity, sim_identity)) {
      stop(f, " was made with different inputs, or another fork, package or R version. ",
           "Move it to rebuild the history.", call. = FALSE)
    }
    return(list(region = region, file = f, reused = TRUE, minutes = h$minutes))
  }
  t0 <- proc.time()[["elapsed"]]
  # As in malariasimulation's own resume tests: build the correlations once,
  # here, and let the later segments restore them from the saved state.
  correlations <- get_correlation_parameters(p)
  invisible(correlations$mvnorm())
  set.seed(SEED_BASE)
  r <- run_resumable_simulation(history_end_day, parameters = p, correlations = correlations)
  minutes <- minutes_since(t0)
  save_whole(list(region = region, parameters = p, seed = SEED_BASE, end_day = history_end_day,
                  fork_version = fork_version, sim_identity = sim_identity, code_md5 = code_md5,
                  conan_id = conan_id, minutes = minutes, state = r$state,
                  data = annotate(r$data, region, "history", NA_integer_, SEED_BASE)), f)
  message(sprintf("History %s: %.1f min", region, minutes))
  list(region = region, file = f, reused = FALSE, minutes = minutes)
}

# Segments 2 and 3 for one region and realisation (job$draw = draws-table row):
# 2022-2024 with the realisation's seed and ITN draw, then every arm resumed
# from that end-2024 state with its random-number state restored. end_day can
# stop the arms early (for tests and timing); the net schedules are unchanged.
run_row <- function(job, end_day = n_steps) {
  region <- job$region; draw <- job$draw
  seed <- SEED_BASE + draw
  fail <- function(arm, e) data.frame(region = region, arm = arm, draw = draw,
                                      message = conditionMessage(e), stringsAsFactors = FALSE)
  minutes <- setNames(rep(NA_real_, length(arms) + 1L), c("2022-2024", arms))
  t0 <- proc.time()[["elapsed"]]
  pre <- tryCatch({
    h <- readRDS(history_file(region))
    p_pre <- build_params(region, "none", draw)
    set.seed(seed)
    run_resumable_simulation(pre_end_day, parameters = p_pre,
                             initial_state = h$state, restore_random_state = FALSE)
  }, error = function(e) fail("2022-2024", e))
  minutes[["2022-2024"]] <- minutes_since(t0)
  if (is.data.frame(pre)) return(list(failures = pre, minutes = minutes))
  out <- lapply(arms, function(arm) {
    t0 <- proc.time()[["elapsed"]]
    x <- tryCatch({
      p <- build_params(region, arm, draw)
      r <- run_resumable_simulation(end_day, parameters = p,
                                    initial_state = pre$state, restore_random_state = TRUE)
      trim_columns(annotate(r$data, region, arm, draw, seed))
    }, error = function(e) fail(arm, e))
    minutes[[arm]] <<- minutes_since(t0)
    x
  })
  failed <- vapply(out, function(x) "message" %in% names(x), logical(1))
  list(pre_projection = trim_columns(annotate(pre$data, region, "2022-2024", draw, seed)),
       results = dplyr::bind_rows(out[!failed]),
       failures = dplyr::bind_rows(out[failed]),
       minutes = minutes)
}

# TRUE if STOP_BY is set and a job started now would not finish by then, taking
# a job to last as long as the slowest one saved so far (any job, before the
# first is saved).
too_late <- function() {
  if (is.null(stop_by)) return(FALSE)
  saved <- list.files(JOB_DIR, "\\.rds$", full.names = TRUE)
  longest <- if (length(saved)) max(vapply(saved, function(g) sum(readRDS(g)$minutes, na.rm = TRUE),
                                           numeric(1))) else 0
  Sys.time() + 60 * longest > stop_by
}

# One job of the sweep: a region and a realisation (job$draw = row of the draws
# table). Skipped if its file exists; otherwise runs the realisation and saves it,
# and returns only a summary, so no results are held in memory. Any error is
# recorded as a failure rather than stopping the sweep. A job in which any
# segment failed is saved in jobs_failed/ instead of jobs/, so the next launch
# runs it again rather than skipping it.
run_job <- function(job, end_day = END_DAY) {
  f <- job_file(job)
  summary <- function(skipped, failures, minutes, started = !skipped, file = f) {
    list(region = job$region, draw = job$draw, file = file, skipped = skipped,
         started = started, failures = failures, minutes = minutes)
  }
  if (file.exists(f)) return(summary(TRUE, NULL, readRDS(f)$minutes))
  if (too_late()) return(summary(FALSE, NULL, NULL, started = FALSE))
  tryCatch({
    x <- run_row(job, end_day)
    out <- if (NROW(x$failures)) job_file(job, FAILED_DIR) else f
    save_whole(c(list(region = job$region, draw = job$draw, seed = SEED_BASE + job$draw,
                      end_day = end_day, fork_version = fork_version,
                      atnet_version = atnet_version, sim_identity = sim_identity,
                      code_md5 = code_md5, site_file_md5 = site_file_md5, conan_id = conan_id,
                      history_md5 = unname(tools::md5sum(history_file(job$region)))), x), out)
    message(sprintf("Job %s realisation %d: %.0f min, %d failure(s)", job$region, job$draw,
                    sum(x$minutes, na.rm = TRUE), NROW(x$failures)))
    summary(FALSE, x$failures, x$minutes, file = out)
  }, error = function(e) {
    summary(FALSE, data.frame(region = job$region, arm = "job", draw = job$draw,
                              message = conditionMessage(e), stringsAsFactors = FALSE), NULL)
  })
}

# ---------------------------------------------------------------------
# 9. Parallel sweep (Windows PSOCK)
#    Skipped when options(country_test_mode = TRUE).
# ---------------------------------------------------------------------
if (isTRUE(getOption("country_test_mode"))) {
  message("country_test_mode = TRUE: setup complete, skipping parallel sweep.")
} else {


# The settings that shape every job's results. The first sweep into OUT_DIR
# records them; a later batch of realisations must match, so one folder never
# mixes runs. On the cluster only the histories stage writes the file.
run_settings <- list(
  country_iso      = COUNTRY_ISO,
  fork_version     = fork_version,
  sim_identity     = sim_identity,
  site_file_md5    = site_file_md5,
  code_md5         = code_md5,
  draws_md5        = unname(tools::md5sum(DRAWS_FILE)),
  itn_draws_source = itn_draws_source()[c("repo", "commit", "md5")],
  human_pop        = human_pop,
  retention_time   = retention_time,
  campaign_cov     = campaign_cov,
  cd_floor_frac    = cd_floor_frac,
  hl_fixed_years   = hl_fixed_years,   # NULL = each realisation's drawn half-life
  arms             = arms,
  n_steps          = n_steps,
  deltaq_use       = deltaq_use,
  seed_base        = SEED_BASE,
  branch_year      = branch_year,
  future_campaign_days = future_campaign_days,
  end_day          = END_DAY,
  keep_columns     = keep_columns,
  irs_corrections  = irs_corrections
)
settings_file <- file.path(OUT_DIR, "settings.rds")
if (file.exists(settings_file)) {
  if (!identical(readRDS(settings_file), run_settings)) {
    stop(OUT_DIR, " holds runs with other settings (see its settings.rds). ",
         "Move the folder, or set OUT_SUFFIX.", call. = FALSE)
  }
} else if (STAGE == "realisations") {
  stop(settings_file, " does not exist: run the histories stage first.", call. = FALSE)
} else {
  save_whole(run_settings, settings_file)
}

# The realisations stage only uses saved histories; building one here could race
# with another job building the same region's history.
if (STAGE == "realisations") {
  no_history <- sweep_regions[!file.exists(vapply(sweep_regions, history_file, character(1)))]
  if (length(no_history)) {
    stop("No saved history for ", paste(no_history, collapse = ", "), " in ", HISTORY_DIR,
         ": run the histories stage first.", call. = FALSE)
  }
}

# One job per region and realisation; region varies fastest, so whole
# realisations finish first. The histories stage runs no job.
jobs <- if (STAGE == "histories") list() else {
  g <- expand.grid(region = sweep_regions, draw = draw_rows, stringsAsFactors = FALSE)
  lapply(seq_len(nrow(g)), function(k) as.list(g[k, ]))
}
if (length(jobs)) {
  message(sprintf("Jobs: %d (%d already saved, skipped)", length(jobs),
                  sum(vapply(jobs, function(j) file.exists(job_file(j)), logical(1)))))
}

# Never more workers than there is work to share between them. The pool is
# stopped explicitly below (before any error too) rather than by on.exit(), which
# does not reliably wait at the top level of a script.
n_workers <- min(N_CORES, max(length(sweep_regions), length(jobs)))
cl <- parallel::makeCluster(n_workers, setup_strategy = "sequential", outfile = "")
stop_pool <- function(...) {
  try(parallel::stopCluster(cl), silent = TRUE)
  if (length(list(...))) stop(..., call. = FALSE)
}

# Each worker loads the same malariasimulation as this session, and this package
# (for expand_interventions()), before the functions and data below are sent to it.
parallel::clusterExport(cl, c("FORK_PATH", "CODE_PATH"))
parallel::clusterEvalQ(cl, {
  suppressMessages({ library(site); library(dplyr); library(tidyr) })
  if (nzchar(FORK_PATH)) {
    pkgload::load_all(FORK_PATH, quiet = TRUE, compile = FALSE)
  } else {
    library(malariasimulation)
  }
  pkgload::load_all(CODE_PATH, quiet = TRUE, compile = FALSE)
})
parallel::clusterExport(cl, c(
  "build_params", "build_past_schedule", "build_future_schedule", "net_draw",
  "row_inputs", "annotate", "trim_columns", "region_tag", "history_file", "job_file",
  "save_whole", "minutes_since", "run_history", "run_row", "run_job", "too_late", "stop_by",
  "site_obj", "regions", "arms", "start_year", "hist_last", "future_yr0",
  "future_start_day", "history_end_day", "pre_end_day", "branch_year", "future_campaign_days",
  "n_steps", "n_future_years", "draws_tab", "itn_draws", "net_of_arm",
  "net_of_site_type", "p_a_of_arm", "hl_fixed_years", "SEED_BASE",
  "cd_cov", "campaign_cov", "cd_interval", "retention_time", "deltaq_use",
  "render_overrides", "form_overrides", "human_pop", "future_interventions",
  "COUNTRY_ISO", "HISTORY_DIR", "JOB_DIR", "FAILED_DIR", "END_DAY", "keep_columns",
  "fork_version", "atnet_version", "sim_identity", "code_md5", "site_file_md5", "conan_id"
))

# Segment 1: the shared histories, one per region (reused when saved; in the
# realisations stage this only checks that each saved history still matches).
t0 <- proc.time()[["elapsed"]]
hist_res <- parallel::parLapplyLB(cl, sweep_regions, function(r)
  tryCatch(run_history(r), error = function(e) list(region = r, error = conditionMessage(e))))
hist_err <- vapply(hist_res, function(x) !is.null(x$error), logical(1))
if (any(hist_err)) {
  stop_pool("History failed for: ", paste(vapply(hist_res[hist_err], function(x)
    paste0(x$region, " (", x$error, ")"), character(1)), collapse = "; "))
}
message(sprintf("Histories: %d built, %d reused, %.0f s",
                sum(!vapply(hist_res, `[[`, logical(1), "reused")),
                sum(vapply(hist_res, `[[`, logical(1), "reused")),
                proc.time()[["elapsed"]] - t0))

if (STAGE == "histories") {
  stop_pool()
  message("STAGE histories: settings and histories in place.")
} else {

# A saved job is skipped only if it was run from the history now in place (for
# example, not from a history since rebuilt after recalibration).
done <- Filter(function(j) file.exists(job_file(j)), jobs)
stale <- Filter(function(j) !identical(readRDS(job_file(j))$history_md5,
                                       unname(tools::md5sum(history_file(j$region)))), done)
if (length(stale)) {
  stop_pool(length(stale), " saved job(s) were run from a different history, e.g. ",
            job_file(stale[[1]]), ". Move the results folder, or set OUT_SUFFIX.")
}

# Segments 2 and 3: one job per region and realisation, its arms run inside the
# job. Each worker saves its own job file.
t0 <- proc.time()[["elapsed"]]
job_res <- parallel::parLapplyLB(cl, jobs, run_job)
stop_pool()
message(sprintf("Done %d region x realisation jobs (%d arms each) in %.0f s; %d not started (STOP_BY)",
                length(jobs), length(arms), proc.time()[["elapsed"]] - t0,
                sum(!vapply(job_res, `[[`, logical(1), "started") &
                    !vapply(job_res, `[[`, logical(1), "skipped"))))

failures <- dplyr::bind_rows(lapply(job_res, `[[`, "failures"))
fail_csv <- sub("\\.rds$", "_failures.csv", SWEEP_FILE)
if (nrow(failures)) {
  message(sprintf(paste("WARNING: %d segment(s) failed. Their jobs are saved in %s and will run",
                        "again at the next launch:"), nrow(failures), FAILED_DIR))
  print(failures)
  write.csv(failures, fail_csv, row.names = FALSE)
} else if (file.exists(fail_csv)) {
  # An earlier launch's failures, now all rerun (its failed job files stay in FAILED_DIR).
  file.remove(fail_csv)
}
# Run time of each segment (minutes), for the timing of larger runs.
job_minutes <- dplyr::bind_rows(lapply(job_res, function(x) {
  if (is.null(x$minutes)) return(NULL)
  data.frame(region = x$region, draw = x$draw, segment = names(x$minutes),
             minutes = unname(x$minutes), skipped = x$skipped, stringsAsFactors = FALSE)
}))
save_whole(list(
  jobs     = data.frame(region = vapply(job_res, `[[`, character(1), "region"),
                        draw   = vapply(job_res, `[[`, integer(1), "draw"),
                        file   = vapply(job_res, `[[`, character(1), "file"),
                        skipped = vapply(job_res, `[[`, logical(1), "skipped"),
                        started = vapply(job_res, `[[`, logical(1), "started"),
                        stringsAsFactors = FALSE),
  failures = failures,
  minutes  = job_minutes,
  history_minutes = vapply(hist_res, function(x) if (is.null(x$minutes)) NA_real_ else x$minutes,
                           numeric(1)),
  shape    = site_obj$shape$level_1,
  meta     = c(run_settings, list(
    param_set      = "draws",
    design         = paste0("history 2000-", branch_year - 1L, " shared per region (seed seed_base); ",
                            branch_year, "-", future_yr0 - 1L, " per realisation i = draws-table row i ",
                            "(seed seed_base + i, ITN draw); ",
                            future_yr0, " on: arms resumed from the realisation's state"),
    history_end_day = history_end_day,
    pre_end_day    = pre_end_day,
    history_files  = vapply(sweep_regions, history_file, character(1)),
    draws_file     = DRAWS_FILE,
    draw_rows      = draw_rows,
    atnet_version  = atnet_version,
    country_name   = country_name,
    start_year     = start_year,
    hist_last      = hist_last,
    future_yr0     = future_yr0,
    future_start_day = future_start_day,
    n_future_years = n_future_years,
    peak_days      = peak_days,
    peak_day       = peak_day,
    campaign_months_before_peak = campaign_months_before_peak,
    campaign_month = campaign_month,
    cd_floor       = cd_floor,
    cd_cov         = cd_cov,
    hl_tag         = hl_tag
  ))
), SWEEP_FILE)
message("Saved -> ", SWEEP_FILE)

} # end if (STAGE != "histories")
} # end if (!country_test_mode)
