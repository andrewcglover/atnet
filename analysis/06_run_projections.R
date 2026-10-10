# =====================================================================
# 06_run_projections.R
#   Country-level admin-1 (no urban/rural split): 7-year forward projections
#   (2025-2031 for Mali) under SEVEN future net arms using Churcher 2024
#   MEDIAN net-efficacy parameters (churcher2024_{only,pbo,cfp}.csv — one row
#   per resistance level, no draws). Mass campaigns every three years (2025,
#   2028, 2031), each at the start of the calendar month two months before the
#   month of peak rainfall (1 June for Mali, whose regions peak 8-13 August),
#   with monthly continuous distribution between:
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
#   Antimalarial decay is fixed at log(2)/(2.64*365) — NOT derived from
#   any read-in net parameter (see ANTIMAL_HL_DAYS below).
#
#   Change COUNTRY_ISO to run any available country. Default: MLI (Mali).
#
#   Moved from the malariasimulation fork's dev/c24med_projection_run.R
#   (atn-dev at 967d2bd). The model code is unchanged; what changed is where
#   things are read from and written to:
#     site file         data_private/site_files/<ISO>.rds
#                       (copied from the fork's dev/site_files/without_split/)
#     net efficacy      data_private/itn_params/churcher2024_{only,pbo,cfp}.csv
#                       (copied from the fork's dev/itn_params/)
#     ATN kernel fits   outputs/fits/<stamp>/{eip_hill,tra}.rds, the newest
#                       full fit unless FIT_STAMP names one
#     results           outputs/projections/
#   and expand_interventions() is now a function of this package.
#
#   Run with the working directory at the root of this package. The
#   malariasimulation fork is loaded from source when ATN_FORK_PATH is set
#   (on the laptop), and otherwise from the installed library (on the
#   cluster, where the fork is installed at a fixed tag):
#       ATN_FORK_PATH=C:/Users/ag4218/Local/GitHub/malariasimulation \
#         Rscript analysis/06_run_projections.R
#   Either way the run stops unless the loaded malariasimulation is the fork.
# =====================================================================

devtools::load_all(quiet = TRUE)   # atnet: expand_interventions()
ATNET_PATH <- normalizePath(".")
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
library(site)
library(dplyr); library(tidyr)
library(parallel)

# ---------------------------------------------------------------------
# 0. Paths & top-level settings
# ---------------------------------------------------------------------
COUNTRY_ISO <- "MLI"   # ISO3 of site file under data_private/site_files/
SITE_FILE   <- file.path("data_private", "site_files", sprintf("%s.rds", COUNTRY_ISO))

# FIT_STAMP: the ATN kernel fit to use, a folder name under outputs/fits/.
# Empty = the newest full fit.
fit_stamp <- Sys.getenv("FIT_STAMP", "")
if (!nzchar(fit_stamp)) {
  stamps <- sort(list.files(file.path("outputs", "fits"), pattern = "^[0-9]{8}_[0-9]{4}$"))
  if (!length(stamps)) {
    stop("No full fit in outputs/fits/. Run analysis/03_fit_models.R first.", call. = FALSE)
  }
  fit_stamp <- stamps[length(stamps)]
}
FIT_DIR <- file.path("outputs", "fits", fit_stamp)

# ANTIMAL_HL_YEARS: antimalarial half-life in years. Override via env var to run
# sensitivity sweeps without editing this file (avoids edit-races during concurrent runs).
#   Default: 2.64 yr (same as before). Sensitivity: 5, 1.
#   Launch as: ANTIMAL_HL_YEARS=5 SWEEP_CORES=10 Rscript analysis/06_run_projections.R
hl_years       <- as.numeric(Sys.getenv("ANTIMAL_HL_YEARS", "2.64"))
hl_tag         <- gsub("\\.", "p", format(hl_years, trim = TRUE))  # 2.64->"2p64", 5->"5", 1->"1"

# ATN_OMEGA: omega_atn, the probability that a repelled mosquito contacts the net (SI omega;
# p_C = s_N + omega * r_N), the same for all net types and ages. Default 0.9; plausible range
# 0.8 to 1.0. Affects every ATN arm. When not 0.9, OUT_FILE gets an _omega{tag} suffix.
# Replaces ATN_CHEM_DOSE (chem_dose_atn, removed 2026-10-09). That variable is refused if set,
# so an old driver script (e.g. run_c24med_overnight.sh) cannot silently run under the new
# contact model with its old meaning.
if (nzchar(Sys.getenv("ATN_CHEM_DOSE", ""))) {
  stop("ATN_CHEM_DOSE was removed on 2026-10-09 (chem_dose_atn replaced by omega_atn); ",
       "set ATN_OMEGA instead")
}
omega_atn_use  <- as.numeric(Sys.getenv("ATN_OMEGA", "0.9"))
omega_tag      <- gsub("\\.", "p", format(omega_atn_use, trim = TRUE))  # 0.8->"0p8", 1->"1"
omega_suffix   <- if (omega_atn_use != 0.9) sprintf("_omega%s", omega_tag) else ""
# OUT_SUFFIX: extra tag appended to the filename (before .rds) so non-default runs
# (e.g. retention sweeps, single-arm add-on runs) never overwrite existing results.
out_extra_suffix <- Sys.getenv("OUT_SUFFIX", "")
OUT_FILE       <- sprintf("outputs/projections/%s_c24med_projection_results_hl%s%s%s.rds",
                          tolower(COUNTRY_ISO), hl_tag, omega_suffix, out_extra_suffix)

# SWEEP_CORES: number of parallel workers. Override via env var for concurrent runs.
sweep_cores_env <- Sys.getenv("SWEEP_CORES", "")
N_CORES <- if (nchar(sweep_cores_env) > 0L) {
  as.integer(sweep_cores_env)
} else {
  min(16L, max(1L, parallel::detectCores() - 1L))
}

human_pop      <- 100000L
n_future_years <- 7L      # 2025-2031: campaigns in 2025, 2028 and 2031
# Mass campaigns are this many calendar months before the month of peak rainfall,
# at the start of the month (so 60-90 days before the peak); see section 5.
campaign_months_before_peak <- 2L
# The seven arms in the paper. pyr_cfp_mc_atn_cd can still be run via SWEEP_ARMS.
arms           <- c("none", "pyr", "pyr_pbo", "pyr_cfp", "atn", "pyr_atn", "pyr_cfp_atn")
# SWEEP_ARMS: comma-separated subset to run (default = the seven above). Used for sensitivity
# sweeps that only need the affected arms (e.g. the ATN arms for ATN_OMEGA).
sweep_arms_env <- Sys.getenv("SWEEP_ARMS", "")
if (nzchar(sweep_arms_env)) {
  arms <- trimws(strsplit(sweep_arms_env, ",")[[1]])
}
# NET_RETENTION_DAYS: override the site mean_retention (days). Empty = use site value.
# Used for the retention sensitivity (site value treated as a half-life -> scale by log(2)).
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

site_obj     <- readRDS(SITE_FILE)
country_name <- unique(site_obj$country)[1]
message(sprintf("Country: %s (%s)", country_name, COUNTRY_ISO))
message(sprintf("ANTIMAL HL: %.2f yr (tag: hl%s) | cores: %d", hl_years, hl_tag, N_CORES))
message(sprintf("omega_atn: %g | arms: %s", omega_atn_use, paste(arms, collapse = ",")))
message(sprintf("ATN kernel fits: %s", FIT_DIR))
message(sprintf("OUT_FILE: %s", OUT_FILE))

site_retention <- unique(site_obj$interventions$mean_retention)
stopifnot(length(site_retention) == 1L)
retention_time <- if (!is.null(retention_override)) retention_override else site_retention
message(sprintf("Net retention = %.1f days (%s)", retention_time,
                if (is.null(retention_override)) "site mean_retention" else "manual override"))

# ---------------------------------------------------------------------
# 1. CD coverage math
# ---------------------------------------------------------------------
campaign_cov <- 0.9
CD_frac      <- 0.155
cd_interval  <- 365 / 12

d_cd     <- exp(-cd_interval / retention_time)
cd_floor <- CD_frac * campaign_cov / (1 - CD_frac * (1 - campaign_cov))
cd_cov   <- cd_floor * (1 - d_cd) / (1 - cd_floor * d_cd)
message(sprintf("CD: campaign_cov=%.2f -> floor=%.4f, monthly_top-up_cov=%.4f",
                campaign_cov, cd_floor, cd_cov))

# ---------------------------------------------------------------------
# 2. Churcher 2024 median net-efficacy params
#    gamman in the CSV files is in years; * 365 -> days.
#    One row per resistance level — med_net() works unchanged.
# ---------------------------------------------------------------------
itn_dir   <- file.path("data_private", "itn_params")
only_pars <- read.csv(file.path(itn_dir, "churcher2024_only.csv")) |> mutate(gamman = gamman * 365)
pbo_pars  <- read.csv(file.path(itn_dir, "churcher2024_pbo.csv"))  |> mutate(gamman = gamman * 365)
cfp_pars  <- read.csv(file.path(itn_dir, "churcher2024_cfp.csv"))  |> mutate(gamman = gamman * 365)

# Pick net params at the nearest resistance-grid level (grid at 0.01 increments).
# gamman from the CSV is already a mean (not a half-life); used directly for set_bednets.
med_net <- function(pars, res) {
  grid    <- sort(unique(pars$resistance))
  res_use <- grid[which.min(abs(grid - res))]
  pars |> filter(resistance == res_use) |>
    summarise(dn0 = median(dn0), rn0 = median(rn0), gamman = median(gamman))
}

# ---------------------------------------------------------------------
# 3. Antimalarial decay — derived from hl_years (env-var overridable, see top).
#    Default half-life = 2.64 years. Sensitivity: 5 yr, 1 yr.
# ---------------------------------------------------------------------
ANTIMAL_HL_DAYS <- hl_years * 365    # antimalarial half-life (days)
gamma_atn       <- log(2) / ANTIMAL_HL_DAYS   # potency decay rate (/day)

# ---------------------------------------------------------------------
# 4. Median ATN kernel parameters (Stan posteriors — unchanged)
# ---------------------------------------------------------------------
eip_samples <- rstan::extract(readRDS(file.path(FIT_DIR, "eip_hill.rds")))
tra_samples <- rstan::extract(readRDS(file.path(FIT_DIR, "tra.rds")))
atn_kern <- list(
  s_half_eip  = median(eip_samples$s_half),
  nH_eip      = median(eip_samples$nH),
  rho_frac    = median(eip_samples$r0_frac),
  s_half_pre  = median(tra_samples$s_half_pre),
  nH_pre      = median(tra_samples$nH_pre),
  B_max_post  = median(tra_samples$b_max),
  s_half_post = median(tra_samples$s_half_post),
  nH_post     = median(tra_samples$nH_post)
)

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
start_year <- min(site_obj$interventions$year)
hist_last  <- max(site_obj$interventions$year)
future_yr0 <- hist_last + 1L
n_years    <- (hist_last - start_year + 1L) + n_future_years
n_steps    <- n_years * 365L
future_start_day <- (future_yr0 - start_year) * 365L

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
build_past_schedule <- function(idf) {
  pd <- idf[!is.na(idf$itn_input_dist) & idf$itn_input_dist > 0, ]
  if (nrow(pd) == 0L) return(NULL)
  list(
    timesteps = (pd$year - start_year) * 365L + pd$itn_distribution_day,
    coverages = pd$itn_input_dist,
    dn0 = pd$dn0,
    rn  = pd$rn0,
    rnm = pd$rnm,
    gam = pd$gamman * 365
  )
}

build_future_schedule <- function(arm, region) {
  if (arm == "none") return(NULL)

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
    gam <- rep(ANTIMAL_HL_DAYS, n)
  } else if (arm == "pyr_cfp_mc_atn_cd") {
    # Mixed-delivery: Pyr-CFP (insecticidal) for mass-campaign rows;
    # non-insecticidal ATN for CD top-up rows.
    p_list_mc <- lapply(res_grid[is_camp], function(r) med_net(cfp_pars, r))
    dn0 <- rn <- gam <- numeric(n)
    rnm <- numeric(n)
    # MC rows: Pyr-CFP efficacy, resistance-projected
    dn0[is_camp]  <- vapply(p_list_mc, `[[`, numeric(1), "dn0")
    rn [is_camp]  <- vapply(p_list_mc, `[[`, numeric(1), "rn0")
    rnm[is_camp]  <- 0.24
    gam[is_camp]  <- vapply(p_list_mc, `[[`, numeric(1), "gamman")
    # CD rows: non-insecticidal ATN
    dn0[!is_camp] <- 0
    rn [!is_camp] <- 0.24
    rnm[!is_camp] <- 0.24 - 1e-9
    gam[!is_camp] <- ANTIMAL_HL_DAYS
  } else {
    # Select the appropriate Churcher 2024 efficacy table for this arm.
    # pyr_atn uses pyr-only params; pyr_cfp_atn uses cfp params.
    pars_use <- switch(arm,
      pyr         = only_pars,
      pyr_atn     = only_pars,
      pyr_pbo     = pbo_pars,
      pyr_cfp     = cfp_pars,
      pyr_cfp_atn = cfp_pars
    )
    p_list <- lapply(res_grid, function(r) med_net(pars_use, r))
    dn0    <- vapply(p_list, `[[`, numeric(1), "dn0")
    rn     <- vapply(p_list, `[[`, numeric(1), "rn0")
    rnm    <- rep(0.24, n)
    gam    <- vapply(p_list, `[[`, numeric(1), "gamman")
  }
  list(timesteps = grid, coverages = cov, dn0 = dn0, rn = rn, rnm = rnm, gam = gam,
       is_camp = is_camp)
}

# ---------------------------------------------------------------------
# 7. Build parameters for one region x arm
# ---------------------------------------------------------------------
build_params <- function(region, arm) {

  site_row <- site_obj$sites[site_obj$sites$name_1 == region, , drop = FALSE]
  ms       <- site::subset_site(site_obj, site_row)

  ms_ext <- expand_interventions(ms, expand_year = n_future_years)

  # WORKAROUND (2026-06-22): the MLI site file has ms_gamma = +0.004 for all DDT
  # IRS rows (2000-2016) — wrong sign, so deterrence rises to 1 instead of decaying
  # to 0 (inflates Z, suppresses biting). actellic rows (2017+) are correct.
  # Flip any positive ms_gamma to negative (magnitude preserved). Raise upstream
  # with the site-file maintainer. -abs() is idempotent (correct rows unaffected).
  ms_ext$interventions$ms_gamma <- -abs(ms_ext$interventions$ms_gamma)

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

  fsch <- build_future_schedule(arm, region)

  # ATN construction-time overrides for antimalarial arms.
  # gamma_atn is fixed (ANTIMAL_HL_DAYS), not derived from net params.
  atn_overrides <- list()
  if (arm %in% c("atn", "pyr_atn", "pyr_cfp_atn")) {
    atn_overrides <- c(list(
      p_atn         = 0.9,
      deltaq        = deltaq_use,
      gamma_atn     = gamma_atn,
      omega_atn     = omega_atn_use,
      Q0_atn        = fsch$coverages,
      t0_atn        = fsch$timesteps,
      n_atn         = length(fsch$timesteps)
    ), atn_kern)
  } else if (arm == "pyr_cfp_mc_atn_cd") {
    # Mixed-delivery: drug events = CD rows only; MC rows are non-drug displacement
    # events that overwrite ATN holders in the IBM but deliver no antimalarial.
    # atn_displace_t0/Q0 lets compute_atn_kernels collapse Q_atn_t at each campaign.
    cd_rows   <- !fsch$is_camp
    mc_rows   <-  fsch$is_camp
    cd_t0  <- fsch$timesteps[cd_rows]
    cd_Q0  <- fsch$coverages[cd_rows]
    mc_t0  <- fsch$timesteps[mc_rows]
    mc_Q0  <- fsch$coverages[mc_rows]
    atn_overrides <- c(list(
      p_atn             = 0.9,
      deltaq            = deltaq_use,
      gamma_atn         = gamma_atn,
      omega_atn         = omega_atn_use,
      t0_atn            = cd_t0,
      Q0_atn            = cd_Q0,
      n_atn             = length(cd_t0),
      atn_displace_t0   = mc_t0,
      atn_displace_Q0   = mc_Q0
    ), atn_kern)
  }

  p <- site::site_parameters(
    interventions = ms_ext$interventions,
    demography    = ms_ext$demography,
    vectors       = ms_ext$vectors$vector_species,
    seasonality   = ms_ext$seasonality$seasonality_parameters,
    eir           = ms$eir$eir,
    overrides     = c(render_overrides, form_overrides, atn_overrides,
                      list(ode_max_steps = 1e8, a_tol = 0.1))
  )

  past <- build_past_schedule(ms$interventions)

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

  p
}

# ---------------------------------------------------------------------
# 8. Run one region x arm row
#    Keeps ALL run_simulation columns (incl. Sv/Ev/Iv_*_count, EIR_<species>).
# ---------------------------------------------------------------------
run_one <- function(row) {
  region <- row$region; arm <- row$arm
  tryCatch({
    p <- build_params(region, arm)
    r <- run_simulation(timesteps = n_steps, parameters = p)
    as.data.frame(r) |>
      mutate(
        region    = region,
        arm       = arm,
        year_rel  = (timestep - future_start_day) / 365,
        pfpr2to10 = n_detect_lm_730_3649 / n_age_730_3649,
        clin_inc  = n_inc_clinical_0_1824 + n_inc_clinical_1825_5474 + n_inc_clinical_5475_36499,
        EIR_total_pp = rowSums(across(starts_with("EIR_"))) / human_pop
      )
  }, error = function(e) {
    list(.__error__ = TRUE, region = region, arm = arm, message = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------
# 9. Parallel sweep (Windows PSOCK)
#    Skipped when options(country_test_mode = TRUE).
# ---------------------------------------------------------------------
if (isTRUE(getOption("country_test_mode"))) {
  message("country_test_mode = TRUE: setup complete, skipping parallel sweep.")
} else {

grid_df <- expand.grid(region = regions, arm = arms, stringsAsFactors = FALSE)
rows    <- split(grid_df, seq_len(nrow(grid_df)))

cl <- parallel::makeCluster(N_CORES, setup_strategy = "sequential", outfile = "")
on.exit(parallel::stopCluster(cl), add = TRUE)

# Each worker loads the same malariasimulation as this session, and this package
# (for expand_interventions), before the functions below are sent to it.
parallel::clusterExport(cl, c("FORK_PATH", "ATNET_PATH"))
parallel::clusterEvalQ(cl, {
  suppressMessages({ library(site); library(dplyr); library(tidyr) })
  if (nzchar(FORK_PATH)) {
    pkgload::load_all(FORK_PATH, quiet = TRUE, compile = FALSE)
  } else {
    library(malariasimulation)
  }
  pkgload::load_all(ATNET_PATH, quiet = TRUE, compile = FALSE)
})
parallel::clusterExport(cl, c(
  "build_params", "build_past_schedule", "build_future_schedule", "med_net",
  "run_one",
  "site_obj", "regions", "start_year", "hist_last", "future_yr0",
  "future_start_day", "future_campaign_days", "n_steps", "n_future_years",
  "cfp_pars", "only_pars", "pbo_pars", "atn_kern", "gamma_atn", "ANTIMAL_HL_DAYS",
  "omega_atn_use",
  "cd_cov", "campaign_cov", "cd_interval", "retention_time", "deltaq_use",
  "render_overrides", "form_overrides", "human_pop",
  "future_interventions"
))

t0 <- proc.time()[["elapsed"]]
res_list <- parallel::parLapplyLB(cl, rows, run_one)
message(sprintf("Done %d runs in %.0f s", length(rows),
                proc.time()[["elapsed"]] - t0))

is_error <- vapply(res_list, function(x) isTRUE(x$.__error__), logical(1))
if (any(is_error)) {
  failures <- dplyr::bind_rows(lapply(res_list[is_error], function(x)
    data.frame(region = x$region, arm = x$arm, message = x$message,
               stringsAsFactors = FALSE)))
  message(sprintf("WARNING: %d job(s) failed and excluded from output:", sum(is_error)))
  print(failures)
  write.csv(failures,
            sub("\\.rds$", "_failures.csv", OUT_FILE),
            row.names = FALSE)
}
df_full <- dplyr::bind_rows(res_list[!is_error])
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  results = df_full,
  shape   = site_obj$shape$level_1,
  meta    = list(
    param_set      = "churcher2024_median",
    country_iso    = COUNTRY_ISO,
    country_name   = country_name,
    fit_dir        = FIT_DIR,
    start_year     = start_year,
    hist_last      = hist_last,
    future_yr0     = future_yr0,
    future_start_day = future_start_day,
    n_future_years = n_future_years,
    peak_days      = peak_days,
    peak_day       = peak_day,
    campaign_months_before_peak = campaign_months_before_peak,
    campaign_month = campaign_month,
    future_campaign_days = future_campaign_days,
    human_pop      = human_pop,
    retention_time = retention_time,
    cd_floor       = cd_floor,
    cd_cov         = cd_cov,
    gamma_atn      = gamma_atn,
    ANTIMAL_HL_DAYS = ANTIMAL_HL_DAYS,
    hl_years       = hl_years,
    hl_tag         = hl_tag,
    omega_atn      = omega_atn_use,
    omega_tag      = omega_tag,
    arms           = arms
  )
), OUT_FILE)
message("Saved -> ", OUT_FILE)

} # end if (!country_test_mode)
