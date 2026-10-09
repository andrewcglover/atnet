# Prior hyperparameters for all four Stan models.
#
# The single source of truth. The models declare the hyperparameters they use
# in their data blocks and receive them from prior_stan_data(), so a prior is
# changed here and nowhere else.
#
# Moved from priors_atn_main.R (malariasimple_ATNs/dev), with the same values.
# All four current posteriors were fitted with these values; the record each
# fit saved is kept as tests/testthat/fixtures/priors_used_20260520.json. That
# record lacks log_nH_eip$lower, which was then fixed inside eip_fit_hill.stan
# at the same value, as was sigma_exp's in both EIP models.

#' Prior hyperparameters for the Stan models
#'
#' Blocking models (TRA and laboratory TBA):
#' * `r_min` ~ Beta(0.5, 0.5), Jeffreys, on the treated arm's relative
#'   oocyst count or prevalence at the peak of the effect;
#' * `log_s_half` ~ N(log 2, 1), the time at which blocking has halved, in
#'   days (median 2, 95% within 0.28 to 14.2);
#' * `log_nH` ~ N(log 5, 1), the Hill exponent (median 5, 95% within 0.70 to
#'   35.5);
#' * `mu_L` ~ N(1, 1.5), the control arm's log mean oocyst count (TRA);
#' * `mu_p` ~ N(2, 1.5), the control arm's logit prevalence (laboratory TBA);
#' * `reciprocal_phi` ~ N+(0, 1), the negative binomial dispersion, as 1/phi.
#'
#' Both blocking models give the pre- and post-infection sides independent
#' copies of the `log_s_half` and `log_nH` priors.
#'
#' Both EIP models:
#' * `eip_floor` = 7 days, a hard minimum on the baseline EIP;
#' * `log_eip_excess` ~ N(log 3 - 0.5, 0.75), the baseline EIP's excess over
#'   the floor (baseline EIP 95% within 7.4 to 14.9 days).
#'
#' Exponential EIP model only:
#' * `r0_frac` ~ Beta(0.5, 0.5), Jeffreys, on the progression rate just after
#'   exposure as a fraction of the baseline rate;
#' * `log_zeta` ~ N(log(log 2 / 1.5), 0.75), the recovery rate, centred on a
#'   1.5-day half-life (95% within 0.34 to 6.5 days).
#'
#' Hill EIP model only, tighter than their counterparts above to keep the
#' sampler away from regions where the Hill fit is numerically unstable:
#' * `r0_frac_eip` ~ Beta(2, 5) (95% within 0.04 to 0.64), keeping prior mass
#'   away from 1, where the suppression vanishes and the kernel cannot be
#'   seen by the likelihood;
#' * `log_s_half_eip` ~ N(log 1.5, 0.4) (95% within 0.69 to 3.3 days);
#' * `log_nH_eip` ~ N(log 5, 0.4) truncated below at `lower` = 0.5, so the
#'   Hill exponent is at least exp(0.5), about 1.65 (95% within 2.3 to 11.0).
#'
#' `sigma_exp` ~ N+(0, 1) is the standard deviation of the per-experiment
#' effect in every model.
#'
#' @return A named list of hyperparameters.
#'
#' @export
atn_priors <- function() {
  list(
    r_min           = list(alpha = 0.5, beta = 0.5),
    log_s_half      = list(mean  = log(2.0), sd = 1.0),
    log_nH          = list(mean  = log(5.0), sd = 1.0),
    mu_L            = list(mean  = 1.0, sd = 1.5),
    mu_p            = list(mean  = 2.0, sd = 1.5),
    sigma_exp       = list(sd    = 1.0),
    reciprocal_phi  = list(sd    = 1.0),
    eip_floor       = 7.0,
    log_eip_excess  = list(mean  = log(3.0) - 0.5, sd = 0.75),
    r0_frac         = list(alpha = 0.5, beta = 0.5),
    log_zeta        = list(mean  = log(log(2.0) / 1.5), sd = 0.75),
    log_s_half_eip  = list(mean  = log(1.5), sd = 0.4),
    log_nH_eip      = list(mean  = log(5.0), sd = 0.4, lower = 0.5),
    r0_frac_eip     = list(alpha = 2.0, beta = 5.0)
  )
}

#' The prior hyperparameters as Stan data
#'
#' Every model is passed the full set; each declares and reads only the subset
#' it uses, and rstan ignores the rest.
#'
#' @param priors A list in the form of [atn_priors()].
#'
#' @return A named list to merge into a model's data list.
#'
#' @export
prior_stan_data <- function(priors = atn_priors()) {
  list(
    r_min_alpha       = priors$r_min$alpha,
    r_min_beta        = priors$r_min$beta,
    log_s_half_mean   = priors$log_s_half$mean,
    log_s_half_sd     = priors$log_s_half$sd,
    log_nH_mean       = priors$log_nH$mean,
    log_nH_sd         = priors$log_nH$sd,
    mu_L_mean         = priors$mu_L$mean,
    mu_L_sd           = priors$mu_L$sd,
    mu_p_mean         = priors$mu_p$mean,
    mu_p_sd           = priors$mu_p$sd,
    sigma_exp_sd      = priors$sigma_exp$sd,
    reciprocal_phi_sd = priors$reciprocal_phi$sd,
    eip_floor             = priors$eip_floor,
    log_eip_excess_mean   = priors$log_eip_excess$mean,
    log_eip_excess_sd     = priors$log_eip_excess$sd,
    r0_frac_alpha         = priors$r0_frac$alpha,
    r0_frac_beta          = priors$r0_frac$beta,
    log_zeta_mean         = priors$log_zeta$mean,
    log_zeta_sd           = priors$log_zeta$sd,
    log_s_half_eip_mean   = priors$log_s_half_eip$mean,
    log_s_half_eip_sd     = priors$log_s_half_eip$sd,
    log_nH_eip_mean       = priors$log_nH_eip$mean,
    log_nH_eip_sd         = priors$log_nH_eip$sd,
    log_nH_eip_lower      = priors$log_nH_eip$lower,
    r0_frac_eip_alpha     = priors$r0_frac_eip$alpha,
    r0_frac_eip_beta      = priors$r0_frac_eip$beta
  )
}

#' Draw blocking-function parameters from their priors
#'
#' For drawing the prior band alongside a posterior.
#'
#' @param n_prior Number of draws.
#' @param variant `"bidir_split"` for the fitted models (separate half-times
#'   and exponents before and after infection); `"bidir_shared"` shares the
#'   exponent; `"uni"` has one of each.
#' @param priors A list in the form of [atn_priors()].
#' @param seed If given, passed to [set.seed()] first.
#'
#' @return A list of vectors of length `n_prior`, including `b_max = 1 - r_min`.
#'
#' @export
sample_kernel_prior <- function(n_prior, variant, priors = atn_priors(),
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  variant <- match.arg(variant, c("uni", "bidir_shared", "bidir_split"))

  out <- list()
  out$r_min <- stats::rbeta(n_prior, priors$r_min$alpha, priors$r_min$beta)
  out$b_max <- 1 - out$r_min

  rln_s <- function() exp(stats::rnorm(n_prior, priors$log_s_half$mean,
                                       priors$log_s_half$sd))
  rln_n <- function() exp(stats::rnorm(n_prior, priors$log_nH$mean,
                                       priors$log_nH$sd))

  if (variant == "uni") {
    out$s_half <- rln_s()
    out$nH     <- rln_n()
  } else if (variant == "bidir_shared") {
    out$s_half_pre  <- rln_s()
    out$s_half_post <- rln_s()
    out$nH          <- rln_n()
  } else {
    out$s_half_pre  <- rln_s()
    out$s_half_post <- rln_s()
    out$nH_pre      <- rln_n()
    out$nH_post     <- rln_n()
  }
  out
}

#' Save the priors alongside a fit
#'
#' @param path Where to write the JSON record.
#' @param priors A list in the form of [atn_priors()].
#'
#' @return `path`, invisibly.
#'
#' @export
write_priors_json <- function(path, priors = atn_priors()) {
  jsonlite::write_json(priors, path = path, pretty = TRUE, auto_unbox = TRUE,
                       digits = 12)
  invisible(path)
}
