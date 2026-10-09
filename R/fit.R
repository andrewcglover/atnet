# Fitting the four Stan models.
#
# The sampler settings and initial values are those of the drivers the current
# posteriors were fitted with (fit_eip_hill.R, fit_eip.R and
# fit_atn_main_bidir_splitnH.R), held here so that every fit runs the same way.

#' The four Stan models and how each is fitted
#'
#' @return A named list, one element per model, giving its Stan file, which
#'   data builder it takes, its sampler control settings, and the parameters
#'   to summarise and to plot in pairs.
#'
#' @export
atn_models <- function() {
  eip_control <- list(adapt_delta = 0.95, max_treedepth = 12)
  blocking_control <- list(adapt_delta = 0.999, max_treedepth = 14)
  blocking_key <- c("r_min", "b_max", "s_half_pre", "log_s_half_pre",
                    "s_half_post", "log_s_half_post", "nH_pre", "log_nH_pre",
                    "nH_post", "log_nH_post", "sigma_exp")
  list(
    eip_hill = list(
      file = "eip_fit_hill.stan", data = "eip", control = eip_control,
      key = c("log_eip_excess", "eip_excess", "eip_baseline", "rho", "rho0",
              "r0_frac", "B", "log_s_half", "s_half", "log_nH", "nH",
              "sigma_exp"),
      pairs = c("log_eip_excess", "r0_frac", "log_s_half", "log_nH",
                "sigma_exp")
    ),
    eip_exponential = list(
      file = "si_comparison/eip_fit.stan", data = "eip", control = eip_control,
      key = c("log_eip_excess", "eip_excess", "eip_baseline", "rho", "rho0",
              "r0_frac", "zeta", "log_zeta", "zeta_halflife", "sigma_exp"),
      pairs = c("log_eip_excess", "r0_frac", "log_zeta", "sigma_exp")
    ),
    tra = list(
      file = "int_fit_nb_rate_hill_bidir_splitnH.stan", data = "tra",
      control = blocking_control,
      key = c("mu_L", blocking_key, "phi", "reciprocal_phi"),
      pairs = c("mu_L", "r_min", "log_s_half_pre", "log_s_half_post",
                "log_nH_pre", "log_nH_post", "reciprocal_phi", "sigma_exp")
    ),
    tba_lab = list(
      file = "si_comparison/int_fit_linear_hill_bidir_splitnH.stan",
      data = "tba_lab", control = blocking_control,
      key = c("mu_p", blocking_key),
      pairs = c("mu_p", "r_min", "log_s_half_pre", "log_s_half_post",
                "log_nH_pre", "log_nH_post", "sigma_exp")
    )
  )
}

#' Initial values for the Hill EIP model
#'
#' Stan's default initial values can place the Hill exponent where the model's
#' integral is very slow to evaluate, which once left a fit stuck at its first
#' iteration for a day. Each chain instead starts near the prior modes, with
#' small differences between chains so that R-hat can still detect a failure
#' to converge.
#'
#' @param priors A list in the form of [atn_priors()].
#' @param K Number of experiments.
#'
#' @return A function of no arguments returning one chain's initial values.
#'
#' @export
eip_hill_inits <- function(priors, K) {
  force(priors)
  force(K)
  function() {
    r0 <- priors$r0_frac_eip
    list(
      log_eip_excess = priors$log_eip_excess$mean + stats::rnorm(1, 0, 0.10),
      r0_frac        = max(0.05, min(0.95, r0$alpha / (r0$alpha + r0$beta) +
                                         stats::runif(1, -0.05, 0.05))),
      log_s_half     = priors$log_s_half_eip$mean + stats::rnorm(1, 0, 0.10),
      log_nH         = priors$log_nH_eip$mean + stats::rnorm(1, 0, 0.10),
      sigma_exp      = max(0.05, 0.3 + stats::runif(1, -0.05, 0.05)),
      z_exp          = stats::rnorm(K, 0, 0.3)
    )
  }
}

#' Fit one of the four models
#'
#' @param model A name from [atn_models()].
#' @param table The sporozoite table for the EIP models, or the oocyst table
#'   for the blocking models.
#' @param priors A list in the form of [atn_priors()].
#' @param chains,iter,warmup,seed Passed to `rstan::sampling()`; `iter`
#'   includes the warmup.
#' @param cores Number of chains run at once.
#' @param ... Further arguments to `rstan::sampling()`.
#'
#' @return A `stanfit` object.
#'
#' @export
fit_atn_model <- function(model, table, priors = atn_priors(), chains = 4L,
                          iter = 3000L, warmup = 1000L, seed = 2025L,
                          cores = min(chains, parallel::detectCores()), ...) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Fitting needs the rstan package.", call. = FALSE)
  }
  spec <- atn_models()[[model]]
  if (is.null(spec)) {
    stop("Unknown model '", model, "'.", call. = FALSE)
  }

  data <- switch(spec$data,
    eip = eip_stan_data(table, priors),
    tra = blocking_stan_data(table, "tra", priors),
    tba_lab = blocking_stan_data(table, "tba_lab", priors)
  )
  init <- if (model == "eip_hill") eip_hill_inits(priors, data$K) else "random"

  compiled <- rstan::stan_model(
    file = system.file("stan", spec$file, package = "atnet", mustWork = TRUE),
    model_name = model
  )
  set.seed(seed)
  rstan::sampling(compiled, data = data, chains = chains, iter = iter,
                  warmup = warmup, seed = seed, init = init,
                  control = spec$control, cores = cores, ...)
}
