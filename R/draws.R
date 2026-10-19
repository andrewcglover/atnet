# The antimalarial columns of the draws table used by the projections.

#' Draw the antimalarial inputs of the projections
#'
#' Each row is one draw of every uncertain antimalarial input:
#' * one whole post-warmup iteration of the EIP fit under the Hill
#'   assumption, so its parameters stay together: `s_half_eip`, `nH_eip` and
#'   `rho_frac` (the fit's `s_half`, `nH` and `r0_frac`);
#' * one whole iteration of the TRA fit, chosen independently of the EIP
#'   iteration: `s_half_pre`, `nH_pre`, `B_max_post` (the fit's `b_max`),
#'   `s_half_post` and `nH_post`, on the laboratory scale;
#' * one value from each prior in [projection_priors()]. Each is drawn as a
#'   uniform quantile, recorded alongside the value (`omega_quantile`,
#'   `p_a_quantile`, `half_life_quantile`), and the value is that quantile of
#'   the prior. `p_a_atn` and `p_a_aitn` are the same quantile, so a draw with
#'   low uptake for an ATN is equally low for an AITN. Keeping the quantiles
#'   means a later revision of a prior can be applied to the same realisations
#'   without disturbing their order: realisation i stays at the same point of
#'   whatever prior is used.
#'
#' Iterations are sampled without replacement, and `eip_chain`, `eip_iter`,
#' `tra_chain` and `tra_iter` record which chain and post-warmup iteration
#' each row came from. Columns are named as the malariasimulation fork's
#' parameters where one exists. Uses the current random number stream: call
#' [set.seed()] first.
#'
#' @param eip_draws,tra_draws Post-warmup draws of the EIP (Hill) and TRA
#'   fits, as returned by `rstan::extract(fit, permuted = FALSE)`: an array of
#'   iterations x chains x parameters.
#' @param n Number of rows.
#' @param priors A list in the form of [projection_priors()].
#'
#' @return A data frame with `n` rows, the first column `draw` = 1, ..., `n`.
#'
#' @export
draw_antimalarial_inputs <- function(eip_draws, tra_draws, n,
                                     priors = projection_priors()) {
  pick <- function(draws, pars) {
    n_iter <- dim(draws)[1]
    n_avail <- n_iter * dim(draws)[2]
    if (n > n_avail) {
      stop(sprintf("Asked for %d draws but the fit has %d.", n, n_avail),
           call. = FALSE)
    }
    missing <- setdiff(pars, dimnames(draws)[[3]])
    if (length(missing)) {
      stop("The fit has no parameter ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    k <- sample.int(n_avail, n)
    iter <- (k - 1L) %% n_iter + 1L
    chain <- (k - 1L) %/% n_iter + 1L
    vals <- lapply(match(pars, dimnames(draws)[[3]]),
                   function(j) unname(draws[cbind(iter, chain, j)]))
    c(list(chain = chain, iter = iter), stats::setNames(vals, names(pars)))
  }

  eip <- pick(eip_draws, c(s_half_eip = "s_half", nH_eip = "nH",
                           rho_frac = "r0_frac"))
  tra <- pick(tra_draws, c(s_half_pre = "s_half_pre", nH_pre = "nH_pre",
                           B_max_post = "b_max", s_half_post = "s_half_post",
                           nH_post = "nH_post"))
  names(eip)[1:2] <- c("eip_chain", "eip_iter")
  names(tra)[1:2] <- c("tra_chain", "tra_iter")

  # Each prior is drawn as a uniform quantile, which is kept, so a revision of
  # any prior can be applied to the same realisations without reordering them.
  # p_a_quantile is drawn first, and before the other two, so that a table made
  # with the same seed keeps the EIP and TRA iterations and the p_A quantiles it
  # had when only the priors changed (2026-10-19).
  p <- priors$p_a
  hl <- priors$half_life_atn_years
  p_a_quantile       <- stats::runif(n)
  omega_quantile     <- stats::runif(n)
  half_life_quantile <- stats::runif(n)
  data.frame(
    draw = seq_len(n),
    eip,
    tra,
    omega_quantile = omega_quantile,
    omega_atn = stats::qbeta(omega_quantile, priors$omega_atn$shape1,
                             priors$omega_atn$shape2),
    p_a_quantile = p_a_quantile,
    p_a_atn = stats::qbeta(p_a_quantile, p$atn$shape1, p$atn$shape2),
    p_a_aitn = stats::qbeta(p_a_quantile, p$aitn$shape1, p$aitn$shape2),
    half_life_quantile = half_life_quantile,
    half_life_atn_years = stats::qlnorm(half_life_quantile,
                                        meanlog = hl$meanlog,
                                        sdlog = hl$sdlog)
  )
}
