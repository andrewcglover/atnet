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
#' * one value from each prior in [projection_priors()]. `p_atn_atn` and
#'   `p_atn_aitn` are the same quantile, `p_atn_u`, of the ATN and the AITN
#'   prior, so a draw with low uptake for an ATN is equally low for an AITN.
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

  p <- priors$p_atn
  p_atn_u <- stats::runif(n)
  data.frame(
    draw = seq_len(n),
    eip,
    tra,
    omega_atn = stats::rbeta(n, priors$omega_atn$shape1,
                             priors$omega_atn$shape2),
    p_atn_u = p_atn_u,
    p_atn_atn = stats::qbeta(p_atn_u, p$atn$shape1, p$atn$shape2),
    p_atn_aitn = stats::qbeta(p_atn_u, p$aitn$shape1, p$aitn$shape2),
    half_life_years = stats::rgamma(n, shape = priors$half_life_years$shape,
                                    scale = priors$half_life_years$scale)
  )
}
