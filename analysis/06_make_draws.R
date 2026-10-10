# 06_make_draws.R
#
# The draws table for the projections: one row per simulation. Simulation i of
# every arm uses row i, so arms are compared on matched inputs.
#
# Antimalarial columns (see ?draw_antimalarial_inputs):
#   eip_chain, eip_iter, s_half_eip, nH_eip, rho_frac
#       one whole post-warmup iteration of the EIP fit under the Hill
#       assumption (eip_hill.rds)
#   tra_chain, tra_iter, s_half_pre, nH_pre, B_max_post, s_half_post, nH_post
#       one whole iteration of the TRA fit (tra.rds), chosen independently of
#       the EIP iteration; on the laboratory scale, since malariasimulation
#       maps them to the field scale itself (use_bompard = TRUE)
#   omega_atn, p_atn_u, p_atn_atn, p_atn_aitn, half_life_years
#       one value from each prior in projection_priors(); p_atn_atn and
#       p_atn_aitn are the same quantile p_atn_u of the ATN and AITN priors
#
# Insecticide-treated net column:
#   itn_draw
#       the draw of the Churcher et al. (2024) ITN parameters (dn0, rn0,
#       gamman) to use for every net type; see ?itn_draws_source. Their 1000
#       draws are each a whole curve across pyrethroid resistance and are
#       linked across net types by draw number, so row i uses draw i of every
#       net type, and none is resampled. Their order is unrelated to their
#       values, so the first N rows use a random subset. The draws are kept,
#       unmodified, in inst/extdata/churcher2024_itn_draws/.
#
# Written to inst/extdata/projection_draws.csv, which is tracked. The fit is
# named below rather than taken as the newest, and the seed is fixed, so
# rerunning this script reproduces the table exactly:
#   Rscript analysis/06_make_draws.R
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages(library(rstan))

fit_stamp <- "20261009_1853"
n_draws   <- 1000L
seed      <- 20261010L
out_file  <- file.path("inst", "extdata", "projection_draws.csv")

fit_dir <- file.path("outputs", "fits", fit_stamp)
eip_draws <- rstan::extract(readRDS(file.path(fit_dir, "eip_hill.rds")),
                            permuted = FALSE)
tra_draws <- rstan::extract(readRDS(file.path(fit_dir, "tra.rds")),
                            permuted = FALSE)

set.seed(seed)
draws <- draw_antimalarial_inputs(eip_draws, tra_draws, n = n_draws)

itn <- read_itn_draws()
n_itn <- vapply(itn, function(x) length(unique(x$draw)), integer(1))
if (any(n_itn < n_draws)) {
  stop("The ITN files hold fewer than ", n_draws, " draws.", call. = FALSE)
}
draws$itn_draw <- seq_len(n_draws)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(draws, out_file, row.names = FALSE)
cat(sprintf("Wrote %d draws from fit %s (seed %d) to %s\n",
            n_draws, fit_stamp, seed, out_file))

# The table against the full posteriors and the priors it was drawn from
summ <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
fitted <- list(
  s_half_eip  = eip_draws[, , "s_half"],
  nH_eip      = eip_draws[, , "nH"],
  rho_frac    = eip_draws[, , "r0_frac"],
  s_half_pre  = tra_draws[, , "s_half_pre"],
  nH_pre      = tra_draws[, , "nH_pre"],
  B_max_post  = tra_draws[, , "b_max"],
  s_half_post = tra_draws[, , "s_half_post"],
  nH_post     = tra_draws[, , "nH_post"]
)
tab <- t(vapply(names(fitted), function(col) {
  c(summ(draws[[col]]), summ(fitted[[col]]))
}, numeric(6)))
colnames(tab) <- c("table 2.5%", "table 50%", "table 97.5%",
                   "fit 2.5%", "fit 50%", "fit 97.5%")
cat("\nPosterior parameters, table against full fit:\n")
print(signif(tab, 4))

pr <- projection_priors()
q <- c(0.025, 0.5, 0.975)
prior_tab <- rbind(
  omega_atn = c(summ(draws$omega_atn),
                stats::qbeta(q, pr$omega_atn$shape1, pr$omega_atn$shape2)),
  p_atn_atn = c(summ(draws$p_atn_atn),
                stats::qbeta(q, pr$p_atn$atn$shape1, pr$p_atn$atn$shape2)),
  p_atn_aitn = c(summ(draws$p_atn_aitn),
                 stats::qbeta(q, pr$p_atn$aitn$shape1, pr$p_atn$aitn$shape2)),
  half_life_years = c(summ(draws$half_life_years),
                      stats::qgamma(q, shape = pr$half_life_years$shape,
                                    scale = pr$half_life_years$scale))
)
colnames(prior_tab) <- c("table 2.5%", "table 50%", "table 97.5%",
                         "prior 2.5%", "prior 50%", "prior 97.5%")
cat("\nPrior inputs, table against prior:\n")
print(signif(prior_tab, 4))

cat("\nITN parameters, median over all draws:\n")
for (net in names(itn)) {
  x <- itn[[net]][itn[[net]]$resistance %in% c(0, 0.5, 0.9), ]
  med <- stats::aggregate(cbind(dn0, rn0, gamman) ~ resistance, x, stats::median)
  cat(net, "\n")
  print(signif(med, 4), row.names = FALSE)
}
