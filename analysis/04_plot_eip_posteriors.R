# Posterior plots for the two EIP fits written by 03_fit_models.R:
#   (i)  eip_exponential: exponential decay in the rate deficit
#   (ii) eip_hill:        Hill decay in the rate deficit (used by the
#                         simulations)
#
# Ported from malariasimple_ATNs/dev/eip_posterior_plots_080626.R. The plots
# are unchanged. What changed: the fits, priors and data are read from this
# package's outputs, the data overlay uses the data each fit was given, and
# code whose output was never plotted is gone (the prior pushforward bands,
# commented out of every panel, and the exploratory analysis at the end).
#
# For each model, written to outputs/fits/<stamp>/figures/:
#   <model>_draws.png, <model>_draws_with_data.png
#       mean EIP against the time of exposure relative to infection: 100
#       posterior draws (orange), with the baseline EIP of the same draws
#       (grey)
#   <model>_cri.png, <model>_cri_with_data.png
#       the same as posterior median and 95% CrI
#   <model>_relrate_draws.png, <model>_relrate_cri.png
#       relative EIP progression rate against time since ATN exposure
# and for both models together:
#   eip_combined_draws.png, eip_combined_cri.png
#       rate (left) and mean EIP with data (right), exponential above Hill
#
# "Mean EIP" throughout is E[T], the mean of the Erlang-threshold completion
# time the model implies -- see the section on it below. Before 2026-10-15
# these panels plotted the time at which a mosquito progressing at the average
# rate would finish, which is a different quantity and produced a spurious
# cliff at late exposure times.
#
# The data overlay gives, for each experiment and arm, the mean EIP implied
# by the observed sporozoite prevalence: a Beta-Jeffreys posterior on the
# prevalence, converted to a mean EIP by rescaling the speed of the model's
# fitted rate profile, at its posterior medians, until the model reproduces
# that prevalence at the dissection time; both ends clipped at the prior's
# EIP floor.
#
# Uses the newest full fit in outputs/fits/, or the one named:
#   Rscript analysis/04_plot_eip_posteriors.R 20261009_1853
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(rstan)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
fit_dir <- if (length(args)) {
  file.path("outputs", "fits", args[[1]])
} else {
  stamps <- sort(list.files(file.path("outputs", "fits"),
                            pattern = "^[0-9]{8}_[0-9]{4}(_eipR[0-9]+)?$"))
  if (!length(stamps)) {
    stop("No full fit in outputs/fits/. Run analysis/03_fit_models.R first.",
         call. = FALSE)
  }
  file.path("outputs", "fits", stamps[length(stamps)])
}
fig_dir <- file.path(fit_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat("Fits from ", fit_dir, "\n", sep = "")

priors <- atn_priors()

# ==== Editable plot configuration =============================================

# Edit the values below to change the plotted x range, the x and y axis
# tick / gridline placement, and (implicitly) the placement of the
# Exposure-before-infection / Exposure-after-infection arrow annotations
# on every output panel.  Both the exp and Hill models inherit this
# configuration.
#
# Arrow / label logic (driven by x_min, x_max via make_ann_layers):
#   - Each arrow's tip lands exactly at the corresponding x-axis bound
#     (left arrow tip at x_min, right arrow tip at x_max), with a small
#     symmetric inset from x = 0 at the tail.
#   - The label belonging to the *shorter* arrow is centred on the
#     midpoint of that arrow; the *longer* arrow's label is mirrored to
#     the same |x|, so the two labels sit at symmetric x-coordinates
#     about x = 0.
#       * |x_min| < |x_max|: "Exposure before infection" centred on the
#         left arrow; "Exposure after infection" placed at the mirrored
#         positive x.
#       * |x_min| > |x_max|: "Exposure after infection" centred on the
#         right arrow; "Exposure before infection" placed at the
#         mirrored negative x.
#       * |x_min| = |x_max|: each label centred on its own arrow.
#
# x_grid_step controls the resolution at which r(x) is evaluated for
# the posterior curves; finer = smoother, slower.
x_min          <- -12
x_max          <-  12
x_grid_step    <-   0.1

# How many posterior draws the mean-EIP curves are computed over. See the note
# in process_model: n_eff is around 3000, so this is thinning redundancy, not
# discarding information. Inf uses every draw and takes roughly four times as
# long.
n_draws_curve  <- 2000L

x_major_breaks <- seq(-14, 14, by = 2)
x_minor_breaks <- seq(-14, 14, by = 1)
y_major_breaks <- seq(0, 50, by = 5)
y_minor_breaks <- seq(0, 50, by = 1)

# ==== Editable rate-plot configuration =======================================
# Controls the "relative EIP progression rate" panels: rho(s)/rho vs
# time s since ATN exposure. The plot is independent of the r(x) plots
# above (different x and y semantics) so it has its own configuration.
#
# For s < 0 the curve is fixed at 1 (pre-exposure baseline). At s = 0 it
# drops vertically from 1 to r0_frac (impaired-rate fraction) and recovers
# back to 1 as s -> infinity, with the kernel shape determined by the
# fitted recovery parameters. The vertical drop at s = 0 is drawn by the
# duplicate-zero trick: s_grid contains 0 twice adjacently, with the first
# 0 carrying y = 1 and the second carrying y = rho(0+)/rho = r0_frac, so
# geom_line connects them with a vertical segment automatically.
rate_s_min          <- -2
rate_s_max          <-  12
rate_s_step         <-   0.05

rate_x_major_breaks <- seq(-2, 12, by = 2)
rate_x_minor_breaks <- seq(-2, 12, by = 1)
rate_y_major_breaks <- seq(0, 1, by = 0.2)
rate_y_minor_breaks <- seq(0, 1, by = 0.1)
rate_y_lim          <- c(0, 1.05)   # ylim leaves a sliver of headroom above 1

# ---- Universal constants (shared by both models) ----
# The Erlang shape (number of EIP stages) is the one the fits were given, read from
# the data each was fitted to rather than assumed: it was 10 up to the fits of
# 9 October 2026 and 40 from 14 October 2026 (ATN_ERLANG_SHAPE in 03_fit_models.R).
# Both EIP models in a folder must share it; process_model() checks each in turn.
d_files <- file.path(fit_dir, paste0(c("eip_hill", "eip_exponential"), "_data.rds"))
d_files <- d_files[file.exists(d_files)]
if (!length(d_files)) {
  stop("No EIP fit data in ", fit_dir, ". Run analysis/03_fit_models.R first.",
       call. = FALSE)
}
d       <- readRDS(d_files[[1]])$Delta_r             # Erlang shape
cat("Erlang shape (EIP stages): ", d, "\n", sep = "")
x_grid  <- seq(x_min, x_max, by = x_grid_step)

# ==== Mean EIP under a time-varying progression rate ==========================

# A mosquito completes the EIP when its accumulated progression G(tau) reaches
# an Erlang(R, 1) threshold (SI eq:trt_erlang_prob; the same statement the Stan
# likelihood evaluates). The threshold is random, so the EIP is a random
# variable, and its mean is
#
#     E[T] = int_0^inf Pr(Gamma(R, 1) > G(t)) dt.
#
# That is the quantity plotted here, and it is the mean of the distribution the
# fitted model implies.
#
# It is NOT the tau solving G(tau) = R. That is the completion time of a
# mosquito progressing at exactly the average rate: the hitting time of the mean
# rather than the mean hitting time. The two coincide only as R -> infinity,
# where the EIP becomes deterministic. At finite R they separate once exposure
# is late enough that some mosquitoes have finished and others have not, because
# the G(tau) = R quantity drops discontinuously to the baseline as soon as
# exposure follows the average completion time -- its single notional mosquito
# has already finished -- whereas at R = 30 some 47% of mosquitoes are still
# mid-sporogony at day 10 and a dose then still delays them.
#
# Until 2026-10-15 these panels plotted the G(tau) = R quantity, which produced
# a near-vertical drop to the baseline at x = the baseline EIP. The fits were
# never affected: the Stan likelihood works with completion probabilities and
# never inverts anything. The two definitions agree to within 0.01 d for
# exposure up to ~6 days after infection, and diverge by at most ~0.7 d beyond
# that. SI section eip_lengthening still states the superseded claim that the
# R -> infinity expressions "describe the mean lengthened-EIP durations" at
# finite R, and its two mean-EIP tables are built on it.
#
# With no suppression G(t) = rho * t and the integral returns R / rho exactly,
# which is baseline_per_draw below.

# Gauss-Legendre nodes and weights mapped to [0, 1], via Golub-Welsch. The same
# scheme eip_fit_hill.stan uses, so the R and Stan integrals agree by
# construction rather than by coincidence.
gauss_legendre_01 <- function(m) {
  i <- seq_len(m - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = (e$vectors[1, o])^2)
}
GL16 <- gauss_legendre_01(16L)
GL32 <- gauss_legendre_01(32L)

# E[T] for a cumulative-progression function G. `rho` is the baseline (maximal)
# rate, used only to locate a safe upper limit. `brk_extra` gives the interior
# times at which the integrand is not analytic or turns sharply; the callers
# supply the moment of exposure (where G' jumps from rho to rho0), the centre of
# the recovery kernel, and the time by which the kernel is spent.
#
# Gauss-Legendre is spectrally accurate on an analytic integrand and merely
# ordinary across a kink, so the panel edges must line up with the places G is
# not smooth. For exposure at x >= 0 the rate is the baseline until t = x and
# suppressed after, and for late exposure that kink falls right inside the
# Erlang transition where the integrand carries all its mass: integrating
# straight through it cost about 1e-3 days. Breaking there, and again at the
# kernel's own bend (which matters only for a near-step kernel, large nH with a
# short s_half), is what makes 16 nodes per panel enough: the worst error over
# 300 draws from the widened posterior is then 3.6e-5 days, about 3 seconds of
# EIP. Checked against adaptive quadrature and a fine Simpson rule; 32 nodes
# would buy another three orders of magnitude and cost twice the runtime, which
# is not a trade worth making for a figure.
mean_eip_from_G <- function(G, rho, d, brk_extra = numeric(0)) {
  # Integrate out to where the survival term has fallen below about 1e-13.
  M    <- d + 10 * sqrt(d)
  t_hi <- M / rho
  # G' <= rho everywhere, so this iteration approaches the root from below and
  # contracts. It replaces a root-find at a fraction of the cost.
  for (k in 1:8) t_hi <- t_hi + (M - G(t_hi)) / rho
  guard <- 0L
  while (G(t_hi) < M && guard < 40L) { t_hi <- t_hi * 1.5; guard <- guard + 1L }
  # The integrand carries all its mass in the Erlang transition, around G(t) = d
  # with width sqrt(d). Beyond the kernel G is near-linear with slope rho, so the
  # times bracketing that transition follow from the t_hi solve for free, and
  # splitting there is what lets 16 nodes per panel suffice.
  trans <- t_hi - (M - d + c(5, -5) * sqrt(d)) / rho
  brk <- sort(unique(pmin(pmax(c(0, brk_extra, trans, t_hi), 0), t_hi)))
  tot <- 0
  for (i in seq_len(length(brk) - 1L)) {
    a <- brk[i]; b <- brk[i + 1L]
    if (b <= a) next
    tt  <- a + (b - a) * GL16$x
    tot <- tot + (b - a) * sum(GL16$w * pgamma(G(tt), shape = d, rate = 1,
                                               lower.tail = FALSE))
  }
  tot
}

# J(x; n) = int_0^x du / (1 + u^n).
#
# This replaces the regularised-incomplete-beta closed form used until
# 2026-10-15, which carries a prefactor pi / (n sin(pi / n)): that form is valid
# only for n > 1 and is singular at n = 1 exactly, so it cannot be used now the
# prior lower bound on nH is 1 rather than exp(0.5).
#
# Two pieces, chosen so both integrands are analytic and bounded:
#   u in [0, min(x, 1)]  the integrand runs between 1/2 and 1; nothing to resolve.
#   u in [1, x]          substituted u = exp(s), so the (1 + u^n)^-1 decay becomes
#                        a clean exponential exp(-(n - 1) s) in s. That piece is
#                        split again where the integrand has died, so the nodes
#                        are not spent on a numerically zero tail at large n.
# eip_fit_hill.stan splits [0, 1] at v* = 1/x instead, which is accurate over the
# small arguments the likelihood sees (x = (tau - Delta)/s_half, at most ~9) but
# loses digits over the much larger arguments these plots reach.
hill_J <- function(x, n, gl = GL32) {
  out <- numeric(length(x))
  ok  <- is.finite(x) & x > 0
  if (!any(ok)) return(out)
  xx  <- x[ok]
  u1  <- pmin(xx, 1)
  res <- u1 * as.vector((1 / (1 + outer(u1, gl$x)^n)) %*% gl$w)
  big <- xx > 1
  if (any(big)) {
    L  <- log(xx[big])
    # exp(s) / (1 + exp(n s)), written so neither term can overflow.
    f  <- function(s) exp((1 - n) * s) / (1 + exp(-n * s))
    q  <- function(a, b) (b - a) * as.vector(f(outer(b - a, gl$x) + a) %*% gl$w)
    sc <- if (n > 1.0001) pmin(L, 40 / (n - 1)) else L
    res[big] <- res[big] + q(0, sc) + q(sc, L)
  }
  out[ok] <- res
  out
}

# Wrapper: matrix of mean-EIP values across a list of par draws ----------------
# Both models go through the same mean_eip_from_G; they differ only in the
# cumulative-progression function, which make_G_exp / make_G_hill supply.
compute_r_matrix_exp <- function(pars, x_grid, d) {
  n   <- length(pars$rho)
  out <- matrix(NA_real_, n, length(x_grid))
  for (k in seq_len(n)) {
    # Exposure, one e-folding time, and spent after about ten.
    off <- c(0, 1, 10) / pars$zeta[k]
    for (j in seq_along(x_grid)) {
      G <- make_G_exp(x_grid[j], pars$rho[k], pars$rho0[k], pars$zeta[k])
      out[k, j] <- mean_eip_from_G(G, pars$rho[k], d,
                                   max(0, x_grid[j]) + off)
    }
  }
  out
}
compute_r_matrix_hill <- function(pars, x_grid, d) {
  n   <- length(pars$rho)
  out <- matrix(NA_real_, n, length(x_grid))
  for (k in seq_len(n)) {
    # Exposure, the half-recovery time, and spent after about six of them.
    off <- c(0, 1, 6) * pars$s_half[k]
    for (j in seq_along(x_grid)) {
      G <- make_G_hill(x_grid[j], pars$rho[k], pars$rho0[k],
                       pars$s_half[k], pars$nH[k])
      out[k, j] <- mean_eip_from_G(G, pars$rho[k], d,
                                   max(0, x_grid[j]) + off)
    }
  }
  out
}

# ==== Relative EIP progression rate at time s since ATN exposure ============
# rho(s) / rho = 1 - (1 - r0_frac) * sigma(s), where sigma is the model's
# suppression kernel: exponential exp(-zeta s) or Hill s_half^n / (s_half^n
# + s^n). Both functions take a vector s_post of *non-negative* s values
# and return a draw-by-row matrix; the s < 0 (pre-exposure) portion is
# trivially 1 and is concatenated outside.
compute_rel_rate_matrix_exp <- function(pars, s_post) {
  n <- length(pars$rho)
  out <- matrix(NA_real_, n, length(s_post))
  r0_post <- pars$rho0 / pars$rho
  for (k in seq_len(n)) {
    out[k, ] <- 1 - (1 - r0_post[k]) * exp(-pars$zeta[k] * s_post)
  }
  out
}
compute_rel_rate_matrix_hill <- function(pars, s_post) {
  n <- length(pars$rho)
  out <- matrix(NA_real_, n, length(s_post))
  r0_post <- pars$rho0 / pars$rho
  for (k in seq_len(n)) {
    nh    <- pars$nH[k]
    sh_n  <- pars$s_half[k]^nh
    out[k, ] <- 1 - (1 - r0_post[k]) * sh_n / (sh_n + s_post^nh)
  }
  out
}

# ==== Empirical-overlay machinery (model-aware) ===============================

# Cumulative rate function G(tau) under the fitted model with exposure
# at delta. Used to invert observed prevalence into implied tau*.
# Both make_G_* functions branch on the sign of delta and match the
# corresponding Stan likelihood exactly:
#   delta >= 0  -> existing closed form (suppression integral starts at
#                  exposure = delta, accumulates to tau).
#   delta <  0  -> exposure preceded infection by |delta|; the
#                  suppression integral from infection (tau = 0) to tau
#                  is the integrated suppression on (|delta|, tau+|delta|),
#                  i.e. integrated over u >= |delta| only.
make_G_exp <- function(delta, rho, rho0, zeta) {
  A <- (rho - rho0) / zeta
  if (delta >= 0) {
    function(tau) {
      ifelse(tau <= delta,
             rho * tau,
             rho * tau - A * (1 - exp(-zeta * (tau - delta))))
    }
  } else {
    function(tau) {
      ifelse(tau <= 0,
             rho * tau,
             rho * tau - A * exp(zeta * delta) * (1 - exp(-zeta * tau)))
    }
  }
}
make_G_hill <- function(delta, rho, rho0, s_half, nH) {
  B <- rho - rho0
  if (delta >= 0) {
    function(tau) {
      out <- rho * tau
      i <- tau > delta
      if (any(i)) {
        out[i] <- out[i] - B * s_half * hill_J((tau[i] - delta) / s_half, nH)
      }
      out
    }
  } else {
    ax <- -delta
    J_lo <- hill_J(ax / s_half, nH)   # constant; recovery already accumulated
                                      # before infection (tau = 0)
    function(tau) {
      out <- rho * tau
      i <- tau > 0
      if (any(i)) {
        out[i] <- out[i] - B * s_half *
                           (hill_J((tau[i] + ax) / s_half, nH) - J_lo)
      }
      out
    }
  }
}

# ---- Turning an observed prevalence into a mean EIP --------------------------
#
# A dissection gives one number: the share of mosquitoes sporozoite-positive at
# tau_obs. Under the model that share is F_Gamma(G(tau_obs); R, 1), so the
# observed prevalence p names the accumulated progression the assay actually
# saw, s_hat = qgamma(p; R, 1), against the G(tau_obs) the fitted model expects.
#
# The overlay reports the mean EIP of the mosquitoes that produced p. We take
# the fitted model's rate PROFILE -- its full time-varying suppression, not a
# constant rate -- and rescale its speed by the single factor
#
#     k = s_hat / G(tau_obs)
#
# so that the rescaled model reproduces the observed prevalence exactly at the
# dissection time. The mean EIP is then E[T] under k * G. When the data match
# the model, k = 1 and the overlay point sits on the curve.
#
# Until 2026-10-15 these functions instead inverted G(tau) = R deterministically
# after an ADDITIVE shift, G(tau) = G(tau_obs) + (R - s_hat). That cannot carry
# over to E[T]: an additive shift leaves the accumulated progression non-zero
# (indeed negative) at the moment of infection, so there is no coherent
# trajectory to take the expectation over. Rescaling the rate keeps G(0) = 0 and
# keeps the suppression shape, and for the control arm it reduces to the exact
# closed form below.
#
# Control arm: no exposure, so the rate really is constant at rho and the
# observed prevalence implies the constant rate s_hat / tau_obs directly. The
# mean EIP is then R divided by that rate, exactly -- no quadrature needed,
# because the EIP is exactly Erlang.
implied_tau_ctl <- function(p, tau_obs, rho, d) {
  if (!is.finite(p)) return(NA_real_)
  p <- pmax(pmin(p, 1 - 1e-12), 1e-12)
  s_hat <- qgamma(p, shape = d, rate = 1)
  pmax(d * tau_obs / s_hat, 0.1)
}

# Treatment arm: rescale the fitted rate profile to match the observation, then
# take E[T] under it. Model dispatch via the spec.
implied_tau_scaled <- function(G, rho, tau_obs, p, d, brk_extra) {
  if (!is.finite(p)) return(NA_real_)
  p     <- max(min(p, 1 - 1e-12), 1e-12)
  s_hat <- qgamma(p, shape = d, rate = 1)
  g_obs <- G(tau_obs)
  if (!is.finite(g_obs) || g_obs <= 0) return(NA_real_)
  k  <- s_hat / g_obs
  Gk <- function(tau) k * G(tau)
  max(mean_eip_from_G(Gk, k * rho, d, brk_extra), 0.1)
}
implied_tau_trt_exp <- function(p, delta_val, tau_obs, med, d) {
  G <- make_G_exp(delta_val, med$rho, med$rho0, med$zeta)
  implied_tau_scaled(G, med$rho, tau_obs, p, d,
                     max(0, delta_val) + c(0, 1, 10) / med$zeta)
}
implied_tau_trt_hill <- function(p, delta_val, tau_obs, med, d) {
  G <- make_G_hill(delta_val, med$rho, med$rho0, med$s_half, med$nH)
  implied_tau_scaled(G, med$rho, tau_obs, p, d,
                     max(0, delta_val) + c(0, 1, 6) * med$s_half)
}

# ==== Model specs =============================================================

model_spec_exp <- list(
  model     = "eip_exponential",
  par_names = c("rho", "rho0", "zeta"),
  compute_r_matrix        = compute_r_matrix_exp,
  compute_rel_rate_matrix = compute_rel_rate_matrix_exp,
  implied_tau_trt_fn      = implied_tau_trt_exp,
  median_pars             = function(post) list(rho = median(post$rho),
                                                rho0 = median(post$rho0),
                                                zeta = median(post$zeta)),
  label     = "exponential decay"
)

model_spec_hill <- list(
  model     = "eip_hill",
  par_names = c("rho", "rho0", "s_half", "nH"),
  compute_r_matrix        = compute_r_matrix_hill,
  compute_rel_rate_matrix = compute_rel_rate_matrix_hill,
  implied_tau_trt_fn      = implied_tau_trt_hill,
  median_pars             = function(post) list(rho    = median(post$rho),
                                                rho0   = median(post$rho0),
                                                s_half = median(post$s_half),
                                                nH     = median(post$nH)),
  label     = "Hill decay"
)

# ==== Shared plotting bits ====================================================

# Annotation arrows + labels (driven by x_min, x_max from the config block).
make_ann_layers <- function(y_top, xmin = x_min, xmax = x_max) {
  arrow_inset <- 0.2
  arrow_head  <- arrow(length = unit(0.18, "inches"), type = "closed")

  left_start  <- -arrow_inset
  left_end    <-  xmin
  right_start <-  arrow_inset
  right_end   <-  xmax

  left_mid  <- (left_start  + left_end ) / 2
  right_mid <- (right_start + right_end) / 2

  if (abs(xmin) < abs(xmax)) {
    before_x <- left_mid
    after_x  <- -before_x
  } else if (abs(xmin) > abs(xmax)) {
    after_x  <- right_mid
    before_x <- -after_x
  } else {
    before_x <- left_mid
    after_x  <- right_mid
  }

  list(
    annotate("segment",
             x = right_start, xend = right_end,
             y = y_top,       yend = y_top,
             arrow = arrow_head, colour = "black"),
    annotate("text", x = after_x,  y = y_top,
             label = "Exposure after infection",  vjust = -0.4),
    annotate("segment",
             x = left_start,  xend = left_end,
             y = y_top,       yend = y_top,
             arrow = arrow_head, colour = "black"),
    annotate("text", x = before_x, y = y_top,
             label = "Exposure before infection", vjust = -0.4)
  )
}

posterior_colour <- "#B85A26"   # orange: EIP under ATN exposure
baseline_colour  <- "#7A7A7A"   # grey: baseline EIP

common_labs <- labs(
  x = "Time of exposure relative to infection (days)",
  y = "Mean EIP (days)"
)
common_scales <- list(
  scale_x_continuous(breaks       = x_major_breaks,
                     minor_breaks = x_minor_breaks),
  scale_y_continuous(breaks       = y_major_breaks,
                     minor_breaks = y_minor_breaks)
)
common_coord <- coord_cartesian(xlim = c(x_min, x_max), clip = "off")

# ==== Per-model processing function ==========================================

process_model <- function(spec, n_draws_to_plot = 100,
                          jitter_w = 0.5, dodge_w = 1.5) {

  fit_path <- file.path(fit_dir, paste0(spec$model, ".rds"))
  if (!file.exists(fit_path)) {
    cat("[", spec$label, "] fit file not found at '", fit_path,
        "' -- skipping.\n", sep = "")
    return(invisible(NULL))
  }

  cat("[", spec$label, "] loading fit from ", fit_path, "\n", sep = "")

  # Local seed so re-running this block (or the other model) is
  # deterministic.
  set.seed(2025)

  # ---- Load posterior draws and the data the model was fitted to -----------
  fit  <- readRDS(fit_path)
  post <- rstan::extract(fit, pars = spec$par_names)
  med  <- spec$median_pars(post)
  dat  <- readRDS(file.path(fit_dir, paste0(spec$model, "_data.rds")))
  stopifnot(dat$Delta_r == d)

  # ---- Thin the draws used for the curves ---------------------------------
  # E[T] costs a quadrature per (draw, x) where the superseded deterministic
  # quantity cost a root-find, so all 8000 draws would take minutes per model
  # rather than seconds. n_eff for the kernel parameters is about 3000, so the
  # chain never carried 8000 independent pieces of information in the first
  # place: a random n_draws_curve of them gives a visually identical band for a
  # quarter of the work. Set n_draws_curve to Inf to use every draw.
  n_draws_all <- length(post$rho)
  idx_curve   <- if (n_draws_curve < n_draws_all) {
    sort(sample(n_draws_all, n_draws_curve))
  } else seq_len(n_draws_all)
  post    <- lapply(post, function(v) v[idx_curve])
  n_draws <- length(idx_curve)

  # ---- Posterior mean-EIP curves ------------------------------------------
  cat("[", spec$label, "] computing mean EIP over ", n_draws, " of ",
      n_draws_all, " posterior draws and ", length(x_grid), " x values...\n",
      sep = "")
  r_all <- spec$compute_r_matrix(post, x_grid, d)

  idx_sub  <- sample(n_draws, n_draws_to_plot)
  df_draws <- do.call(rbind, lapply(idx_sub, function(k) {
    data.frame(x = x_grid, r = r_all[k, ], draw = k)
  }))
  df_cri <- data.frame(
    x      = x_grid,
    median = apply(r_all, 2, median),
    lower  = apply(r_all, 2, quantile, probs = 0.025),
    upper  = apply(r_all, 2, quantile, probs = 0.975)
  )

  # ---- Baseline-EIP per posterior draw ------------------------------------
  # Same thinned draws as the curves, so the grey and orange bands are built
  # from one set of mosquitoes rather than two.
  baseline_per_draw <- d / post$rho

  df_draws_baseline <- do.call(rbind, lapply(idx_sub, function(k) {
    data.frame(x = x_grid, r = baseline_per_draw[k], draw = k)
  }))
  df_baseline_cri <- data.frame(
    x      = x_grid,
    median = median(baseline_per_draw),
    lower  = unname(quantile(baseline_per_draw, 0.025)),
    upper  = unname(quantile(baseline_per_draw, 0.975))
  )

  # ---- Empirical EIP overlay (uses model medians) -------------------------
  days <- sort(unique(dat$tau_obs))
  make_row <- function(pos, n, tau_obs, delta_val, arm_label,
                       eip_floor = priors$eip_floor) {
    p_hat <- (pos + 0.5) / (n + 1)
    p_lo  <- qbeta(0.025, pos + 0.5, n - pos + 0.5)
    p_hi  <- qbeta(0.975, pos + 0.5, n - pos + 0.5)

    if (arm_label == "Control") {
      est <- implied_tau_ctl(p_hat, tau_obs, med$rho, d)
      lo  <- implied_tau_ctl(p_hi,  tau_obs, med$rho, d)
      hi  <- implied_tau_ctl(p_lo,  tau_obs, med$rho, d)
    } else {
      est <- spec$implied_tau_trt_fn(p_hat, delta_val, tau_obs, med, d)
      lo  <- spec$implied_tau_trt_fn(p_hi,  delta_val, tau_obs, med, d)
      hi  <- spec$implied_tau_trt_fn(p_lo,  delta_val, tau_obs, med, d)
    }

    est <- pmax(est, eip_floor)
    lo  <- pmax(lo,  eip_floor)
    hi  <- pmax(hi,  eip_floor)

    data.frame(x = delta_val, arm = arm_label,
               dpi = factor(tau_obs, levels = days),
               mean_eip = est, eip_lo = lo, eip_hi = hi)
  }

  data_overlay <- bind_rows(
    lapply(seq_len(dat$N), function(i) {
      bind_rows(
        make_row(dat$pos_ctl[i], dat$n_ctl[i], dat$tau_obs[i], dat$delta[i],
                 "Control"),
        make_row(dat$pos_trt[i], dat$n_trt[i], dat$tau_obs[i], dat$delta[i],
                 "200 mg")
      )
    })
  )

  # Jitter (matches 05_plot_blocking_posteriors.R convention).
  add_jitter <- function(df, seed = 2025) {
    set.seed(seed)
    df$x <- df$x + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
    df
  }
  data_overlay <- add_jitter(data_overlay)

  # ---- Layer lists --------------------------------------------------------
  data_layers <- list(
    geom_errorbar(
      data = data_overlay,
      aes(x = x, ymin = eip_lo, ymax = eip_hi,
          colour = arm, group = interaction(arm, dpi)),
      width = 0.25, alpha = 0.4, linewidth = 0.4,
      position = position_dodge(width = dodge_w),
      inherit.aes = FALSE
    ),
    geom_point(
      data = data_overlay,
      aes(x = x, y = mean_eip, colour = arm, shape = dpi,
          group = interaction(arm, dpi)),
      size = 2.2, alpha = 0.65, stroke = 0.8,
      position = position_dodge(width = dodge_w),
      inherit.aes = FALSE
    ),
    scale_colour_manual(values = c("Control" = baseline_colour,
                                   "200 mg" = posterior_colour),
                        labels = c("Control" = "Control",
                                   "200 mg" = "ATN")),
    scale_shape_manual(values = c(16, 17, 15, 18)[seq_along(days)]),
    labs(colour = "Arm", shape = "Days post infection")
  )

  baseline_cri_layers <- list(
    geom_ribbon(
      data        = df_baseline_cri,
      aes(x = x, ymin = lower, ymax = upper),
      fill        = baseline_colour, alpha = 0.18,
      inherit.aes = FALSE
    ),
    geom_line(
      data        = df_baseline_cri,
      aes(x = x, y = median),
      colour = baseline_colour, linewidth = 0.8,
      inherit.aes = FALSE
    )
  )
  baseline_draws_layers <- list(
    geom_line(
      data        = df_draws_baseline,
      aes(x = x, y = r, group = draw),
      colour = baseline_colour, alpha = 0.15,
      inherit.aes = FALSE
    )
  )

  y_top_base <- max(c(df_cri$upper, df_baseline_cri$upper),
                    na.rm = TRUE) * 1.05
  y_top_data <- max(c(df_cri$upper, df_baseline_cri$upper,
                      data_overlay$eip_hi), na.rm = TRUE) * 1.05

  # ---- Four plots --------------------------------------------------------
  p_draws_base <- ggplot(df_draws, aes(x = x, y = r, group = draw)) +
    baseline_draws_layers +
    geom_line(alpha = 0.18, colour = posterior_colour) +
    make_ann_layers(y_top_base) +
    common_labs +
    common_scales +
    common_coord +
    theme_minimal()

  p_draws_data <- ggplot(df_draws, aes(x = x, y = r, group = draw)) +
    baseline_draws_layers +
    geom_line(alpha = 0.18, colour = posterior_colour) +
    data_layers +
    make_ann_layers(y_top_data) +
    common_labs +
    common_scales +
    common_coord +
    theme_minimal()

  p_cri_base <- ggplot(df_cri, aes(x = x)) +
    baseline_cri_layers +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = posterior_colour) +
    geom_line(aes(y = median), colour = posterior_colour, linewidth = 0.8) +
    make_ann_layers(y_top_base) +
    common_labs +
    common_scales +
    common_coord +
    theme_minimal()

  p_cri_data <- ggplot(df_cri, aes(x = x)) +
    baseline_cri_layers +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = posterior_colour) +
    geom_line(aes(y = median), colour = posterior_colour, linewidth = 0.8) +
    data_layers +
    make_ann_layers(y_top_data) +
    common_labs +
    common_scales +
    common_coord +
    theme_minimal()

  # ---- Save with the model name as prefix ---------------------------------
  out_path <- function(type) {
    file.path(fig_dir, paste0(spec$model, "_", type, ".png"))
  }

  ggsave(out_path("draws"),           p_draws_base,
         width = 8, height = 5, dpi = 200, bg = "white")
  ggsave(out_path("draws_with_data"), p_draws_data,
         width = 8, height = 5, dpi = 200, bg = "white")
  ggsave(out_path("cri"),             p_cri_base,
         width = 8, height = 5, dpi = 200, bg = "white")
  ggsave(out_path("cri_with_data"),   p_cri_data,
         width = 8, height = 5, dpi = 200, bg = "white")

  # ---- Relative EIP progression rate panels ------------------------------
  # rho(s) / rho vs time s since ATN exposure. The s_grid contains 0
  # twice in adjacent positions: the first 0 carries y = 1 (baseline,
  # pre-exposure) and the second carries y = r0_frac (impaired rate at
  # s = 0+); geom_line then draws the vertical drop at s = 0
  # automatically.
  s_grid_pre  <- seq(rate_s_min, 0, by = rate_s_step)   # ends at 0
  s_grid_post <- seq(0, rate_s_max, by = rate_s_step)   # starts at 0
  s_grid_rate <- c(s_grid_pre, s_grid_post)
  n_pre       <- length(s_grid_pre)

  rate_post_mat <- cbind(
    matrix(1, n_draws, n_pre),
    spec$compute_rel_rate_matrix(post, s_grid_post)
  )

  df_rate_draws <- do.call(rbind, lapply(idx_sub, function(k) {
    data.frame(s = s_grid_rate, rate = rate_post_mat[k, ], draw = k)
  }))
  df_rate_cri <- data.frame(
    s      = s_grid_rate,
    median = apply(rate_post_mat, 2, median),
    lower  = apply(rate_post_mat, 2, quantile, probs = 0.025),
    upper  = apply(rate_post_mat, 2, quantile, probs = 0.975)
  )

  rate_scales <- list(
    scale_x_continuous(breaks       = rate_x_major_breaks,
                       minor_breaks = rate_x_minor_breaks),
    scale_y_continuous(breaks       = rate_y_major_breaks,
                       minor_breaks = rate_y_minor_breaks)
  )
  rate_coord <- coord_cartesian(xlim = c(rate_s_min, rate_s_max),
                                ylim = rate_y_lim,
                                clip = "off")
  rate_labs <- labs(
    x = "Time since ATN exposure (days)",
    y = "Relative EIP progression rate"
  )

  p_rate_draws <- ggplot(df_rate_draws, aes(x = s, y = rate, group = draw)) +
    geom_line(alpha = 0.18, colour = posterior_colour) +
    rate_labs +
    rate_scales +
    rate_coord +
    theme_minimal()

  p_rate_cri <- ggplot(df_rate_cri, aes(x = s)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = posterior_colour) +
    geom_line(aes(y = median), colour = posterior_colour, linewidth = 0.8) +
    rate_labs +
    rate_scales +
    rate_coord +
    theme_minimal()

  ggsave(out_path("relrate_draws"), p_rate_draws,
         width = 8, height = 5, dpi = 200, bg = "white")
  ggsave(out_path("relrate_cri"),   p_rate_cri,
         width = 8, height = 5, dpi = 200, bg = "white")

  cat("[", spec$label, "] saved ", spec$model, "_{draws,draws_with_data,",
      "cri,cri_with_data,relrate_draws,relrate_cri}.png\n", sep = "")

  invisible(list(p_draws_data = p_draws_data, p_cri_data = p_cri_data,
                 p_rate_draws = p_rate_draws, p_rate_cri = p_rate_cri))
}

# ==== Run both models =========================================================

list_exp  <- process_model(model_spec_exp)
list_hill <- process_model(model_spec_hill)

patch_layout <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 3),
  area(t = 2, l = 1, b = 2, r = 1),
  area(t = 2, l = 2, b = 2, r = 3)
)

if (!is.null(list_exp) && !is.null(list_hill)) {
  p_draws_combined <-
    list_exp$p_rate_draws +
    list_exp$p_draws_data +
    list_hill$p_rate_draws +
    list_hill$p_draws_data +
    plot_layout(design = patch_layout)

  p_cri_combined <-
    list_exp$p_rate_cri +
    list_exp$p_cri_data +
    list_hill$p_rate_cri +
    list_hill$p_cri_data +
    plot_layout(design = patch_layout)

  ggsave(file.path(fig_dir, "eip_combined_draws.png"), p_draws_combined,
         width = 12, height = 8, dpi = 450, bg = "white")
  ggsave(file.path(fig_dir, "eip_combined_cri.png"),   p_cri_combined,
         width = 12, height = 8, dpi = 450, bg = "white")
}

cat("\nDone. Figures in ", fig_dir, "\n", sep = "")
