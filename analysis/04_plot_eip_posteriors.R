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
# The data overlay gives, for each experiment and arm, the mean EIP implied
# by the observed sporozoite prevalence: a Beta-Jeffreys posterior on the
# prevalence, inverted through the cumulative-rate function of the model at
# its posterior medians; both ends clipped at the prior's EIP floor.
#
# Uses the newest full fit in outputs/fits/, or the one named:
#   Rscript analysis/04_plot_eip_posteriors.R 20261009_1853
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(rstan)
  library(lamW)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
fit_dir <- if (length(args)) {
  file.path("outputs", "fits", args[[1]])
} else {
  stamps <- sort(list.files(file.path("outputs", "fits"),
                            pattern = "^[0-9]{8}_[0-9]{4}$"))
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
d       <- 10                                        # Erlang shape
x_grid  <- seq(x_min, x_max, by = x_grid_step)

# ==== Model-specific r(x) computation =========================================

# Exponential-decay model: exact Lambert-W inversion of the cumulative
# integrated rate (matches the closed form used by eip_fit.stan).
compute_r_vec_exp <- function(rho, rho0, zeta, x, d = 10) {
  A        <- (rho - rho0) / zeta
  K        <- (rho - rho0) / rho
  baseline <- d / rho
  r        <- numeric(length(x))
  W_floor  <- -1 / exp(1) + 1e-12

  for (i in seq_along(x)) {
    if (x[i] >= 0) {
      d0 <- rho * x[i]
      if (d0 >= d) { r[i] <- baseline; next }
      r0 <- x[i]
      d1 <- d - d0
      B  <- (d1 + A) / rho
      r1 <- B + (1 / zeta) * lambertW0(max(-K * exp(-zeta * B), W_floor))
      r[i] <- r0 + r1
    } else {
      Ax <- A * exp(-zeta * abs(x[i]))
      Kx <- K * exp(-zeta * abs(x[i]))
      Bx <- (d + Ax) / rho
      r[i] <- Bx + (1 / zeta) * lambertW0(max(-Kx * exp(-zeta * Bx), W_floor))
    }
  }
  r[!is.finite(r) | r < baseline] <- baseline
  r
}

# Hill-decay model: numerical inversion (uniroot) of the cumulative
# integrated rate, since the Hill J(.) integral does not admit an
# elementary inverse for general non-integer nH.
hill_J <- function(x, n) {
  out <- numeric(length(x))
  pos <- is.finite(x) & x > 0
  if (any(pos)) {
    z   <- x[pos]^n / (1 + x[pos]^n)
    out[pos] <- pi / (n * sin(pi / n)) * pbeta(z, 1/n, 1 - 1/n)
  }
  out
}

# Cumulative rate S(tau; x) under the Hill model with exposure at x.
# Handles x >= 0 and x < 0 analogously to compute_r_vec_exp.
hill_S <- function(tau, x_exp, rho, rho0, s_half, nH) {
  B <- rho - rho0
  out <- rho * tau
  if (x_exp >= 0) {
    if (tau > x_exp) {
      out <- out - B * s_half * hill_J((tau - x_exp) / s_half, nH)
    }
  } else {
    ax <- -x_exp
    out <- out - B * s_half * (
      hill_J((tau + ax) / s_half, nH) - hill_J(ax / s_half, nH)
    )
  }
  out
}

compute_r_vec_hill <- function(rho, rho0, s_half, nH, x, d = 10) {
  baseline <- d / rho
  r <- numeric(length(x))
  for (i in seq_along(x)) {
    # Quick exit: exposure happens after the EIP would have completed at
    # baseline rate (only for x >= 0). EIP = baseline (no extension).
    if (x[i] >= 0 && rho * x[i] >= d) {
      r[i] <- baseline
      next
    }
    # Numerical inversion of S(tau; x[i]) = d.
    f <- function(tau) hill_S(tau, x[i], rho, rho0, s_half, nH) - d
    upper <- max(baseline, abs(x[i])) * 5 + 10
    sol <- tryCatch(
      uniroot(f, interval = c(1e-3, upper), extendInt = "upX")$root,
      error = function(e) NA_real_
    )
    r[i] <- if (is.finite(sol) && sol >= baseline) sol else baseline
  }
  r
}

# Wrapper: matrix of r values across a list of par draws ----------------
compute_r_matrix_exp <- function(pars, x_grid, d) {
  n <- length(pars$rho)
  out <- matrix(NA_real_, n, length(x_grid))
  for (k in seq_len(n)) {
    out[k, ] <- compute_r_vec_exp(pars$rho[k], pars$rho0[k], pars$zeta[k],
                                  x_grid, d = d)
  }
  out
}
compute_r_matrix_hill <- function(pars, x_grid, d) {
  n <- length(pars$rho)
  out <- matrix(NA_real_, n, length(x_grid))
  for (k in seq_len(n)) {
    out[k, ] <- compute_r_vec_hill(pars$rho[k], pars$rho0[k],
                                   pars$s_half[k], pars$nH[k],
                                   x_grid, d = d)
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
      sapply(tau, function(t) {
        if (t <= delta) rho * t
        else rho * t - B * s_half * hill_J((t - delta) / s_half, nH)
      })
    }
  } else {
    ax <- -delta
    J_lo <- hill_J(ax / s_half, nH)   # constant; recovery already accumulated
                                      # before infection (tau = 0)
    function(tau) {
      sapply(tau, function(t) {
        if (t <= 0) rho * t
        else rho * t - B * s_half *
                       (hill_J((t + ax) / s_half, nH) - J_lo)
      })
    }
  }
}

# Control arm: rate is rho throughout (no exposure). Closed form,
# identical for both models.
implied_tau_ctl <- function(p, tau_obs, rho, d = 10) {
  if (!is.finite(p)) return(NA_real_)
  p <- pmax(pmin(p, 1 - 1e-12), 1e-12)
  s_hat <- qgamma(p, shape = d, rate = 1)
  pmax(tau_obs + (d - s_hat) / rho, 0.1)
}

# Treatment arm: numerical inversion of G. Model dispatch via the spec.
implied_tau_trt_exp <- function(p, delta_val, tau_obs, med, d = 10) {
  if (!is.finite(p)) return(NA_real_)
  p <- max(min(p, 1 - 1e-12), 1e-12)
  s_hat  <- qgamma(p, shape = d, rate = 1)
  G      <- make_G_exp(delta_val, med$rho, med$rho0, med$zeta)
  target <- G(tau_obs) + (d - s_hat)
  if (target <= 0) return(0.1)
  f <- function(tau) G(tau) - target
  res <- tryCatch(uniroot(f, interval = c(0.01, 200), extendInt = "upX")$root,
                  error = function(e) NA_real_)
  max(res, 0.1)
}
implied_tau_trt_hill <- function(p, delta_val, tau_obs, med, d = 10) {
  if (!is.finite(p)) return(NA_real_)
  p <- max(min(p, 1 - 1e-12), 1e-12)
  s_hat  <- qgamma(p, shape = d, rate = 1)
  G      <- make_G_hill(delta_val, med$rho, med$rho0, med$s_half, med$nH)
  target <- G(tau_obs) + (d - s_hat)
  if (target <= 0) return(0.1)
  f <- function(tau) G(tau) - target
  res <- tryCatch(uniroot(f, interval = c(0.01, 200), extendInt = "upX")$root,
                  error = function(e) NA_real_)
  max(res, 0.1)
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

  n_draws <- length(post$rho)

  # ---- Posterior r(x) curves ----------------------------------------------
  cat("[", spec$label, "] computing r(x) over ", n_draws,
      " posterior draws and ", length(x_grid), " x values...\n", sep = "")
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
