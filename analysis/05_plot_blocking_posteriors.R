# Posterior plots for the two blocking fits written by 03_fit_models.R (tra,
# tba_lab), and for field TBA, which is not fitted but obtained from the TRA
# posterior through Bompard's closed-form relationship, Challenger et al.
# 2023 JID 228:212-23 (SI p.4):
#
#   TBA(TRA; m, k) = [1 / (1 - (k/(k+m))^k)] *
#                    [ (k/(k+m*(1-TRA)))^k - (k/(k+m))^k ]
#
# with m = 0.000157, k = 0.00000495 (Challenger 2023 SI, Burkina Faso wild
# mosquitoes, Bompard 2020), the m_bompard and k_bompard defaults of the
# malariasimulation fork. Sensitivity to m: m/4, m, 2m (k held fixed).
#
# Ported from two scripts in malariasimple_ATNs/dev/:
#   Part A  plot_atn_main_bidir_splitnH.R               (TRA, lab TBA)
#   Part B  convert_tra_to_field_tba_bidir_splitnH_edit.R (field TBA)
# The plots are unchanged. What changed: the fits and priors are read from
# this package, the data overlay uses the mosquitoes the TRA fit was given,
# the two scripts' shared code appears once, and unused code is gone.
#
# Written to outputs/fits/<stamp>/figures/:
#   Part A, for <model> in tra, tba_lab:
#     <model>_draws.png, <model>_draws_with_data.png    100 posterior draws
#     <model>_cri.png, <model>_cri_with_data.png        median and 95% CrI
#     (the _cri_with_data plot also shows the prior 95% band, dashed)
#   Part B:
#     field_tba_cri.png, field_tba_cri_with_data.png    median and 95% CrI
#     field_tba_sensitivity.png, field_tba_sensitivity_with_data.png
#     blocking_three_curves.png                         TRA, field TBA, lab TBA
#     blocking_three_curves_with_data_4days.png         the same, +-4 days
#     blocking_three_panels_with_data.png               one panel each
#     blocking_three_panels_with_data_4days.png         the same, +-4 days
#     field_tba_summary.csv                             values at key times
#
# Data overlay, per experiment: TRA from a bootstrap of the control and
# treated oocyst counts; lab TBA from Beta-Jeffreys posteriors on the share
# of mosquitoes with any oocyst; field TBA by putting the TRA bootstrap
# through the Bompard map.
#
# Uses the newest full fit in outputs/fits/, or the one named:
#   Rscript analysis/05_plot_blocking_posteriors.R 20261009_1853
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(rstan)
  library(ggplot2)
  library(dplyr)
  library(ggokabeito)
  library(ggnewscale)
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
fig <- function(name) file.path(fig_dir, name)
cat("Fits from ", fit_dir, "\n", sep = "")

set.seed(2025)

# ---- Bompard parameters ----
m_central <- 0.000157
k_bompard <- 0.00000495
m_low     <- m_central / 4
m_high    <- m_central * 2

# ---- Curves grid (single uninterrupted grid; continuous at 0) ----
x_full <- seq(-10, 10, by = 0.02)

# ---- Model curves (Hill kernel, side-dependent s_half AND nH) ----
kernel_hill <- function(s_half, nH, s) s_half^nH / (s_half^nH + s^nH)

# b(s) for each draw. x = 0 is routed via the pre branch (kernel = 1
# either way), so the curve is continuous at the origin.
b_matrix <- function(b_max, s_half_pre, s_half_post, nH_pre, nH_post) {
  abs_x   <- abs(x_full)
  is_post <- x_full > 0
  b_all   <- matrix(NA_real_, nrow = length(b_max), ncol = length(x_full))
  for (k in seq_along(b_max)) {
    kern <- ifelse(is_post,
                   kernel_hill(s_half_post[k], nH_post[k], abs_x),
                   kernel_hill(s_half_pre[k],  nH_pre[k],  abs_x))
    b_all[k, ] <- b_max[k] * kern
  }
  b_all
}

compute_b_all <- function(fit) {
  post <- rstan::extract(fit, pars = c("r_min",
                                       "s_half_pre", "s_half_post",
                                       "nH_pre", "nH_post"))
  b_matrix(1 - post$r_min, post$s_half_pre, post$s_half_post,
           post$nH_pre, post$nH_post)
}

# ---- Numerically stable Bompard TRA -> TBA ----
tba_from_tra <- function(tra, m, k) {
  tra_safe <- pmin(pmax(tra, 0), 1 - 1e-12)
  log_a    <- k * (log(k) - log(k + m))
  log_b    <- k * (log(k) - log(k + m * (1 - tra_safe)))
  num      <- exp(log_a) * expm1(log_b - log_a)
  denom    <- -expm1(log_a)
  out      <- num / denom
  out[tra_safe <= 0] <- 0
  pmin(pmax(out, 0), 1)
}

convert_to_field_tba <- function(tra_mat, m, k) {
  matrix(tba_from_tra(as.vector(tra_mat), m, k),
         nrow = nrow(tra_mat), ncol = ncol(tra_mat))
}

# ---- Posterior summary helper (output in percent, 0-100) ----
summarise_b <- function(b_mat) {
  data.frame(
    x      = x_full,
    median = 100 * apply(b_mat, 2, median),
    lower  = 100 * apply(b_mat, 2, quantile, probs = 0.025),
    upper  = 100 * apply(b_mat, 2, quantile, probs = 0.975)
  )
}

quant_band <- function(b_mat) {
  data.frame(
    x     = x_full,
    lower = 100 * apply(b_mat, 2, quantile, probs = 0.025),
    upper = 100 * apply(b_mat, 2, quantile, probs = 0.975)
  )
}

# ---- Prior 95% band on b(s) (bidir split-nH kernel) ----------------------
# Draws (r_min, s_half_pre, s_half_post, nH_pre, nH_post) from atn_priors().
# log_nH_pre and log_nH_post are independent draws. OLRE / mu priors do not
# affect b(s) and are excluded. The same lab-scale prior on b(s) applies to
# the TRA and lab TBA fits; the field-scale band is its Bompard image.
pp <- sample_kernel_prior(10000, variant = "bidir_split", seed = 2025)
prior_b_lab    <- b_matrix(pp$b_max, pp$s_half_pre, pp$s_half_post,
                           pp$nH_pre, pp$nH_post)
prior_b_field  <- convert_to_field_tba(prior_b_lab, m_central, k_bompard)
df_prior_lab   <- quant_band(prior_b_lab)
df_prior_field <- quant_band(prior_b_field)

# ---- Load fits ----
fit_tra   <- readRDS(file.path(fit_dir, "tra.rds"))
b_tra_all <- compute_b_all(fit_tra)
cat(sprintf("TRA fit loaded: %d posterior draws.\n", nrow(b_tra_all)))

fit_tba_lab <- readRDS(file.path(fit_dir, "tba_lab.rds"))
b_tba_lab   <- compute_b_all(fit_tba_lab)
cat(sprintf("Lab TBA fit loaded: %d posterior draws.\n", nrow(b_tba_lab)))

b_tba_field      <- convert_to_field_tba(b_tra_all, m_central, k_bompard)
b_tba_field_low  <- convert_to_field_tba(b_tra_all, m_low,     k_bompard)
b_tba_field_high <- convert_to_field_tba(b_tra_all, m_high,    k_bompard)

df_tra       <- summarise_b(b_tra_all)
df_tba_lab   <- summarise_b(b_tba_lab)
df_tba_field <- summarise_b(b_tba_field)
df_tba_low   <- summarise_b(b_tba_field_low)
df_tba_high  <- summarise_b(b_tba_field_high)

# ---- Empirical overlays ---------------------------------------------------
# The mosquitoes the TRA fit was given, one row per mosquito. Exposure time
# relative to infection in days: positive after infection (side 1).
mos  <- readRDS(file.path(fit_dir, "tra_data.rds"))
x_of <- ifelse(mos$side == 1L, mos$s, -mos$s)

B   <- 4000
eps <- 1e-8

bootstrap_tra_row <- function(y_ctl, y_trt) {
  b_draws <- numeric(B)
  for (b in seq_len(B)) {
    yc <- sample(y_ctl, length(y_ctl), replace = TRUE)
    yt <- sample(y_trt, length(y_trt), replace = TRUE)
    b_draws[b] <- 1 - mean(yt) / max(mean(yc), eps)
  }
  b_draws
}

jeffreys_tba_row <- function(y_ctl, y_trt) {
  pos_c <- sum(y_ctl > 0); n_c <- length(y_ctl)
  pos_t <- sum(y_trt > 0); n_t <- length(y_trt)
  p_c   <- rbeta(B, pos_c + 0.5, n_c - pos_c + 0.5)
  p_t   <- rbeta(B, pos_t + 0.5, n_t - pos_t + 0.5)
  1 - p_t / pmax(p_c, eps)
}

field_tba_draws_row <- function(y_ctl, y_trt) {
  tba_from_tra(bootstrap_tra_row(y_ctl, y_trt), m_central, k_bompard)
}

build_overlay <- function(draw_fn) {
  bind_rows(lapply(seq_len(mos$N), function(i) {
    in_row <- mos$row_id == i
    y_ctl  <- mos$y[in_row & mos$arm == 0L]
    y_trt  <- mos$y[in_row & mos$arm == 1L]
    d      <- draw_fn(y_ctl, y_trt)
    data.frame(
      x     = x_of[i],
      b_est = max(median(d), 0),
      b_lo  = max(quantile(d, 0.025, names = FALSE), 0),
      b_hi  = min(quantile(d, 0.975, names = FALSE), 1)
    )
  }))
}

overlay_tra       <- build_overlay(bootstrap_tra_row)
overlay_field_tba <- build_overlay(field_tba_draws_row)
overlay_tba_lab   <- build_overlay(jeffreys_tba_row)

add_jitter <- function(df, jitter_w, seed) {
  set.seed(seed)
  df$x <- df$x + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
  df
}

# ---- Shared plot components ----
y_floor_pct <- -12.0

common_x <- labs(x = "Time of exposure relative to infection (days)")
common_y <- coord_cartesian(ylim = c(y_floor_pct, 100))

theme_transparent <- theme(
  plot.background  = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA)
)

point_colour <- "grey20"
prior_alpha  <- 0.3
prior_lw     <- 0.5

data_layers <- function(overlay_df) {
  list(
    geom_errorbar(
      data = overlay_df,
      aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi),
      colour = point_colour, width = 0.1, alpha = 0.5, linewidth = 0.4,
      inherit.aes = FALSE
    ),
    geom_point(
      data = overlay_df,
      aes(x = x, y = 100 * b_est),
      colour = point_colour, size = 2.2, alpha = 0.85, stroke = 0.8,
      inherit.aes = FALSE
    )
  )
}

prior_layers <- function(df_prior, colour) {
  list(
    geom_line(
      data = df_prior, aes(x = x, y = lower),
      colour = colour, linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    ),
    geom_line(
      data = df_prior, aes(x = x, y = upper),
      colour = colour, linetype = "dashed",
      linewidth = prior_lw, alpha = prior_alpha,
      inherit.aes = FALSE
    )
  )
}

# Arrows to +-10 days (Part A) or +-4 days (Part B).
make_ann_layers <- function(reach, label_x) {
  list(
    annotate("segment", x = 0.2, xend = reach, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = label_x, y = -10.0,
             label = "Exposure after infection", vjust = 0),
    annotate("segment", x = -0.2, xend = -reach, y = -5.0, yend = -5.0,
             arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
             colour = "black"),
    annotate("text",    x = -label_x, y = -10.0,
             label = "Exposure before infection", vjust = 0)
  )
}

# ==== Part A: TRA and lab TBA posteriors ==================================

ann_a <- make_ann_layers(reach = 10, label_x = 4.6)

plot_variant <- function(model, y_lab, b_all, overlay_df) {
  n_draws   <- nrow(b_all)
  b_all_pct <- 100 * b_all

  idx_sub <- sample(n_draws, min(100, n_draws))
  df_draws <- do.call(rbind, lapply(idx_sub, function(k) {
    data.frame(x = x_full, b = b_all_pct[k, ], draw = k)
  }))
  df_cri <- summarise_b(b_all)

  y_lab_layer <- labs(y = y_lab)

  p_draws_base <- ggplot(df_draws, aes(x = x, y = b, group = draw)) +
    geom_line(alpha = 0.18, colour = "dodgerblue") +
    ann_a +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_draws_data <- ggplot(df_draws, aes(x = x, y = b, group = draw)) +
    geom_line(alpha = 0.18, colour = "dodgerblue") +
    data_layers(overlay_df) +
    ann_a +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_cri_base <- ggplot(df_cri, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = "dodgerblue") +
    geom_line(aes(y = median), colour = "dodgerblue", linewidth = 0.8) +
    ann_a +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  p_cri_data <- ggplot(df_cri, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.3, fill = "dodgerblue") +
    geom_line(aes(y = median), colour = "dodgerblue", linewidth = 0.8) +
    prior_layers(df_prior_lab, "dodgerblue") +
    data_layers(overlay_df) +
    ann_a +
    common_x + y_lab_layer + common_y +
    theme_minimal() + theme_transparent

  ggsave(fig(paste0(model, "_draws.png")), p_draws_base,
         width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(fig(paste0(model, "_draws_with_data.png")), p_draws_data,
         width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(fig(paste0(model, "_cri.png")), p_cri_base,
         width = 8, height = 5, dpi = 450, bg = "transparent")
  ggsave(fig(paste0(model, "_cri_with_data.png")), p_cri_data,
         width = 8, height = 5, dpi = 450, bg = "transparent")
  cat(sprintf("Saved 4 plots for '%s'.\n", model))
}

plot_variant("tba_lab", "TBA (%)", b_tba_lab,
             add_jitter(overlay_tba_lab, jitter_w = 0.15, seed = 2025))
plot_variant("tra", "TRA (%)", b_tra_all,
             add_jitter(overlay_tra, jitter_w = 0.15, seed = 2025))

# ==== Part B: field TBA ===================================================

ann_b <- make_ann_layers(reach = 4, label_x = 2.1)

overlay_tra       <- add_jitter(overlay_tra,       jitter_w = 0.3, seed = 123)
overlay_field_tba <- add_jitter(overlay_field_tba, jitter_w = 0.3, seed = 123)
overlay_tba_lab   <- add_jitter(overlay_tba_lab,   jitter_w = 0.3, seed = 123)

# ---- Plot 1: field TBA posterior median + 95% CrI ----
field_TBA_col <- ggokabeito::palette_okabe_ito()[2]
p_field <- ggplot(df_tba_field, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.3, fill = field_TBA_col) +
  geom_line(aes(y = median), colour = field_TBA_col, linewidth = 0.8) +
  ann_b +
  common_x + labs(y = "Field TBA (%)") + common_y +
  theme_minimal() + theme_transparent

ggsave(fig("field_tba_cri.png"), p_field,
       width = 8, height = 5, dpi = 450, bg = "transparent")

# ---- Plot 2: sensitivity to m (m/4, m, 2m) ----
sens_disc <- c("slateblue2", "#56B4E9", "seagreen3")

df_sens <- bind_rows(
  df_tba_low   %>% mutate(label = sprintf("m = %.4g", m_low)),
  df_tba_field %>% mutate(label = sprintf("m = %.4g", m_central)),
  df_tba_high  %>% mutate(label = sprintf("m = %.4g", m_high))
)
df_sens$label <- factor(df_sens$label, levels = unique(df_sens$label))

p_sens <- ggplot(df_sens, aes(x = x, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  scale_colour_manual(values = sens_disc) +
  scale_fill_manual(values = sens_disc) +
  ann_b +
  common_x + labs(y = "Field TBA (%)", colour = NULL, fill = NULL) +
  common_y +
  theme_minimal() + theme_transparent +
  theme(legend.position = "bottom")

ggsave(fig("field_tba_sensitivity.png"), p_sens,
       width = 8, height = 5.5, dpi = 450, bg = "transparent")

# ---- Plot 3: three-curve comparison (legend mode) ----
df_compare <- bind_rows(
  df_tra       %>% mutate(curve = "Lab TRA"),
  df_tba_field %>% mutate(curve = "Field TBA"),
  df_tba_lab   %>% mutate(curve = "Lab TBA")
)
df_compare$curve <- factor(df_compare$curve, levels = unique(df_compare$curve))

p_compare <- ggplot(df_compare, aes(x = x, colour = curve, fill = curve)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.8) +
  ann_b +
  common_x + labs(y = "Blocking probability (%)",
                  colour = NULL, fill = NULL) +
  common_y +
  theme_minimal() + theme_transparent +
  scale_colour_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  scale_fill_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  theme(legend.position = "bottom")

ggsave(fig("blocking_three_curves.png"), p_compare,
       width = 8, height = 5.5, dpi = 450, bg = "white")

# ---- Plot 1b: field TBA CrI with overlay + prior band ----
p_field_data <- p_field + prior_layers(df_prior_field, "darkorange") +
  data_layers(overlay_field_tba)
ggsave(fig("field_tba_cri_with_data.png"), p_field_data,
       width = 8, height = 5, dpi = 450, bg = "transparent")

# ---- Plot 2b: sensitivity with central-m overlay (no prior band) ----
p_sens_data <- p_sens + data_layers(overlay_field_tba)
ggsave(fig("field_tba_sensitivity_with_data.png"), p_sens_data,
       width = 8, height = 5.5, dpi = 450, bg = "transparent")

# ---- Plot 3b: three-curve comparison faceted with overlays + prior bands ----
# Each panel shows one model curve, its model-matched empirical overlay,
# and the panel's prior 95% band as dashed lines (lab-scale prior for
# Lab TRA / Lab TBA, field-scale prior at central m for Field TBA).
# Panel y-axis labels (e.g. "Lab TRA (%)") are rendered via strip.position
# = "left" + a labeller; no shared y-axis title.
panel_levels <- c("Lab TRA", "Field TBA", "Lab TBA")
panel_labels <- c("Lab TRA"   = "Lab Transmission Reduction Activity (%)",
                  "Field TBA" = "Field Transmission Blocking Activity (%)",
                  "Lab TBA"   = "Lab Transmission Blocking Activity (%)")

as_panels <- function(blocks) {
  out <- bind_rows(Map(function(df, p) mutate(df, panel = p),
                       blocks, panel_levels))
  out$panel <- factor(out$panel, levels = panel_levels)
  out
}
df_facet_curve   <- as_panels(list(df_tra, df_tba_field, df_tba_lab))
df_facet_overlay <- as_panels(list(overlay_tra, overlay_field_tba,
                                   overlay_tba_lab))
df_facet_prior   <- as_panels(list(df_prior_lab, df_prior_field,
                                   df_prior_lab))

p_facet <- ggplot(df_facet_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = panel),
              alpha = 0.25, colour = NA) +
  geom_line(aes(y = median, colour = panel), linewidth = 0.8) +
  geom_line(
    data = df_facet_prior, aes(x = x, y = lower, colour = panel),
    linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_line(
    data = df_facet_prior, aes(x = x, y = upper, colour = panel),
    linetype = "dashed", linewidth = prior_lw, alpha = prior_alpha,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi),
    colour = point_colour, width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    aes(x = x, y = 100 * b_est),
    colour = point_colour, size = 2.0, alpha = 0.85, stroke = NA,
    inherit.aes = FALSE
  ) +
  scale_color_okabe_ito() +
  scale_fill_okabe_ito() +
  facet_wrap(~ panel, ncol = 1, strip.position = "left",
             labeller = labeller(panel = panel_labels)) +
  ann_b +
  common_x + labs(y = NULL) + common_y +
  theme_minimal() + theme_transparent +
  theme(legend.position   = "none",
        strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 90, face = "bold"))

ggsave(fig("blocking_three_panels_with_data.png"), p_facet,
       width = 8, height = 11, dpi = 450, bg = "white")

# ---- Plot 3c: as 3b over +-4 days, curve colours, no prior bands ----
p_facet_4days <- ggplot(df_facet_curve, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = panel),
              alpha = 0.25, colour = NA) +
  geom_line(aes(y = median, colour = panel), linewidth = 0.8) +
  scale_colour_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +
  scale_fill_manual(values = c("#C89B2B", "#EE7733", "#9A4D5E")) +

  ggnewscale::new_scale_colour() +

  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi, colour = panel),
    width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    aes(x = x, y = 100 * b_est, colour = panel),
    size = 2.0, alpha = 0.85, stroke = NA,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("#A77816", "#B85A26", "#7A2F3E")) +

  facet_wrap(~ panel, ncol = 1, strip.position = "left",
             labeller = labeller(panel = panel_labels)) +
  ann_b +
  common_x + labs(y = NULL) + common_y +
  theme_minimal() + theme_transparent +
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    limits = c(-4, 4)
  ) +
  theme(
    legend.position   = "none",
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold")
  )

ggsave(fig("blocking_three_panels_with_data_4days.png"), p_facet_4days,
       width = 8, height = 11, dpi = 450, bg = "white")

# ---- Plot 3d: three-curve comparison with overlays, +-4 days ----
p_compare_data <- p_compare +
  geom_errorbar(
    data = df_facet_overlay,
    aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi, colour = panel),
    width = 0.1, alpha = 0.5, linewidth = 0.4,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_point(
    data = df_facet_overlay,
    aes(x = x, y = 100 * b_est, colour = panel),
    size = 2.0, alpha = 0.85, stroke = NA,
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1),
    limits = c(-4, 4)
  ) +
  scale_colour_manual(values = c("#A77816", "#B85A26", "#7A2F3E"))

ggsave(fig("blocking_three_curves_with_data_4days.png"), p_compare_data,
       width = 8, height = 5, dpi = 450, bg = "white")

# ---- Summary table at key |s| values, both sides ----
key_s <- c(0, 0.5, 1, 2, 3, 5)

get_at_x <- function(b_mat, x_val) {
  idx  <- which.min(abs(x_full - x_val))
  vals <- b_mat[, idx]
  c(median = median(vals),
    lower  = quantile(vals, 0.025, names = FALSE),
    upper  = quantile(vals, 0.975, names = FALSE))
}

build_summary_side <- function(side_label, sign) {
  do.call(rbind, lapply(key_s, function(s_val) {
    x_val  <- sign * s_val
    q_tra  <- get_at_x(b_tra_all,        x_val)
    q_fld  <- get_at_x(b_tba_field,      x_val)
    q_lo   <- get_at_x(b_tba_field_low,  x_val)
    q_hi   <- get_at_x(b_tba_field_high, x_val)
    data.frame(
      side                 = side_label,
      s_days               = s_val,
      lab_TRA_med_pct      = round(100 * q_tra["median"], 2),
      lab_TRA_lo_pct       = round(100 * q_tra["lower"],  2),
      lab_TRA_hi_pct       = round(100 * q_tra["upper"],  2),
      field_TBA_med_pct    = round(100 * q_fld["median"], 2),
      field_TBA_lo_pct     = round(100 * q_fld["lower"],  2),
      field_TBA_hi_pct     = round(100 * q_fld["upper"],  2),
      field_TBA_low_m_med  = round(100 * q_lo["median"],  2),
      field_TBA_high_m_med = round(100 * q_hi["median"],  2)
    )
  }))
}

summary_tbl <- rbind(
  build_summary_side("pre",  -1),
  build_summary_side("post", +1)
)
rownames(summary_tbl) <- NULL

write.csv(summary_tbl, fig("field_tba_summary.csv"), row.names = FALSE)
cat("\nSummary table:\n")
print(summary_tbl, row.names = FALSE)
cat("\nDone. Figures in ", fig_dir, "\n", sep = "")
