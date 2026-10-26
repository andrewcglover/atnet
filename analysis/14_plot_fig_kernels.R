# Two-panel kernel figure for the main paper.
#
#   A  the three blocking curves on one pair of axes, distinguished by colour:
#      lab TRA, field TBA, lab TBA. Lab TRA and lab TBA carry their empirical
#      overlays; field TBA deliberately does not, because it is not fitted to
#      anything -- it is the lab TRA posterior pushed through the Bompard map,
#      so points drawn under it would be the TRA points a second time.
#   B  mean EIP against time of exposure, treated and untreated netting, with
#      the empirical overlay and the spread of individual mosquitoes.
#
# Both panels share one x axis -- the same quantity, time of exposure relative
# to infection -- so the axis is drawn once, under B, and one pair of
# "Exposure before / after infection" arrows sits along the foot of B and
# serves both.
#
# The axis runs -4 to +16 days, asymmetric on purpose, and panel A is cut off
# at +6: beyond there its curves are flat on zero and say nothing, so plotting
# and tick-marking them would only assert emptiness at length. What the cut
# buys is the whole right-hand side of A for the figure's single legend, and
# the shape of the argument stays visible -- read down the figure and blocking
# is spent within a day of infection while EIP lengthening is still unwinding
# at +16. Blocking matters for exposure BEFORE infection, EIP lengthening for
# exposure after it.
#
# Two versions of panel B are written. The plain one shows the mean EIP and its
# credible interval. The "_mosqrange" one adds dashed lines at the central 95%
# of INDIVIDUAL mosquitoes -- posterior-predictive percentiles, which carry
# parameter uncertainty inside them, so there is no second band to be mistaken
# for the first. Readers have taken the CrI on the mean to be the range for
# mosquitoes; these lines are what they were looking for, and at x = 0 they are
# 7.4 times wider (8.55-15.76 d against 11.37-12.35 d).
#
# The relative-EIP-progression-rate panel is NOT here: it says the same thing
# as B in units a reader has to be taught, so it belongs in the supplement.
#
# Nothing is recomputed. Both panels are drawn from the plot-ready frames that
# 04_plot_eip_posteriors.R and 05_plot_blocking_posteriors.R write beside their
# own figures, so this figure cannot drift from the panels those scripts
# publish, and a figure does not cost a fresh pass of quadrature over 2000
# posterior draws. Re-run those two scripts first if the fits change.
#
# The two fit folders are the pair 06_make_draws.R samples the projection
# draws table from; keep them in step with it.
#
#   Rscript analysis/14_plot_fig_kernels.R [x_lo x_hi]
#
# Run with the working directory at the root of the package.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
})

eip_stamp <- "20261015_1306_eipR30"   # EIP kernels, Hill model, 30 Erlang stages
tra_stamp <- "20261009_1853"          # blocking kernels

x_lo_default <- -6
x_hi_default <- 14
a_x_max      <- 6     # panel A is drawn and tick-marked only this far

# Panel B's vertical frame. Only 2 of the 48 empirical intervals reach past 16
# days, and letting them set the scale wasted the top third of the panel on two
# of its least informative points, so they are capped and the top tick reads
# ">16". The floor is 4 rather than the lowest interval (5.25) to clear a band
# for the exposure arrows; it is tick-marked like the rest of the axis, since an
# axis that stopped labelling at 6 while the panel ran to 4 would read as a
# truncation artefact.
b_y_lo     <- 4
b_y_cap    <- 16
b_y_breaks <- seq(4, 16, by = 2)
b_arrow_y  <- 4.5     # between the 4 and 5 gridlines, clear of every interval

ribbon_alpha <- 0.3
pt_size      <- 2.0   # panel A's size, now used in both
bar_width    <- 0.22

out_dir <- file.path("outputs", "figures", "kernels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

blk_file <- file.path("outputs", "fits", tra_stamp, "figures",
                      "blocking_plotdata.rds")
eip_file <- file.path("outputs", "fits", eip_stamp, "figures",
                      "eip_hill_plotdata.rds")
for (f in c(blk_file, eip_file)) {
  if (!file.exists(f)) {
    stop("missing ", f, "\nRun analysis/04_plot_eip_posteriors.R and ",
         "analysis/05_plot_blocking_posteriors.R first.", call. = FALSE)
  }
}
blk <- readRDS(blk_file)
eip <- readRDS(eip_file)
has_range <- !is.null(eip$pred_range)

# ---- Colours ---------------------------------------------------------------
# Paul Tol's "medium contrast" scheme, which ships as dark-line / light-fill
# PAIRS of one hue -- built for this structure, so three shaded credible
# intervals can overlap and each still reads as its own curve.
#
# Deliberately NOT the projection palette (10_plot_fig5.R: Tol "vibrant",
# #BBBBBB none, #EE7733 Pyr-only, #CC3311 Pyr-PBO, #EE3377 Pyr-CFP, #009988
# ATN, #0077BB Pyr-ATN, #33BBEE Pyr-CFP-ATN). Every colour there already means
# a net type, and these kernels are properties of the DRUG -- they apply to
# every A(I)TN -- so an arm colour would tell a reader who has learned fig 5
# something false.
#
# Blue is the quantity the transmission model actually consumes: field TBA in
# panel A, and the treated arm in panel B. One colour, one meaning, across the
# figure. The two warm hues are the laboratory endpoints behind it.
# Paul Tol's "dark" scheme (sronpersonalpages.nl/~pault/, 11 named colours).
#
# Two orders, deliberately different. draw_order decides what sits on top --
# field TBA last, because it is the curve the transmission model consumes and
# the one a reader should be able to trace unbroken where the three cross.
# legend_order is what the key reads, which follows the vertical order of the
# curves on their rise rather than the drawing order.
draw_order   <- c("Lab TRA", "Lab TBA", "Field TBA")
legend_order <- c("Lab TRA", "Field TBA", "Lab TBA")
curve_labels <- c("Lab TRA"   = "Laboratory transmission-reducing activity",
                  "Field TBA" = "Field transmission-blocking activity",
                  "Lab TBA"   = "Laboratory transmission-blocking activity")
line_cols    <- c("Lab TRA"   = "#2255AA",   # sapphire
                  "Field TBA" = "#006644",   # pine
                  "Lab TBA"   = "#993366")   # maroon
fill_cols    <- line_cols   # one hue per curve; the ribbon alpha lightens it

atn_colour      <- "#117788"   # teal
atn_fill        <- "#117788"
baseline_colour <- "#444444"   # grey
baseline_fill   <- "#444444"
net_line <- c("200 mg" = atn_colour,  "Control" = baseline_colour)
net_fill <- c("200 mg" = atn_fill,    "Control" = baseline_fill)

# Levels as stored in the fit data, and what the reader should see instead.
# 200 mg/m2 is the concentration every treated arm was built from
# (build_fitting_tables.R); the compound is ELQ-453.
arm_breaks <- c("200 mg", "Control")
arm_labels <- expression("ELQ-453 200 mg m"^-2, "Untreated")

# Dissection day. One shape scale spans both panels: panel A's oocyst assays
# were all read at day 7 and panel B's sporozoite assays at 10 and 13. The day
# 7 is from the assay protocol, not from the data -- the oocyst table carries
# no dissection-day column, because every row of it is day 7.
dpi_levels <- c("7", "10", "13")
dpi_shapes <- c("7" = 15, "10" = 16, "13" = 17)

# Dashed reference lines: ONE linetype and ONE linewidth for both arms. Dash
# length in grid scales with linewidth, so an otherwise identical "22" drawn
# thinner comes out looking dotted rather than dashed.
range_lty <- "22"
range_lwd <- 0.55

# ---- Panel A data ----------------------------------------------------------
curves <- blk$curves
curves$panel <- factor(as.character(curves$panel), levels = draw_order)

# Field TBA is a transformation of lab TRA, not an independent measurement.
overlay_a <- blk$overlay
overlay_a <- overlay_a[as.character(overlay_a$panel) != "Field TBA", ]
overlay_a$panel <- factor(as.character(overlay_a$panel), levels = draw_order)
overlay_a$dpi   <- factor("7", levels = dpi_levels)

# ---- Panel B data ----------------------------------------------------------
eip_cri      <- eip$cri
eip_baseline <- eip$baseline_cri
overlay_b    <- eip$overlay
overlay_b$arm <- factor(as.character(overlay_b$arm), levels = arm_breaks)
overlay_b$eip_hi <- pmin(overlay_b$eip_hi, b_y_cap)
overlay_b$dpi <- factor(as.character(overlay_b$dpi), levels = dpi_levels)

# ==== Figure ================================================================

# Arrows along the foot of panel B. They label the shared axis, so they belong
# to the figure rather than to either panel. Each arrow starts a fixed step
# either side of zero and its caption sits a fixed distance out, so the two
# read as a matched pair even though the right arrow runs four times further.
arrow_layers <- function(x_lo, x_hi, y, gap = 0.4, label_x = 3,
                         text_size = 2.9) {
  head <- arrow(length = unit(0.11, "inches"), type = "closed")
  list(
    annotate("segment", x =  gap, xend = x_hi, y = y, yend = y,
             arrow = head, colour = "black", linewidth = 0.45),
    annotate("segment", x = -gap, xend = x_lo, y = y, yend = y,
             arrow = head, colour = "black", linewidth = 0.45),
    annotate("text", x =  label_x, y = y, vjust = -0.55, size = text_size,
             label = "Exposure after infection"),
    annotate("text", x = -label_x, y = y, vjust = -0.55, size = text_size,
             label = "Exposure before infection")
  )
}

build_figure <- function(x_lo, x_hi, show_range = FALSE) {

  trim_a <- function(df) df[df$x >= x_lo & df$x <= a_x_max, , drop = FALSE]
  grid_end <- a_x_max + 0.05 * (a_x_max - x_lo)   # ggplot's own 5% expansion
  trim_b <- function(df) df[df$x >= x_lo & df$x <= x_hi,   , drop = FALSE]

  # ---- A -------------------------------------------------------------------
  # The legend lives in the space the +6 cut frees. Its three blocks are built
  # here, including the two that describe panel B, so the figure carries one
  # legend rather than three scattered across two panels. The panel-B blocks
  # come from dummy layers parked below the y range and clipped away: that way
  # their keys inherit the real ribbon and line styling instead of needing
  # override.aes, so "Netting" reads as a band with a line through it exactly
  # as "Activity" does.
  dummy_net <- data.frame(
    x    = rep(c(x_lo, x_lo + 0.01), times = 2),
    y    = -20, ymin = -21, ymax = -19,
    net  = factor(rep(arm_breaks, each = 2), levels = arm_breaks))
  dummy_dpi <- data.frame(
    x = x_lo, y = -20,
    dpi = factor(dpi_levels, levels = dpi_levels))

  p_a <- ggplot(trim_a(curves), aes(x = x)) +
    # Horizontal gridlines drawn by hand so they stop at the cut instead of
    # running under the legend. Shortening the coord is not an option:
    # patchwork aligns panel AREAS, so A's +6 would sit above B's +16 and the
    # shared axis would break. Added first, so they stay behind everything.
    annotate("segment", x = -Inf, xend = grid_end,
             y = seq(0, 100, 25), yend = seq(0, 100, 25),
             colour = "grey92", linewidth = 0.5) +
    annotate("segment", x = -Inf, xend = grid_end,
             y = seq(12.5, 87.5, 25), yend = seq(12.5, 87.5, 25),
             colour = "grey92", linewidth = 0.25) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = panel),
                alpha = ribbon_alpha, colour = NA) +
    geom_line(aes(y = median, colour = panel), linewidth = 0.8) +
    geom_errorbar(
      data = trim_a(overlay_a),
      aes(x = x, ymin = 100 * b_lo, ymax = 100 * b_hi, colour = panel),
      width = bar_width, alpha = 0.5, linewidth = 0.4,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_point(
      data = trim_a(overlay_a),
      aes(x = x, y = 100 * b_est, colour = panel, shape = dpi),
      size = pt_size, alpha = 0.85, stroke = NA,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    # Breaks and labels must match on both scales or the legend splits in two.
    scale_colour_manual(values = line_cols, breaks = legend_order,
                        labels = curve_labels, name = "Activity",
                        guide = guide_legend(order = 1)) +
    scale_fill_manual(values = fill_cols, breaks = legend_order,
                      labels = curve_labels, name = "Activity",
                      guide = guide_legend(order = 1)) +
    new_scale_colour() + new_scale_fill() +
    geom_ribbon(data = dummy_net,
                aes(x = x, ymin = ymin, ymax = ymax, fill = net),
                alpha = ribbon_alpha, colour = NA, inherit.aes = FALSE) +
    geom_line(data = dummy_net, aes(x = x, y = y, colour = net),
              linewidth = 0.8, inherit.aes = FALSE) +
    scale_colour_manual(values = net_line, breaks = arm_breaks,
                        labels = arm_labels, name = "Netting",
                        guide = guide_legend(order = 2)) +
    scale_fill_manual(values = net_fill, breaks = arm_breaks,
                      labels = arm_labels, name = "Netting",
                      guide = guide_legend(order = 2)) +
    geom_point(data = dummy_dpi, aes(x = x, y = y, shape = dpi),
               colour = "grey30", size = pt_size, alpha = 0.85,
               inherit.aes = FALSE) +
    scale_shape_manual(values = dpi_shapes, breaks = dpi_levels,
                       drop = FALSE, name = "Dissection day",
                       guide = guide_legend(order = 3)) +
    scale_x_continuous(breaks       = seq(x_lo, a_x_max, by = 2),
                       minor_breaks = seq(x_lo, a_x_max, by = 1)) +
    scale_y_continuous(breaks       = seq(0, 100, 25),
                       minor_breaks = seq(0, 100, 12.5)) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0, 100), clip = "on") +
    # Two of the three curves are reducing activity, not blocking, so the old
    # "Blocking probability (%)" was wrong for them; the legend now names each.
    labs(x = NULL, y = "Activity (%)") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x            = element_blank(),
      axis.title.x           = element_blank(),
      panel.grid.major.y     = element_blank(),
      panel.grid.minor.y     = element_blank(),
      plot.margin            = margin(t = 6, r = 6, b = 2, l = 6),
      legend.position        = "inside",
      legend.position.inside = c(0.995, 0.5),
      legend.justification.inside = c(1, 0.5),
      legend.background      = element_blank(),
      legend.key             = element_blank(),
      legend.title           = element_text(size = 10, face = "bold"),
      legend.text            = element_text(size = 9),
      legend.margin          = margin(0, 0, 0, 0),
      legend.spacing.y       = unit(3, "pt"),
      legend.key.size        = unit(11, "pt")
    )

  # ---- B -------------------------------------------------------------------
  draw_range <- show_range && has_range
  range_atn  <- if (draw_range) trim_b(eip$pred_range)     else NULL
  range_base <- if (draw_range) trim_b(eip$baseline_range) else NULL

  p_b <- ggplot(trim_b(eip_cri), aes(x = x)) +
    geom_ribbon(data = trim_b(eip_baseline),
                aes(x = x, ymin = lower, ymax = upper),
                fill = baseline_fill, alpha = ribbon_alpha,
                inherit.aes = FALSE) +
    geom_line(data = trim_b(eip_baseline), aes(x = x, y = median),
              colour = baseline_colour, linewidth = 0.8, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = atn_fill, alpha = ribbon_alpha) +
    geom_line(aes(y = median), colour = atn_colour, linewidth = 0.8)

  if (draw_range) {
    p_b <- p_b +
      geom_line(data = range_base, aes(x = x, y = value, group = bound),
                colour = baseline_colour, linetype = range_lty,
                linewidth = range_lwd, inherit.aes = FALSE) +
      geom_line(data = range_atn, aes(x = x, y = value, group = bound),
                colour = atn_colour, linetype = range_lty,
                linewidth = range_lwd, inherit.aes = FALSE)
  }

  p_b <- p_b +
    geom_errorbar(
      data = trim_b(overlay_b),
      aes(x = x, ymin = eip_lo, ymax = eip_hi, colour = arm),
      width = bar_width, alpha = 0.4, linewidth = 0.4, inherit.aes = FALSE
    ) +
    geom_point(
      data = trim_b(overlay_b),
      aes(x = x, y = mean_eip, colour = arm, shape = dpi),
      size = pt_size, alpha = 0.65, stroke = 0.8, inherit.aes = FALSE
    ) +
    arrow_layers(x_lo, x_hi, y = b_arrow_y) +
    scale_colour_manual(values = net_line, guide = "none") +
    scale_shape_manual(values = dpi_shapes, drop = FALSE,
                       guide = "none") +
    scale_x_continuous(breaks       = seq(x_lo, x_hi, by = 2),
                       minor_breaks = seq(x_lo, x_hi, by = 1)) +
    scale_y_continuous(breaks       = b_y_breaks,
                       labels       = c(head(b_y_breaks, -1),
                                        paste0(">", b_y_cap)),
                       minor_breaks = seq(b_y_lo, b_y_cap, by = 1)) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(b_y_lo, b_y_cap),
                    clip = "on") +
    labs(x = "Time of exposure relative to infection (days)",
         y = "EIP (days)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.margin     = margin(t = 4, r = 6, b = 4, l = 6))

  (p_a / p_b) +
    plot_layout(heights = c(0.45, 0.55)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 13))
}

# ==== Write =================================================================
args <- commandArgs(trailingOnly = TRUE)
rng  <- if (length(args) >= 2) as.numeric(args[1:2]) else
  c(x_lo_default, x_hi_default)

variants <- list(list(suffix = "", range = FALSE))
if (has_range) {
  variants <- c(variants, list(list(suffix = "_mosqrange", range = TRUE)))
}

for (v in variants) {
  p   <- build_figure(rng[1], rng[2], show_range = v$range)
  out <- file.path(out_dir, sprintf("fig_kernels_x%g_%g%s.png",
                                    rng[1], rng[2], v$suffix))
  ggsave(out, p, width = 8, height = 8, dpi = 400, bg = "white")
  cat("wrote ", out, "\n", sep = "")
}

cat(sprintf("\nPoint jitter, both panels: uniform +-%.2f days.\n", 0.3 / 2))
cat(sprintf("Empirical intervals capped at %g d: %d of %d.\n",
            b_y_cap, sum(eip$overlay$eip_hi > b_y_cap), nrow(eip$overlay)))
