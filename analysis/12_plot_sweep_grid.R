# The two-way sensitivity grid: clinical cases averted at each
# combination of the antimalarial's half-life and the nets' retention.
#
# Four figures, one per comparison, each a grid of 20 panels -- net
# retention across, half-life down. In every panel the three A(I)TN arms
# are a point and a 95% uncertainty interval, and the three ITNs are a
# dashed line at their median over a faint band holding their interval.
# The ITNs carry no antimalarial, so their backdrop does not change down
# a column; it changes across one, retention being theirs too.
#
#   sweep_grid_vs_nonets_<form>   averted against no future nets
#   sweep_grid_vs_pyr_<form>      additional, against Pyr
#   sweep_grid_vs_pyr_pbo_<form>  additional, against Pyr-PBO
#   sweep_grid_vs_pyr_cfp_<form>  additional, against Pyr-CFP
#
# National values only, and everything is paired within a draws-table
# row before any median is taken. Every panel is held to the rows its
# own A(I)TN cells hold, so the ITN lines are drawn from the same rows
# as the points they are compared with.
#
# Two things to keep in the caption. The cells hold 160 rows where the
# main run holds 640, so these are the noisier figures of the two. And
# no cell is the main result: the main run draws the half-life instead
# of fixing it, so it lies inside this grid rather than on it.
#
# Reads outputs/summaries/projections_v2/sweep_values.rds, written by
# analysis/11_summarise_sweeps.R, and runs no simulations.
#   Rscript analysis/12_plot_sweep_grid.R                 smooth-compact
#   Rscript analysis/12_plot_sweep_grid.R exp             exponential
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
form <- if (length(args)) args[1] else "smoothcompact"
sweeps <- file.path("outputs", "summaries", "projections_v2",
                    "sweep_values.rds")
if (!file.exists(sweeps)) {
  stop(sweeps, " does not exist. Run analysis/11_summarise_sweeps.R.",
       call. = FALSE)
}
out_dir <- file.path("outputs", "figures", "sweeps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

arms <- projection_arms()
arms$colour <- c("#BBBBBB", "#EE7733", "#CC3311", "#EE3377", "#009988",
                 "#0077BB", "#33BBEE")
atn_arms <- c("atn", "pyr_atn", "pyr_cfp_atn")
itn_arms <- c("pyr", "pyr_pbo", "pyr_cfp")
arm_colours <- setNames(arms$colour, arms$arm)
arm_labels <- setNames(arms$label, arms$arm)

all <- readRDS(sweeps)
all <- all[all$region == "National", ]
all$retention_years <- round(all$retention_days / 365, 1)

# The grid's own cells: a fixed half-life, the arms' own mechanism, and
# nothing else held away from its prior. That last pair of conditions is
# what keeps the min and max runs out, which fix the half-life as well
# but fix the contact probability and p_A with it.
cells <- all[all$decay_form == form & !is.na(all$hl_years) &
               all$mechanism == "both" & is.na(all$omega_q) &
               is.na(all$p_a_q) & all$arm %in% atn_arms, ]
if (!nrow(cells)) {
  stop("No grid cells for decay form ", form, ".", call. = FALSE)
}

# The comparators come from the seven-arm runs, which hold no
# antimalarial and so appear once per retention whatever the grid cell.
itns <- all[all$arm %in% itn_arms, ]

quantities <- c(vs_nonets = "averted_per1000", vs_pyr = "vs_pyr",
                vs_pyr_pbo = "vs_pyr_pbo", vs_pyr_cfp = "vs_pyr_cfp")
against <- c(vs_nonets = "no future nets", vs_pyr = "Pyr",
             vs_pyr_pbo = "Pyr-PBO", vs_pyr_cfp = "Pyr-CFP")

# Median and 95% interval of one column, by whatever else is given.
band <- function(d, column, by) {
  d$value <- d[[column]]
  d %>%
    group_by(across(all_of(by))) %>%
    summarise(median = median(value),
              lower = quantile(value, 0.025, names = FALSE),
              upper = quantile(value, 0.975, names = FALSE),
              .groups = "drop")
}

grid_figure <- function(key) {
  column <- quantities[[key]]

  # A contrast is missing only when the seven-arm run at that retention
  # has not been summarised; once 11 has run over everything, none is.
  have <- cells[!is.na(cells[[column]]), ]
  if (nrow(have) < nrow(cells)) {
    message("  ", nrow(cells) - nrow(have), " of ", nrow(cells),
            " cell lines have no comparator at their retention and are ",
            "left out")
  }
  if (!nrow(have)) stop("No cell has a comparator for ", key, ".",
                        call. = FALSE)

  # Each panel is held to the rows its own cells hold, so a line and the
  # point beside it are drawn from the same rows.
  rows <- have %>%
    group_by(hl_years, retention_years) %>%
    summarise(draws = list(sort(unique(draw))), .groups = "drop")

  points <- band(have, column, c("hl_years", "retention_years", "arm"))
  points$arm <- factor(points$arm, levels = atn_arms)

  lines <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
    d <- itns[itns$retention_years == rows$retention_years[i] &
                itns$draw %in% rows$draws[[i]], ]
    d <- d[!is.na(d[[column]]), ]
    if (!nrow(d)) return(NULL)
    b <- band(d, column, "arm")
    b$hl_years <- rows$hl_years[i]
    b$retention_years <- rows$retention_years[i]
    b
  }))
  if (is.null(lines)) {
    message("  no comparator runs at these retentions yet")
    lines <- points[0, ]
  }
  lines$arm <- factor(lines$arm, levels = itn_arms)
  # The arm a comparison is against sits at nothing, with nothing either
  # side of it: one line says that, three on top of each other do not.
  ends <- lines[!(lines$median == 0 & lines$lower == 0 & lines$upper == 0), ]

  missing <- nrow(rows) - nrow(unique(lines[c("hl_years", "retention_years")]))
  if (missing > 0) {
    message("  ", missing, " panel(s) have no comparator yet: the ",
            "seven-arm run at that retention is not summarised")
  }

  ggplot(points, aes(arm, median, colour = arm)) +
    # The comparators' intervals are a backdrop, so they are shaded
    # rather than ruled: three pairs of dotted lines crowded the panels
    # and read as content. They are drawn first and faintly, the three
    # of them overlapping over much of the range.
    geom_rect(data = ends, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper,
                  fill = arm), alpha = 0.10) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_hline(data = lines, aes(yintercept = median, colour = arm),
               linetype = "dashed", linewidth = 0.5) +
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.9) +
    geom_point(size = 2.4) +
    facet_grid(hl_years ~ retention_years,
               labeller = labeller(
                 hl_years = function(x) paste0("Half-life ", x, " y"),
                 retention_years = function(x) paste0("Retention ", x, " y"))) +
    scale_colour_manual(values = arm_colours,
                        breaks = c(itn_arms, atn_arms),
                        labels = arm_labels[c(itn_arms, atn_arms)],
                        name = NULL) +
    scale_fill_manual(values = arm_colours, guide = "none") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    labs(x = NULL,
         y = sprintf("Clinical cases averted per 1,000 people per year, vs %s",
                     against[[key]])) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 11),
          legend.position = "bottom",
          axis.ticks.x = element_blank())
}

for (key in names(quantities)) {
  message("Drawing ", key)
  fig <- grid_figure(key)
  stem <- file.path(out_dir, sprintf("sweep_grid_%s_%s", key, form))
  ggsave(paste0(stem, ".png"), fig, width = 13, height = 9, dpi = 300,
         device = ragg::agg_png, bg = "white")
  ggsave(paste0(stem, ".pdf"), fig, width = 13, height = 9,
         device = grDevices::cairo_pdf)
  message("Saved: ", stem, ".(png|pdf)")
}
