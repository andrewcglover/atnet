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
  library(patchwork)
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
against <- c(vs_nonets = "no future nets", vs_pyr = "Pyr-only",
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

# What one comparison's grid is drawn from. Separated from the drawing
# so that the four can be worked out first and share a y scale.
grid_data <- function(key) {
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
  list(points = points, lines = lines, ends = ends)
}

# Everything a grid draws, for setting a scale the four can share.
grid_range <- function(d) {
  range(c(d$points$lower, d$points$upper, d$lines$lower, d$lines$upper),
        na.rm = TRUE)
}

grid_figure <- function(key, d = grid_data(key), ylim = NULL,
                        legend = TRUE, title = NULL, ylab = NULL,
                        row_strips = TRUE, base = 12, yticks = TRUE,
                        tag = NULL, border = FALSE) {
  # Dropping the tick labels says the scale is somebody else's to read;
  # that is only true when one has been handed in.
  stopifnot(yticks || !is.null(ylim))
  points <- d$points; lines <- d$lines; ends <- d$ends

  ggplot(points, aes(arm, median, colour = arm)) +
    # The comparators are a backdrop in three weights: their interval as
    # a tint, its ends as faint dashed rules, and the median as a solid
    # line held back from full strength, so the A(I)TN points stay the
    # darkest thing in the panel. The tints are faint: they are wide,
    # overlapping, and cover much of a panel, so anything stronger
    # reads as a background colour rather than as an interval.
    geom_rect(data = ends, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper,
                  fill = arm), alpha = 0.10) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_hline(data = ends, aes(yintercept = lower, colour = arm),
               linetype = "dashed", linewidth = 0.4, alpha = 0.35) +
    geom_hline(data = ends, aes(yintercept = upper, colour = arm),
               linetype = "dashed", linewidth = 0.4, alpha = 0.35) +
    # `ends` rather than `lines`: the arm a comparison is against is at
    # nothing by construction, and the grey rule at zero says so.
    geom_hline(data = ends, aes(yintercept = median, colour = arm),
               linewidth = 0.5, alpha = 0.65) +
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.9) +
    geom_point(size = 2.4) +
    # Retention down, half-life across: the comparators depend on the
    # retention alone, so laid out this way their lines are the same
    # the length of a row and read as one backdrop to it.
    facet_grid(retention_years ~ hl_years,
               labeller = labeller(
                 hl_years = function(x) paste0("Half-life ", x, " y"),
                 retention_years = function(x) paste0("Retention ", x, " y"))) +
    scale_colour_manual(values = arm_colours,
                        breaks = c(itn_arms, atn_arms),
                        labels = arm_labels[c(itn_arms, atn_arms)],
                        name = NULL,
                        # Four grids in one figure want one legend. They
                        # are drawn from different data, so patchwork
                        # does not take them for copies of each other
                        # and would lay all four side by side.
                        guide = if (legend) guide_legend(nrow = 1)
                                else "none") +
    scale_fill_manual(values = arm_colours, guide = "none") +
    scale_x_discrete(labels = NULL, breaks = NULL) +
    # Labelled every 200 with a line every 100 between. The sequence
    # runs well past any panel; the scale shows only what it reaches.
    scale_y_continuous(breaks = seq(-2000, 2000, 200),
                       minor_breaks = seq(-2000, 2000, 100)) +
    coord_cartesian(ylim = ylim) +
    labs(x = NULL, title = title, tag = tag,
         y = if (is.null(ylab))
           sprintf("Clinical cases averted per 1,000 people per year vs %s",
                   against[[key]])
         else if (isFALSE(ylab)) NULL else ylab) +
    theme_minimal(base_size = base) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text = element_text(size = base - 1),
          strip.text.y = if (row_strips) element_text(size = base - 1)
                         else element_blank(),
          # The title carries the letter beside it, so the room below
          # it is the room between both of them and the strips.
          plot.title = element_text(size = base + 1, hjust = 0.5,
                                    margin = margin(b = 12)),
          # A grid with no scale of its own gives up the room for one,
          # so its panels start at its rule.
          axis.text.y = if (yticks) element_text() else element_blank(),
          axis.title.y = if (isFALSE(ylab)) element_blank()
                         else element_text(angle = 90),
          axis.ticks.length.y = if (yticks) unit(2.75, "pt") else unit(0, "pt"),
          legend.position = "bottom",
          axis.ticks.x = element_blank(),
          # The rule belongs to a grid, so it is set here and not over
          # the whole figure, which would draw one round everything and
          # another round the legend.
          plot.background = if (border)
            element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
          else element_blank(),
          # The columns of a row carry one backdrop between them, so
          # they run together rather than sitting apart. The rows do not,
          # and are given a little more room than the default.
          panel.spacing.x = unit(0, "pt"),
          panel.spacing.y = unit(10, "pt"))
}

save_figure <- function(fig, stem, width, height) {
  ggsave(paste0(stem, ".png"), fig, width = width, height = height, dpi = 300,
         device = ragg::agg_png, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".pdf"), fig, width = width, height = height,
         device = grDevices::cairo_pdf, limitsize = FALSE)
  message("Saved: ", stem, ".(png|pdf)")
}

parts <- lapply(names(quantities), function(key) {
  message("Reading ", key)
  grid_data(key)
})
names(parts) <- names(quantities)

# One figure per comparison, each on the scale its own numbers need.
for (key in names(quantities)) {
  save_figure(grid_figure(key, parts[[key]]),
              file.path(out_dir, sprintf("sweep_grid_%s_%s", key, form)),
              11, 12)
}

# And the four together, A to D, on one scale so that the size of a
# difference can be carried from one comparison to the next. The scale
# has to hold everything from the deepest the weakest arm falls below
# Pyr-CFP to the most the strongest averts against no nets, so each
# grid uses less of its height than it does on its own.
needed <- range(vapply(parts, grid_range, numeric(2)))
shared <- c(-600, 1000)
message(sprintf("Shared y scale: %.0f to %.0f (holding %.0f to %.0f)",
                shared[1], shared[2], needed[1], needed[2]))
if (needed[1] < shared[1] || needed[2] > shared[2]) {
  warning("The shared scale cuts off ", sprintf("%.0f to %.0f", needed[1],
          needed[2]), call. = FALSE)
}
grids <- lapply(seq_along(quantities), function(i) {
  key <- names(quantities)[i]
  # One y label, on A, naming the quantity without naming a comparison;
  # each grid says which comparison it is in its title. The retention
  # strips are on D alone: the rows line up across all four, so one set
  # of labels does for them.
  grid_figure(key, parts[[key]], ylim = shared, legend = i == 1L,
              title = paste0("Comparator: ", against[[key]]),
              ylab = if (i == 1L)
                "Clinical cases averted per 1,000 people per year vs comparator"
              else FALSE,
              row_strips = i == length(quantities), base = 10,
              yticks = i == 1L, tag = LETTERS[i], border = TRUE)
})
# The legend gets a row of its own, below the grids, so that their rules
# close under them instead of reaching down around it. The letters are
# each grid's own, patchwork's tagging having nothing to tag the legend's
# slot with.
combined <- patchwork::wrap_plots(
  c(grids, list(patchwork::guide_area())),
  design = "ABCD\nEEEE", heights = c(1, 0.05), guides = "collect") &
  # The letter is placed on the grid rather than beside it. Left where
  # it defaults, it takes a column of its own at the grid's left edge
  # and holds the panels off their rule by its own width.
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 20, face = "bold",
                                                  hjust = 0, vjust = 1),
                 plot.tag.position = c(0.004, 0.997),
                 # The grids meet: with no y labels on B to D and no
                 # retention strips on A to C there is nothing between
                 # them to collide.
                 plot.margin = ggplot2::margin(4, 0, 4, 0))
# Four grids in a row, each about two thirds the width it has on its own.
save_figure(combined, file.path(out_dir, sprintf("sweep_grid_all_%s", form)),
            14, 10.5)
