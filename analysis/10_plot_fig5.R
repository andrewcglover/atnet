# Fig. 5: projected PfPR2-10 in each region of Mali (A), and national
# clinical cases averted by each net type over the quoted campaign cycle (B).
#
# A follows the June figure (analysis/legacy/c24med_projection_plots_select.R,
# section 7): the regions laid out roughly as on the map, 2024-2031, with a
# faint line for the daily median across draws-table rows and a bold line for
# the median of each row's centred 365-day rolling mean. Until 2025 every arm
# of a row shares one simulation, so that stretch is drawn once, in dark grey;
# the bold lines part half a year before 2025, where the rolling window first
# reaches the arms.
# B gives clinical cases averted per 1,000 people per year against no future
# nets, from the 2028 campaign to the day before the 2031 campaign, for the
# country as a whole: each region's cases are scaled to its population before
# they are added up, arms are compared within each row, and the dot and bar
# are the median and 95% interval across rows (see
# analysis/09_summarise_projections.R).
#
# Two versions, written to outputs/figures/<first run folder's name>/, with
# the antimalarial's decay form (exp or smoothcompact) in the file name:
#   fig5_medians_<form>.(png|pdf)  medians for all seven arms (main text)
#   fig5_bands_<form>.(png|pdf)    medians and 95% intervals for Pyr-CFP,
#                            ATN, Pyr-ATN and Pyr-CFP-ATN (appendix); each
#                            band is for one arm on its own, although the
#                            arms are paired within rows
#
# NB for the caption: the bold lines are a CENTRED 365-day mean, so they
# begin to separate half a year before the arms do. The arms themselves part
# on the dashed line, the last day of 2024; at that point half of the
# window is already arm days. The faint daily lines part on the day itself.
#
# Needs every region of the site file. Reads the outputs only; it runs no
# simulations. Name one run folder, or several when the arms are spread
# over more than one run: each arm is taken from the first folder named
# that holds it, so a folder named later only fills gaps. The exponential
# main run, split by row:
#   Rscript analysis/10_plot_fig5.R \
#     outputs/projections_v2/mli_hldrawn_v2_main \
#     outputs/projections_v2/mli_hldrawn_v2_main_161_640
# The smooth-compact run, which holds the ATN arms alone, with the
# comparator arms from the exponential run (they do not depend on the decay
# form):
#   Rscript analysis/10_plot_fig5.R \
#     outputs/projections_v2/mli_hldrawn_v2_scmain \
#     outputs/projections_v2/mli_hldrawn_v2_main \
#     outputs/projections_v2/mli_hldrawn_v2_main_161_640
#
# Reading 640 rows for nine regions takes some minutes, so what is read is
# cached in the output folder and reused when the same folders are named
# again. Add the word "refresh" to the command to read the runs afresh.
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(sf)
})

args <- commandArgs(trailingOnly = TRUE)
refresh <- "refresh" %in% args
run_dirs <- setdiff(args, "refresh")
if (!length(run_dirs)) stop("Name at least one run folder.", call. = FALSE)
settings <- readRDS(file.path(run_dirs[1], "settings.rds"))
out_dir <- file.path("outputs", "figures", basename(run_dirs[1]))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_file <- file.path(out_dir, "fig5_data.rds")

site <- readRDS(file.path("data_private", "site_files",
                          paste0(settings$country_iso, ".rds")))
country <- unique(site$eir$country)

# The June layout: regions placed roughly as on the map.
geo_grid <- geofacet::grid_auto(site$shape$level_1, names = "name_1", seed = 1)
names(geo_grid)[names(geo_grid) == "name_name_1"] <- "region"
taken <- paste(geo_grid$row, geo_grid$col)
corners <- c(paste(rep(1:2, 2), rep(1:2, each = 2)),   # legend
             paste(rep(3:4, 2), rep(4:5, each = 2)))   # panel B
if (any(corners %in% taken)) {
  stop("The legend and panel B corners are not empty in this grid.",
       call. = FALSE)
}

# Paul Tol's "vibrant" scheme, as in June, plus dark grey for 2022-24.
arms <- projection_arms()
arms$colour <- c("#BBBBBB", "#EE7733", "#CC3311", "#EE3377", "#009988",
                 "#0077BB", "#33BBEE")
key_colours <- c(setNames(arms$colour, arms$arm), shared = "#555555")
key_labels <- c(setNames(arms$label, arms$arm), shared = "Before the switch (all arms)")
band_arms <- c("pyr_cfp", "atn", "pyr_atn", "pyr_cfp_atn")

model_year <- function(t) 2000 + (t - 1) / 365
half_window <- 182L

# Panel A starts just before 2024: a year of the shared history is enough
# to place the arms against, and the earlier years crowd the panels. The
# few days before 2024 are there so the distribution on 1 January 2024 is
# drawn as a line rather than along the edge of the panel.
x_from <- 2023.95
x_to <- 2032


# Daily and rolling medians (and the rolling 95% interval) across rows for
# one region, each arm's series being 2022-24 followed by the arm.
region_series <- function(region) {
  message("  reading ", region)
  x <- read_projection_jobs(run_dirs, regions = region, columns = "pfpr2to10")
  pre <- x[x$arm == "2022-2024", ]
  x <- rbind(x[x$arm %in% arms$arm, ],
             do.call(rbind, lapply(arms$arm, function(a) transform(pre, arm = a))))
  x %>%
    arrange(draw, arm, timestep) %>%
    group_by(draw, arm) %>%
    mutate(roll = zoo::rollmean(100 * pfpr2to10, 365L, fill = NA,
                                align = "center")) %>%
    group_by(arm, timestep) %>%
    summarise(daily = median(100 * pfpr2to10),
              roll_med = median(roll),
              roll_lo = quantile(roll, 0.025, na.rm = TRUE, names = FALSE),
              roll_hi = quantile(roll, 0.975, na.rm = TRUE, names = FALSE),
              .groups = "drop") %>%
    mutate(region = region, year = model_year(timestep))
}

# Everything read from the runs, in one place: the regional series for A,
# the per-row totals for B, and the day the arms part. Reading 640 rows for
# nine regions takes some minutes, so it is kept and reused.
cycle <- campaign_cycle(settings$future_campaign_days)
cached <- if (!refresh && file.exists(cache_file)) readRDS(cache_file)
if (is.null(cached) || !identical(cached$run_dirs, run_dirs) ||
    !identical(cached$cycle, cycle)) {
  message("Reading ", paste(basename(run_dirs), collapse = " + "), " ...")
  cached <- list(
    run_dirs = run_dirs,
    cycle = cycle,
    series = bind_rows(lapply(geo_grid$region, region_series)),
    values = projection_row_values(run_dirs, cycle["from"], cycle["to"]),
    t0 = min(readRDS(projection_job_files(run_dirs)$file[1])$results$timestep)
  )
  saveRDS(cached, cache_file)
} else {
  message("Reusing ", cache_file, " (add \"refresh\" to read the runs again).")
}
# The cache holds every day read; panel A shows from x_from on. Dropping
# the rest here rather than at the scale keeps the rolling means intact
# (they were worked out over the whole series) and saves ggplot discarding
# the same rows, with a warning, once per panel.
series <- cached$series %>% filter(year >= x_from)
values <- cached$values
t0 <- cached$t0

# Every mass distribution the panels show. The site file's own ITN
# distributions are yearly, on 1 January (07), and the one in 2024 is the
# last before the arms part: it is what they all start from. The future
# campaigns are read from the run, the start of the month two months
# before the month of peak rainfall, so 1 June of 2025, 2028 and 2031.
hist_years <- seq(ceiling(x_from), floor(model_year(t0)))
campaign_years <- c(model_year((hist_years - 2000L) * 365L + 1L),
                    model_year(sort(settings$future_campaign_days)))

# Everything before the arms part is drawn in grey, the rolling lines
# included. A centred 365-day mean reaches the arms half a year early, so
# over the second half of 2024 the grey rolling lines do fan out a little;
# that half year still belongs to the old nets, which hold the majority of
# the window until the last day of 2024 (182 arm days of 365 on the 30th,
# 183 on the 31st). So the colours take over on the day the arms part, and
# the fan before it reads as what it is: the window filling up, not a
# change of net. The two stretches share their end day, so the lines join.
split_shared <- function(s, value, shared_until) {
  s$value <- s[[value]]
  bind_rows(
    s %>% filter(timestep <= shared_until) %>% mutate(key = "shared"),
    s %>% filter(timestep >= shared_until) %>% mutate(key = arm)
  ) %>%
    mutate(key = factor(key, levels = c("shared", arms$arm)),
           grp = paste(key, arm))
}
daily <- split_shared(series, "daily", t0 - 1L)
rolling <- split_shared(series, "roll_med", t0 - 1L)

# One grey band over the shared stretch, not one per arm: translucent bands
# laid on top of each other would read as a darker block where they are
# merely identical. It is the union of the arms' bands, so it is an
# envelope rather than an interval -- which is immaterial here, since the
# arms differ by almost nothing until they part.
shared_band <- function(display_arms) {
  rolling %>%
    filter(key == "shared", arm %in% display_arms, !is.na(roll_lo)) %>%
    group_by(region, timestep, year) %>%
    summarise(roll_lo = min(roll_lo), roll_hi = max(roll_hi), .groups = "drop") %>%
    mutate(key = factor("shared", levels = c("shared", arms$arm)),
           grp = "shared")
}
band_data <- function(display_arms) {
  bind_rows(
    shared_band(display_arms),
    rolling %>% filter(key != "shared", arm %in% display_arms, !is.na(roll_lo))
  )
}

# The lines of one figure, for the arms it shows. Over the shared stretch
# the daily series is one line, the arms being the same simulation there,
# while the rolling series keeps one line per arm, since a centred window
# has started to tell them apart.
line_data <- function(d, display_arms, shared_one) {
  shared <- d %>% filter(key == "shared")
  shared <- if (shared_one) shared %>% filter(arm == arms$arm[1])
            else shared %>% filter(arm %in% display_arms)
  bind_rows(shared, d %>% filter(key != "shared", arm %in% display_arms))
}
# The y scale is shared by every panel, and set by what is on show.
y_max <- 10 * ceiling(max(c(series$daily, series$roll_hi), na.rm = TRUE) / 10)

region_panel <- function(region, keys, show_daily, show_bands) {
  pos <- geo_grid[geo_grid$region == region, ]
  show_x <- !any(geo_grid$row == pos$row + 1 & geo_grid$col == pos$col)
  show_y <- !any(geo_grid$row == pos$row & geo_grid$col == pos$col - 1)
  # The grey stretch is kept only for the arms this figure shows, so the
  # fan before the arms part has as many lines as the colours after it.
  display_arms <- setdiff(keys, "shared")
  pick <- function(d) d[d$region == region & d$key %in% keys, ]
  p <- ggplot(mapping = aes(year, value, colour = key, fill = key, group = grp)) +
    geom_vline(xintercept = campaign_years, linetype = "dotted",
               colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = model_year(t0), linetype = "dashed",
               colour = "grey40", linewidth = 0.4)
  if (show_bands) {
    p <- p + geom_ribbon(data = pick(band_data(display_arms)),
                         aes(ymin = roll_lo, ymax = roll_hi),
                         colour = NA, alpha = 0.2)
  }
  if (show_daily) {
    p <- p + geom_line(data = pick(line_data(daily, display_arms, TRUE)),
                       linewidth = 0.4, alpha = 0.2)
  }
  p +
    geom_line(data = pick(line_data(rolling, display_arms, FALSE)),
              linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ region) +
    scale_colour_manual(values = key_colours, guide = "none") +
    scale_fill_manual(values = key_colours, guide = "none") +
    scale_x_continuous(breaks = seq(2024, 2030, 2), limits = c(x_from, x_to)) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, y_max)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = if (show_x) element_text() else element_blank(),
          axis.text.y = if (show_y) element_text() else element_blank(),
          strip.text = element_text(size = 13))
}

# A plot holding only the legend (no data are drawn), which also carries
# the "A" tag.
legend_panel <- function(keys, title, show_bands) {
  d <- data.frame(key = factor(keys, levels = keys), year = NA_real_,
                  value = NA_real_)
  p <- ggplot(d, aes(year, value, colour = key, fill = key))
  if (show_bands) {
    p <- p + geom_ribbon(aes(ymin = value, ymax = value), colour = NA,
                         alpha = 0.2, na.rm = TRUE)
  }
  p +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = key_colours, limits = keys,
                        labels = key_labels[keys], name = title) +
    scale_fill_manual(values = key_colours, limits = keys,
                      labels = key_labels[keys], name = title) +
    labs(tag = "A") +
    theme_void(base_size = 13) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.5, 0.5),
          legend.key.width = unit(2, "lines"),
          legend.text = element_text(size = 13),
          plot.tag = element_text(size = 18, face = "bold"))
}

# Panel B: national cases averted per 1,000 per year over the quoted cycle.
national <- national_values(values, region_populations(site), name = country)
averted <- summarise_rows(arm_contrasts(national, comparator = "none"),
                          "cases_averted_per1000", digits = 0)
averted$label <- factor(arms$label[match(averted$arm, arms$arm)],
                        levels = rev(arms$label[-1]))
model_month <- function(t) {
  format(as.Date("2001-01-01") + (t - 1) %% 365, "%b")
}
cycle_text <- sprintf("%s, %s %d to %s %d", country,
                      model_month(cycle["from"]), floor(model_year(cycle["from"])),
                      model_month(cycle["to"]), floor(model_year(cycle["to"])))
panel_b <- ggplot(averted, aes(median, label, colour = arm)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = key_colours, guide = "none") +
  scale_x_continuous(limits = c(min(0, averted$lower), NA),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Clinical cases averted per 1,000\npeople per year, vs no future nets",
       y = NULL, title = cycle_text, tag = "B") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 13),
        plot.title.position = "plot",
        plot.tag = element_text(size = 18, face = "bold"))

y_title <- wrap_elements(full = grid::textGrob(
  expression(italic(Pf) * PR[2 - 10] * " (%)"), rot = 90,
  gp = grid::gpar(fontsize = 15)))
x_title <- wrap_elements(full = grid::textGrob("Year",
                                               gp = grid::gpar(fontsize = 15)))

build_figure <- function(keys, legend_title, show_daily, show_bands) {
  design <- c(
    area(1, 1, 4, 1),   # y title
    area(5, 2, 5, 6),   # x title
    area(1, 2, 2, 3),   # legend, tag A
    area(3, 5, 4, 6),   # panel B
    do.call(c, lapply(seq_len(nrow(geo_grid)), function(i) {
      area(geo_grid$row[i], geo_grid$col[i] + 1)
    }))
  )
  panels <- lapply(geo_grid$region, region_panel, keys = keys,
                   show_daily = show_daily, show_bands = show_bands)
  plots <- c(list(y_title, x_title,
                  wrap_elements(full = legend_panel(keys, legend_title,
                                                    show_bands)),
                  wrap_elements(full = panel_b)),
             panels)
  wrap_plots(plots, design = design, widths = c(0.12, rep(1, 5)),
             heights = c(rep(1, 4), 0.1))
}

# The antimalarial's decay form goes in the file name: the two forms give
# two versions of this figure and they must not be told apart by their
# folder alone. A run built before the parameter existed used the
# exponential.
decay_form <- if (is.null(settings$atn_decay_form)) "exp" else settings$atn_decay_form

save_figure <- function(fig, name) {
  stem <- file.path(out_dir, paste0(name, "_", decay_form))
  ggsave(paste0(stem, ".png"), fig, width = 13, height = 9, dpi = 300,
         device = ragg::agg_png, bg = "white")
  ggsave(paste0(stem, ".pdf"), fig, width = 13, height = 9,
         device = grDevices::cairo_pdf)
  message("Saved: ", stem, ".(png|pdf)")
}

save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE),
            "fig5_medians")
save_figure(build_figure(c(band_arms, "shared"),
                         "Median and 95% uncertainty\ninterval across draws",
                         show_daily = FALSE, show_bands = TRUE),
            "fig5_bands")
