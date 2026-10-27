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
# Two layouts, each in two versions, written to
# outputs/figures/<first run folder's name>/ with the antimalarial's decay
# form (exp or smoothcompact) in the file name:
#   fig5_medians_<form>.(png|pdf)  medians for all seven arms (main text)
#   fig5_bands_<form>.(png|pdf)    medians and 95% intervals for Pyr-CFP,
#                            ATN, Pyr-ATN and Pyr-CFP-ATN (appendix); each
#                            band is for one arm on its own, although the
#                            arms are paired within rows
#   fig5_national_medians_<form>, fig5_national_bands_<form>
#                            the same two, with the legend moved to the
#                            empty cell in the bottom left and the top
#                            left corner given to the national series
#
# The national series is not an average of the regional ones: each
# region's daily counts are scaled to its real population (the site
# file's) before they are added up and the rate is formed, exactly as
# national_values() does for panel B, so the two national panels weight
# the regions identically.
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
# "complete" plots a run that is still filling in, over the realisations every region
# holds. Without it a part-finished run stops in read_all_series (regions must line up).
only_complete <- "complete" %in% args
# "boundsonly" writes ONLY the bounds figures, for iterating on their design without
# re-rendering the ten others at 300 dpi each time. Never use it for a final set.
only_bounds <- "boundsonly" %in% args
# pessimistic=<dir> and optimistic=<dir> switch on the bounds version of the A-F figure.
# Each names a folder under outputs/summaries -- a bounds run ALREADY SUMMARISED by 09 --
# whose National lines for the three antimalarial arms are drawn beside the primary run's
# as extra rows in C-F. The summaries are read rather than the runs themselves because the
# National totals there come from the same national_values() call this script makes, and
# reading them costs nothing next to a second pass over the sweep files. Only the
# antimalarial arms are taken: the others carry no antimalarial, so the bounds cannot move
# them, and their rows stay single.
#   Rscript analysis/10_plot_fig5.R outputs/projections_v4/mli_hldrawn_v4_scmain \
#     pessimistic=mli_hl2_v4_scmin optimistic=mli_hl5_v4_scmax
bound_args <- grep("^(pessimistic|optimistic)=", args, value = TRUE)
BOUNDS <- stats::setNames(sub("^[a-z]+=", "", bound_args), sub("=.*", "", bound_args))
if (length(BOUNDS) == 1L) {
  stop("Name BOTH pessimistic= and optimistic=, or neither.", call. = FALSE)
}
run_dirs <- setdiff(args, c("refresh", "complete", "boundsonly", bound_args))
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
# "Pyr-only" rather than "Pyr" here, so the arm reads as a choice of net alongside Pyr-PBO
# and Pyr-CFP rather than as a prefix (user, 2026-10-26). Overridden in this script rather
# than in projection_arms(), which 09 and 11-13 also use: changing it there would quietly
# relabel every other figure and summary table too.
arms$label[arms$arm == "pyr"] <- "Pyr-only"
arms$colour <- c("#BBBBBB", "#EE7733", "#CC3311", "#EE3377", "#009988",
                 "#0077BB", "#33BBEE")
key_colours <- c(setNames(arms$colour, arms$arm), shared = "#555555")
key_labels <- c(setNames(arms$label, arms$arm), shared = "Historic")
band_arms <- c("pyr_cfp", "atn", "pyr_atn", "pyr_cfp_atn")

# Drawing order, so that where lines and bands overlap the later arm lies
# on top: the shared stretch first, then the arms in the order of
# projection_arms(), ending with Pyr-CFP-ATN. ggplot draws groups in the
# order of the grouping factor's levels, so the order has to be set here
# and not left to the alphabet.
grp_levels <- c("shared", paste("shared", arms$arm), paste(arms$arm, arms$arm))
order_groups <- function(d) {
  d$grp <- factor(d$grp, levels = grp_levels)
  d[order(d$grp), , drop = FALSE]
}

model_year <- function(t) 2000 + (t - 1) / 365
half_window <- 182L

# Panel A starts just before 2024: a year of the shared history is enough
# to place the arms against, and the earlier years crowd the panels. The
# few days before 2024 are there so the distribution on 1 January 2024 is
# drawn as a line rather than along the edge of the panel.
x_from <- 2023.95
x_to <- 2032


# Daily and rolling medians, and the rolling 95% interval, across rows.
# Each arm's series is the shared 2022-24 stretch followed by the arm.
summarise_series <- function(x, label) {
  x %>%
    group_by(draw, arm) %>%
    mutate(roll = zoo::rollmean(100 * pfpr2to10, 365L, fill = NA,
                                align = "center")) %>%
    group_by(arm, timestep) %>%
    summarise(daily = median(100 * pfpr2to10),
              roll_med = median(roll),
              roll_lo = quantile(roll, 0.025, na.rm = TRUE, names = FALSE),
              roll_hi = quantile(roll, 0.975, na.rm = TRUE, names = FALSE),
              .groups = "drop") %>%
    mutate(region = label, year = model_year(timestep))
}

# Realisations a still-filling run holds in EVERY region. The national series adds
# population-scaled counts across regions, so a realisation missing one region would make
# the country look smaller rather than merely noisier -- read_all_series' stopifnot catches
# the mismatch and stops. With "complete" the run is plotted over the realisations all nine
# regions hold: the honest partial picture rather than no picture at all.
complete_draws <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      f <- unlist(lapply(run_dirs, function(d) list.files(d, "^sweep_.*_rows.*\\.rds$")))
      blocks <- sub("\\.rds$", "", sub("^.*_rows", "", f))
      n <- table(blocks)
      whole <- names(n)[n == length(geo_grid$region)]
      cache <<- sort(unlist(lapply(strsplit(whole, "-"), function(p)
        seq.int(as.integer(p[[1]]), as.integer(p[[2]])))))
    }
    cache
  }
})

# One region's daily counts, the 2022-24 stretch repeated for every arm,
# sorted so that regions line up with each other row for row.
region_counts <- function(region) {
  x <- read_projection_jobs(run_dirs, regions = region,
                            columns = c("n_detect_lm_730_3649",
                                        "n_age_730_3649"))
  x <- x[c("timestep", "arm", "draw", "n_detect_lm_730_3649",
           "n_age_730_3649")]
  if (only_complete) x <- x[x$draw %in% complete_draws(), , drop = FALSE]
  pre <- x[x$arm == "2022-2024", ]
  x <- rbind(x[x$arm %in% arms$arm, ],
             do.call(rbind, lapply(arms$arm, function(a) transform(pre, arm = a))))
  x[order(x$draw, x$arm, x$timestep), ]
}

# Every region's own series, and the national one, from a single pass over
# the files: PfPR2-10 is n_detect_lm_730_3649 / n_age_730_3649 exactly, so
# the counts give both. The national series adds up population-scaled
# counts before forming the rate, as national_values() does for the
# cycle, and takes the median across rows afterwards -- a median of
# regional medians would be a different quantity and not the country's.
# The scale is each region's real population over its simulated one, which
# is read from the cycle totals so that panel A's national line and panel
# B's national bars weight the regions identically.
read_all_series <- function(values) {
  pop <- region_populations(site)
  sim_pop <- 365 * values$person_years / values$days
  sim_key <- paste(region_tag(values$region), values$draw, values$arm)
  regional <- list()
  nat_detect <- NULL; nat_age <- NULL; shape <- NULL
  for (rg in geo_grid$region) {
    message("  reading ", rg)
    x <- region_counts(rg)
    regional[[rg]] <- summarise_series(
      transform(x, pfpr2to10 = n_detect_lm_730_3649 / n_age_730_3649), rg)

    # Blocks of one (row, arm) each, in the order x is sorted in, so the
    # per-block scale can simply be repeated down the column.
    draws <- sort(unique(x$draw)); use_arms <- sort(unique(x$arm))
    days <- nrow(x) / (length(draws) * length(use_arms))
    this_shape <- list(draws = draws, arms = use_arms, days = days,
                       timesteps = sort(unique(x$timestep)))
    stopifnot(days == length(this_shape$timesteps))
    if (is.null(shape)) shape <- this_shape
    stopifnot(identical(shape, this_shape))   # regions must line up

    blocks <- expand.grid(arm = use_arms, draw = draws,
                          stringsAsFactors = FALSE)
    scale <- pop$pop[match(region_tag(rg), region_tag(pop$region))] /
      sim_pop[match(paste(region_tag(rg), blocks$draw, blocks$arm), sim_key)]
    stopifnot(!anyNA(scale))
    scale <- rep(scale, each = days)
    if (is.null(nat_detect)) {
      nat_detect <- scale * x$n_detect_lm_730_3649
      nat_age <- scale * x$n_age_730_3649
    } else {
      nat_detect <- nat_detect + scale * x$n_detect_lm_730_3649
      nat_age <- nat_age + scale * x$n_age_730_3649
    }
    rm(x)
  }
  message("  national")
  nat <- data.frame(
    draw = rep(shape$draws, each = length(shape$arms) * shape$days),
    arm = rep(rep(shape$arms, each = shape$days), times = length(shape$draws)),
    timestep = rep(shape$timesteps,
                   times = length(shape$draws) * length(shape$arms)),
    pfpr2to10 = nat_detect / nat_age
  )
  list(series = bind_rows(regional), nat_series = summarise_series(nat, country))
}

# Everything read from the runs, in one place: the regional and national
# series for A, the per-row totals for B, and the day the arms part.
# Reading 640 rows for nine regions takes some minutes, so it is kept and
# reused. `made` is bumped when what is kept changes shape, so an older
# cache is read afresh rather than half-used.
cycle <- campaign_cycle(settings$future_campaign_days)
cached <- if (!refresh && file.exists(cache_file)) readRDS(cache_file)
if (is.null(cached) || !identical(cached$made, 2L) ||
    !identical(cached$run_dirs, run_dirs) || !identical(cached$cycle, cycle)) {
  message("Reading ", paste(basename(run_dirs), collapse = " + "), " ...")
  vals <- projection_row_values(run_dirs, cycle["from"], cycle["to"])
  both <- read_all_series(vals)
  cached <- list(
    made = 2L,
    run_dirs = run_dirs,
    cycle = cycle,
    series = both$series,
    nat_series = both$nat_series,
    values = vals,
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
nat_series <- cached$nat_series %>% filter(year >= x_from)
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
nat_daily <- split_shared(nat_series, "daily", t0 - 1L)
nat_rolling <- split_shared(nat_series, "roll_med", t0 - 1L)

# One grey band over the shared stretch, not one per arm: translucent bands
# laid on top of each other would read as a darker block where they are
# merely identical. It is the union of the arms' bands, so it is an
# envelope rather than an interval -- which is immaterial here, since the
# arms differ by almost nothing until they part.
shared_band <- function(roll, display_arms) {
  roll %>%
    filter(key == "shared", arm %in% display_arms, !is.na(roll_lo)) %>%
    group_by(region, timestep, year) %>%
    summarise(roll_lo = min(roll_lo), roll_hi = max(roll_hi), .groups = "drop") %>%
    mutate(key = factor("shared", levels = c("shared", arms$arm)),
           grp = "shared")
}
band_data <- function(roll, display_arms) {
  order_groups(bind_rows(
    shared_band(roll, display_arms),
    roll %>% filter(key != "shared", arm %in% display_arms, !is.na(roll_lo))
  ))
}

# The lines of one figure, for the arms it shows. Over the shared stretch
# the daily series is one line, the arms being the same simulation there,
# while the rolling series keeps one line per arm, since a centred window
# has started to tell them apart.
line_data <- function(d, display_arms, shared_one) {
  shared <- d %>% filter(key == "shared")
  shared <- if (shared_one) shared %>% filter(arm == arms$arm[1])
            else shared %>% filter(arm %in% display_arms)
  order_groups(bind_rows(shared, d %>% filter(key != "shared",
                                              arm %in% display_arms)))
}
# The y scale is shared by every panel, and set by what is on show.
y_max <- 10 * ceiling(max(c(series$daily, series$roll_hi), na.rm = TRUE) / 10)

region_panel <- function(region, keys, show_daily, show_bands, tag = NULL) {
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
    p <- p + geom_ribbon(data = pick(band_data(rolling, display_arms)),
                         aes(ymin = roll_lo, ymax = roll_hi),
                         colour = NA, alpha = 0.2)
  }
  if (show_daily) {
    p <- p + geom_line(data = pick(line_data(daily, display_arms, TRUE)),
                       linewidth = 0.4, alpha = 0.4)
  }
  p +
    geom_line(data = pick(line_data(rolling, display_arms, FALSE)),
              linewidth = 0.8, na.rm = TRUE) +
    facet_wrap(~ region) +
    scale_colour_manual(values = key_colours, guide = "none") +
    scale_fill_manual(values = key_colours, guide = "none") +
    scale_x_continuous(breaks = seq(2024, 2030, 2), limits = c(x_from, x_to)) +
    scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, y_max)) +
    labs(x = NULL, y = NULL, tag = tag) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = if (show_x) element_text() else element_blank(),
          axis.text.y = if (show_y) element_text() else element_blank(),
          strip.text = element_text(size = 13),
          plot.tag = element_text(size = 18, face = "bold"))
}

# The national series, as an inset drawn like a region's panel but on its
# own y scale, since a population-weighted national rate is well inside
# the range the regions cover.
national_panel <- function(keys, show_daily, show_bands) {
  display_arms <- setdiff(keys, "shared")
  # Headed the way a region is, by a centred strip rather than a title,
  # so the inset reads as one more panel of the same kind.
  pick <- function(d) {
    d <- d[d$key %in% keys, ]
    d$panel <- "National"
    d
  }
  # Its own scale, set by what the inset actually draws: the daily lines
  # reach well above the bands, so an inset showing bands alone does not
  # want the room they need.
  shown <- nat_series[nat_series$arm %in% display_arms, ]
  top <- 10 * ceiling(max(c(shown$roll_med,
                            if (show_daily) shown$daily,
                            if (show_bands) shown$roll_hi),
                          na.rm = TRUE) / 10)
  p <- ggplot(mapping = aes(year, value, colour = key, fill = key, group = grp)) +
    geom_vline(xintercept = campaign_years, linetype = "dotted",
               colour = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = model_year(t0), linetype = "dashed",
               colour = "grey40", linewidth = 0.4)
  if (show_bands) {
    p <- p + geom_ribbon(data = pick(band_data(nat_rolling, display_arms)),
                         aes(ymin = roll_lo, ymax = roll_hi),
                         colour = NA, alpha = 0.2)
  }
  if (show_daily) {
    p <- p + geom_line(data = pick(line_data(nat_daily, display_arms, TRUE)),
                       linewidth = 0.4, alpha = 0.4)
  }
  p +
    geom_line(data = pick(line_data(nat_rolling, display_arms, FALSE)),
              linewidth = 0.8, na.rm = TRUE) +
    scale_colour_manual(values = key_colours, guide = "none") +
    scale_fill_manual(values = key_colours, guide = "none") +
    # Minor gridlines every year and every 5 PfPR points, so the national panel can be read
    # off between its labelled ticks (user, 2026-10-26). The minor grid must not be blanked
    # in the theme below or these are silently discarded.
    scale_x_continuous(breaks = seq(2024, 2030, 2), minor_breaks = seq(2024, 2032, 1),
                       limits = c(x_from, x_to)) +
    facet_wrap(~ panel) +
    scale_y_continuous(breaks = seq(0, 100, 10), minor_breaks = seq(0, 100, 5),
                       limits = c(0, top)) +
    labs(x = NULL, y = NULL, tag = "A") +
    theme_minimal(base_size = 13) +
    theme(
          strip.text = element_text(size = 13),
          plot.tag = element_text(size = 18, face = "bold"),
          plot.background = element_rect(colour = "grey40", fill = NA,
                                         linewidth = 0.4),
          plot.margin = margin(12, 12, 12, 12))
}

# A plot holding only the legend (no data are drawn). It carries the "A"
# tag in the layout where it sits in the top left corner, and none in the
# one where the national inset has taken that corner.
# `scale` shrinks the key and its text together (0.82 = about 18% smaller), for layouts
# where the legend shares a row with map panels rather than having a corner to itself.
# 0.92 after 0.82 proved too small to read comfortably (user, 2026-10-26).
legend_panel <- function(keys, title, show_bands, tag = "A", scale = 1) {
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
    labs(tag = tag) +
    theme_void(base_size = 13 * scale) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.5, 0.5),
          legend.key.width = unit(2 * scale, "lines"),
          legend.key.height = unit(1.2 * scale, "lines"),
          legend.text = element_text(size = 13 * scale),
          legend.title = element_text(size = 13 * scale),
          plot.tag = element_text(size = 18, face = "bold"))
}

# Panel B: national cases averted per 1,000 per year over the quoted cycle.
national <- national_values(values, region_populations(site), name = country)
averted <- summarise_rows(arm_contrasts(national, comparator = "none"),
                          "cases_averted_per1000", digits = 0)
averted$label <- factor(arms$label[match(averted$arm, arms$arm)],
                        levels = rev(arms$label[-1]))
# The period panel B covers is left to the caption rather than titled, so
# it is printed here for whoever writes it.
model_month <- function(t) {
  format(as.Date("2001-01-01") + (t - 1) %% 365, "%b")
}
message(sprintf("Panel B covers %s, %s %d to %s %d", country,
                model_month(cycle["from"]), floor(model_year(cycle["from"])),
                model_month(cycle["to"]), floor(model_year(cycle["to"]))))

# The same panel against any comparator, so the figure can show what each net adds over no
# future nets AND over the nets a programme might otherwise deploy. Contrasts stay paired
# within the realisation, and the comparator arm drops out (its contrast with itself is 0).
averted_vs <- function(cmp) {
  raw <- arm_contrasts(national, comparator = cmp)
  a <- summarise_rows(raw, "cases_averted_per1000", digits = 0)
  # "No future nets" is dropped -- panel C already covers it, and it would stretch the
  # shared axis past -400 and squash the contrasts these panels exist to show. The
  # COMPARATOR ITSELF IS KEPT: its contrast with itself is exactly 0, so it shows as a dot
  # on the line, which keeps every panel in the same order and makes the reference explicit
  # rather than leaving a gap to work out (user, 2026-10-26).
  keep <- setdiff(arms$arm, "none")
  a <- a[a$arm %in% keep, , drop = FALSE]
  # arm_contrasts() returns no row for the comparator against itself, so add the zero it
  # stands for: without it the row is an empty slot rather than a dot on the line.
  if (!cmp %in% a$arm) {
    z <- a[1, , drop = FALSE]
    z[] <- lapply(z, function(col) if (is.numeric(col)) 0 else col[[1]])
    z$arm <- cmp
    a <- rbind(a, z)
  }
  a$label <- factor(arms$label[match(a$arm, arms$arm)],
                    levels = rev(arms$label[match(keep, arms$arm)]))
  # Posterior probability that this arm averts MORE cases than the comparator: the share of
  # realisations whose paired contrast is positive. Exactly 0 for the comparator against
  # itself, which is meaningless, so that row is left blank.
  p <- tapply(raw$cases_averted_per1000 > 0, raw$arm, mean)
  a$p_better <- ifelse(a$arm == cmp, NA_real_, p[a$arm])
  a
}

# ---- the bounds version of C-F ------------------------------------------------------
# Each antimalarial arm gains two more rows: the same quantity with the antimalarial's
# half-life, omega and p_A all held at the ends of their 95% prior ranges rather than drawn.
# They are NOT a wider interval round the primary -- each carries its own interval, and the
# primary already propagates the uncertainty these two replace with a fixed assumption.
VAR_LEV    <- c("pessimistic", "primary", "optimistic")
VAR_LABS   <- c(pessimistic = "Minimum plausible", primary = "Primary analysis",
                optimistic = "Maximum plausible")
VAR_SHAPES <- c(pessimistic = 25, primary = 16, optimistic = 17)  # solid down, dot, up
VAR_SIZES  <- c(pessimistic = 2.0, primary = 2.35, optimistic = 2.0)
# One whisker width throughout, the bounds' simply lighter. The MARKERS stay solid: they
# carry the median under each assumption, and a ~6 pt triangle at 45% on white washes out
# where a long whisker at the same alpha still reads (user, 2026-10-27, after seeing both).
VAR_ALPHA  <- c(pessimistic = 0.5, primary = 1, optimistic = 0.5)
VAR_OFFSET <- 0.29
atn_arm_keys <- c("atn", "pyr_atn", "pyr_cfp_atn")

if (length(BOUNDS)) {
  read_bounds <- function(dir) {
    f <- file.path("outputs", "summaries", dir, "row_values.csv")
    if (!file.exists(f)) {
      stop("No summary for ", dir, ": run analysis/09_summarise_projections.R on that ",
           "run first (it writes ", f, ").", call. = FALSE)
    }
    v <- utils::read.csv(f, stringsAsFactors = FALSE)
    v <- v[v$region == "National" & v$arm %in% atn_arm_keys, , drop = FALSE]
    if (!nrow(v)) stop("No National antimalarial lines in ", f, call. = FALSE)
    v$region <- country          # 09 calls it "National", this script calls it the country
    v[c("region", "draw", "arm", "clin_cases", "person_years")]
  }
  nat_b <- national[c("region", "draw", "arm", "clin_cases", "person_years")]
  nat_b$variant <- "primary"
  for (vn in names(BOUNDS)) {
    d <- read_bounds(BOUNDS[[vn]])
    # Contrasts are paired WITHIN a realisation, so a bounds run summarised over a
    # different set of rows would pair the wrong simulations together.
    if (!setequal(d$draw, national$draw)) {
      stop(BOUNDS[[vn]], " covers ", length(unique(d$draw)), " realisation(s) against ",
           length(unique(national$draw)), " here; summarise both over the same rows.",
           call. = FALSE)
    }
    d$variant <- vn
    nat_b <- rbind(nat_b, d)
  }
  nat_b$base <- nat_b$arm
  nat_b$arm  <- ifelse(nat_b$variant == "primary", nat_b$base,
                       paste0(nat_b$base, "__", nat_b$variant))
  bkey <- unique(nat_b[c("arm", "base", "variant")])

  averted_vs_b <- function(cmp) {
    raw <- arm_contrasts(nat_b, comparator = cmp)
    a   <- summarise_rows(raw, "cases_averted_per1000", digits = 0)
    a   <- merge(a, bkey, by = "arm")
    a   <- a[a$base %in% setdiff(arms$arm, "none"), , drop = FALSE]
    # the comparator's contrast with itself is exactly 0; keep the row so every panel holds
    # the same arms. Not for "none", which these panels drop.
    if (cmp %in% setdiff(arms$arm, "none") && !cmp %in% a$arm) {
      z <- a[1, , drop = FALSE]
      z[] <- lapply(z, function(col) if (is.numeric(col)) 0 else col[[1]])
      z$arm <- cmp; z$base <- cmp; z$variant <- "primary"
      a <- rbind(a, z)
    }
    pb <- tapply(raw$cases_averted_per1000 > 0, raw$arm, mean)
    a$p_better <- ifelse(a$arm == cmp, NA_real_, pb[a$arm])
    a$label   <- factor(arms$label[match(a$base, arms$arm)], levels = levels(averted$label))
    a$variant <- factor(a$variant, levels = VAR_LEV)
    # y set by hand rather than by position_dodge, so a one-variant arm sits exactly on its
    # label and a three-variant arm splits evenly, however dodge happens to count groups
    nv <- table(a$label)[as.character(a$label)]
    a$y <- as.integer(a$label) +
      ifelse(nv == 1, 0, (as.integer(a$variant) - 2L) * VAR_OFFSET)
    a
  }

  # The same panel as averted_panel(), on a continuous y so the three rows can be offset.
  averted_panel_b <- function(a, xlab, tag, xlim = NULL, y_labels = TRUE,
                              x_text = NULL, x_breaks = waiver(), x_minor = waiver()) {
    lev <- levels(a$label)
    p <- ggplot(a, aes(median, y, colour = base)) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
      geom_linerange(aes(xmin = lower, xmax = upper, alpha = variant),
                     linewidth = 0.85) +
      geom_point(aes(shape = variant, size = variant, fill = base), alpha = 1) +
      scale_colour_manual(values = key_colours, guide = "none") +
      scale_fill_manual(values = key_colours, guide = "none") +
      scale_shape_manual(values = VAR_SHAPES, guide = "none") +
      scale_size_manual(values = VAR_SIZES, guide = "none") +
      scale_alpha_manual(values = VAR_ALPHA, guide = "none") +
      scale_y_continuous(breaks = seq_along(lev),
                         labels = if (y_labels) lev else NULL,
                         limits = c(0.4, length(lev) + 1.1)) +
      scale_x_continuous(limits = if (is.null(xlim)) c(min(0, a$lower), NA) else xlim,
                         breaks = x_breaks, minor_breaks = x_minor,
                         expand = expansion(mult = c(0, 0.05))) +
      labs(x = xlab, y = NULL, tag = tag)
    if (!is.null(x_text)) {
      # ONE line per arm, on the primary's row: the primary's probability with the two
      # bounds in brackets. Stacking the three put them about 8 pt apart against ~7.7 pt
      # text, so they collided, and shrinking them far enough to fit left the bounds at
      # roughly 3 pt once the figure is reduced to a journal column (user, 2026-10-27).
      # A single fixed-colour layer, not a colour or size aesthetic: either scale added
      # here REPLACES the arm colour scale, and the dots and whiskers then try to read
      # "atn" as a colour name.
      tp <- a[a$variant == "primary", ]
      key <- paste(a$base, a$variant)
      lo  <- a$p_better[match(paste(tp$base, "pessimistic"), key)]
      hi  <- a$p_better[match(paste(tp$base, "optimistic"), key)]
      tp$shown <- ifelse(is.na(tp$p_better), "\u2013",
                  ifelse(is.na(lo) | is.na(hi),
                         sprintf("%.0f%%", 100 * tp$p_better),
                         sprintf("%.0f%% (%.0f, %.0f)",
                                 100 * tp$p_better, 100 * lo, 100 * hi)))
      p <- p +
        geom_text(data = tp, aes(x = x_text, y = y, label = shown), inherit.aes = FALSE,
                  hjust = 0.5, colour = "grey20", size = 3.9) +
        annotate("text", x = x_text, y = length(lev) + 0.85, label = "P(averts more)",
                 size = 3.9, colour = "grey25", hjust = 0.5) +
        coord_cartesian(clip = "off")
    }
    p + theme_minimal(base_size = 13) +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
            axis.text.y = element_text(size = 13, colour = "black"),
            plot.tag = element_text(size = 18, face = "bold"),
            plot.background = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
            plot.margin = margin(12, 12, 12, 12)) +
      (if (is.null(x_text)) NULL else theme(axis.title.x = element_text(size = 11.6)))
  }

  # A key for the three markers alone, in neutral grey: the arm is already keyed by colour,
  # so repeating it here would say the shapes belong to one arm.
  shape_legend_panel <- function(scale = 0.92) {
    d <- data.frame(v = factor(VAR_LEV, levels = VAR_LEV),
                    x = NA_real_, y = NA_real_)
    ggplot(d, aes(x, y, shape = v, fill = v)) +
      geom_point(size = 3.0 * scale, colour = "grey25", na.rm = TRUE) +
      # Reversed, so the key reads maximum, primary, minimum down the page exactly as the
      # three markers sit in each arm's row (user, 2026-10-27). Only the KEY is reversed:
      # VAR_LEV still drives the offsets, so reordering it would flip the panels too.
      scale_shape_manual(values = VAR_SHAPES, labels = VAR_LABS,
                         name = "ATN/AITN impact",
                         guide = guide_legend(reverse = TRUE)) +
      scale_fill_manual(values = stats::setNames(rep("grey25", 3), VAR_LEV),
                        labels = VAR_LABS, name = "ATN/AITN impact",
                        guide = guide_legend(reverse = TRUE)) +
      theme_void(base_size = 13 * scale) +
      theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5),
            legend.key.width = unit(1.3 * scale, "lines"),
            legend.key.height = unit(1.25 * scale, "lines"),
            legend.text = element_text(size = 13 * scale),
            legend.title = element_text(size = 13 * scale, face = "bold"))
  }
}

averted_panel <- function(a, xlab, tag, xlim = NULL, y_labels = TRUE,
                          x_text = NULL, x_breaks = waiver(), x_minor = waiver()) {
  # The probability column is drawn INSIDE the panel at a fixed x rather than on a secondary
  # axis: an axis title can only sit rotated beside its axis, and the column wants a plain
  # horizontal heading above it (user, 2026-10-26). Placing the figures at a y given by the
  # discrete level keeps each one on its own dot and whisker.
  lev <- levels(a$label)
  p <- ggplot(a, aes(median, label, colour = arm)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 1) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = key_colours, guide = "none") +
  scale_y_discrete(limits = lev, labels = if (y_labels) waiver() else NULL,
                   expand = expansion(add = c(0.6, 1.1))) +
  scale_x_continuous(limits = if (is.null(xlim)) c(min(0, a$lower), NA) else xlim,
                     breaks = x_breaks, minor_breaks = x_minor,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = xlab, y = NULL, tag = tag)
  if (!is.null(x_text)) {
    txt <- data.frame(label = factor(lev, levels = lev),
                      p = a$p_better[match(lev, a$label)])
    # A dash for the comparator against itself, rather than a blank: it marks the row as
    # deliberately not applicable instead of looking like a value that failed to print.
    # Written as an escape so the file's encoding cannot mangle it.
    txt$shown <- ifelse(is.na(txt$p), "–", sprintf("%.0f%%", 100 * txt$p))
    p <- p +
      # Both centred on x_text, which sits far enough inside the axis end that the heading
      # -- much wider than "100%" -- still clears the panel edge. Centring it at the
      # old 500 against a 600 axis was what clipped F's closing bracket.
      geom_text(data = txt, aes(x = x_text, y = label, label = shown),
                inherit.aes = FALSE, size = 4.2, colour = "grey25", hjust = 0.5) +
      annotate("text", x = x_text, y = length(lev) + 0.85, label = "P(averts more)",
               size = 3.9, colour = "grey25", fontface = "plain", hjust = 0.5) +
      # Let the heading overflow into the plot margin rather than be cut by the panel. The
      # rightmost of the three panels has slightly less usable width than the other two --
      # it abuts the figure edge -- so with clipping on it lost its closing bracket at every
      # column position tried, however much clear space sat to its right.
      coord_cartesian(clip = "off")
  }
  p +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        # Only the Y minor grid is blanked. Blanking the whole minor grid silently threw
        # away the minor_breaks set on the x scale, so the fine ticks never appeared.
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 13, colour = "black"),
        plot.tag = element_text(size = 18, face = "bold"),
        # A rule round the inset, set in from the edge of its cell so the
        # inset reads as a panel of its own rather than as the corner of
        # the figure.
        plot.background = element_rect(colour = "grey40", fill = NA,
                                       linewidth = 0.4),
        plot.margin = margin(12, 12, 12, 12)) +
  # D-F carry their label on one line across a panel a third of the figure wide, so it
  # needs a point less than panel C's two-line version to sit inside its panel.
  (if (is.null(x_text)) NULL else theme(axis.title.x = element_text(size = 11.6)))
}

panel_b <- averted_panel(averted,
  "Annual clinical cases averted per\n1,000 people vs no future nets", "B",
  x_minor = seq(0, 900, by = 50))

# D-F: the same against the nets a programme might otherwise deploy. One x range across the
# three so the benefit shrinking as the comparator strengthens is read off the page rather
# than from the axis labels.
extra_cmp <- c(pyr = "Pyr-only", pyr_pbo = "Pyr-PBO", pyr_cfp = "Pyr-CFP")
extra_averted <- lapply(names(extra_cmp), averted_vs)
names(extra_averted) <- names(extra_cmp)
# A fixed range for D-F, shared so the three are directly comparable, with the right-hand
# stretch beyond 400 left clear of ticks so the probability column does not collide with
# axis labels or gridlines (user, 2026-10-26). Ticks therefore stop at 400 while the axis
# runs to 600, and the figures sit at 500.
extra_xlim  <- c(-600, 1150)
extra_breaks <- seq(-600, 600, by = 200)
extra_minor <- seq(-600, 650, by = 50)
extra_xtext <- 950
local({
  r <- range(unlist(lapply(extra_averted, function(a) c(a$lower, a$upper))))
  if (r[[1]] < extra_xlim[[1]] || r[[2]] > 550) {
    message(sprintf("NOTE: contrasts span %.0f to %.0f; the fixed D-F axis is %.0f to %.0f ",
                    r[[1]], r[[2]], extra_xlim[[1]], extra_xlim[[2]]),
            "with ticks to 400, so some intervals may run under the probability column.")
  }
})
# No y tick labels in D-F: the arms are already named down the side of panel C and keyed by
# colour in the legend, and repeating seven names three more times crowds panels that are a
# third of the width (user, 2026-10-26).
extra_panels <- Map(function(nm, tag) {
  averted_panel(extra_averted[[nm]],
                sprintf("Annual cases averted per 1,000 vs %s", extra_cmp[[nm]]),
                tag, xlim = extra_xlim, y_labels = FALSE, x_text = extra_xtext,
                x_breaks = extra_breaks, x_minor = extra_minor)
}, names(extra_cmp), c("D", "E", "F"))

# The same two, with the bounds rows. Built only when the bounds folders were named.
if (length(BOUNDS)) {
  panel_b_b <- averted_panel_b(averted_vs_b("none"),
    "Annual clinical cases averted per\n1,000 people vs no future nets", "B",
    x_breaks = seq(0, 800, by = 200), x_minor = seq(0, 900, by = 50))
  # D-F's own data region, narrower than the one the figures without the bounds use
  # (user, 2026-10-27). Their own, not a change to extra_xlim: that is shared with those
  # figures, and this must leave them exactly as they were. Shape is the same -- data
  # region, then a clear strip of ticks holding the probability column.
  bounds_xlim   <- c(-450, 1150)
  bounds_breaks <- seq(-400, 600, by = 200)
  bounds_minor  <- seq(-450, 600, by = 50)
  bounds_dat    <- lapply(stats::setNames(names(extra_cmp), names(extra_cmp)), averted_vs_b)
  # scale_x_continuous LIMITS DROP a row rather than zoom to it, so a whisker reaching past
  # the end would take its point with it and simply not be drawn. At the realisations this
  # was built on the widest is -442 against -450, eight units of headroom, so say so loudly
  # rather than let a later run lose an arm in silence.
  local({
    r <- range(unlist(lapply(bounds_dat, function(a) c(a$lower, a$upper))))
    if (r[[1]] < bounds_xlim[[1]] || r[[2]] > 600) {
      message(sprintf("NOTE: bounds contrasts span %.0f to %.0f, outside the D-F data ",
                      r[[1]], r[[2]]),
              sprintf("region of %.0f to 600: whatever falls outside is DROPPED, not ",
                      bounds_xlim[[1]]),
              "clipped. Widen bounds_xlim.")
    }
  })
  extra_panels_b <- Map(function(nm, tag) {
    averted_panel_b(bounds_dat[[nm]],
                    sprintf("Annual cases averted per 1,000 vs %s", extra_cmp[[nm]]),
                    tag, xlim = bounds_xlim, y_labels = FALSE, x_text = extra_xtext,
                    x_breaks = bounds_breaks, x_minor = bounds_minor)
  }, names(extra_cmp), c("D", "E", "F"))
}

# G-I: the same again against the three ATN arms, for the figure that carries all of them.
# Its own x range, not D-F's: these contrasts reach -546 (Pyr-only against Pyr-CFP-ATN),
# which D-F's -450 would drop entirely rather than truncate. Same shape though -- data
# region, then a clear strip of the same width holding the probability column.
atn_cmp <- c(atn = "ATN", pyr_atn = "Pyr-ATN", pyr_cfp_atn = "Pyr-CFP-ATN")
atn_averted <- lapply(names(atn_cmp), averted_vs)
names(atn_averted) <- names(atn_cmp)
atn_xlim   <- extra_xlim
atn_breaks <- extra_breaks
atn_minor  <- extra_minor
atn_xtext  <- extra_xtext
local({
  r <- range(unlist(lapply(atn_averted, function(a) c(a$lower, a$upper))))
  if (r[[1]] < atn_xlim[[1]] || r[[2]] > 450) {
    message(sprintf("NOTE: ATN-comparator contrasts span %.0f to %.0f, outside the G-I ",
                    r[[1]], r[[2]]), "data region of -600 to 450.")
  }
})
atn_panels <- Map(function(nm, tag) {
  averted_panel(atn_averted[[nm]],
                sprintf("Annual cases averted per 1,000 vs %s", atn_cmp[[nm]]),
                tag, xlim = atn_xlim, y_labels = FALSE, x_text = atn_xtext,
                x_breaks = atn_breaks, x_minor = atn_minor)
}, names(atn_cmp), c("G", "H", "I"))

y_title <- wrap_elements(full = grid::textGrob(
  expression(italic(Pf) * PR[2 - 10] * " (%)"), rot = 90,
  gp = grid::gpar(fontsize = 15)))
# There is no x title: the years label themselves, and the caption says
# what the panels are.

# A plot laid into a share of the space it is given, so panel B sits a
# little inside its corner of the figure rather than filling it.
scaled <- function(p, s = 0.9) {
  grid::grobTree(ggplot2::ggplotGrob(p),
                 vp = grid::viewport(width = s, height = s))
}

# Two layouts. In both the regions keep their places on the map and panel
# B sits in the bottom right. In the first the legend has the top left
# corner; in the second the national series takes that corner and the
# legend drops to the empty cell south of Kayes and west of Sikasso.
# `nat_bands` lets the national inset show its bands while the regions
# show their medians: the country has room for them where nine small
# panels do not. `nat_keys` lets it show fewer arms than they do, since
# seven bands laid over each other are a muddle however much room there
# is.
#
# `comparators` adds a bottom row of three more cases-averted panels (D-F), against the
# nets a programme might otherwise deploy rather than against no future nets. For that row
# the column grid is refined threefold -- five map columns do not divide into thirds -- so
# each region spans three narrow columns and each of D, E and F spans five, exactly a third
# of the width. Its row is given twice the height of a map row, which brings A, C, D, E and
# F to roughly one size. The legend then moves to the cell north of Gao and east of Kidal,
# which the map leaves empty (user, 2026-10-26).
build_figure <- function(keys, legend_title, show_daily, show_bands,
                         national = FALSE, nat_bands = show_bands,
                         nat_keys = keys, comparators = FALSE, atn_comparators = FALSE,
                         bounds = FALSE) {
  if (bounds && !(national && comparators)) {
    stop("bounds = TRUE is only defined for the national A-F layout.", call. = FALSE)
  }
  pb <- if (bounds) panel_b_b else panel_b
  ep <- if (bounds) extra_panels_b else extra_panels
  w <- if (comparators) 3L else 1L     # narrow columns per map column
  map_col <- function(c) (c - 1L) * w + 2L            # first narrow column of map column c
  region_areas <- do.call(c, lapply(seq_len(nrow(geo_grid)), function(i) {
    area(geo_grid$row[i], map_col(geo_grid$col[i]),
         geo_grid$row[i], map_col(geo_grid$col[i]) + w - 1L)
  }))
  # In the layout with the national inset the array of regions is panel B,
  # and its letter goes on the first panel read, top left of the top row
  # and right against the inset's corner. In the other layout the legend
  # holds the top left corner and carries A for the figure as a whole.
  first_region <- geo_grid$region[order(geo_grid$row, geo_grid$col)][1]
  panels <- lapply(geo_grid$region, function(rg) {
    region_panel(rg, keys = keys, show_daily = show_daily,
                 show_bands = show_bands,
                 tag = if (national && identical(rg, first_region)) "B")
  })
  n_row <- if (atn_comparators) 6L else if (comparators) 5L else 4L
  last_col <- map_col(max(geo_grid$col)) + w - 1L
  fixed <- list(y_title)
  if (national) {
    # With the comparator row the legend sits north of Gao and east of Kidal, the cell the
    # map leaves empty; without it, in the bottom-left cell south of Kayes as before.
    # With the comparator row the legend sat in the north-east cell until 2026-10-27, when
    # it moved to the one south of Kayes and west of Sikasso, next to the panels that use
    # it rather than diagonally opposite them (user). Without that row, the bottom-left
    # cell south of Kayes as before.
    legend_area <- if (comparators)
      area(4, map_col(1), 4, map_col(1) + w - 1L) else area(4, 2)
    # The y title spans the FOUR map rows, not every row: it reads PfPR2-10, which is what
    # A and B show, while C-I are cases averted and carry their own x labels. Spanning
    # n_row centred it over the comparator rows as well, so with six rows it sat beside
    # Sikasso rather than beside Kayes (user, 2026-10-26).
    design <- c(area(1, 1, 4, 1),                        # y title
                area(1, 2, 2, map_col(2) + w - 1L),      # national inset, tag A
                legend_area,
                area(3, map_col(4), 4, last_col),        # panel C
                region_areas)
    plots <- c(fixed,
               list(wrap_elements(full = scaled(national_panel(
                      nat_keys, show_daily = show_daily && !nat_bands,
                      show_bands = nat_bands))),
                    wrap_elements(full = legend_panel(keys, legend_title,
                                                      show_bands, tag = NULL,
                                                      scale = if (comparators) 0.92 else 1)),
                    wrap_elements(full = scaled(pb + labs(tag = "C")))),
               panels)
    if (comparators) {
      # D, E, F across the bottom as ONE cell holding a row of three, not three cells of
      # five columns each. Five columns each is equal in the grid but not on the page: the
      # map panels above add fixed-width axis decorations into the same columns, so the
      # three spans came out unequal and the last was visibly the narrowest. A nested
      # wrap_plots divides its own width into exact thirds, whatever sits above it.
      design <- c(design, area(5, 2, 5, last_col))
      plots <- c(plots, list(wrap_elements(full = wrap_plots(
        lapply(ep, scaled), nrow = 1))))
      # The marker key goes in the cell south of Segou and east of Sikasso, which the map
      # leaves empty (user, 2026-10-27). Appended to design and plots together, so the two
      # stay matched position by position.
      if (bounds) {
        design <- c(design, area(4, map_col(3), 4, map_col(3) + w - 1L))
        plots  <- c(plots, list(wrap_elements(full = shape_legend_panel())))
      }
    }
    if (atn_comparators) {
      # G, H, I in a sixth row, nested the same way so their thirds match D-F's.
      design <- c(design, area(6, 2, 6, last_col))
      plots <- c(plots, list(wrap_elements(full = wrap_plots(
        lapply(atn_panels, scaled), nrow = 1))))
    }
  } else {
    design <- c(area(1, 1, 4, 1),   # y title
                area(1, 2, 2, 3),   # legend, tag A
                area(3, 5, 4, 6),   # panel B
                region_areas)
    plots <- c(fixed,
               list(wrap_elements(full = legend_panel(keys, legend_title,
                                                      show_bands)),
                    wrap_elements(full = scaled(pb))),
               panels)
  }
  # The y title keeps the same share of the width whatever the column grid, and the
  # comparator row is twice a map row's height so A, C and D-F come out a similar size.
  wrap_plots(plots, design = design,
             widths = c(0.12 * w, rep(1, last_col - 1L)),
             heights = if (atn_comparators) c(rep(1, 4), 2, 2) else
                       if (comparators) c(rep(1, 4), 2) else
                       rep(1, 4))
}

# BOTH forms go in the file name, as <drug>drug_<ret>ret -- scdrug_scret, scdrug_expret,
# expdrug_expret (user, 2026-10-26). The antimalarial's decay form alone is not enough:
# since v4 the net RETENTION can be smooth-compact too, and a run with both would otherwise
# carry exactly the same name as one with smooth-compact decay and exponential retention,
# told apart only by its folder -- too easy to confuse once the files are out of it. A run
# built before either parameter existed used the exponential.
short_form <- function(x) if (is.null(x) || identical(x, "exp")) "exp" else
  if (identical(x, "smoothcompact")) "sc" else x
form_tag <- paste0(short_form(settings$atn_decay_form), "drug_",
                   short_form(settings$retention_form), "ret")

# A figure drawn from a run that is still filling in says so in its name, as _n<used>, so a
# draft is never mistaken for the final version once it is out of this folder (user,
# 2026-10-26). Nothing is appended when every realisation the run holds was used, so the
# finished figures keep clean names.
# Counted from the CACHE, not from the sweep files on disk. Those two disagree whenever a
# run finishes between one build and the next: the cache is only re-read when the folders,
# the cycle or "refresh" change, so on 2026-10-26 a figure drawn from 152 cached
# realisations was named as though it held all 160, which is the marker asserting exactly
# the thing it exists to deny. What the figure is drawn from is what its name must report.
local({
  v <- cached$values
  in_cache <- unique(v$draw[v$region != country])
  nine <- table(v$draw[v$region != country & v$arm == "none"])
  used <- as.integer(names(nine)[nine == length(geo_grid$region)])
  if (length(used) < length(in_cache)) {
    form_tag <<- paste0(form_tag, "_n", length(used))
    message(sprintf("PARTIAL run: %d of %d realisation(s) in the cache are complete in all ",
                    length(used), length(in_cache)),
            sprintf("%d regions; the file names carry _n%d.",
                    length(geo_grid$region), length(used)))
  }
})

save_figure <- function(fig, name, height = 9) {
  if (only_bounds && !grepl("_bounds$", name)) return(invisible(NULL))
  stem <- file.path(out_dir, paste0(name, "_", form_tag))
  ggsave(paste0(stem, ".png"), fig, width = 13, height = height, dpi = 300,
         device = ragg::agg_png, bg = "white")
  ggsave(paste0(stem, ".pdf"), fig, width = 13, height = height,
         device = grDevices::cairo_pdf)
  message("Saved: ", stem, ".(png|pdf)")
}

# The two layouts, each in its two versions. Neither legend is titled:
# what the bands are is left to the caption.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE),
            "fig5_medians")
save_figure(build_figure(c(band_arms, "shared"), NULL,
                         show_daily = FALSE, show_bands = TRUE),
            "fig5_bands")
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE),
            "fig5_natmed")
save_figure(build_figure(c(band_arms, "shared"), NULL,
                         show_daily = FALSE, show_bands = TRUE,
                         national = TRUE),
            "fig5_natbands")
# The medians figure, but with the national inset banded: the regions by
# their medians, all seven arms, and the country with its 95% intervals
# as well. The inset keeps five of the arms, the four the banded figure
# shows and no future nets, whose band sits clear above the rest; all
# seven would lay their bands over each other. The legend is the medians
# one, so what the shading is belongs in the caption.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE,
                         nat_bands = TRUE,
                         nat_keys = c("none", band_arms, "shared")),
            "fig5_natmed_bandinset")
# The same again without Pyr-CFP-ATN in the inset, whose band runs under
# the others for the whole period.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE,
                         nat_bands = TRUE,
                         nat_keys = c("none",
                                      setdiff(band_arms, "pyr_cfp_atn"),
                                      "shared")),
            "fig5_natmed_bandinset_v2")

# The same again with the comparator row: A-C as before, plus D-F giving cases averted
# against pyrethroid-only, pyrethroid-PBO and pyrethroid-CFP nets rather than against no
# future nets. Saved as its own figure so the six above are unchanged.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE,
                         comparators = TRUE),
            "fig5_natmed_vscomp", height = 13.5)

# The same layout with panel A showing its 95% intervals rather than medians alone, as the
# bandinset figures do: the country has room for bands where nine small panels do not.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE, comparators = TRUE,
                         nat_bands = TRUE,
                         nat_keys = c("none",
                                      setdiff(band_arms, "pyr_cfp_atn"),
                                      "shared")),
            "fig5_natband_vscomp", height = 13.5)

# A-F again with the pessimistic and optimistic rows in C-F. Written only when the bounds
# folders were named, under their own names, so every figure above is untouched.
if (length(BOUNDS) == 2L) {
  save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                           show_bands = FALSE, national = TRUE,
                           comparators = TRUE, bounds = TRUE),
              "fig5_natmed_vscomp_bounds", height = 13.5)
  save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                           show_bands = FALSE, national = TRUE, comparators = TRUE,
                           bounds = TRUE, nat_bands = TRUE,
                           nat_keys = c("none",
                                        setdiff(band_arms, "pyr_cfp_atn"),
                                        "shared")),
              "fig5_natband_vscomp_bounds", height = 13.5)
}

# A-I: the same again with a sixth row, G-I, giving cases averted against each of the three
# ATN arms. Saved separately so the A-F figures above are unchanged. Taller by one 2-unit
# row, so 18 inches against 13.5.
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE,
                         comparators = TRUE, atn_comparators = TRUE),
            "fig5_natmed_vsall", height = 18)
save_figure(build_figure(c(arms$arm, "shared"), NULL, show_daily = TRUE,
                         show_bands = FALSE, national = TRUE,
                         comparators = TRUE, atn_comparators = TRUE,
                         nat_bands = TRUE,
                         nat_keys = c("none",
                                      setdiff(band_arms, "pyr_cfp_atn"),
                                      "shared")),
            "fig5_natband_vsall", height = 18)
for (nm in names(atn_cmp)) {
  a <- atn_averted[[nm]]
  message(sprintf("  vs %-16s %s", atn_cmp[[nm]],
                  paste(sprintf("%s %+.0f [%+.0f, %+.0f]", a$arm, a$median, a$lower, a$upper),
                        collapse = " | ")))
}
for (nm in names(extra_cmp)) {
  a <- extra_averted[[nm]]
  message(sprintf("  vs %-16s %s", extra_cmp[[nm]],
                  paste(sprintf("%s %+.0f [%+.0f, %+.0f]", a$arm, a$median, a$lower, a$upper),
                        collapse = " | ")))
}
