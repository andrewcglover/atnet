# Indicative 3x3 grid: clinical cases averted by each net type, for nine
# hypothetical settings -- EIR 5/50/500 x pyrethroid resistance 50/80/95% --
# over the quoted campaign cycle (three years from the second arm campaign).
#
# Each panel is drawn exactly as panels C-F of fig. 5: one row per arm, the
# dot and bar being the median and 95% interval of the paired within-
# realisation contrast against the comparator, and the three antimalarial
# arms carrying three markers rather than one -- the primary analysis, and
# the same quantity with the antimalarial's half-life, omega and p_A all
# held at the ends of their 95% prior ranges. Those two are NOT a wider
# interval round the primary: each carries its own interval, and the primary
# already propagates the uncertainty they replace with a fixed assumption.
#
# Four figures, one per comparator:
#   grid_vs_none.(png|pdf)     against no future nets
#   grid_vs_pyr.(png|pdf)      against Pyr-only
#   grid_vs_pyr_pbo.(png|pdf)  against Pyr-PBO
#   grid_vs_pyr_cfp.(png|pdf)  against Pyr-CFP
#
# LAYOUT. Columns are transmission and rows are resistance, not the other way
# round. The x range is set mostly by transmission, so with transmission down
# the columns each column spans a similar range and the shared x axis does not
# squash the low-EIR cells against the zero line (user, 2026-11-26).
#
# The x axis is shared by all nine panels of a figure, so the transmission and
# resistance gradients can be read straight off the page. Within a figure the
# cells really do differ several-fold, and that is the point of the grid.
#
# P(averts more) is the share of realisations whose paired contrast is
# positive, printed above its own marker. It is exactly 0 for a comparator
# against itself, which is meaningless, so that row is left blank.
#
# NO NATIONAL ROW anywhere: the nine cells are independent hypothetical
# settings, not the regions of a country. See analysis/09g_summarise_generic.R.
#
# Reads the summaries only; it runs no simulations, and runs in seconds.
# Needs analysis/09g_summarise_generic.R to have been run on all three
# configurations first. Run with the working directory at the root of the
# package:
#   Rscript analysis/15_plot_grid_contrasts.R

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

PRIMARY <- "gen_hldrawn_gprimary"
BOUNDS  <- c(pessimistic = "gen_hl2_gmin", optimistic = "gen_hl5_gmax")
OUT_DIR <- file.path("outputs", "figures", "gen_grid")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- arms, colours and the three variants, as in analysis/10_plot_fig5.R ----
arms <- projection_arms()
arms$label[arms$arm == "pyr"] <- "Pyr-only"
arms$colour <- c("#BBBBBB", "#EE7733", "#CC3311", "#EE3377", "#009988",
                 "#0077BB", "#33BBEE")
key_colours <- setNames(arms$colour, arms$arm)

VAR_LEV    <- c("pessimistic", "primary", "optimistic")
VAR_LABS   <- c(pessimistic = "Minimum plausible", primary = "Primary analysis",
                optimistic = "Maximum plausible")
VAR_SHAPES <- c(pessimistic = 25, primary = 16, optimistic = 17)  # down, dot, up
VAR_SIZES  <- c(pessimistic = 1.7, primary = 2.0, optimistic = 1.7)
VAR_ALPHA  <- c(pessimistic = 0.5, primary = 1, optimistic = 0.5)
# Wider than fig. 5's 0.29: the probability sits above each marker here, so the
# three variants need room for a line of text between them (user, 2026-11-26).
VAR_OFFSET <- 0.34
ATN_ARMS   <- c("atn", "pyr_atn", "pyr_cfp_atn")

# ---- read the three configurations ------------------------------------------
read_values <- function(dir) {
  f <- file.path("outputs", "summaries", dir, "row_values.csv")
  if (!file.exists(f)) {
    stop("No summary for ", dir, ": run analysis/09g_summarise_generic.R on that ",
         "run first (it writes ", f, ").", call. = FALSE)
  }
  utils::read.csv(f, stringsAsFactors = FALSE)
}

primary <- read_values(PRIMARY)[c("region", "draw", "arm", "clin_cases",
                                  "person_years")]
primary$variant <- "primary"

v <- primary
for (vn in names(BOUNDS)) {
  d <- read_values(BOUNDS[[vn]])
  # The bounds runs were summarised WITH the primary folder appended, so their
  # row_values.csv also holds the four pyrethroid arms, copied from primary.
  # Keep only the three antimalarial arms: fixing the antimalarial's half-life,
  # omega and p_A cannot move an arm that carries no antimalarial, so those
  # copies would be the primary's own numbers wearing a bounds label.
  d <- d[d$arm %in% ATN_ARMS, c("region", "draw", "arm", "clin_cases",
                                "person_years")]
  if (!nrow(d)) stop("No antimalarial arms in ", BOUNDS[[vn]], call. = FALSE)
  # Contrasts are paired WITHIN a realisation, so a bounds run summarised over
  # a different set of rows would pair the wrong simulations together.
  if (!setequal(d$draw, primary$draw)) {
    stop(BOUNDS[[vn]], " covers ", length(unique(d$draw)), " realisation(s) ",
         "against ", length(unique(primary$draw)), " in ", PRIMARY,
         "; summarise both over the same rows.", call. = FALSE)
  }
  if (!setequal(d$region, primary$region)) {
    stop(BOUNDS[[vn]], " covers different cells from ", PRIMARY, ".", call. = FALSE)
  }
  d$variant <- vn
  v <- rbind(v, d)
}

# Each variant becomes its own "arm", so that arm_contrasts() pairs it against
# the comparator within the realisation. The comparator itself exists only in
# the primary, which is what every variant should be measured against.
v$base <- v$arm
v$arm  <- ifelse(v$variant == "primary", v$base, paste0(v$base, "__", v$variant))
bkey   <- unique(v[c("arm", "base", "variant")])

# ---- the nine cells ----------------------------------------------------------
# Region names are built by analysis/07g_build_generic_site.R as E<eir>R<res>,
# e.g. E050R80 = EIR 50, 80% pyrethroid resistance.
cells <- unique(v$region)
stopifnot(length(cells) == 9L, all(grepl("^E[0-9]{3}R[0-9]{2}$", cells)))
cell_eir <- as.integer(substr(cells, 2, 4))
cell_res <- as.integer(substr(cells, 6, 7))

EIR_LEV <- sort(unique(cell_eir))
RES_LEV <- sort(unique(cell_res))
EIR_LAB <- setNames(sprintf("%s transmission  (EIR %d)",
                            c("Low", "Moderate", "High"), EIR_LEV), EIR_LEV)
RES_LAB <- setNames(sprintf("%d%% pyrethroid resistance", RES_LEV), RES_LEV)

cell_facets <- data.frame(
  region = cells,
  eir = factor(EIR_LAB[as.character(cell_eir)], levels = EIR_LAB),
  res = factor(RES_LAB[as.character(cell_res)], levels = RES_LAB),
  stringsAsFactors = FALSE)

# ---- one comparator's data ---------------------------------------------------
averted_vs <- function(cmp) {
  raw <- arm_contrasts(v, comparator = cmp)
  a   <- summarise_rows(raw, "cases_averted_per1000", digits = 0)
  a   <- merge(a, bkey, by = "arm")
  # "No future nets" never appears as an arm: against it the figure IS the
  # whole comparison, and against any net it would stretch the shared axis far
  # enough to squash the contrasts these panels exist to show.
  a   <- a[a$base %in% setdiff(arms$arm, "none"), , drop = FALSE]

  # arm_contrasts() returns no line for the comparator against itself, so add
  # the zero it stands for, in every cell: without it the row is an empty slot
  # rather than a dot on the line, and the reference has to be worked out.
  if (cmp %in% setdiff(arms$arm, "none")) {
    z <- do.call(rbind, lapply(cells, function(rg) {
      data.frame(arm = cmp, region = rg, quantity = "cases_averted_per1000",
                 n_rows = length(unique(v$draw)), median = 0, lower = 0,
                 upper = 0, text = "0", base = cmp, variant = "primary",
                 stringsAsFactors = FALSE)
    }))
    a <- rbind(a, z[names(a)])
  }

  # Posterior probability that this arm averts MORE cases than the comparator.
  # PER CELL AND ARM: fig. 5 pools with tapply() over the arm alone, which is
  # right there because it has a single National line, but here it would
  # silently average the nine independent settings into one number.
  raw$better <- raw$cases_averted_per1000 > 0
  pb <- stats::aggregate(better ~ region + arm, data = raw, FUN = mean)
  a$p_better <- pb$better[match(paste(a$region, a$arm),
                                paste(pb$region, pb$arm))]
  a$p_better[a$arm == cmp] <- NA_real_

  a$label <- factor(arms$label[match(a$base, arms$arm)],
                    levels = rev(arms$label[match(setdiff(arms$arm, "none"),
                                                  arms$arm)]))
  a$variant <- factor(a$variant, levels = VAR_LEV)
  # y set by hand rather than by position_dodge, so a one-variant arm sits
  # exactly on its label and a three-variant arm splits evenly, however dodge
  # happens to count the groups present.
  nv <- stats::ave(rep(1L, nrow(a)), paste(a$region, a$label), FUN = length)
  a$y <- as.integer(a$label) +
    ifelse(nv == 1L, 0, (as.integer(a$variant) - 2L) * VAR_OFFSET)
  merge(a, cell_facets, by = "region")
}

# ---- the figure --------------------------------------------------------------
grid_figure <- function(cmp) {
  a <- averted_vs(cmp)
  lev <- levels(a$label)
  cmp_label <- if (cmp == "none") "no future nets" else
    arms$label[match(cmp, arms$arm)]

  a$shown <- ifelse(is.na(a$p_better), NA_character_,
                    sprintf("%.0f%%", 100 * a$p_better))

  ggplot(a, aes(median, y, colour = base)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_linerange(aes(xmin = lower, xmax = upper, alpha = variant),
                   linewidth = 0.7) +
    geom_point(aes(shape = variant, size = variant, fill = base), alpha = 1) +
    # The percentage sits above its own marker rather than in a column down the
    # right-hand side: at 3x3 a column would be 54 numbers competing with the
    # panels (user, 2026-11-26).
    geom_text(aes(label = shown), nudge_y = 0.155, size = 2.15, colour = "grey25",
              na.rm = TRUE, show.legend = FALSE) +
    facet_grid(res ~ eir) +
    scale_colour_manual(values = key_colours, guide = "none") +
    scale_fill_manual(values = key_colours, guide = "none") +
    scale_shape_manual(values = VAR_SHAPES, labels = VAR_LABS,
                       breaks = VAR_LEV, name = NULL) +
    scale_size_manual(values = VAR_SIZES, guide = "none") +
    scale_alpha_manual(values = VAR_ALPHA, guide = "none") +
    scale_y_continuous(breaks = seq_along(lev), labels = lev,
                       limits = c(0.35, length(lev) + 0.85)) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
    guides(shape = guide_legend(override.aes = list(
      colour = "grey25", fill = "grey25", size = 2.4, alpha = 1))) +
    labs(x = sprintf("Clinical cases averted per 1,000 people per year, against %s",
                     cmp_label),
         y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(10, "pt"),
          panel.spacing.y = unit(8, "pt"),
          axis.text.y = element_text(size = 8.5, colour = "black"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10, margin = margin(t = 8)),
          strip.text = element_text(size = 9.5, face = "bold"),
          strip.background = element_rect(fill = "grey94", colour = NA),
          legend.position = "bottom",
          legend.text = element_text(size = 9),
          plot.margin = margin(10, 14, 8, 8))
}

for (cmp in c("none", "pyr", "pyr_pbo", "pyr_cfp")) {
  p <- grid_figure(cmp)
  stem <- file.path(OUT_DIR, paste0("grid_vs_", cmp))
  ggsave(paste0(stem, ".png"), p, width = 11.5, height = 10.6, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 11.5, height = 10.6)
  message("Written ", stem, ".png and .pdf")
}

message("\nAll four figures in ", OUT_DIR)
