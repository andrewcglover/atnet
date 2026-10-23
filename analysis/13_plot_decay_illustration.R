# Two illustrative figures, for explaining the antimalarial's decay:
#
#   decay_forms.(png|pdf)    the antimalarial effect P(a) left on a net
#                            of age a, under each decay form, at four
#                            effective half-lives
#   half_life_prior.(png|pdf) the half-life the projections draw from
#
# P(a) is the share of the antimalarial's effect that is left: every
# one of its three mechanisms is multiplied by it, so P = 1 is a net
# giving its full effect and P = 0 one giving none. It is not an
# efficacy against malaria. What is drawn is the half-life of that
# effect -- the "effective antimalarial half-life" -- and not the
# half-life of the chemical on the net. The fork calls P potency, so
# the function below keeps that name; the figures do not.
#
# The two forms are the fork's `atn_potency()` (R/biting_process.R),
# repeated here so the figure can be drawn without the fork loaded, and
# checked below against the one thing both must do: hand back a half at
# the half-life.
#
#   Rscript analysis/13_plot_decay_illustration.R
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages(library(ggplot2))

out_dir <- file.path("outputs", "figures", "illustrative")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zero_years <- 20   # atn_decay_zero_days, the fork's default

potency <- function(age, half_life, form, l = zero_years) {
  if (form == "exp") return(exp(-log(2) / half_life * age))
  k <- log(2) * ((l / half_life)^2 - 1)
  x <- pmin(age / l, 1)
  ifelse(age >= l, 0, exp(-k * x^2 / (1 - x^2)))
}

# Both forms take the same half-life, which is what lets one prior serve
# for either: at that age each is left with half of its effect.
for (f in c("exp", "smoothcompact")) {
  for (h in c(2, 3, 4, 5)) {
    stopifnot(isTRUE(all.equal(potency(h, h, f), 0.5)))
  }
}

half_lives <- c(2, 3, 4, 5)
forms <- c(exp = "Exponential", smoothcompact = "Smooth-compact")
age <- seq(0, 6, by = 0.01)
curves <- do.call(rbind, lapply(names(forms), function(f) {
  do.call(rbind, lapply(half_lives, function(h) {
    data.frame(age = age, potency = potency(age, h, f),
               half_life = factor(paste(h, "y"),
                                  levels = paste(half_lives, "y")),
               form = factor(forms[[f]], levels = unname(forms)))
  }))
}))

ggsave(
  file.path(out_dir, "decay_forms.png"),
  ggplot(curves, aes(age, potency, colour = half_life)) +
    geom_hline(yintercept = 0.5, colour = "grey60", linetype = "dashed",
               linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ form) +
    scale_colour_viridis_d(end = 0.9, name = "Effective half-life") +
    scale_x_continuous(breaks = 0:6) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Net age (years)", y = "Antimalarial effect remaining") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(size = 13)),
  width = 9, height = 4.2, dpi = 300, device = ragg::agg_png, bg = "white")

ggsave(file.path(out_dir, "decay_forms.pdf"), last_plot(),
       width = 9, height = 4.2, device = grDevices::cairo_pdf)
message("Saved: ", file.path(out_dir, "decay_forms.(png|pdf)"))

# The half-life the projections draw from, from the package's own prior
# so the figure cannot drift from the runs.
hl <- projection_priors()$half_life_atn_years
q <- stats::qlnorm(c(0.025, 0.5, 0.975), hl$meanlog, hl$sdlog)
message(sprintf("Effective antimalarial half-life: median %.2f y, 95%% %.2f to %.2f y", q[2], q[1], q[3]))

x <- seq(1, 8, by = 0.005)
dens <- data.frame(x = x, d = stats::dlnorm(x, hl$meanlog, hl$sdlog))

ggsave(
  file.path(out_dir, "half_life_prior.png"),
  ggplot(dens, aes(x, d)) +
    geom_area(data = subset(dens, x >= q[1] & x <= q[3]), fill = "#0077BB",
              alpha = 0.2) +
    geom_line(linewidth = 0.9, colour = "#0077BB") +
    geom_vline(xintercept = q[2], colour = "#0077BB", linetype = "dashed",
               linewidth = 0.5) +
    scale_x_continuous(breaks = 1:8) +
    labs(x = "Effective antimalarial half-life (years)", y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank()),
  width = 7, height = 3.6, dpi = 300, device = ragg::agg_png, bg = "white")

ggsave(file.path(out_dir, "half_life_prior.pdf"), last_plot(),
       width = 7, height = 3.6, device = grDevices::cairo_pdf)
message("Saved: ", file.path(out_dir, "half_life_prior.(png|pdf)"))
