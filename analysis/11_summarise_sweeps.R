# Cycle values for every run of the sensitivity programme, in one table:
# the input the sensitivity figures are drawn from.
#
# For each run folder it works out, over the quoted campaign cycle, the
# totals and rates for each region, draws-table row and arm, adds the
# national values (each region's counts scaled to its real population
# before they are added up), and tags every line with the run's design
# point read from its own settings.rds:
#   decay_form        exp or smoothcompact
#   hl_years          the antimalarial's half-life, NA when it is drawn
#   retention_days    the arms' net retention (the baseline when the run
#                     does not change it)
#   mechanism         both, blocking or eip
#   omega_q, p_a_q    the quantile the contact probability and p_A were
#                     held at, NA when they are drawn
#
# It does NOT work out cases averted. Most runs hold the three ATN arms
# alone, so their comparator comes from the matching seven-arm run and the
# arms have to be paired by draws-table row across runs; that pairing
# belongs with the figure that needs it, not here.
#
# Written to outputs/summaries/<root folder's name>/:
#   runs/<run>.rds    one run's lines, kept so the pass can be resumed
#   sweep_values.rds  every run's lines, bound together
#
# Reading 109,000 job files takes hours, so each run is saved as it is
# finished and a run already saved is skipped. Deleting a file in runs/ is
# enough to have it done again.
#
# Reads the outputs only; it runs no simulations.
#   Rscript analysis/11_summarise_sweeps.R                  every run
#   Rscript analysis/11_summarise_sweeps.R tw_ret           runs matching
#   Rscript analysis/11_summarise_sweeps.R . outputs/projections_v2
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
pattern <- if (length(args) >= 1) args[1] else "."
root <- if (length(args) >= 2) args[2] else "outputs/projections_v2"

out_dir <- file.path("outputs", "summaries", basename(root))
run_dir <- file.path(out_dir, "runs")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

runs <- list.dirs(root, recursive = FALSE)
runs <- runs[file.exists(file.path(runs, "settings.rds"))]
runs <- runs[grepl(pattern, basename(runs))]
if (!length(runs)) stop("No run folder under ", root, " matches ", pattern,
                        call. = FALSE)

site <- NULL   # read once, from the first run's country
pop <- NULL

# The design point of one run, as the figures need to read it.
design_of <- function(s, run) {
  one <- function(x, empty = NA) if (is.null(x)) empty else unname(x[[1]])
  data.frame(
    run = run,
    decay_form = one(s$atn_decay_form, "exp"),
    hl_years = as.numeric(one(s$hl_fixed_years)),
    # The arms keep the baseline retention unless the run changes it.
    retention_days = as.numeric(one(s$future_retention, one(s$retention_time))),
    mechanism = one(s$atn_mechanism, "both"),
    omega_q = as.numeric(one(s$omega_quantile)),
    p_a_q = as.numeric(one(s$p_a_quantile)),
    stringsAsFactors = FALSE
  )
}

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...))
started <- Sys.time()
done <- 0L; skipped <- 0L

for (r in runs) {
  run <- basename(r)
  file <- file.path(run_dir, paste0(run, ".rds"))
  if (file.exists(file)) { skipped <- skipped + 1L; next }

  s <- readRDS(file.path(r, "settings.rds"))
  if (is.null(site)) {
    site <- readRDS(file.path("data_private", "site_files",
                              paste0(s$country_iso, ".rds")))
    pop <- region_populations(site)
  }
  cycle <- campaign_cycle(s$future_campaign_days)
  if (s$end_day < cycle["to"]) {
    warning(run, " stops on day ", s$end_day, ", before the cycle ends on ",
            cycle["to"], "; left out.", call. = FALSE)
    next
  }

  v <- projection_row_values(r, cycle["from"], cycle["to"])
  v <- v[v$arm %in% projection_arms()$arm, ]
  if (all(region_tag(pop$region) %in% region_tag(v$region))) {
    v <- rbind(v, national_values(v, pop))
  } else {
    message("  ", run, ": not every region, so no national values")
  }
  saveRDS(cbind(design_of(s, run), v[setdiff(names(v), "run")]), file)

  done <- done + 1L
  left <- length(runs) - skipped - done
  per <- as.numeric(difftime(Sys.time(), started, units = "mins")) / done
  say("%-32s %5d rows   %2d done, %2d left, ~%.0f min to go",
      run, nrow(v), done, left, per * left)
}

say("%d run(s) summarised, %d already had a file", done, skipped)

files <- list.files(run_dir, pattern = "\\.rds$", full.names = TRUE)
all <- do.call(rbind, lapply(files, readRDS))
rownames(all) <- NULL
saveRDS(all, file.path(out_dir, "sweep_values.rds"))
say("%d lines from %d run(s) -> %s", nrow(all), length(files),
    file.path(out_dir, "sweep_values.rds"))

cat("\nDesign points held:\n")
pts <- unique(all[c("run", "decay_form", "hl_years", "retention_days",
                    "mechanism", "omega_q", "p_a_q")])
print(pts[order(pts$decay_form, pts$run), ], row.names = FALSE)
