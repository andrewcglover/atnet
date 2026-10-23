# Summaries of a projection run, for the figures and the manuscript's
# numbers: totals and rates over the quoted campaign cycle (from the 2028
# campaign to the day before the 2031 campaign) for each region, draws-table
# row and arm; clinical cases averted against no future nets, paired within
# rows; and the median [95% interval] across rows.
#
# National values are added when the run holds every region of the site
# file. Each region's totals are scaled to its real population (the site
# file's mean for 2029 and 2030) before they are added up, so national
# percentages are never averages of regional ones.
#
# Written to outputs/summaries/<run folder name>/:
#   row_values.csv      totals and rates, one line per region, row and arm
#   row_contrasts.csv   cases averted against no future nets, per row
#   summary.csv         median [95% interval] across rows of both
#
# Reads the outputs only; it runs no simulations. Uses the dummy run unless
# another run folder is named:
#   Rscript analysis/09_summarise_projections.R outputs/projections/mli_hldrawn
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args)) args[1] else "outputs/projections/mli_hldrawn_dummy"
settings <- readRDS(file.path(run_dir, "settings.rds"))
out_dir <- file.path("outputs", "summaries", basename(run_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The model counts 365-day years from 1 Jan 2000 (timestep 1).
model_date <- function(t) {
  year <- 2000 + (t - 1) %/% 365
  day <- as.Date("2001-01-01") + (t - 1) %% 365
  paste(format(day, "%d %b"), year)
}

cycle <- campaign_cycle(settings$future_campaign_days)
arms <- projection_arms()
values <- projection_row_values(run_dir, cycle["from"], cycle["to"])
values <- values[values$arm %in% arms$arm, ]

site <- readRDS(file.path("data_private", "site_files",
                          paste0(settings$country_iso, ".rds")))
pop <- region_populations(site)
if (all(region_tag(pop$region) %in% region_tag(values$region))) {
  values <- rbind(values, national_values(values, pop))
} else {
  message("Not every region is in this run, so there are no national values.")
}
contrasts <- arm_contrasts(values, comparator = "none")

digits <- c(pfpr2to10_pct = 1, clin_inc_per1000 = 0, severe_inc_per1000 = 1,
            net_use_pct = 1, annual_eir = 0, atn_exposed_pct = 1,
            cases_averted_per1000 = 0, pct_cases_averted = 1)
summary <- rbind(
  summarise_rows(values, c("pfpr2to10_pct", "clin_inc_per1000",
                           "severe_inc_per1000", "net_use_pct", "annual_eir",
                           "atn_exposed_pct"), digits = digits),
  summarise_rows(contrasts, c("cases_averted_per1000", "pct_cases_averted"),
                 digits = digits)
)
summary$label <- arms$label[match(summary$arm, arms$arm)]
summary <- summary[order(summary$region == "National", summary$region,
                         match(summary$arm, arms$arm),
                         match(summary$quantity, names(digits))), ]

write.csv(values, file.path(out_dir, "row_values.csv"), row.names = FALSE)
write.csv(contrasts, file.path(out_dir, "row_contrasts.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)

cat(sprintf("%s: %d row(s); cycle = days %d-%d (%s to %s)\n", run_dir,
            length(unique(values$draw)), cycle["from"], cycle["to"],
            model_date(cycle["from"]), model_date(cycle["to"])))
cat("Median [95% interval] across rows; cases averted are against no future nets.\n")
for (r in unique(summary$region)) {
  s <- summary[summary$region == r, c("label", "quantity", "text")]
  wide <- reshape(s, idvar = "label", timevar = "quantity", direction = "wide")
  names(wide) <- sub("^text\\.", "", names(wide))
  cat("\n", r, "\n", sep = "")
  print(wide, row.names = FALSE, right = FALSE)
}
cat("\nWritten to", out_dir, "\n")
