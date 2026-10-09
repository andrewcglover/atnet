# Build the two tables the fits consume, and write them to outputs/.
#
# Each row carries the workbook, worksheet and columns it came from, so any
# number can be traced back to a cell. The corrections applied are those listed
# by assay_corrections(), and nothing else; data_private/DATA_SOURCES.md gives
# the evidence for each.
#
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

dir.create("outputs", showWarnings = FALSE)
stamp <- format(Sys.Date(), "%Y%m%d")

spz <- build_sporozoite_table()
ooc <- build_oocyst_table()

# The pooled workbook is an independent record of the same experiments, so it
# must agree on how many mosquitoes exist once the corrections are applied.
pooled <- find_workbook("data_private/fitting", "pooled")
for (d in c(10, 13)) {
  sheet <- if (d == 10) "#7 72h post spz10" else "#8 72h post spz13"
  from_pooled <- sum(tabulate_assay_workbook(pooled, sheet)$n)
  ours <- with(spz[spz$dissection_day == d & spz$hours_after_infection == 72 &
                     spz$worksheet != "24h, 72h post 9.2.25 spz13", ],
               sum(n_control + n_treated))
  if (from_pooled != ours) {
    stop(sprintf("72h exposures read at day %d: the pooled workbook has %d ",
                 d, from_pooled),
         sprintf("mosquitoes but this table has %d. Re-check against ", ours),
         "data_private/DATA_SOURCES.md before fitting.", call. = FALSE)
  }
  message(sprintf("72h exposures, day %2d: %d mosquitoes, agrees with the pooled workbook",
                  d, ours))
}

spz_path <- file.path("outputs", sprintf("%s_sporozoite_experiments.csv", stamp))
ooc_path <- file.path("outputs", sprintf("%s_oocyst_mosquitoes.csv", stamp))
write.csv(spz, spz_path, row.names = FALSE)
write.csv(ooc, ooc_path, row.names = FALSE)

message(sprintf("\nsporozoite: %d experiments, %d mosquitoes -> %s",
                nrow(spz), sum(spz$n_control + spz$n_treated), spz_path))
message(sprintf("oocyst    : %d experiments, %d mosquitoes -> %s",
                length(unique(ooc$experiment_id)), nrow(ooc), ooc_path))
message("\ncorrections applied:")
for (i in seq_len(nrow(assay_corrections()))) {
  cr <- assay_corrections()[i, ]
  message(sprintf("  %-22s (%s)  %s", cr$id, cr$applies_to, cr$evidence))
}
