# Provenance check for the oocyst data used by the blocking fits.
#
# The blocking fits were previously given mos_level_int_summary_pre_post_full.csv,
# one row per mosquito. This script compares that file with the oocyst table
# rebuilt from the workbooks by build_oocyst_table(), and reports:
#   - whether the same mosquitoes, with the same oocyst counts, appear at each
#     exposure time in each arm;
#   - which experiment of the rebuilt table each old experiment reproduces.
# It writes nothing and changes nothing, and prints no oocyst counts.
#
# The old file coded the six-minute rest on the net as an exposure time of
# 0.1 hours before infection; it was concurrent exposure, so it is read as 0.
#
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

old_path <- "data_private/archive/raw/2026-05_exp_data/mos_level_int_summary_pre_post_full.csv"
old <- read.csv(old_path)
old$hours_before_infection <- ifelse(abs(old$pre_time - 0.1) < 1e-9, 0,
                                     old$pre_time)
old$arm <- ifelse(old$group == "control", "control", "treated")
old$experiment <- paste(old$exp_id, old$pre_time)

new <- build_oocyst_table()

# ---- The same mosquitoes at each exposure time and arm ----
group_key <- function(x) paste(x$hours_before_infection, x$arm)
burdens <- function(x) lapply(split(as.numeric(x$oocysts), group_key(x)), sort)
old_b <- burdens(data.frame(oocysts = old$int,
                            hours_before_infection = old$hours_before_infection,
                            arm = old$arm))
new_b <- burdens(new)

message(sprintf("mosquitoes: old %d, rebuilt %d", nrow(old), nrow(new)))
groups <- union(names(old_b), names(new_b))
same <- vapply(groups, function(g) identical(old_b[[g]], new_b[[g]]), logical(1))
message(sprintf("exposure time x arm groups with identical oocyst counts: %d of %d",
                sum(same), length(same)))
if (!all(same)) {
  message("  differ: ", paste(groups[!same], collapse = "; "))
}

# ---- Which rebuilt experiment each old experiment reproduces ----
# An experiment is identified by the sorted oocyst counts of its two arms.
signature <- function(oocysts, arm) {
  paste(paste(sort(oocysts[arm == "control"]), collapse = ","),
        paste(sort(oocysts[arm == "treated"]), collapse = ","), sep = " | ")
}
old_sig <- vapply(split(old, old$experiment),
                  function(e) signature(e$int, e$arm), character(1))
new_sig <- vapply(split(new, new$experiment_id),
                  function(e) signature(e$oocysts, e$arm), character(1))

matched <- old_sig %in% new_sig
message(sprintf("\nold experiments: %d; rebuilt experiments: %d",
                length(old_sig), length(new_sig)))
message(sprintf("old experiments reproduced exactly by one rebuilt experiment: %d",
                sum(matched)))

# The rest should combine into rebuilt experiments, by the corrections listed in
# assay_corrections().
spare_new <- setdiff(names(new_sig), names(new_sig)[new_sig %in% old_sig])
spare_old <- names(old_sig)[!matched]
for (id in spare_new) {
  e <- new[new$experiment_id == id, ]
  parts <- old[old$experiment %in% spare_old &
                 old$hours_before_infection == e$hours_before_infection[1], ]
  combined <- identical(signature(parts$int, parts$arm), new_sig[[id]])
  message(sprintf(paste0("rebuilt experiment %s (worksheet '%s', %g h after ",
                         "infection) = old experiments %s combined: %s"),
                  id, e$worksheet[1], e$hours_after_infection[1],
                  paste(unique(parts$experiment), collapse = " + "), combined))
}
