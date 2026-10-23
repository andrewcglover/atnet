# FAKE all-region projections, made from the dummy Mopti run, for developing
# maps and regional plots while the real runs are on the cluster. They are
# NOT model output: never use them for any number or interpretation, and
# delete them (and this script) once the real all-region runs arrive.
#
# Every region is a copy of Mopti's job files, row by row, scaled by the
# region's baseline EIR relative to Mopti's (the site file's `$eir`):
#   - EIR and mosquito numbers by the ratio itself;
#   - prevalence counts and clinical and severe cases by the ratio to the
#     power 0.3, a damped ratio, since prevalence rises far less than
#     linearly with EIR; the power gives roughly the ranges of the June
#     figure.
# Population, net use and the share of mosquitoes exposed to ATN are left as
# they are. Every region therefore has Mopti's relative effects (the same %
# of cases averted), and only absolute levels differ.
#
# Written to outputs/fake/FAKE_mli_hldrawn_dummy/: settings.rds and
# jobs/FAKE_mli_<region>_row<i>.rds, all carrying fake = TRUE. They sit
# outside outputs/projections and are named FAKE so they cannot be mistaken
# for real runs, and the loader warns whenever it reads them.
#
# Reads and writes files only; it runs no simulations.
#   Rscript analysis/dev_make_FAKE_all_regions.R
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)

src_dir <- "outputs/projections/mli_hldrawn_dummy"
src_region <- "Mopti"
out_dir <- "outputs/fake/FAKE_mli_hldrawn_dummy"
damping <- 0.3

settings <- readRDS(file.path(src_dir, "settings.rds"))
site <- readRDS(file.path("data_private", "site_files",
                          paste0(settings$country_iso, ".rds")))
eir <- site$eir[site$eir$sp == "pf", c("name_1", "eir")]
eir$ratio <- eir$eir / eir$eir[eir$name_1 == src_region]
eir$damped <- eir$ratio^damping

fake_region <- function(x, region, ratio, damped) {
  by_ratio <- grep("^EIR_|^total_M_|^M_exposed_", names(x), value = TRUE)
  by_damped <- grep("^n_detect_lm_|^n_inc_clinical_|^n_inc_severe_",
                    names(x), value = TRUE)
  x[by_ratio] <- x[by_ratio] * ratio
  x[by_damped] <- x[by_damped] * damped
  x$pfpr2to10 <- x$n_detect_lm_730_3649 / x$n_age_730_3649
  x$clin_inc <- rowSums(x[grep("^n_inc_clinical_", names(x))])
  x$region <- region
  x
}

dir.create(file.path(out_dir, "jobs"), recursive = TRUE, showWarnings = FALSE)
jobs <- projection_job_files(src_dir, regions = src_region)
for (i in seq_len(nrow(eir))) {
  region <- eir$name_1[i]
  for (f in jobs$file) {
    job <- readRDS(f)
    job$results <- fake_region(job$results, region, eir$ratio[i], eir$damped[i])
    job$pre_projection <- fake_region(job$pre_projection, region, eir$ratio[i],
                                      eir$damped[i])
    job$region <- region
    job$fake <- TRUE
    job$fake_note <- sprintf(
      "FAKE: %s row %d scaled to %s (EIR ratio %.4f; prevalence and cases x %.4f)",
      src_region, job$draw, region, eir$ratio[i], eir$damped[i])
    saveRDS(job, file.path(out_dir, "jobs",
                           sprintf("FAKE_%s_%s_row%04d.rds",
                                   tolower(settings$country_iso),
                                   region_tag(region), job$draw)))
  }
}
settings$fake <- TRUE
settings$fake_note <- sprintf(
  "FAKE: every region is %s from %s, scaled by EIR ratio (EIR, mosquitoes) and its power %.1f (prevalence, cases).",
  src_region, src_dir, damping)
saveRDS(settings, file.path(out_dir, "settings.rds"))

eir$file_tag <- region_tag(eir$name_1)
print(eir, row.names = FALSE, digits = 3)
cat(nrow(eir) * nrow(jobs), "FAKE job files written to", out_dir, "\n")
