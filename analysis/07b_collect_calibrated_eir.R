# Collect the regional calibrations of 07a_calibrate_eir.R into one table the projections read.
#   Rscript analysis/07b_collect_calibrated_eir.R 20261016_dq30_sl30
# The argument is a folder of outputs/calibration/. Writes inst/extdata/calibrated_eir.csv, which
# 07_run_projections.R uses in place of the site file's baseline EIRs. The table is tracked and
# lives next to projection_draws.csv, not under outputs/, because the cluster runs a frozen copy
# of the code at a commit (08_submit_cluster.R) and so can only read what has been committed.
#
# One row per region. Besides the EIR, each row records everything the calibration depends on
# (compartment counts, net retention, fitting window, and the fork, site and cali versions):
# an EIR is only valid for the model it was fitted with, and 07 refuses to use the table when
# its compartment counts or retention differ from the run's.
#
# Regions are written without accents ("Segou"); 07 matches them by their ASCII form, as it
# does for SWEEP_REGIONS. The script refuses a calibration that did not end in cali's own
# "Success" (the other exit, at the attempt limit, returns the mean of the last two proposals
# rather than an EIR that brackets the target), that mixes configurations, or that misses a
# region of the site file. Run with the working directory at the package root.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Give the calibration folder, e.g. 20261016_dq30_sl30")
dir <- file.path("outputs", "calibration", args[[1]])
fs  <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(fs)) stop("No regional results in ", dir, call. = FALSE)

ascii <- function(x) iconv(x, to = "ASCII//TRANSLIT")
rows <- lapply(fs, function(f) {
  x   <- readRDS(f)
  log <- file.path(dir, paste0("log_", sub("\\.rds$", "", basename(f)), ".log"))
  lg  <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  ended <- if (any(grepl("^Success", lg))) "Success" else
           if (any(grepl("Terminating as max attempts", lg))) "MAX ATTEMPTS" else "unknown"
  data.frame(
    region          = ascii(x$region),
    calibrated_eir  = x$eir,
    site_file_eir   = x$site_file_eir,
    target_mean_pfpr = mean(x$target),
    fit_first_year  = x$window[[1]],
    fit_last_year   = x$window[[2]],
    deltaq          = x$deltaq,
    spor_len        = x$spor_len,
    retention_days  = x$retention_days,
    populations     = paste(x$populations, collapse = ";"),
    attempts        = sum(grepl("^Attempt", lg)),
    ended           = ended,
    minutes         = round(x$minutes, 1),
    cali_version    = x$cali,
    fork_version    = x$fork_version,
    site_version    = unname(x$sim_identity$packages[["site"]]),
    source          = args[[1]],
    stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$region), ]

if (any(tab$ended != "Success")) {
  stop("Not every calibration ended in Success: ",
       paste(tab$region[tab$ended != "Success"], tab$ended[tab$ended != "Success"], sep = "=",
             collapse = ", "), ". Rerun those regions.", call. = FALSE)
}
same <- function(col) length(unique(tab[[col]])) == 1L
for (col in c("fit_first_year", "fit_last_year", "deltaq", "spor_len", "retention_days",
              "cali_version", "fork_version", "site_version")) {
  if (!same(col)) stop("The regions disagree on ", col, ": ", toString(unique(tab[[col]])),
                       call. = FALSE)
}
if (anyDuplicated(tab$region)) stop("A region appears twice.", call. = FALSE)

site_file <- file.path("data_private", "site_files", "MLI.rds")
if (file.exists(site_file)) {
  want <- ascii(readRDS(site_file)$sites$name_1)
  if (!setequal(want, tab$region)) {
    stop("Regions differ from the site file's. Missing: ", toString(setdiff(want, tab$region)),
         "; extra: ", toString(setdiff(tab$region, want)), call. = FALSE)
  }
}

out <- file.path("inst", "extdata", "calibrated_eir.csv")
write.csv(tab, out, row.names = FALSE)
message(sprintf("Wrote %d regions to %s (deltaq %d, spor_len %d, retention %g d, fitted %d-%d)",
                nrow(tab), out, tab$deltaq[1], tab$spor_len[1], tab$retention_days[1],
                tab$fit_first_year[1], tab$fit_last_year[1]))
print(tab[, c("region", "site_file_eir", "calibrated_eir", "attempts", "ended")], row.names = FALSE)
