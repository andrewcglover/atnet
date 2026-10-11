# =====================================================================
# 00_download_site_files.R
# Download all sub-Saharan African site files at admin level 1, in both
# the urban/rural-split and non-split versions, into data_private/site_files/.
# Run with the working directory at the root of this package. The projections
# read data_private/site_files/<ISO>.rds, a copy of the country's file from
# without_split/ (this project uses no urban/rural split).
# Moved from the malariasimulation fork's dev/download_ssa_site_files.R.
#
# Prereqs (already done in your session): on the malariaverse GitHub team,
# `orderly` installed, and authenticated once via a manual fetch_site().
#
# Versions are PER-product/PER-country (the split and non-split sets are
# different releases), so we resolve each from available_sites() rather than
# hard-coding a tag. We take the latest packet per combo (max id = newest).
# =====================================================================

library(site)

dir_split   <- "data_private/site_files/with_split"     # urban_rural = TRUE
dir_nosplit <- "data_private/site_files/without_split"  # urban_rural = FALSE
dir.create(dir_split,   recursive = TRUE, showWarnings = FALSE)
dir.create(dir_nosplit, recursive = TRUE, showWarnings = FALSE)

# FAILSAFE: with overwrite = FALSE (default), existing site files are KEPT (not
# re-downloaded) — a re-run only fills in any missing ones and never clobbers what
# you already have. Set overwrite = TRUE only to deliberately refresh to newer versions.
overwrite <- FALSE

# Sub-Saharan African ISO3c codes. Any not on the server are skipped.
ssa_iso3c <- c(
  # West Africa
  "BEN","BFA","CPV","CIV","GMB","GHA","GIN","GNB","LBR","MLI","MRT","NER","NGA","SEN","SLE","TGO",
  # Central Africa
  "AGO","CMR","CAF","TCD","COG","COD","GNQ","GAB","STP",
  # East Africa
  "BDI","COM","DJI","ERI","ETH","KEN","MDG","MWI","MOZ","RWA","SOM","SSD","SDN","TZA","UGA","ZMB","ZWE",
  # Southern Africa
  "BWA","SWZ","LSO","NAM","ZAF"
)

# ---- one metadata pull; resolve the latest version per combo ----
av <- site::available_sites()

# version of the most recent packet for (iso3c, admin_level, urban_rural).
# version is a list-column; id is a chronological string, so max(id) = newest.
resolve_version <- function(av, iso, adm, ur) {
  rows <- subset(av, iso3c == iso & admin_level == adm & urban_rural == ur)
  if (nrow(rows) == 0) return(NA_character_)
  vers <- unlist(rows$version)
  ids  <- as.character(rows$id)
  vers[match(max(ids), ids)]
}

fetch_one <- function(iso, ur, outdir) {
  tag <- if (ur) "split" else "no-split"
  ver <- resolve_version(av, iso, 1L, ur)
  if (is.na(ver)) {
    message(sprintf("  --   %s (%s): no site file on server", iso, tag))
    return(data.frame(iso3c = iso, urban_rural = ur, version = NA_character_, status = "absent"))
  }
  out <- file.path(outdir, paste0(iso, ".rds"))
  if (file.exists(out) && !overwrite) {        # FAILSAFE: keep existing, don't clobber
    message(sprintf("  KEEP %s (%s)  [already present; set overwrite = TRUE to refresh]", iso, tag))
    return(data.frame(iso3c = iso, urban_rural = ur, version = ver, status = "kept"))
  }
  status <- tryCatch({
    s <- site::fetch_site(iso, admin_level = 1, urban_rural = ur, version = ver)
    saveRDS(s, out)            # NB assumes fetch_site() returns the site object
    "ok"
  }, error = function(e) paste0("error: ", conditionMessage(e)))
  message(sprintf("  %-4s %s (%s)  [%s]", if (status == "ok") "OK" else "SKIP", iso, tag, ver))
  data.frame(iso3c = iso, urban_rural = ur, version = ver, status = status)
}

# ---- loop countries x {split, non-split}, collect a provenance manifest ----
existing <- length(list.files(dir_split, "\\.rds$")) + length(list.files(dir_nosplit, "\\.rds$"))
if (existing > 0 && !overwrite)
  message(sprintf("Note: %d site file(s) already present — existing files will be KEPT (overwrite = FALSE).\n",
                  existing))

manifest <- do.call(rbind, lapply(ssa_iso3c, function(iso) {
  message(iso, " ...")
  rbind(
    fetch_one(iso, TRUE,  dir_split),
    fetch_one(iso, FALSE, dir_nosplit)
  )
}))

# The manifest records exactly which version of each file was downloaded. It
# sits beside the files in data_private/, which git ignores.
write.csv(manifest, "data_private/site_files/site_files_manifest.csv", row.names = FALSE)

n_ok   <- sum(manifest$status == "ok")
n_kept <- sum(manifest$status == "kept")
message(sprintf("\nDone: %d downloaded, %d kept (already present), %d rows total. Manifest: data_private/site_files/site_files_manifest.csv",
                n_ok, n_kept, nrow(manifest)))
message("Files in:\n  ", dir_split, "\n  ", dir_nosplit)
