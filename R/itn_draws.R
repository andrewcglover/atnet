# The insecticide-treated net (ITN) parameter draws of Churcher et al. (2024)
# used by the projections.

#' Where the ITN parameter draws come from
#'
#' The posterior draws of the ITN parameters of Churcher et al. (2024) for
#' pyrethroid-only, pyrethroid-PBO and pyrethroid-pyrrole (chlorfenapyr) nets,
#' from the `parameters/` folder of the repository below (MIT licence), pinned
#' to one commit. Each file holds 1000 draws at each of 101 levels of
#' pyrethroid resistance, 0 to 1 in steps of 0.01. The MD5 checksums are of
#' the files as downloaded, which matched the repository's own (git blob)
#' hashes.
#'
#' @return A list: `repo`, `commit`, and named vectors `files` and `md5`, the
#'   names being the net types `only`, `pbo` and `cfp`.
#'
#' @export
itn_draws_source <- function() {
  list(
    repo   = "EllieSherrardSmith/Mosq-Net-Efficacy",
    commit = "057ac811879b2a58722c542d535cd829328eafb3",
    files  = c(
      only = "pyrethroid_uncertainty_LancetGH2024",
      pbo  = "pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024",
      cfp  = "pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024"
    ),
    md5 = c(
      only = "912b291e27ad07c35b325912e5059971",
      pbo  = "f8b28649337ed5783bfba7ce1e7af36e",
      cfp  = "baf1ab0fbad49f12e53b67c5278d115c"
    )
  )
}

#' Read the ITN parameter draws, downloading them first if needed
#'
#' Any of the three files in [itn_draws_source()] missing from `dir` is
#' downloaded from the pinned commit; all three are then checked against
#' their checksums before being read.
#'
#' @param dir The folder holding the files.
#'
#' @return A named list (`only`, `pbo`, `cfp`) of data frames in the form
#'   returned by [standardise_itn_draws()].
#'
#' @export
read_itn_draws <- function(dir) {
  src <- itn_draws_source()
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  paths <- stats::setNames(file.path(dir, src$files), names(src$files))
  for (net in names(paths)) {
    if (!file.exists(paths[[net]])) {
      url <- sprintf("https://raw.githubusercontent.com/%s/%s/parameters/%s",
                     src$repo, src$commit, src$files[[net]])
      utils::download.file(url, paths[[net]], mode = "wb", quiet = TRUE)
    }
  }
  md5 <- tools::md5sum(paths)
  bad <- names(paths)[unname(md5) != src$md5[names(paths)]]
  if (length(bad)) {
    stop("Checksum mismatch for ", paste(paths[bad], collapse = ", "),
         "; delete the file to download it again.", call. = FALSE)
  }
  lapply(stats::setNames(names(paths), names(paths)),
         function(net) standardise_itn_draws(readRDS(paths[[net]]), net))
}

#' Put one file of ITN parameter draws into a common form
#'
#' Extracts the parameters as Part 6 of the source repository does: `dn0`;
#' `rn0` from `rn_pyr` (pyrethroid-only) or `rn_pbo` (both other nets); and
#' `gamman` as the file's durability (`mean_duration` or `mn_dur`) divided by
#' log(2). Draw `k` is taken from the row id, ids `(k - 1) * 101 + 1` to
#' `k * 101` being draw `k` at resistance 0 to 1, which is also the k-th row
#' at each level in file order, the order Part 6 relies on. Both are checked.
#'
#' @param x One file of draws, as read by [readRDS()].
#' @param net `"only"`, `"pbo"` or `"cfp"`.
#'
#' @return A data frame sorted by draw then resistance, with columns `draw`,
#'   `resistance` (rounded to 2 decimal places), `dn0`, `rn0` and `gamman`
#'   (in years; multiply by 365 for malariasimulation).
#'
#' @export
standardise_itn_draws <- function(x, net) {
  net <- match.arg(net, c("only", "pbo", "cfp"))
  if (net == "only") {
    id <- x$id
    res <- x$resistance
    rn0 <- x$rn_pyr
    durability <- x$mean_duration
  } else {
    id <- x$id.x
    res <- x$resistance.x
    rn0 <- x$rn_pbo
    durability <- x$mn_dur
  }
  n_res <- length(unique(round(res, 2)))
  level <- (id - 1L) %% n_res
  if (anyDuplicated(id) || any(abs(res - level / (n_res - 1L)) > 1e-9)) {
    stop("Row ids of the ", net, " draws do not follow the expected layout.",
         call. = FALSE)
  }
  out <- data.frame(draw = (id - 1L) %/% n_res + 1L,
                    resistance = round(res, 2),
                    dn0 = x$dn0, rn0 = rn0, gamman = durability / log(2))
  file_order <- unlist(lapply(split(out$draw, out$resistance), function(d) {
    identical(as.numeric(d), as.numeric(seq_along(d)))
  }))
  if (!all(file_order)) {
    stop("The ", net, " draws are not in draw order within each resistance ",
         "level.", call. = FALSE)
  }
  out <- out[order(out$draw, out$resistance), ]
  rownames(out) <- NULL
  out
}
