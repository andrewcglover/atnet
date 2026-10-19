# =====================================================================
# 08_submit_cluster.R
#   Runs analysis/07_run_projections.R on the DIDE cluster with hipercow, in two
#   stages:
#     1. one "histories" job: the results folder's settings.rds and every
#        region's shared 2000-2021 history, one region per core;
#     2. "realisations" jobs, each one region x a block of k consecutive
#        realisations on k cores, one realisation per core (realisation i =
#        row i of the draws table: its 2022-2024 segment, then every arm).
#   Every job writes its results to disk and returns only a short status.
#
#   Sourcing this file defines functions and submits nothing. In an R session
#   in the atnet folder on M: (the hipercow root):
#
#     source("analysis/08_submit_cluster.R")
#     b <- new_batch(c(OUT_SUFFIX = "_test", HUMAN_POP = "1000"), regions = c("Mopti", "Kidal"))
#     submit_histories(b)                    # shows the job; submit = TRUE submits it
#     batch_status(b)                        # until the histories job has finished
#     submit_realisations(b, 1, 10, k = 10)  # shows the jobs; submit = TRUE submits them
#     batch_status(b)
#     x <- collect_batch(b)
#
#   new_batch() checks that atnet's and the fork's code are committed and that
#   the fork's code matches the copy installed on the cluster; freezes atnet's
#   code at the last commit that changed it, in
#   outputs/projections/code/<commit>/, which the jobs run, so later edits
#   cannot reach a queued job; and runs a pre-flight on this computer: 07's
#   setup and every region's parameters (the history, and every arm at the
#   first and last realisations of the draws table), with no simulation. Its
#   settings are 07's for this results folder, as name = value strings (see
#   BATCH_SETTINGS); STAGE, SWEEP_CORES, SWEEP_REGIONS and DRAW_ROWS are set
#   per job. A later batch of realisations into the same folder must use
#   the same settings and code (07 checks settings.rds).
# =====================================================================

CLUSTER_MAX_CORES <- 32L        # cores per node
CLUSTER_R         <- "4.6.1"    # as set by hipercow::hipercow_configure()
CLUSTER_LIB       <- file.path("hipercow", "lib", "windows", CLUSTER_R)
# This computer's copy of the fork, and a library holding the site package
# version installed on the cluster (for the pre-flight).
FORK_LOCAL     <- Sys.getenv("ATN_FORK_PATH", "C:/Users/ag4218/Local/github/malariasimulation")
SITE_LIB_LOCAL <- Sys.getenv("ATNET_SITE_LIB", "C:/Users/ag4218/AppData/Local/R/atnet-lib/4.4")
# Where the frozen copies of atnet's code go (tests may point elsewhere).
SNAPSHOT_ROOT  <- Sys.getenv("ATNET_SNAPSHOT_ROOT", file.path("outputs", "projections", "code"))
# The code each repository's results depend on.
ATNET_CODE <- c("analysis/07_run_projections.R", "R", "DESCRIPTION", "NAMESPACE", "inst")
FORK_CODE  <- c("R", "src", "DESCRIPTION", "NAMESPACE", "inst")
# 07's settings that a batch may set (see 07's header).
BATCH_SETTINGS <- c("OUT_SUFFIX", "HUMAN_POP", "END_DAY", "HISTORY_DIR", "ANTIMAL_HL_YEARS",
                    "SWEEP_ARMS", "NET_RETENTION_DAYS", "ATN_DELTAQ", "ATN_SPOR_LEN",
                    "OMEGA_QUANTILE", "P_A_QUANTILE", "ATN_MECHANISM",
                    "FUTURE_RETENTION_DAYS")

# Nodes to keep work off, comma separated (HIPERCOW_EXCLUDE_NODES; "" for none). A node that
# cannot start work fails in seconds and so is free again at once, which makes the scheduler
# prefer it: one broken node will otherwise take a job from every submission. 2026-10-17:
# wpia-029 failed two jobs with "the trust relationship between this workstation and the
# primary domain failed". Remove it from this list once DIDE have repaired the node.
EXCLUDE_NODES <- trimws(strsplit(Sys.getenv("HIPERCOW_EXCLUDE_NODES", "wpia-029"), ",")[[1]])
EXCLUDE_NODES <- EXCLUDE_NODES[nzchar(EXCLUDE_NODES)]

# The nodes a task may run on: every node of the cluster but those. NULL = let hipercow choose
# (no exclusions wanted, or the node list could not be read).
usable_nodes <- function() {
  if (!length(EXCLUDE_NODES)) return(NULL)
  all <- tryCatch(hipercow::hipercow_cluster_info()$resources$nodes, error = function(e) NULL)
  if (is.null(all)) {
    message("Note: could not read the cluster's node list; not excluding ", toString(EXCLUDE_NODES))
    return(NULL)
  }
  keep <- setdiff(all, EXCLUDE_NODES)
  if (!length(keep)) stop("Excluding ", toString(EXCLUDE_NODES), " would leave no nodes.", call. = FALSE)
  keep
}

git <- function(dir, ...) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(dir), ...), stdout = TRUE, stderr = TRUE))
  list(out = out, ok = is.null(attr(out, "status")))
}
ascii      <- function(x) iconv(x, to = "ASCII//TRANSLIT")
region_tag <- function(region) gsub("[^A-Za-z0-9]", "", ascii(region))   # as in 07

# ---------------------------------------------------------------------
# Checks before anything is frozen or submitted
# ---------------------------------------------------------------------
check_clean <- function() {
  if (!dir.exists("hipercow") || !file.exists("analysis/07_run_projections.R")) {
    stop("Run this from the atnet folder on M: (the hipercow root).", call. = FALSE)
  }
  # atnet: the code the jobs run must be committed.
  st <- git(".", "status", "--porcelain", "--untracked-files=no", "--", ATNET_CODE)
  if (length(st$out)) {
    stop("atnet has uncommitted changes to its code; commit them first:\n",
         paste(st$out, collapse = "\n"), call. = FALSE)
  }
  # The last commit that changed that code: the frozen copy and the results are
  # keyed on it, so commits elsewhere in atnet (such as to this script) reuse it.
  atnet_sha <- git(".", "log", "-1", "--format=%H", "--", ATNET_CODE)$out
  if (git(".", "fetch", "-q")$ok) {
    n <- as.integer(strsplit(git(".", "rev-list", "--left-right", "--count", "HEAD...origin/main")$out,
                             "\\s+")[[1]])
    if (n[1] > 0) message("Note: atnet has ", n[1], " commit(s) not yet pushed to GitHub.")
    if (n[2] > 0) message("Note: GitHub has ", n[2], " atnet commit(s) not yet pulled here.")
  }

  # The fork: committed, and its code identical to the copy installed on the cluster.
  st <- git(FORK_LOCAL, "status", "--porcelain", "--untracked-files=no", "--", FORK_CODE)
  if (length(st$out)) {
    stop("The fork (", FORK_LOCAL, ") has uncommitted changes to its code:\n",
         paste(st$out, collapse = "\n"), call. = FALSE)
  }
  cl <- read.dcf(file.path(CLUSTER_LIB, "malariasimulation", "DESCRIPTION"),
                 fields = c("Version", "RemoteRef", "RemoteSha"))
  if (!git(FORK_LOCAL, "cat-file", "-e", cl[, "RemoteSha"])$ok) {
    stop("The cluster's fork commit ", cl[, "RemoteSha"], " is not in ", FORK_LOCAL,
         "; fetch it there first.", call. = FALSE)
  }
  if (!git(FORK_LOCAL, "diff", "--quiet", cl[, "RemoteSha"], "HEAD", "--", FORK_CODE)$ok) {
    stop("The fork's code here differs from the cluster's installed copy (", cl[, "RemoteRef"], ", ",
         substr(cl[, "RemoteSha"], 1, 7), "). Check that version out here, or re-provision ",
         "the cluster with the new one.", call. = FALSE)
  }

  # The pre-flight must use the site package version the cluster has.
  site_cluster <- unname(read.dcf(file.path(CLUSTER_LIB, "site", "DESCRIPTION"), fields = "Version")[1, 1])
  site_local   <- unname(read.dcf(file.path(SITE_LIB_LOCAL, "site", "DESCRIPTION"), fields = "Version")[1, 1])
  if (!identical(site_cluster, site_local)) {
    stop("The site package differs: ", site_cluster, " on the cluster, ", site_local, " in ",
         SITE_LIB_LOCAL, ".", call. = FALSE)
  }
  list(atnet_sha = atnet_sha, fork = paste(cl[, "RemoteRef"], cl[, "RemoteSha"]),
       site_version = site_cluster)
}

# ---------------------------------------------------------------------
# A frozen copy of atnet's code at a commit (made once per commit, never changed)
# ---------------------------------------------------------------------
snapshot_code <- function(sha) {
  dest <- file.path(SNAPSHOT_ROOT, substr(sha, 1, 12))
  if (dir.exists(dest)) {
    man <- read.csv(file.path(dest, "MANIFEST.csv"), stringsAsFactors = FALSE)
    if (!identical(unname(tools::md5sum(file.path(dest, man$file))), man$md5)) {
      stop(dest, " has changed since it was made. Move it to make a fresh copy.", call. = FALSE)
    }
    return(normalizePath(dest, winslash = "/"))
  }
  tmp <- paste0(dest, "_partial")
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE)
  dir.create(tmp, recursive = TRUE)
  zip <- tempfile(fileext = ".zip")
  if (!git(".", "archive", "--format=zip", "-o", shQuote(zip), sha, ATNET_CODE)$ok) {
    stop("git archive failed for atnet commit ", sha, call. = FALSE)
  }
  utils::unzip(zip, exdir = tmp)
  files <- list.files(tmp, recursive = TRUE)
  write.csv(data.frame(file = files, md5 = unname(tools::md5sum(file.path(tmp, files)))),
            file.path(tmp, "MANIFEST.csv"), row.names = FALSE)
  if (!file.rename(tmp, dest)) stop("Could not rename ", tmp, " to ", dest, call. = FALSE)
  message("Froze atnet's code at ", substr(sha, 1, 12), " in ", dest)
  normalizePath(dest, winslash = "/")
}

# ---------------------------------------------------------------------
# The pre-flight: 07's setup and parameters on this computer, no simulation
# ---------------------------------------------------------------------
# Run in a separate R process (Rscript) with the fork from source, the cluster's
# site package version and the frozen code. rows: realisations, "last" = the
# draws table's last.
preflight_main <- function(script, out, rows, history) {
  options(country_test_mode = TRUE)
  source(script, local = FALSE)
  rows <- strsplit(rows, ",")[[1]]
  rows <- unique(ifelse(rows == "last", nrow(draws_tab), as.integer(rows)))
  n <- 0L
  for (r in sweep_regions) {
    if (isTRUE(as.logical(history))) { invisible(build_params(r, "history", NULL)); n <- n + 1L }
    for (i in rows) for (a in arms) { invisible(build_params(r, a, i)); n <- n + 1L }
  }
  saveRDS(list(built = n, rows = rows, regions = sweep_regions,
               all_regions = setequal(sweep_regions, regions), arms = arms,
               n_draws = nrow(draws_tab), iso = tolower(COUNTRY_ISO),
               out_dir = normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE),
               history_files = vapply(sweep_regions, function(r)
                 normalizePath(history_file(r), winslash = "/", mustWork = FALSE), character(1)),
               irs_corrections = irs_corrections, code_md5 = code_md5), out)
}

preflight <- function(b, rows = c("1", "last"), history = TRUE) {
  fn <- tempfile(fileext = ".rds"); out <- tempfile(fileext = ".rds"); log <- tempfile(fileext = ".log")
  f <- preflight_main
  environment(f) <- globalenv()
  saveRDS(f, fn)
  env <- c(b$settings, ATN_FORK_PATH = FORK_LOCAL, ATNET_CODE_DIR = b$snapshot,
           ATNET_VERSION = b$atnet_version, R_LIBS = SITE_LIB_LOCAL, DRAW_ROWS = "1",
           SWEEP_REGIONS = if (is.null(b$regions_env)) NA else b$regions_env,
           STAGE = NA, SWEEP_CORES = NA, STOP_BY = NA)
  t0 <- Sys.time()
  status <- withr::with_envvar(env, system2(
    file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote("do.call(readRDS(commandArgs(TRUE)[1]), as.list(commandArgs(TRUE)[-1]))"),
      shQuote(fn), shQuote(file.path(b$snapshot, "analysis", "07_run_projections.R")), shQuote(out),
      paste(rows, collapse = ","), history),
    stdout = log, stderr = log))
  if (status != 0 || !file.exists(out)) {
    stop("Pre-flight failed. End of its log (", log, "):\n",
         paste(utils::tail(readLines(log, warn = FALSE), 25), collapse = "\n"), call. = FALSE)
  }
  pf <- readRDS(out)
  message(sprintf("Pre-flight passed: %d parameter sets (%d region(s): %s; %severy arm at realisation(s) %s) in %.0f s",
                  pf$built, length(pf$regions), paste(pf$regions, collapse = ", "),
                  if (history) "history and " else "", paste(pf$rows, collapse = " and "),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  pf
}

# ---------------------------------------------------------------------
# A batch: one results folder, its settings, and the frozen code
# ---------------------------------------------------------------------
new_batch <- function(settings = character(), regions = NULL) {
  if (length(settings) && (is.null(names(settings)) || !all(names(settings) %in% BATCH_SETTINGS))) {
    stop("settings must be named, from: ", paste(BATCH_SETTINGS, collapse = ", "), call. = FALSE)
  }
  cc <- check_clean()
  b <- list(settings = settings,
            regions_env = if (is.null(regions)) NULL else paste(ascii(regions), collapse = ","),
            snapshot = snapshot_code(cc$atnet_sha), atnet_version = cc$atnet_sha,
            fork = cc$fork, site_version = cc$site_version)
  pf <- preflight(b)
  b[c("regions", "all_regions", "arms", "n_draws", "iso", "out_dir", "history_files", "code_md5")] <-
    pf[c("regions", "all_regions", "arms", "n_draws", "iso", "out_dir", "history_files", "code_md5")]
  message("Results folder: ", b$out_dir)
  b
}

# ---------------------------------------------------------------------
# Submitting
# ---------------------------------------------------------------------
# What every job runs on its node: 07 from the frozen code, in a fresh R process
# whose output goes to a log in the results folder. run_env and script are the
# batch's; stage, region, realisations, cores and log_file are the job's.
CLUSTER_TASK <- quote({
  Sys.unsetenv(c("ATN_FORK_PATH", "SWEEP_REGIONS", "DRAW_ROWS", "STOP_BY"))
  do.call(Sys.setenv, as.list(run_env))
  Sys.setenv(STAGE = stage, SWEEP_CORES = as.character(cores),
             R_LIBS = paste(.libPaths(), collapse = ";"))
  if (nzchar(region)) Sys.setenv(SWEEP_REGIONS = region)
  if (nzchar(realisations)) Sys.setenv(DRAW_ROWS = realisations)
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                    stdout = log_file, stderr = log_file)
  if (status != 0) {
    stop("07 stopped (status ", status, "). End of ", log_file, ":\n",
         paste(utils::tail(readLines(log_file, warn = FALSE), 20), collapse = "\n"), call. = FALSE)
  }
  list(status = status, node = Sys.info()[["nodename"]], log = log_file)
})

submissions_file <- function(out_dir) file.path(out_dir, "submissions.rds")

read_submissions <- function(out_dir) {
  f <- submissions_file(out_dir)
  if (file.exists(f)) readRDS(f) else NULL
}

log_path <- function(b, stage, region = "", realisations = "") {
  file.path(b$out_dir, "logs", paste0(stage, if (nzchar(region)) paste0("_", region_tag(region)),
                                      if (nzchar(realisations)) paste0("_rows", realisations),
                                      format(Sys.time(), "_%Y%m%d-%H%M%S"), ".log"))
}

# Submits the jobs in tasks (one per row: stage, region, realisations, cores,
# log_file), one bulk submission per number of cores, and records them in the
# results folder's submissions.rds.
submit_tasks <- function(b, tasks, name) {
  run_env <- c(b$settings, ATNET_CODE_DIR = b$snapshot, ATNET_VERSION = b$atnet_version)
  script  <- file.path(b$snapshot, "analysis", "07_run_projections.R")
  dir.create(b$out_dir, recursive = TRUE, showWarnings = FALSE)
  nodes <- usable_nodes()
  if (!is.null(nodes)) message("Keeping work off: ", toString(EXCLUDE_NODES))
  for (k in unique(tasks$cores)) {
    data <- tasks[tasks$cores == k, , drop = FALSE]
    bundle <- eval(bquote(hipercow::task_create_bulk_expr(
      .(CLUSTER_TASK), data = data, bundle_name = .(paste0(name, "_", k, "cores")),
      resources = hipercow::hipercow_resources(cores = .(k), requested_nodes = .(nodes)))))
    rec <- cbind(task = bundle$ids, bundle = bundle$name, data, submitted = Sys.time(),
                 snapshot = b$snapshot, atnet_version = b$atnet_version, stringsAsFactors = FALSE)
    saveRDS(rbind(read_submissions(b$out_dir), rec), submissions_file(b$out_dir))
    message(sprintf("Submitted %d job(s) of %d core(s) as bundle %s", nrow(data), k, bundle$name))
  }
  invisible(read_submissions(b$out_dir))
}

# Jobs of this folder still queued or running, with their stage, region and realisations.
active_tasks <- function(out_dir) {
  s <- read_submissions(out_dir)
  if (is.null(s)) return(s)
  s$status <- hipercow::task_status(s$task)
  s[s$status %in% c("created", "submitted", "running"), , drop = FALSE]
}

submit_histories <- function(b, submit = FALSE) {
  cores <- min(length(b$regions), CLUSTER_MAX_CORES)
  tasks <- data.frame(stage = "histories", region = if (is.null(b$regions_env)) "" else b$regions_env,
                      realisations = "", cores = cores, stringsAsFactors = FALSE)
  tasks$log_file <- log_path(b, "histories")
  message(sprintf("1 histories job: %d region(s) (%s) on %d core(s); results folder %s",
                  length(b$regions), paste(b$regions, collapse = ", "), cores, b$out_dir))
  if (!submit) {
    message("Not submitted: run again with submit = TRUE to submit.")
    return(invisible(tasks))
  }
  if (NROW(a <- active_tasks(b$out_dir)) && any(a$stage == "histories")) {
    stop("A histories job for this folder is already queued or running.", call. = FALSE)
  }
  submit_tasks(b, tasks, paste0(basename(b$out_dir), "_histories"))
}

# Is a sweep summary file complete (no failed or unstarted jobs)? As 07 judges it.
sweep_complete <- function(f) {
  if (!file.exists(f)) return(FALSE)
  s <- readRDS(f)
  !NROW(s$failures) && all(s$jobs$started | s$jobs$skipped)
}

# One job per region x block of k consecutive realisations (first to last).
submit_realisations <- function(b, first, last, k = 10L, regions = b$regions, submit = FALSE) {
  first <- as.integer(first); last <- as.integer(last); k <- as.integer(k)
  if (first < 1L || last > b$n_draws || last < first) {
    stop("Realisations must lie within 1-", b$n_draws, " (rows of the draws table).", call. = FALSE)
  }
  if (k < 1L || k > CLUSTER_MAX_CORES) stop("k must be 1 to ", CLUSTER_MAX_CORES, ".", call. = FALSE)
  if (!all(ascii(regions) %in% ascii(b$regions))) {
    stop("This batch's regions are ", paste(b$regions, collapse = ", "), call. = FALSE)
  }
  starts <- seq.int(first, last, by = k)
  ends   <- pmin(starts + k - 1L, last)
  # Every region's job for a block before the next block, so whole realisations finish first.
  tasks  <- expand.grid(region = ascii(regions), block = seq_along(starts), stringsAsFactors = FALSE)
  tasks  <- data.frame(stage = "realisations", region = tasks$region,
                       realisations = sprintf("%d-%d", starts[tasks$block], ends[tasks$block]),
                       cores = ends[tasks$block] - starts[tasks$block] + 1L, stringsAsFactors = FALSE)

  # Leave out region x block jobs already complete, or queued or running.
  sweep <- file.path(b$out_dir, sprintf("sweep_%srows%s.rds",
    if (b$all_regions && length(b$regions) == 1L) "" else paste0(region_tag(tasks$region), "_"),
    tasks$realisations))
  done <- vapply(sweep, sweep_complete, logical(1))
  a    <- active_tasks(b$out_dir)
  busy <- paste(tasks$region, tasks$realisations) %in% paste(a$region, a$realisations)
  if (any(done | busy)) {
    message(sprintf("Leaving out %d job(s) already complete and %d already queued or running.",
                    sum(done), sum(busy & !done)))
  }
  tasks <- tasks[!(done | busy), , drop = FALSE]
  if (!nrow(tasks)) return(invisible(tasks))
  tasks$log_file <- unname(mapply(function(r, x) log_path(b, "realisations", r, x),
                                  tasks$region, tasks$realisations))
  message(sprintf("%d realisations job(s): %d region(s) x blocks of up to %d realisations (%d-%d), %d cores in all",
                  nrow(tasks), length(unique(tasks$region)), k, first, last, sum(tasks$cores)))

  # The histories stage must have finished: settings.rds and every history in place.
  missing <- c(if (!file.exists(file.path(b$out_dir, "settings.rds"))) "settings.rds",
               sprintf("history of %s", names(b$history_files)[!file.exists(b$history_files)]))
  if (!submit) {
    if (length(missing)) {
      message("Not yet possible: the histories job has not finished (missing: ",
              paste(missing, collapse = ", "), ").")
    }
    message("Not submitted: run again with submit = TRUE to submit.")
    return(invisible(tasks))
  }
  if (length(missing)) {
    stop("The histories job has not finished (missing: ", paste(missing, collapse = ", "), ").",
         call. = FALSE)
  }
  # Parameters of the first and last realisation of this range, every region and arm.
  preflight(b, rows = as.character(unique(c(first, last))), history = FALSE)
  submit_tasks(b, tasks, sprintf("%s_rows%d-%d", basename(b$out_dir), first, last))
}

# ---------------------------------------------------------------------
# Following and collecting
# ---------------------------------------------------------------------
# Region tag and realisation of each job file in a folder.
job_files <- function(dir) {
  f <- list.files(dir, "_row[0-9]+\\.rds$")
  m <- regmatches(f, regexec("^[a-z]+_(.+)_row([0-9]+)\\.rds$", f))
  data.frame(tag = vapply(m, `[`, "", 2), draw = as.integer(vapply(m, `[`, "", 3)),
             stringsAsFactors = FALSE)
}

# The state of each job submitted for a results folder (a batch, or the folder's
# path), from hipercow and from the files on disk.
batch_status <- function(b) {
  out_dir <- if (is.character(b)) b else b$out_dir
  s <- read_submissions(out_dir)
  if (is.null(s)) {
    message("Nothing submitted yet for ", out_dir)
    return(invisible(NULL))
  }
  s$status <- hipercow::task_status(s$task)
  saved  <- job_files(file.path(out_dir, "jobs"))
  failed <- job_files(file.path(out_dir, "jobs_failed"))
  count <- function(x, region, realisations) {
    if (!nzchar(realisations)) return(NA_integer_)
    r <- as.integer(strsplit(realisations, "-")[[1]])
    sum(x$tag == region_tag(region) & x$draw >= r[1] & x$draw <= r[2])
  }
  s$saved  <- mapply(count, list(saved), s$region, s$realisations)
  s$failed <- mapply(function(region, realisations) {
    if (!nzchar(realisations)) return(NA_integer_)
    r <- as.integer(strsplit(realisations, "-")[[1]])
    f <- failed[failed$tag == region_tag(region) & failed$draw >= r[1] & failed$draw <= r[2], ]
    sum(!paste(f$tag, f$draw) %in% paste(saved$tag, saved$draw))   # not since rerun
  }, s$region, s$realisations)
  if (!is.character(b)) {
    message(sprintf("Histories saved: %d of %d", sum(file.exists(b$history_files)), length(b$history_files)))
  }
  print(s[, c("stage", "region", "realisations", "cores", "status", "saved", "failed", "submitted")],
        row.names = FALSE)
  message("Logs: ", file.path(out_dir, "logs"), "  (or hipercow::task_log_show(<task id>))")
  invisible(s)
}

# The sweep summaries of a results folder (a batch, or the folder's path): jobs,
# failed segments and run times, read from disk.
collect_batch <- function(b) {
  out_dir <- if (is.character(b)) b else b$out_dir
  files <- list.files(out_dir, "^sweep_.*\\.rds$", full.names = TRUE)
  if (!length(files)) {
    message("No sweep summaries yet in ", out_dir)
    return(invisible(NULL))
  }
  sw <- lapply(files, readRDS)
  jobs <- do.call(rbind, lapply(seq_along(sw), function(i) cbind(sweep = basename(files[i]), sw[[i]]$jobs)))
  failures <- dplyr::bind_rows(lapply(sw, `[[`, "failures"))
  minutes  <- dplyr::bind_rows(lapply(sw, `[[`, "minutes"))
  saved    <- job_files(file.path(out_dir, "jobs"))
  message(sprintf("%d sweep summaries; %d region x realisation jobs saved; %d failed segment(s)",
                  length(files), nrow(saved), NROW(failures)))
  if (NROW(minutes)) {
    m <- minutes[!minutes$skipped, , drop = FALSE]
    if (nrow(m)) {
      message("Median minutes per segment (jobs run, not skipped):")
      print(round(tapply(m$minutes, factor(m$segment, unique(m$segment)), stats::median), 1))
    }
  }
  invisible(list(settings = readRDS(file.path(out_dir, "settings.rds")), jobs = jobs,
                 saved = saved, failures = failures, minutes = minutes))
}
