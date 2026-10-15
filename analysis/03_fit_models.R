# Fit the Stan models to the tables written by 02_build_fitting_data.R.
#
# Fits all four models by default, or those named on the command line:
#   Rscript analysis/03_fit_models.R eip_hill tra
# The names are those of atn_models(): eip_hill and tra feed the simulations;
# eip_exponential and tba_lab are fitted for comparison in the SI.
#
# --smoke runs one chain of 200 iterations (100 warmup) per model instead, to
# check the pipeline end to end and time it. Its draws mean nothing; its
# diagnostics file estimates how long the full fit would take.
#   Rscript analysis/03_fit_models.R --smoke
#
# The EIP models fit an Erlang EIP with ATN_ERLANG_SHAPE stages (default 30, which
# must match the spor_len of 07_run_projections.R, since the kernels fitted here are
# applied to that Erlang; the note there says why 30). Earlier fits used other
# shapes: outputs/fits/20261009_1853 used 10 and 20261014_1554_eipR40 used 40.
# The blocking models do not use it.
#   ATN_ERLANG_SHAPE=10 Rscript analysis/03_fit_models.R eip_hill
#
# Each run writes to outputs/fits/<YYYYMMDD_HHMM>/, with "_eipR<shape>" appended
# when an EIP model is fitted and "_smoke" for a smoke run, for each model:
#   <model>.rds              the stanfit object (every draw)
#   <model>_data.rds         the exact data list the model was given
#   <model>_priors.json      the priors it was fitted with
#   <model>_diagnostics.txt  sampler diagnostics and parameter summaries
#   <model>_trace.png, <model>_pairs.png, <model>_ppc.png
# and inputs.txt, naming the tables the fits were given and the Erlang shape.
#
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(bayesplot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% args
models <- setdiff(args, "--smoke")
full <- list(chains = 4L, iter = 3000L, warmup = 1000L)
run <- if (smoke) list(chains = 1L, iter = 200L, warmup = 100L) else full
if (!length(models)) {
  models <- names(atn_models())
}
unknown <- setdiff(models, names(atn_models()))
if (length(unknown)) {
  stop("Unknown model(s): ", toString(unknown), call. = FALSE)
}
erlang_shape <- as.integer(Sys.getenv("ATN_ERLANG_SHAPE", "30"))
if (is.na(erlang_shape) || erlang_shape < 1L) {
  stop("ATN_ERLANG_SHAPE must be a positive whole number.", call. = FALSE)
}
fits_eip <- any(vapply(atn_models()[models], function(s) s$data == "eip", logical(1)))

latest <- function(pattern) {
  f <- sort(list.files("outputs", pattern = pattern, full.names = TRUE))
  if (!length(f)) {
    stop("No table matching '", pattern, "' in outputs/. Run ",
         "analysis/02_build_fitting_data.R first.", call. = FALSE)
  }
  f[length(f)]
}
spz_path <- latest("^[0-9]{8}_sporozoite_experiments\\.csv$")
ooc_path <- latest("^[0-9]{8}_oocyst_mosquitoes\\.csv$")
spz <- read.csv(spz_path, check.names = FALSE)
ooc <- read.csv(ooc_path, check.names = FALSE)

out_dir <- file.path("outputs", "fits", paste0(format(Sys.time(), "%Y%m%d_%H%M"),
                                              if (fits_eip) paste0("_eipR", erlang_shape),
                                              if (smoke) "_smoke"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("sporozoite table:", spz_path),
             paste("oocyst table:    ", ooc_path),
             paste("models:          ", toString(models)),
             if (fits_eip) sprintf("erlang shape:     %d (EIP stages; EIP models only)",
                                   erlang_shape),
             sprintf("run:              %d chain(s), %d iterations, %d warmup%s",
                     run$chains, run$iter, run$warmup,
                     if (smoke) " (smoke run)" else ""),
             paste("atnet commit:    ",
                   tryCatch(system2("git", c("rev-parse", "--short", "HEAD"),
                                    stdout = TRUE), error = function(e) NA))),
           file.path(out_dir, "inputs.txt"))

# Posterior predictive check: observed against replicated counts per
# experiment and arm. For TRA the oocyst counts are summed within each.
ppc_plot <- function(fit, model, data) {
  if (model == "tra") {
    yrep <- rstan::extract(fit, pars = "y_rep")$y_rep
    cell <- interaction(data$row_id, data$arm, drop = TRUE, lex.order = TRUE)
    idx <- split(seq_len(data$M), cell)
    obs <- vapply(idx, function(ix) sum(data$y[ix]), numeric(1))
    rep <- vapply(idx, function(ix) rowSums(yrep[, ix, drop = FALSE]),
                  numeric(nrow(yrep)))
    arm <- ifelse(grepl("\\.1$", names(idx)), "treated", "control")
    experiment <- as.integer(sub("\\..*$", "", names(idx)))
    ylab <- "total oocysts"
  } else {
    obs <- c(data$pos_ctl, data$pos_trt)
    rep <- cbind(rstan::extract(fit, pars = "pos_ctl_rep")[[1]],
                 rstan::extract(fit, pars = "pos_trt_rep")[[1]])
    arm <- rep(c("control", "treated"), each = data$N)
    experiment <- rep(seq_len(data$N), 2L)
    ylab <- if (model == "tba_lab") "oocyst-positive" else "sporozoite-positive"
  }
  ppc_intervals_grouped(y = obs, yrep = rep, group = arm, x = experiment) +
    labs(title = sprintf("%s: posterior predictive check", model),
         x = "experiment", y = ylab)
}

for (model in models) {
  spec <- atn_models()[[model]]
  table <- if (spec$data == "eip") spz else ooc
  message(sprintf("\n==== %s: %s ====", model, spec$file))
  started <- Sys.time()
  # Warnings are kept for the diagnostics file rather than lost at the end.
  warned <- character(0)
  fit <- withCallingHandlers(
    fit_atn_model(model, table, chains = run$chains, iter = run$iter,
                  warmup = run$warmup, erlang_shape = erlang_shape),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  # Sampling time per chain, apart from compiling. Chains run side by side, so
  # the slowest chain sets the time; scaled up to the full run's iterations.
  elapsed <- rstan::get_elapsed_time(fit)
  slowest <- elapsed[which.max(rowSums(elapsed)), ]
  per_warmup <- slowest[["warmup"]] / run$warmup
  per_sample <- slowest[["sample"]] / (run$iter - run$warmup)
  full_minutes <- (per_warmup * full$warmup +
                     per_sample * (full$iter - full$warmup)) / 60

  path <- function(suffix) file.path(out_dir, paste0(model, suffix))
  saveRDS(fit, path(".rds"))
  data <- switch(spec$data, eip = eip_stan_data(table, erlang_shape = erlang_shape),
                 tra = blocking_stan_data(table, "tra"),
                 tba_lab = blocking_stan_data(table, "tba_lab"))
  saveRDS(data, path("_data.rds"))
  write_priors_json(path("_priors.json"))

  summ <- rstan::summary(fit, pars = spec$key)$summary
  re <- rstan::summary(fit, pars = "beta_exp")$summary
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  # R-hat and effective size over everything but the posterior predictive
  # draws, of which there are many and whose R-hat strays by chance.
  checked <- grep("_rep$", fit@model_pars, value = TRUE, invert = TRUE)
  all_summ <- rstan::summary(fit, pars = checked)$summary
  hmc <- utils::capture.output(rstan::check_hmc_diagnostics(fit),
                               type = "message")
  sink(path("_diagnostics.txt"))
  cat(sprintf("%s fitted in %.1f minutes, including compiling\n", model, minutes))
  if (spec$data == "eip") {
    cat(sprintf("Erlang shape (EIP stages): %d\n", data$Delta_r))
  }
  cat(sprintf("slowest chain: %.0f s warmup, %.0f s sampling\n",
              slowest[["warmup"]], slowest[["sample"]]))
  if (smoke) {
    cat(sprintf(paste0("full run (%d warmup + %d sampling, chains side by side) ",
                       "would take about %.0f minutes plus compiling\n"),
                full$warmup, full$iter - full$warmup, full_minutes))
  }
  cat("\n")
  cat(sprintf("divergent transitions:   %d\n",
              sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))))
  cat(sprintf("at maximum tree depth:   %d\n",
              sum(vapply(sp, function(x) {
                sum(x[, "treedepth__"] >= spec$control$max_treedepth)
              }, numeric(1)))))
  cat(sprintf("largest R-hat:           %.4f\n", max(all_summ[, "Rhat"], na.rm = TRUE)))
  cat(sprintf("smallest effective size: %.0f\n\n", min(all_summ[, "n_eff"], na.rm = TRUE)))
  writeLines(hmc)
  cat(sprintf("\nWarnings while fitting: %d\n", length(warned)))
  if (length(warned)) {
    writeLines(paste("-", trimws(warned)))
  }
  cat("\nKey parameters\n")
  print(round(summ, 4))
  cat("\nPer-experiment effects (beta_exp)\n")
  print(round(re, 4))
  sink()

  arr <- as.array(fit, pars = spec$pairs)
  ggsave(path("_trace.png"), mcmc_trace(arr), width = 11, height = 7, dpi = 200)
  ggsave(path("_pairs.png"),
         mcmc_pairs(arr, off_diag_args = list(size = 0.4, alpha = 0.25)),
         width = 10, height = 10, dpi = 200)
  ggsave(path("_ppc.png"), ppc_plot(fit, model, data),
         width = 11, height = 5, dpi = 200)

  message(sprintf("%s done in %.1f minutes, %d warning(s) -> %s", model,
                  minutes, length(warned), out_dir))
  if (smoke) {
    message(sprintf("%s full run estimate: about %.0f minutes plus compiling",
                    model, full_minutes))
  }
}
