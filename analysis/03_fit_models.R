# Fit the Stan models to the tables written by 02_build_fitting_data.R.
#
# Fits all four models by default, or those named on the command line:
#   Rscript analysis/03_fit_models.R eip_hill tra
# The names are those of atn_models(): eip_hill and tra feed the simulations;
# eip_exponential and tba_lab are fitted for comparison in the SI.
#
# Each run writes to outputs/fits/<YYYYMMDD_HHMM>/, for each model:
#   <model>.rds              the stanfit object
#   <model>_priors.json      the priors it was fitted with
#   <model>_diagnostics.txt  sampler diagnostics and parameter summaries
#   <model>_trace.png, <model>_pairs.png, <model>_ppc.png
# and inputs.txt, naming the tables the fits were given.
#
# Run with the working directory at the root of the package.

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(bayesplot)
  library(ggplot2)
})

models <- commandArgs(trailingOnly = TRUE)
if (!length(models)) {
  models <- names(atn_models())
}
unknown <- setdiff(models, names(atn_models()))
if (length(unknown)) {
  stop("Unknown model(s): ", toString(unknown), call. = FALSE)
}

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

out_dir <- file.path("outputs", "fits", format(Sys.time(), "%Y%m%d_%H%M"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("sporozoite table:", spz_path),
             paste("oocyst table:    ", ooc_path),
             paste("models:          ", toString(models)),
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
    ylab <- "total oocysts"
  } else {
    obs <- c(data$pos_ctl, data$pos_trt)
    rep <- cbind(rstan::extract(fit, pars = "pos_ctl_rep")[[1]],
                 rstan::extract(fit, pars = "pos_trt_rep")[[1]])
    arm <- rep(c("control", "treated"), each = data$N)
    ylab <- if (model == "tba_lab") "oocyst-positive" else "sporozoite-positive"
  }
  ppc_intervals_grouped(y = obs, yrep = rep, group = arm,
                        x = seq_along(obs)) +
    labs(title = sprintf("%s: posterior predictive check", model),
         x = "experiment and arm", y = ylab)
}

for (model in models) {
  spec <- atn_models()[[model]]
  table <- if (spec$data == "eip") spz else ooc
  message(sprintf("\n==== %s: %s ====", model, spec$file))
  started <- Sys.time()
  fit <- fit_atn_model(model, table)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

  path <- function(suffix) file.path(out_dir, paste0(model, suffix))
  saveRDS(fit, path(".rds"))
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
  cat(sprintf("%s fitted in %.1f minutes\n\n", model, minutes))
  cat(sprintf("divergent transitions:   %d\n",
              sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))))
  cat(sprintf("at maximum tree depth:   %d\n",
              sum(vapply(sp, function(x) {
                sum(x[, "treedepth__"] >= spec$control$max_treedepth)
              }, numeric(1)))))
  cat(sprintf("largest R-hat:           %.4f\n", max(all_summ[, "Rhat"], na.rm = TRUE)))
  cat(sprintf("smallest effective size: %.0f\n\n", min(all_summ[, "n_eff"], na.rm = TRUE)))
  writeLines(hmc)
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
  data <- switch(spec$data, eip = eip_stan_data(table),
                 tra = blocking_stan_data(table, "tra"),
                 tba_lab = blocking_stan_data(table, "tba_lab"))
  ggsave(path("_ppc.png"), ppc_plot(fit, model, data),
         width = 11, height = 5, dpi = 200)

  message(sprintf("%s done in %.1f minutes -> %s", model, minutes, out_dir))
}
