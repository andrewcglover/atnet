test_that("the priors are those the current posteriors were fitted with", {
  # Saved by each of the four fits in May 2026; all four records are identical.
  # The record lacks the truncation point, then fixed inside the Hill model
  # (`real<lower=0.5> log_nH`), as was sigma_exp's standard deviation of 1 in
  # both EIP models.
  used <- jsonlite::read_json(test_path("fixtures", "priors_used_20260520.json"))
  now <- atn_priors()
  expect_equal(now$log_nH_eip$lower, 0.5)
  now$log_nH_eip$lower <- NULL
  expect_equal(now, used, tolerance = 1e-12)
  expect_equal(atn_priors()$sigma_exp$sd, 1)
})

stan_files <- c("eip_fit_hill.stan", "int_fit_nb_rate_hill_bidir_splitnH.stan",
                "si_comparison/eip_fit.stan",
                "si_comparison/int_fit_linear_hill_bidir_splitnH.stan")

stan_code <- function(file) {
  sub("//.*$", "", readLines(system.file("stan", file, package = "atnet")))
}

# Prior values written into a model rather than passed in. A half-normal is
# written normal(0, sd), so a leading zero is allowed; bounds of 0 and 1 are
# the natural limits of a scale or a probability.
fixed_priors <- function(code) {
  draws <- unlist(regmatches(code, gregexpr("~\\s*\\w+\\([^)]*\\)", code)))
  args <- strsplit(sub("^~\\s*\\w+\\((.*)\\)$", "\\1", draws), ",")
  literal <- function(a) grepl("^[0-9.]+$", trimws(a))
  in_draws <- draws[vapply(args, function(a) {
    any(literal(a[-1])) || (length(a) > 0L && literal(a[1]) && trimws(a[1]) != "0")
  }, logical(1))]
  bounds <- unlist(regmatches(code, gregexpr("(lower|upper)=[0-9.]+", code)))
  c(in_draws, bounds[!sub(".*=", "", bounds) %in% c("0", "1")])
}

test_that("no prior value is fixed inside a model", {
  for (file in stan_files) {
    expect_equal(fixed_priors(stan_code(file)), character(0), info = file)
  }
  # The check would have caught both values that used to be fixed.
  expect_length(fixed_priors(c("real<lower=0.5> log_nH;",
                               "sigma_exp ~ normal(0, 1);")), 2L)
})

test_that("every hyperparameter a model declares is supplied", {
  # The scalar reals in a model's data block are its prior hyperparameters;
  # everything else there is an integer size or a data array.
  declared <- function(file) {
    lines <- stan_code(file)
    start <- grep("^data\\s*\\{", lines)
    end <- start + grep("^\\}", lines[-seq_len(start)])[1]
    block <- lines[start:end]
    hits <- regmatches(block, regexec("^\\s*real(<[^>]*>)?\\s+(\\w+)\\s*;", block))
    unlist(lapply(hits, function(h) if (length(h)) h[[3]]))
  }
  supplied <- names(prior_stan_data())
  for (file in stan_files) {
    wanted <- declared(file)
    expect_gt(length(wanted), 0L)
    expect_true(all(wanted %in% supplied),
                info = paste(file, "lacks", toString(setdiff(wanted, supplied))))
  }
})

test_that("the prior draws have the requested shape", {
  d <- sample_kernel_prior(50, "bidir_split", seed = 1)
  expect_named(d, c("r_min", "b_max", "s_half_pre", "s_half_post",
                    "nH_pre", "nH_post"))
  expect_true(all(lengths(d) == 50L))
  expect_equal(d$b_max, 1 - d$r_min)
})

test_that("the priors round-trip through the saved record", {
  path <- withr::local_tempfile(fileext = ".json")
  write_priors_json(path)
  expect_equal(jsonlite::read_json(path), atn_priors(), tolerance = 1e-12)
})
