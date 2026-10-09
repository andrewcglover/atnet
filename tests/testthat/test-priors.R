test_that("the priors are those the current posteriors were fitted with", {
  # Saved by each of the four fits in May 2026; all four records are identical.
  used <- jsonlite::read_json(test_path("fixtures", "priors_used_20260520.json"))
  expect_equal(atn_priors(), used, tolerance = 1e-12)
})

test_that("every hyperparameter a model declares is supplied", {
  # The scalar reals in a model's data block are its prior hyperparameters;
  # everything else there is an integer size or a data array.
  declared <- function(file) {
    lines <- readLines(system.file("stan", file, package = "atnet"))
    start <- grep("^data\\s*\\{", lines)
    end <- start + grep("^\\}", lines[-seq_len(start)])[1]
    block <- sub("//.*$", "", lines[start:end])
    hits <- regmatches(block, regexec("^\\s*real(<[^>]*>)?\\s+(\\w+)\\s*;", block))
    unlist(lapply(hits, function(h) if (length(h)) h[[3]]))
  }
  supplied <- names(prior_stan_data())
  for (file in c("eip_fit_hill.stan", "int_fit_nb_rate_hill_bidir_splitnH.stan",
                 "si_comparison/eip_fit.stan",
                 "si_comparison/int_fit_linear_hill_bidir_splitnH.stan")) {
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
