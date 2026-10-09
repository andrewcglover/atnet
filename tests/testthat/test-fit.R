test_that("every model's file exists and names only quantities it declares", {
  for (name in names(atn_models())) {
    spec <- atn_models()[[name]]
    file <- system.file("stan", spec$file, package = "atnet")
    expect_true(nzchar(file), info = name)
    expect_true(spec$data %in% c("eip", "tra", "tba_lab"), info = name)
    code <- paste(sub("//.*$", "", readLines(file)), collapse = "\n")
    for (par in c(spec$key, spec$pairs, "beta_exp")) {
      expect_true(grepl(paste0("\\b", par, "\\b"), code),
                  info = paste(name, "lacks", par))
    }
  }
})

test_that("the models use the sampler settings of the original drivers", {
  m <- atn_models()
  expect_equal(m$eip_hill$control, list(adapt_delta = 0.95, max_treedepth = 12))
  expect_equal(m$eip_exponential$control, m$eip_hill$control)
  expect_equal(m$tra$control, list(adapt_delta = 0.999, max_treedepth = 14))
  expect_equal(m$tba_lab$control, m$tra$control)
})

test_that("the Hill EIP initial values sit inside the model's bounds", {
  priors <- atn_priors()
  init <- eip_hill_inits(priors, K = 7L)
  set.seed(1)
  for (i in 1:20) {
    v <- init()
    expect_named(v, c("log_eip_excess", "r0_frac", "log_s_half", "log_nH",
                      "sigma_exp", "z_exp"))
    expect_length(v$z_exp, 7L)
    expect_true(v$r0_frac > 0 && v$r0_frac < 1)
    expect_true(v$sigma_exp > 0)
    expect_true(v$log_nH > priors$log_nH_eip$lower)
  }
})

test_that("an unknown model is refused before anything is compiled", {
  skip_if_not_installed("rstan")
  expect_error(fit_atn_model("eip_quadratic", data.frame()), "Unknown model")
})
