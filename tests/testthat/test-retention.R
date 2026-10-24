test_that("sc_retention_params hits the mean it was asked for", {
  for (m in seq(1.2, 3.6, 0.6) * 365) {
    p <- sc_retention_params(m)
    expect_equal(p$mean, m, tolerance = 1e-6)
  }
})

test_that("the returned pair implies the zero day it was given", {
  # set_bednets() rebuilds l from (half_life, k); it must come back to 20 y.
  for (m in c(400, 657, 1300)) {
    p <- sc_retention_params(m)
    l <- p$half_life / sqrt(1 - p$k / (p$k - log(0.5)))
    expect_equal(l, sc_zero_days(), tolerance = 1e-6)
  }
})

test_that("the half-life is where the curve reaches a half", {
  p <- sc_retention_params(1.8 * 365)
  x <- p$half_life / p$zero_days
  expect_equal(exp(-p$k * x^2 / (1 - x^2)), 0.5, tolerance = 1e-9)
})

test_that("retention and the antimalarial decay are the same function", {
  # The decay form fixes l and derives k from the half-life; retention takes
  # (half-life, k) and derives l. Both must be k = log2 ((l/h)^2 - 1).
  for (h in c(300, 657, 1500)) {
    for (k in c(0.25, 1, 95.47)) {
      l <- h / sqrt(1 - k / (k - log(0.5)))
      expect_equal(k, log(2) * ((l / h)^2 - 1))
    }
  }
})

test_that("a smooth-compact curve keeps more of its mean inside a campaign cycle", {
  # Both have a mean of 1.8 years, but the exponential spends part of it past
  # the next campaign, where nets are replaced anyway.
  m <- 1.8 * 365
  p <- sc_retention_params(m)
  cycle <- 3 * 365
  sc <- stats::integrate(function(a) {
    x <- a / p$zero_days
    exp(-p$k * x^2 / (1 - x^2))
  }, 0, cycle)$value
  expon <- m * (1 - exp(-cycle / m))
  expect_gt(sc, expon)
  expect_equal(sc / expon, 1.19, tolerance = 0.01)
})

test_that("impossible means are refused", {
  expect_error(sc_retention_params(21 * 365), "impossible")
  expect_error(sc_retention_params(-1), "one positive number")
  expect_error(sc_retention_params(c(100, 200)), "one positive number")
})
