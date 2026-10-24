test_that("net_usage checks its arguments", {
  S <- function(a) exp(-a / 657)
  expect_error(net_usage(1, c(0.3, 0.3), 1, S), "same length")
  expect_error(net_usage(1, c(0.3, 0.3), c(366, 1), S), "in order")
  expect_error(net_usage(1, 1.2, 1, S), "between 0 and 1")
})

test_that("no distribution yet means no use", {
  S <- function(a) exp(-a / 657)
  expect_identical(net_usage(c(0, 100), 0.5, 200, S), c(0, 0))
})

test_that("one distribution is just the retention curve", {
  S <- function(a) exp(-a / 657)
  t <- c(0, 100, 657, 2000)
  expect_equal(net_usage(t, 0.4, 0, S), 0.4 * S(t))
})

test_that("a second round replaces the first, it does not add to it", {
  S <- function(a) 1   # nothing is ever lost, so only replacement is at work
  # 0.5 then 0.5 of everyone at random: 75% have been reached at least once.
  expect_equal(net_usage(10, c(0.5, 0.5), c(0, 5), S), 0.75)
  # full coverage leaves nobody without
  expect_equal(net_usage(10, c(0.3, 1), c(0, 5), S), 1)
})

test_that("use never exceeds one, however many rounds", {
  S <- function(a) exp(-a / 657)
  u <- net_usage(seq(0, 7300, 50), rep(0.9, 20), seq(0, 19) * 365, S)
  expect_true(all(u >= 0 & u <= 1))
})

test_that("net_usage matches netz exactly for the exponential", {
  skip_if_not_installed("netz")
  # netz's recursion ages the standing stock from the last round, which is
  # exact when the curve is memoryless. This is that case, so the two must
  # agree; the next test is the case where they do not.
  set.seed(1)
  dt <- seq(0, 21) * 365
  cov <- runif(length(dt), 0, 0.6)
  ut <- dt + 183
  S <- function(a) netz::net_loss_exp(a, 657)
  mine <- net_usage(ut, cov, dt, S)
  theirs <- netz::model_distribution_to_usage(
    ut, cov, dt, net_loss_function = netz::net_loss_exp, mean_retention = 657)
  expect_equal(mine, theirs, tolerance = 1e-12)
})

test_that("net_usage parts company with netz under smooth-compact", {
  skip_if_not_installed("netz")
  sc <- sc_retention_params(657)
  dt <- seq(0, 21) * 365
  cov <- rep(0.5, length(dt))
  ut <- dt + 183
  S <- function(a) netz::net_loss_map(a, sc$half_life, sc$k)
  mine <- net_usage(ut, cov, dt, S)
  theirs <- netz::model_distribution_to_usage(
    ut, cov, dt, net_loss_function = netz::net_loss_map,
    half_life = sc$half_life, k = sc$k)
  # netz re-ages every net to zero each round, so it keeps nets that are gone.
  expect_gt(mean(theirs[-1] - mine[-1]), 0.02)
})

test_that("net_usage matches a Monte Carlo of the model's own bookkeeping", {
  sc <- sc_retention_params(657)
  l <- sc$half_life / sqrt(1 - sc$k / (sc$k - log(0.5)))
  dt <- seq(0, 15) * 365
  cov <- rep(0.45, length(dt))
  ut <- dt + 183
  S <- function(a) {
    x <- pmin(a / l, 1)
    ifelse(a >= l, 0, exp(-sc$k * x^2 / (1 - x^2)))
  }
  set.seed(42)
  n <- 200000L
  last <- rep(-Inf, n); life <- rep(0, n); mc <- numeric(length(ut)); j <- 1L
  for (g in seq_along(ut)) {
    while (j <= length(dt) && dt[[j]] <= ut[[g]]) {
      got <- stats::runif(n) < cov[[j]]
      last[got] <- dt[[j]]
      a <- -log(stats::runif(sum(got))) / sc$k        # the model's sampler
      life[got] <- pmin(l * sqrt(a / (1 + a)), l)
      j <- j + 1L
    }
    mc[[g]] <- mean(ut[[g]] - last < life)
  }
  expect_equal(net_usage(ut, cov, dt, S), mc, tolerance = 0.005)
})

test_that("usage_to_distribution inverts net_usage", {
  sc <- sc_retention_params(657)
  S <- function(a) netz::net_loss_map(a, sc$half_life, sc$k)
  dt <- seq(0, 21) * 365
  ut <- dt + 183
  want <- seq(0.05, 0.6, length.out = length(dt))
  cov <- usage_to_distribution(want, ut, dt, S)
  expect_equal(net_usage(ut, cov, dt, S), want, tolerance = 1e-10)
})

test_that("usage_to_distribution respects its bounds", {
  S <- function(a) exp(-a / 657)
  dt <- seq(0, 5) * 365
  ut <- dt + 183
  # unreachable in one round: capped, and the use falls short
  cov <- usage_to_distribution(rep(0.99, length(dt)), ut, dt, S)
  expect_true(all(cov <= 1))
  expect_true(any(cov > 1 - 1e-9))
  expect_true(all(net_usage(ut, cov, dt, S) <= 0.99 + 1e-12))
  # a target already met takes no nets away
  cov0 <- usage_to_distribution(c(0.8, 0.01), c(183, 548), c(0, 365), S)
  expect_equal(cov0[[2]], 0)
})

test_that("usage_to_distribution reproduces netz's answer for the exponential", {
  skip_if_not_installed("netz")
  dt <- seq(0, 21) * 365
  ut <- dt + 183
  set.seed(3)
  want <- stats::runif(length(dt), 0.05, 0.7)
  S <- function(a) netz::net_loss_exp(a, 657)
  mine <- usage_to_distribution(want, ut, dt, S)
  theirs <- netz::usage_to_model_distribution(
    usage = want, usage_timesteps = ut, distribution_timesteps = dt,
    distribution_lower = rep(0, length(dt)), distribution_upper = rep(1, length(dt)),
    net_loss_function = netz::net_loss_exp, mean_retention = 657)
  expect_equal(mine, theirs, tolerance = 1e-6)
})

test_that("usage_to_distribution refuses points it cannot solve in turn", {
  S <- function(a) exp(-a / 657)
  expect_error(usage_to_distribution(c(0.3, 0.3), c(183, 548), c(0, 365, 730), S),
               "same length")
  expect_error(usage_to_distribution(c(0.3, 0.3), c(400, 548), c(0, 365), S),
               "before the next distribution")
  expect_error(usage_to_distribution(c(0.3, 0.3), c(-10, 548), c(0, 365), S),
               "on or after its own distribution")
})
