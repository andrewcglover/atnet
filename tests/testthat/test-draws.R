# A fake set of post-warmup draws in which every parameter of iteration i of
# chain c equals 10 * i + c, so a row taken from one whole iteration has all
# its parameters equal and recoverable from its chain and iteration.
fake_draws <- function(pars, n_iter = 50L, n_chain = 4L) {
  vals <- outer(seq_len(n_iter), seq_len(n_chain), function(i, c) 10 * i + c)
  array(rep(vals, length(pars)), dim = c(n_iter, n_chain, length(pars)),
        dimnames = list(iterations = NULL, chains = NULL, parameters = pars))
}
eip_fake <- fake_draws(c("s_half", "nH", "r0_frac", "lp__"))
tra_fake <- fake_draws(c("s_half_pre", "nH_pre", "b_max", "s_half_post",
                         "nH_post", "lp__"))

test_that("the projection priors give the summaries in the SI priors table", {
  pr <- projection_priors()
  q <- c(0.025, 0.5, 0.975)
  expect_equal(round(qbeta(q, pr$omega_atn$shape1, pr$omega_atn$shape2), 2),
               c(0.80, 0.91, 0.97))
  expect_equal(round(qbeta(q, pr$p_atn$atn$shape1, pr$p_atn$atn$shape2), 3),
               c(0.800, 0.950, 0.996))
  expect_equal(round(qbeta(q, pr$p_atn$aitn$shape1, pr$p_atn$aitn$shape2), 2),
               c(0.70, 0.90, 0.99))
  expect_equal(round(qgamma(q, shape = pr$half_life_years$shape,
                            scale = pr$half_life_years$scale), 2),
               c(1.00, 2.50, 5.04))
})

test_that("each row keeps whole iterations, sampled without replacement", {
  d <- withr::with_seed(1, draw_antimalarial_inputs(eip_fake, tra_fake, 200))
  expect_equal(nrow(d), 200)
  expect_equal(d$draw, 1:200)
  eip_cols <- c("s_half_eip", "nH_eip", "rho_frac")
  tra_cols <- c("s_half_pre", "nH_pre", "B_max_post", "s_half_post", "nH_post")
  for (col in eip_cols) {
    expect_equal(d[[col]], 10 * d$eip_iter + d$eip_chain, info = col)
  }
  for (col in tra_cols) {
    expect_equal(d[[col]], 10 * d$tra_iter + d$tra_chain, info = col)
  }
  expect_false(anyDuplicated(d[c("eip_chain", "eip_iter")]) > 0)
  expect_false(anyDuplicated(d[c("tra_chain", "tra_iter")]) > 0)
  # The two fits are sampled independently, not at the same iterations.
  expect_false(identical(d$eip_iter, d$tra_iter))
})

test_that("the ATN and AITN values of p_A are the same quantile", {
  d <- withr::with_seed(2, draw_antimalarial_inputs(eip_fake, tra_fake, 100))
  pr <- projection_priors()$p_atn
  expect_equal(pbeta(d$p_atn_atn, pr$atn$shape1, pr$atn$shape2), d$p_atn_u)
  expect_equal(pbeta(d$p_atn_aitn, pr$aitn$shape1, pr$aitn$shape2), d$p_atn_u)
  expect_true(all(d$p_atn_aitn < d$p_atn_atn))
})

test_that("the same seed gives the same table", {
  a <- withr::with_seed(3, draw_antimalarial_inputs(eip_fake, tra_fake, 50))
  b <- withr::with_seed(3, draw_antimalarial_inputs(eip_fake, tra_fake, 50))
  expect_identical(a, b)
})

test_that("asking for more draws than the fit has, or a missing parameter, errors", {
  expect_error(draw_antimalarial_inputs(eip_fake, tra_fake, 201),
               "Asked for 201 draws")
  expect_error(draw_antimalarial_inputs(eip_fake[, , -2], tra_fake, 10),
               "no parameter nH")
})
