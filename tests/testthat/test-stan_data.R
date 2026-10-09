# Made-up tables in the form the build functions write; no laboratory data.

spz_example <- function() {
  data.frame(
    experiment_id = 1:3,
    hours_after_infection = c(6, 24, 144),
    hours_before_infection = -c(6, 24, 144),
    dissection_day = c(13L, 10L, 13L),
    n_control = c(20, 15, 12), n_positive_control = c(10, 8, 11),
    n_treated = c(18, 16, 12), n_positive_treated = c(2, 5, 9)
  )
}

ooc_example <- function() {
  # Experiment 1 exposed 24 h before infection, experiment 2 72 h after.
  data.frame(
    experiment_id = c(1, 1, 1, 1, 2, 2, 2),
    arm = c("control", "control", "treated", "treated",
            "control", "treated", "treated"),
    oocysts = c(5, 0, 0, 0, 12, 3, 0),
    hours_after_infection = c(-24, -24, -24, -24, 72, 72, 72),
    hours_before_infection = c(24, 24, 24, 24, -72, -72, -72)
  )
}

test_that("the EIP data carry each experiment's counts and timing", {
  d <- eip_stan_data(spz_example())
  expect_equal(d$N, 3L)
  expect_equal(d$K, 3L)
  expect_equal(d$exp_id, 1:3)
  expect_equal(d$delta, c(0.25, 1, 6))
  expect_equal(d$tau_obs, c(13, 10, 13))
  expect_identical(d$pos_trt, c(2L, 5L, 9L))
  expect_identical(d$Delta_r, 10L)
  expect_equal(d$gl_M, 32L)
  expect_equal(d$log_nH_eip_lower, 0.5)
})

test_that("the EIP data refuse impossible counts", {
  x <- spz_example()
  x$n_positive_treated[1] <- 30
  expect_error(eip_stan_data(x))
  x <- spz_example()
  x$n_control[1] <- 20.5
  expect_error(eip_stan_data(x), "whole numbers")
})

test_that("the quadrature integrates polynomials exactly", {
  gl <- gauss_legendre_01(32L)
  expect_equal(sum(gl$gl_weights_01), 1, tolerance = 1e-12)
  expect_equal(sum(gl$gl_weights_01 * gl$gl_nodes_01^5), 1 / 6, tolerance = 1e-12)
})

test_that("the blocking data put each exposure on the right side", {
  d <- blocking_stan_data(ooc_example(), "tra")
  expect_equal(d$N, 2L)
  expect_equal(d$s, c(1, 3))
  # Before infection is side 0, after infection side 1.
  expect_identical(d$side, c(0L, 1L))
  expect_equal(d$M, 7L)
  expect_identical(d$row_id, c(1L, 1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(d$arm, c(0L, 0L, 1L, 1L, 0L, 1L, 1L))
  expect_identical(d$y, c(5L, 0L, 0L, 0L, 12L, 3L, 0L))
})

test_that("the laboratory TBA data count mosquitoes with any oocyst", {
  d <- blocking_stan_data(ooc_example(), "tba_lab")
  expect_identical(d$n_ctl, c(2L, 1L))
  expect_identical(d$pos_ctl, c(1L, 1L))
  expect_identical(d$n_trt, c(2L, 2L))
  expect_identical(d$pos_trt, c(0L, 1L))
  expect_null(d$y)
})

test_that("the blocking data refuse malformed experiments", {
  x <- ooc_example()
  x$hours_before_infection[1] <- 12
  x$hours_after_infection[1] <- -12
  expect_error(blocking_stan_data(x), "more than one exposure time")

  x <- ooc_example()
  x <- x[!(x$experiment_id == 2 & x$arm == "control"), ]
  expect_error(blocking_stan_data(x), "lack a control or treated arm")

  x <- ooc_example()
  x$hours_before_infection[1] <- -24
  expect_error(blocking_stan_data(x))
})
