test_that("the last year of interventions is carried forward unchanged", {
  site_data <- list(interventions = data.frame(
    year = c(2023, 2024),
    itn_input_dist = c(0.5, 0.1),
    tx_cov = c(0.4, 0.45)
  ))
  out <- expand_interventions(site_data, expand_year = 2)$interventions
  expect_equal(out$year, c(2023, 2024, 2025, 2026))
  expect_equal(out$itn_input_dist, c(0.5, 0.1, 0.1, 0.1))
  expect_equal(out$tx_cov, c(0.4, 0.45, 0.45, 0.45))
})

test_that("unimplemented options warn rather than act", {
  site_data <- list(interventions = data.frame(year = 2024, tx_cov = 0.4))
  expect_warning(
    out <- expand_interventions(site_data, expand_year = 1, delay = 1),
    "not implemented"
  )
  expect_equal(out$interventions$year, c(2024, 2025))
})
