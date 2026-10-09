test_that("the dissection day is read from the worksheet name", {
  expect_equal(dissection_day("6h 72h 2BF post 09.19.25 spz10"), 10L)
  expect_equal(dissection_day("#8 72h post spz13"), 13L)
  expect_true(is.na(dissection_day("#5 6h 72h post int")))
})

test_that("timing is counted forward from the infectious blood meal", {
  expect_equal(exposure_sign("#5 6h 72h post int"), 1)
  expect_equal(exposure_sign("#3 12h pre int"), -1)
  expect_true(is.na(exposure_sign("#1 0h int")))
})

test_that("arms are paired within a worksheet into one row per experiment", {
  tab <- data.frame(
    sheet = rep("w", 4),
    column = c("CTL 72h 1BF", "200mg 72h 1BF", "CTL 72h 2BF", "200mg 72h 2BF"),
    arm = c("control", "treated", "control", "treated"),
    concentration = c(NA, 200, NA, 200),
    exposure_h = c(72, 72, 72, 72),
    duration_min = NA_real_,
    n_bloodmeals = c(1L, 1L, 2L, 2L),
    subgroup = c(1L, 1L, 1L, 1L),
    n = c(10L, 11L, 12L, 13L),
    n_positive = c(5L, 1L, 6L, 2L),
    stringsAsFactors = FALSE
  )
  out <- pair_arms(tab)

  expect_equal(nrow(out), 2L)
  expect_equal(out$n_bloodmeals, c(1L, 2L))
  expect_equal(out$n_control, c(10L, 12L))
  expect_equal(out$n_treated, c(11L, 13L))
})

test_that("a treated arm at another concentration is not paired", {
  tab <- data.frame(
    sheet = rep("w", 3), column = c("CTL", "200mg", "100mg"),
    arm = c("control", "treated", "treated"),
    concentration = c(NA, 200, 100), exposure_h = 0, duration_min = NA_real_,
    n_bloodmeals = 1L, subgroup = 1L, n = c(20L, 20L, 20L),
    n_positive = c(18L, 0L, 6L), stringsAsFactors = FALSE
  )
  out <- pair_arms(tab, concentration = 200)

  expect_equal(nrow(out), 1L)
  expect_equal(out$n_positive_treated, 0L)
})

test_that("an unpaired arm yields no row rather than a half row", {
  tab <- data.frame(
    sheet = "w", column = "CTL 72h 1BF", arm = "control",
    concentration = NA_real_, exposure_h = 72,
    duration_min = NA_real_, n_bloodmeals = 1L, subgroup = 1L, n = 10L,
    n_positive = 5L, stringsAsFactors = FALSE
  )
  expect_equal(nrow(pair_arms(tab)), 0L)
})

test_that("a sub-group repeating an earlier one is reported", {
  expect_equal(repeated_subgroups(list(c("1", "2", "3"), c("1", "2", "3"))), 2L)
  expect_equal(repeated_subgroups(list(c("1", "2"), c("3", "4"))), integer(0))
})

test_that("columns taking a single value are not called repeats", {
  # Two arms in which nothing was infected match each other trivially.
  expect_equal(repeated_subgroups(list(c("0", "0", "0"), c("0", "0", "0"))),
               integer(0))
})

test_that("every correction carries its evidence", {
  cr <- assay_corrections()
  expect_true(all(nzchar(cr$evidence)))
  expect_true(all(cr$applies_to %in% c("sporozoite", "oocyst", "both")))
})

test_that("the workbooks are named for when they were received", {
  wb <- assay_workbooks()
  expect_true(all(grepl("^[0-9]{12}_", wb$file)))
})
