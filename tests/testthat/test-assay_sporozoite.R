# These tests use small sheets built here, so they need none of the laboratory
# data. The checks against the real workbooks live in analysis/.

test_that("headers give the arm, exposure time and blood meals", {
  expect_equal(parse_assay_header("CTL 72h 1BF"),
               list(arm = "control", exposure_h = 72, n_bloodmeals = 1L))
  expect_equal(parse_assay_header("200mg 72h 2BF"),
               list(arm = "treated", exposure_h = 72, n_bloodmeals = 2L))
})

test_that("a missing blood meal count means one blood meal", {
  expect_equal(parse_assay_header("CTL 24h")$n_bloodmeals, 1L)
})

test_that("exposure in days is converted to hours", {
  expect_equal(parse_assay_header("200mg 6d")$exposure_h, 144)
})

test_that("line breaks within a header are tolerated", {
  expect_equal(parse_assay_header("CTL\n72h\n1BF"),
               list(arm = "control", exposure_h = 72, n_bloodmeals = 1L))
})

test_that("a column that is not an assay group is ignored", {
  expect_true(is.na(parse_assay_header("NOTE: see methods")$arm))
  expect_true(is.na(parse_assay_header("")$arm))
})

test_that("blanks within a column separate sub-groups", {
  expect_equal(split_at_blanks(c("1", "2", NA, "3")), list(c("1", "2"), "3"))
})

test_that("a run of blanks separates one sub-group, not several", {
  expect_equal(split_at_blanks(c("1", NA, NA, "2")), list("1", "2"))
})

test_that("trailing blanks are discarded", {
  expect_equal(split_at_blanks(c("1", "2", NA, NA)), list(c("1", "2")))
  expect_equal(split_at_blanks(c(NA, NA)), list())
})

test_that("a dissected mosquito is a cell holding a number", {
  expect_equal(count_group(c("0", "250", "1000")),
               c(n = 3L, n_positive = 2L))
})

test_that("an annotated zero counts as dissected but not infected", {
  expect_equal(count_group(c("0*", "250")), c(n = 2L, n_positive = 1L))
})

test_that("a sheet is tabulated one row per sub-group", {
  sheet <- data.frame(
    a = c("CTL 72h 1BF", "0", "250", NA, "500"),
    b = c("200mg 72h 1BF", "0", "0", "0", NA),
    stringsAsFactors = FALSE
  )
  out <- tabulate_assay_sheet(sheet, sheet_name = "example")

  expect_equal(nrow(out), 3L)
  expect_equal(out$arm, c("control", "control", "treated"))
  expect_equal(out$subgroup, c(1L, 2L, 1L))
  expect_equal(out$n, c(2L, 1L, 3L))
  expect_equal(out$n_positive, c(1L, 1L, 0L))
  expect_true(all(out$sheet == "example"))
})

test_that("the columns are separated independently of one another", {
  # The blank in the first column must not break the second.
  sheet <- data.frame(
    a = c("CTL 24h", "1", NA, "2"),
    b = c("200mg 24h", "1", "2", "3"),
    stringsAsFactors = FALSE
  )
  out <- tabulate_assay_sheet(sheet)

  expect_equal(out$n[out$arm == "control"], c(1L, 1L))
  expect_equal(out$n[out$arm == "treated"], 3L)
})

test_that("a sheet with no assay columns gives no rows", {
  sheet <- data.frame(a = c("NOTE", "text"), stringsAsFactors = FALSE)
  expect_equal(nrow(tabulate_assay_sheet(sheet)), 0L)
})
