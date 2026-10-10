# Fake files in the layout of the source repository: 3 draws at 3 resistance
# levels (0, 0.5, 1), ids (k - 1) * 3 + 1 to k * 3 for draw k. The
# pyrethroid-only file is in id order; the other two are merged files whose
# rows are in draw order within each level but not in id order overall.
fake_only <- function() {
  id <- 1:9
  data.frame(id = id, resistance = ((id - 1) %% 3) / 2,
             dn0 = 0.1 * id, rn_pyr = 0.01 * id, mean_duration = id)
}
fake_merged <- function() {
  id <- c(1, 4, 7, 2, 5, 8, 3, 6, 9)
  data.frame(id.x = id, resistance.x = ((id - 1) %% 3) / 2,
             dn0 = 0.1 * id, rn_pbo = 0.01 * id, mean_duration.x = 0,
             id.y = NA, resistance.y = NA, rn_pyr = NA,
             mean_duration.y = NA, mn_dur = id)
}

test_that("draws are extracted as in Part 6 of the source repository", {
  for (x in list(only = fake_only(), pbo = fake_merged())) {
    net <- if ("id" %in% names(x)) "only" else "pbo"
    d <- standardise_itn_draws(x, net)
    id <- (d$draw - 1) * 3 + d$resistance * 2 + 1
    expect_equal(d$draw, rep(1:3, each = 3), info = net)
    expect_equal(d$resistance, rep(c(0, 0.5, 1), 3), info = net)
    expect_equal(d$dn0, 0.1 * id, info = net)
    expect_equal(d$rn0, 0.01 * id, info = net)
    expect_equal(d$gamman, id / log(2), info = net)
  }
})

test_that("a file that breaks the layout is refused", {
  x <- fake_only()
  x$resistance[2] <- 1
  expect_error(standardise_itn_draws(x, "only"), "expected layout")
  x <- fake_merged()[c(2, 1, 3:9), ]
  expect_error(standardise_itn_draws(x, "pbo"), "not in draw order")
})

test_that("the source is pinned to one commit with a checksum per file", {
  src <- itn_draws_source()
  expect_match(src$commit, "^[0-9a-f]{40}$")
  expect_named(src$files, c("only", "pbo", "cfp"))
  expect_named(src$md5, c("only", "pbo", "cfp"))
  expect_true(all(grepl("^[0-9a-f]{32}$", src$md5)))
})
