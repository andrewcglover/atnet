# Turning the clean tables into the data each Stan model expects.
#
# These follow the data preparation in the drivers the current posteriors were
# fitted with (fit_eip_hill.R, fit_eip.R and fit_atn_main_bidir_splitnH.R), so
# a fit here differs from those only in the data it is given.

#' Stan data for the EIP models
#'
#' Each row of the sporozoite table is one experiment and receives its own
#' random effect. The same list serves both EIP models; the exponential model
#' ignores the quadrature nodes and the Hill-only priors.
#'
#' @param spz The sporozoite table, from [build_sporozoite_table()] or the CSV
#'   it is written to.
#' @param priors A list in the form of [atn_priors()].
#' @param erlang_shape Number of stages in the extrinsic incubation period.
#' @param gl_nodes Number of Gauss-Legendre nodes for the Hill model's
#'   integral.
#'
#' @return A named list to pass to `rstan::sampling()`.
#'
#' @export
eip_stan_data <- function(spz, priors = atn_priors(), erlang_shape = 10L,
                          gl_nodes = 32L) {
  check_columns(spz, c("experiment_id", "hours_after_infection",
                       "dissection_day", "n_control", "n_positive_control",
                       "n_treated", "n_positive_treated"))
  stopifnot(!anyDuplicated(spz$experiment_id),
            all(spz$n_positive_control <= spz$n_control),
            all(spz$n_positive_treated <= spz$n_treated))

  exp_id <- as.integer(factor(spz$experiment_id))
  gl <- gauss_legendre_01(gl_nodes)
  c(
    list(
      N       = nrow(spz),
      K       = length(unique(exp_id)),
      exp_id  = exp_id,
      delta   = spz$hours_after_infection / 24,
      tau_obs = spz$dissection_day,
      n_ctl   = whole(spz$n_control),
      pos_ctl = whole(spz$n_positive_control),
      n_trt   = whole(spz$n_treated),
      pos_trt = whole(spz$n_positive_treated),
      Delta_r = as.integer(erlang_shape),
      gl_M          = length(gl$gl_nodes_01),
      gl_nodes_01   = gl$gl_nodes_01,
      gl_weights_01 = gl$gl_weights_01
    ),
    prior_stan_data(priors)
  )
}

#' Stan data for the blocking models
#'
#' Each experiment is one row of the model, at a single exposure time, with
#' its own random effect. Exposure before infection is side 0 and exposure
#' after it side 1, with `s` the time between them in days.
#'
#' @param ooc The oocyst table, one row per mosquito, from
#'   [build_oocyst_table()] or the CSV it is written to.
#' @param measure `"tra"` for the model of oocyst counts
#'   (`int_fit_nb_rate_hill_bidir_splitnH.stan`), or `"tba_lab"` for the
#'   model of the proportion with any oocyst
#'   (`int_fit_linear_hill_bidir_splitnH.stan`).
#' @param priors A list in the form of [atn_priors()].
#'
#' @return A named list to pass to `rstan::sampling()`.
#'
#' @export
blocking_stan_data <- function(ooc, measure = c("tra", "tba_lab"),
                               priors = atn_priors()) {
  measure <- match.arg(measure)
  check_columns(ooc, c("experiment_id", "arm", "oocysts",
                       "hours_after_infection", "hours_before_infection"))
  stopifnot(all(ooc$arm %in% c("control", "treated")),
            all(ooc$oocysts >= 0),
            all(abs(ooc$hours_after_infection + ooc$hours_before_infection) < 1e-6))

  ids <- sort(unique(ooc$experiment_id))
  row <- match(ooc$experiment_id, ids)
  by_row <- function(x) split(x, factor(row, levels = seq_along(ids)))

  pre_time <- unname(vapply(by_row(ooc$hours_before_infection), function(v) {
    if (length(unique(v)) != 1L) {
      stop("An experiment has more than one exposure time.", call. = FALSE)
    }
    v[[1L]]
  }, numeric(1)))
  arms_ok <- vapply(by_row(ooc$arm), function(a) {
    all(c("control", "treated") %in% a)
  }, logical(1))
  if (!all(arms_ok)) {
    stop("Experiment(s) ", toString(ids[!arms_ok]), " lack a control or ",
         "treated arm.", call. = FALSE)
  }

  rows <- list(
    N      = length(ids),
    K      = length(ids),
    exp_id = seq_along(ids),
    s      = abs(pre_time) / 24,
    side   = as.integer(pre_time < 0)
  )
  y <- whole(ooc$oocysts)
  treated <- ooc$arm == "treated"

  per_mosquito <- if (measure == "tra") {
    list(
      M      = nrow(ooc),
      row_id = row,
      arm    = as.integer(treated),
      y      = y
    )
  } else {
    n <- length(ids)
    list(
      n_ctl   = tabulate(row[!treated], n),
      pos_ctl = tabulate(row[!treated & y > 0L], n),
      n_trt   = tabulate(row[treated], n),
      pos_trt = tabulate(row[treated & y > 0L], n)
    )
  }
  c(rows, per_mosquito, prior_stan_data(priors))
}

#' Gauss-Legendre nodes and weights on \[0, 1\]
#'
#' By the Golub-Welsch algorithm, mapped from \[-1, 1\] so the weights sum to
#' one and `sum(w * f(x))` approximates the integral of `f` over \[0, 1\].
#'
#' @param M Number of nodes.
#'
#' @return A list of `gl_nodes_01` and `gl_weights_01`.
#'
#' @keywords internal
gauss_legendre_01 <- function(M) {
  i        <- seq_len(M - 1)
  off_diag <- i / sqrt(4 * i^2 - 1)
  Tmat     <- matrix(0, M, M)
  Tmat[cbind(i,     i + 1)] <- off_diag
  Tmat[cbind(i + 1, i    )] <- off_diag
  eig <- eigen(Tmat, symmetric = TRUE)
  ord <- order(eig$values)
  out <- list(
    gl_nodes_01   = (eig$values[ord] + 1) / 2,
    gl_weights_01 = (2 * eig$vectors[1, ord]^2) / 2
  )
  stopifnot(isTRUE(all.equal(sum(out$gl_weights_01), 1, tolerance = 1e-12)))
  out
}

check_columns <- function(x, need) {
  missing <- setdiff(need, names(x))
  if (length(missing)) {
    stop("Missing column(s): ", toString(missing), call. = FALSE)
  }
}

whole <- function(x) {
  if (any(is.na(x)) || any(x != round(x))) {
    stop("Counts must be whole numbers.", call. = FALSE)
  }
  as.integer(x)
}
