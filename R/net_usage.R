# Turning a net-use time series into the distributions a run should make.
#
# 06b fits, for each region and year, the share of people handed a net on the
# distribution day so that simulated net USE matches the site file's estimates.
# netz does this with model_distribution_to_usage(), whose recursion carries the
# standing stock forward by
#
#   remaining <- stock * S(time since the PREVIOUS distribution)
#
# i.e. it ages every net still in use as though it had been handed out at the
# last round. That is exact when S is memoryless, S(a + b) = S(a) S(b), which
# is to say for the exponential and for nothing else. Under the smooth-compact
# curve, which is nearly flat at age zero (S(183 days) = 0.94 at a mean of 657),
# the stock is re-aged to zero every year and so never dies: netz then reports
# the target use from roughly half the nets that actually deliver it. Fitted
# that way, a smooth-compact world was handed 0.41-0.49 times the exponential's
# distributions and reached 0.395 use against a 0.565 target.
#
# What follows replaces that recursion with the bookkeeping malariasimulation
# itself does. A distribution gives a net to a random share of people whatever
# they already hold, and a net received replaces one held, so a person is using
# a net exactly when the net from their MOST RECENT receipt is still alive:
#
#   U(t) = sum_j c_j prod_{i > j} (1 - c_i) S(t - t_j)
#
# over the distributions at t_j <= t, with the product over the rounds since.
# Each net is held at its own age, so the result is right for any S; for the
# exponential it reduces to netz's recursion, and the two agree to rounding
# (tested). Confirmed against a Monte Carlo of the individual-based model:
# both give 0.395 where netz gives 0.560.

#' Net use implied by a series of distributions
#'
#' The share of people holding a live net at each of `times`, under
#' malariasimulation's own accounting: a distribution reaches a random share of
#' people regardless of what they hold, a net received replaces one held, and
#' nets are kept for an age drawn from `retained`.
#'
#' @param times Days at which use is wanted.
#' @param coverages Share of people given a net at each distribution.
#' @param distribution_times Day of each distribution, same length as
#'   `coverages`.
#' @param retained A function of age in days giving the share of nets still
#'   held, e.g. `function(a) netz::net_loss_exp(a, 657)`.
#' @return Numeric vector of use, one per element of `times`.
#' @export
net_usage <- function(times, coverages, distribution_times, retained) {
  if (length(coverages) != length(distribution_times)) {
    stop("coverages and distribution_times must be the same length.", call. = FALSE)
  }
  if (is.unsorted(distribution_times)) {
    stop("distribution_times must be in order.", call. = FALSE)
  }
  if (any(coverages < 0 | coverages > 1)) {
    stop("coverages must each be between 0 and 1.", call. = FALSE)
  }
  vapply(times, function(t) {
    j <- which(distribution_times <= t)
    if (!length(j)) return(0)
    # Later rounds replace earlier ones, so a cohort survives into the total
    # only for the people no later round reached.
    not_since <- rev(cumprod(c(1, rev(1 - coverages[j])[-length(j)])))
    sum(coverages[j] * not_since * retained(t - distribution_times[j]))
  }, numeric(1))
}

#' Distributions that deliver a wanted net-use series
#'
#' The inverse of [net_usage()]: the share to hand out at each distribution so
#' that use matches `usage` at `usage_times`. Each distribution is solved in
#' turn, which is exact as long as the use each one is matched to falls before
#' the next distribution -- one use point per round, in order.
#'
#' Writing `R` for the use left at the use point from earlier rounds and `s` for
#' the share of this round's nets still held by then, use is
#' `c s + (1 - c) R`, so `c = (u - R) / (s - R)`: no search is needed. A round
#' asked for more than it can deliver is capped at `upper` (and the shortfall
#' shows up in the use it returns); one whose target is already met is held at
#' `lower` rather than taking nets away.
#'
#' @param usage Wanted net use, one value per distribution.
#' @param usage_times Day each `usage` value is measured on.
#' @param distribution_times Day of each distribution.
#' @param retained A function of age in days giving the share of nets still held.
#' @param lower,upper Bounds on each share, recycled to length.
#' @return Numeric vector of shares, one per distribution.
#' @export
usage_to_distribution <- function(usage, usage_times, distribution_times,
                                  retained, lower = 0, upper = 1) {
  n <- length(distribution_times)
  if (length(usage) != n || length(usage_times) != n) {
    stop("usage, usage_times and distribution_times must be the same length.",
         call. = FALSE)
  }
  if (is.unsorted(distribution_times) || is.unsorted(usage_times)) {
    stop("distribution_times and usage_times must be in order.", call. = FALSE)
  }
  if (any(usage_times < distribution_times)) {
    stop("Each use point must fall on or after its own distribution.", call. = FALSE)
  }
  if (n > 1L && any(usage_times[-n] >= distribution_times[-1])) {
    stop("Each use point must fall before the next distribution, or the ",
         "rounds cannot be solved one at a time.", call. = FALSE)
  }
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)
  cov <- numeric(n)
  for (j in seq_len(n)) {
    t <- usage_times[[j]]
    # Use at t from the rounds before this one, after their own replacements.
    before <- if (j == 1L) 0 else
      net_usage(t, cov[seq_len(j - 1L)], distribution_times[seq_len(j - 1L)], retained)
    s <- retained(t - distribution_times[[j]])
    cov[[j]] <- if (s <= before) upper[[j]] else (usage[[j]] - before) / (s - before)
    cov[[j]] <- min(max(cov[[j]], lower[[j]]), upper[[j]])
  }
  cov
}
