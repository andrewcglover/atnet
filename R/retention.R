# Net retention under the smooth-compact loss form.
#
# malariasimulation offers two forms for how long a person keeps a net:
#
#   exponential     S(a) = exp(-a / m), mean m, set by set_bednets(retention = m)
#   smooth-compact  S(a) = exp(-k x^2 / (1 - x^2)), x = a / l, zero beyond l,
#                   set by set_bednets(logistic_half_life =, logistic_k =)
#
# The second is the same function as the antimalarial's smooth-compact decay
# (the fork's atn_potency()), written the other way round. Both obey
#
#   k = log2 ((l / h)^2 - 1)
#
# where h is the half-life and l the age by which everything is gone; the decay
# code fixes l and derives k from h, while set_bednets() takes (h, k) and
# derives l. Verified against the fork's own logistic_net_retention_time():
# 4e6 draws match the formula above to 2e-4 on survival.
#
# Retention is estimated as a MEAN (the exponential's parameter), so a
# smooth-compact run is matched to the same mean rather than the same
# half-life: that is the quantity the evidence reports. Fixing l at the decay
# form's own 20 years then leaves one free parameter, found numerically below
# because the mean has no closed form.

#' The age by which a smooth-compact curve has reached zero
#'
#' Shared by the antimalarial's decay and by net retention, so the two
#' conventions cannot drift apart. Matches the fork's `atn_decay_zero_days`
#' default.
#'
#' @return Days.
#'
#' @export
sc_zero_days <- function() 20 * 365

# Mean of the smooth-compact survival curve, l * int_0^1 exp(-k x^2/(1-x^2)) dx.
sc_mean <- function(k, zero_days = sc_zero_days()) {
  zero_days * stats::integrate(
    function(x) exp(-k * x^2 / (1 - x^2)), 0, 1, rel.tol = 1e-11)$value
}

#' Smooth-compact retention parameters for a given mean
#'
#' Finds the smooth-compact net-loss curve that has the same mean retention as
#' an exponential curve with that mean, holding the age at which every net has
#' gone at `zero_days`. The mean has no closed form, so `k` is found by
#' [stats::uniroot()]; the half-life then follows from
#' `h = l sqrt(log2 / (k + log2))`.
#'
#' Pass both returned values to [malariasimulation::set_bednets()]:
#' it recomputes `l` from the pair, so supplying only one, or a pair that does
#' not belong together, silently gives a different curve.
#'
#' @param mean_days The mean retention to match, in days: the value that would
#'   be given to `set_bednets(retention = )` under the exponential form.
#' @param zero_days The age by which every net has gone; by default
#'   [sc_zero_days()].
#'
#' @return A list with `half_life` and `k` (`logistic_half_life` and
#'   `logistic_k`), plus the `mean` achieved and the `zero_days` used.
#'
#' @examples
#' sc_retention_params(1.8 * 365)
#'
#' @export
sc_retention_params <- function(mean_days, zero_days = sc_zero_days()) {
  if (length(mean_days) != 1L || !is.finite(mean_days) || mean_days <= 0) {
    stop("mean_days must be one positive number, not ", deparse(mean_days),
         call. = FALSE)
  }
  if (mean_days >= zero_days) {
    stop("A mean of ", mean_days, " days is impossible with everything gone by ",
         zero_days, " days.", call. = FALSE)
  }
  # sc_mean() falls from zero_days at k = 0 towards 0 as k grows, so the root
  # is unique. The upper bracket covers means down to a few days.
  k <- stats::uniroot(function(k) sc_mean(k, zero_days) - mean_days,
                      c(1e-8, 1e9), tol = .Machine$double.eps^0.5)$root
  list(half_life = zero_days * sqrt(log(2) / (k + log(2))),
       k         = k,
       mean      = sc_mean(k, zero_days),
       zero_days = zero_days)
}
