#' Eclipse-phase target-cell-limited viral kinetics
#'
#' The within-host model used in the SARS-CoV-2 application:
#' \deqn{\dot U = -\beta U V, \quad \dot I_1 = \beta U V - \kappa I_1, \quad
#'       \dot I_2 = \kappa I_1 - \delta I_2, \quad \dot V = p I_2 - cV,}
#' with \eqn{U(0) = T_0}, \eqn{I_1(0) = I_2(0) = 0} and \eqn{V(0) = V_0}.
#' The eclipse compartment \eqn{I_1} is what separates this from the basic
#' target-cell-limited model: without it a newly infected cell produces virus
#' immediately, which the observed rise is too slow to support.
#'
#' @param times Numeric vector of times (days) at which to return viral load.
#' @param beta,delta,p,V0 Free parameters: infection rate constant,
#'   infected-cell loss rate, virion production rate, inoculum.
#' @param kappa,c_clear,T0 Fixed constants. Defaults are the values adopted by
#'   Ke et al. (2021) for the upper respiratory tract. `T0` is a scale
#'   convention rather than a literature value: it enters only through the
#'   products \eqn{\beta T_0} and \eqn{p T_0}, so rescaling it and
#'   compensating in `beta` and `p` leaves the trajectory unchanged.
#' @param grid Internal solver grid.
#'
#' @return A numeric vector of viral load at `times`.
#' @references
#' Baccam P et al. (2006) Kinetics of influenza A virus infection in humans.
#' *J Virol* 80, 7590-7599.
#'
#' Ke R et al. (2021) In vivo kinetics of SARS-CoV-2 infection and its
#' relationship with a person's infectiousness. *PNAS* 118, e2111477118.
#' @seealso [eclipse_design()]
#' @export
#' @examples
#' if (requireNamespace("deSolve", quietly = TRUE)) {
#'   v <- eclipse_trajectory(0:20, beta = 3e-8, delta = 1.7, p = 45, V0 = 10)
#'   round(log10(v[c(1, 6, 11)]), 2)
#' }
eclipse_trajectory <- function(times, beta, delta, p, V0,
                               kappa = 4, c_clear = 10, T0 = 1e7,
                               grid = seq(0, max(times) + 1, by = 0.1)) {
  if (!requireNamespace("deSolve", quietly = TRUE))
    stop("Package 'deSolve' is required for eclipse_trajectory().")
  rhs <- function(t, y, parms) {
    with(as.list(y), {
      dU  <- -beta * U * V
      dI1 <-  beta * U * V - kappa * I1
      dI2 <-  kappa * I1 - delta * I2
      dV  <-  p * I2 - c_clear * V
      list(c(dU, dI1, dI2, dV))
    })
  }
  grid <- sort(unique(c(0, grid, times)))
  out <- deSolve::ode(c(U = T0, I1 = 0, I2 = 0, V = V0), grid, rhs, NULL)
  stats::approx(out[, "time"], pmax(out[, "V"], 0), xout = times, rule = 2)$y
}

#' Candidate design on identifiable viral coordinates
#'
#' Draws candidate trajectories on the coordinates the data determine --
#' doubling time, peak capacity, clearance rate and inoculum -- and derives
#' the ODE rate constants from them. A box placed directly on
#' \eqn{(\beta, \delta, p)} forces \eqn{\log_{10} K \ge 5.65} for every
#' candidate, so an established infection cannot peak low; the fitted
#' population then acquires a spurious non-establishing family in its place.
#'
#' Because the doubling time is sampled positive, every candidate has
#' \eqn{r > 0} and hence \eqn{R_0 > 1}: the class-1 support consists of
#' establishing kinetics **by construction** rather than by filtering fitted
#' output. That is a property of the design, not a finding.
#'
#' @param n Number of candidates.
#' @param log10_K Range of the peak capacity \eqn{K = p T_0 / c}, log10 scale.
#' @param Td Range of the doubling time in days, sampled log-uniformly.
#' @param delta Range of the infected-cell loss rate, per day.
#' @param log10_V0 Range of the inoculum, log10 scale.
#' @param kappa,c_clear,T0 Fixed constants, as in [eclipse_trajectory()].
#' @param seed Optional seed, so a design can be reproduced exactly.
#'
#' @return A data frame of candidates with the identifiable coordinates and
#'   the derived `beta` and `p`.
#' @export
#' @examples
#' head(eclipse_design(5, seed = 1))
eclipse_design <- function(n, log10_K = c(3, 9.5), Td = c(0.3, 3),
                           delta = c(0.05, 5), log10_V0 = c(0, 3),
                           kappa = 4, c_clear = 10, T0 = 1e7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  u <- function(a, b) stats::runif(n, a, b)
  K   <- 10^u(log10_K[1], log10_K[2])
  Tdv <- 10^u(log10(Td[1]), log10(Td[2]))
  dl  <- u(delta[1], delta[2])
  V0  <- 10^u(log10_V0[1], log10_V0[2])
  r   <- log(2) / Tdv
  p   <- K * c_clear / T0
  beta <- (r + kappa) * (r + dl) * (r + c_clear) / (kappa * K * c_clear)
  data.frame(Td = Tdv, K = K, delta = dl, V0 = V0, beta = beta, p = p,
             kappa = kappa, c_clear = c_clear, T0 = T0)
}
