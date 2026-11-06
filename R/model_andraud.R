#' Reduced prime-boost antibody kinetics
#'
#' A reduced form of the plasma-cell and antibody model of Andraud et al.
#' (2012), in which antibody is produced by a short-lived and a long-lived
#' source and cleared at rate \eqn{m_A}:
#' \deqn{\dot S = -m_s S, \qquad \dot A = W_l + S - m_A A,}
#' with \eqn{S(0) = W_s} and \eqn{A(0) = A_0}. Solving gives a sum of three
#' exponentials, and at the boost the two fluxes receive increments, so the
#' trajectory is continuous there by construction.
#'
#' The three levels the classes are built on are read straight off the model.
#' \eqn{A_0} is the initial condition; setting \eqn{\dot A = 0} once the
#' short-lived source is spent gives the prime-only plateau
#' \eqn{C_1 = W_l / m_A}; the same argument after the boost gives
#' \eqn{D_1 = (W_l + \Delta W_l)/m_A}. The transient peak four weeks after the
#' boost reflects \eqn{\Delta W_s} and decays, so it is not the durable level.
#'
#' Holding \eqn{W_l} constant is what makes the plateaus well defined. It is
#' reasonable over a study of this length: antibody clearance is fast, yet
#' antibody to vaccine antigens persists for decades, which given that
#' clearance rate requires continuing production.
#'
#' @param t Times in months since the prime.
#' @param A0,C1,C2,D1,D2 Levels on the titre scale: baseline, prime-only
#'   plateau and its short-lived amplitude, post-boost plateau and its
#'   short-lived amplitude.
#' @param ms Decay rate of the short-lived source, per month.
#' @param mA Antibody clearance rate, per month. The default corresponds to an
#'   IgG half-life of about 26 days.
#' @param t_boost Time of the boost, months since the prime.
#' @param floor_at Value at which non-positive predictions are floored. A sum
#'   of exponentials with mixed signs is not constrained positive; pass `-Inf`
#'   to inspect the raw solution, which is what a positivity screen needs.
#'
#' @return A numeric vector of antibody concentration at `t`.
#' @references Andraud M et al. (2012) Living on three time scales: the
#'   dynamics of plasma cell and antibody populations. *PLoS Comput Biol* 8,
#'   e1002418.
#' @seealso [amp_to_flux()], [andraud_levels()]
#' @export
#' @examples
#' round(andraud_trajectory(c(0, 1, 2.8, 3.7, 14), A0 = 400, C1 = 1200,
#'                          C2 = 0, D1 = 3000, D2 = 2500, ms = 0.3))
andraud_trajectory <- function(t, A0, C1, C2, D1, D2, ms,
                               mA = log(2) / (26 / 30), t_boost = 84 / 30,
                               floor_at = 0.01) {
  if (abs(mA - ms) < 1e-10) ms <- ms + 1e-8
  f <- amp_to_flux(A0, C1, C2, D1, D2, ms, mA, t_boost)
  C3 <- A0 - C1 - C2
  AB <- C1 + C2 * exp(-ms * t_boost) + C3 * exp(-mA * t_boost)
  D3 <- AB - D1 - D2
  out <- ifelse(
    t < t_boost,
    C1 + C2 * exp(-ms * t) + C3 * exp(-mA * t),
    D1 + D2 * exp(-ms * (t - t_boost)) + D3 * exp(-mA * (t - t_boost)))
  if (is.finite(floor_at)) pmax(out, floor_at) else out
}

#' Convert between amplitudes and production fluxes
#'
#' The implementation samples the levels the data speak to --
#' \eqn{(A_0, C_1, C_2, D_1, D_2)}, all in titre units -- and derives the
#' fluxes. The fluxes are not in titre units: their scale depends on
#' \eqn{m_A} and \eqn{m_s}, so a box on them induces an uncontrolled range of
#' attainable curves. The map is a bijection for fixed \eqn{(m_A, m_s)}, so
#' this is a reparameterisation and not a different model.
#'
#' @inheritParams andraud_trajectory
#' @param Ws,Wl,dWs,dWl Production fluxes and their increments at the boost.
#' @return A list of the converted parameters.
#' @export
#' @examples
#' f <- amp_to_flux(400, 1200, 0, 3000, 2500, ms = 0.3)
#' str(flux_to_amp(f$A0, f$Ws, f$Wl, f$dWs, f$dWl, f$ms))
amp_to_flux <- function(A0, C1, C2, D1, D2, ms,
                        mA = log(2) / (26 / 30), t_boost = 84 / 30) {
  list(A0 = A0, Ws = C2 * (mA - ms), Wl = C1 * mA,
       dWs = (mA - ms) * (D2 - C2 * exp(-ms * t_boost)),
       dWl = (D1 - C1) * mA, ms = ms)
}

#' @rdname amp_to_flux
#' @export
flux_to_amp <- function(A0, Ws, Wl, dWs, dWl, ms,
                        mA = log(2) / (26 / 30), t_boost = 84 / 30) {
  Ws2 <- Ws * exp(-ms * t_boost) + dWs
  list(A0 = A0, C1 = Wl / mA, C2 = Ws / (mA - ms),
       D1 = (Wl + dWl) / mA, D2 = Ws2 / (mA - ms), ms = ms)
}

#' Class coordinate for a prime-boost question
#'
#' Returns the log-ratio that defines a class region. With `ref = "C1"` the
#' coordinate is \eqn{g = \log_{10}(D_1/C_1)}, asking whether the *boost*
#' added a durable gain; with `ref = "A0"` it is \eqn{g = \log_{10}(D_1/A_0)},
#' asking whether the *regimen* left one. In both cases class 1 is
#' \eqn{\{g > 0\}}, so zero is the only boundary and no threshold is imposed.
#'
#' The reference level belongs to the candidate design and must be fixed
#' before fitting. Two questions therefore need two fitted populations; they
#' are not two regions carved from one.
#'
#' @param A0,C1,D1 Baseline, prime-only plateau, post-boost plateau.
#' @param ref Either `"C1"` or `"A0"`.
#' @return A numeric vector of \eqn{g}.
#' @export
#' @examples
#' andraud_levels(A0 = 400, C1 = 1200, D1 = 3000, ref = "C1")
#' andraud_levels(A0 = 400, C1 = 1200, D1 = 3000, ref = "A0")
andraud_levels <- function(A0, C1, D1, ref = c("C1", "A0")) {
  ref <- match.arg(ref)
  denom <- if (ref == "C1") C1 else A0
  log10(pmax(D1, 1e-12) / pmax(denom, 1e-12))
}
