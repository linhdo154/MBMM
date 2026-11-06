#' Integrate a nuisance pool out of a candidate likelihood
#'
#' A mechanistic design often splits into coordinates the question is about and
#' coordinates it is not. In the prime-boost application the class coordinate
#' \eqn{g} and the decay rate \eqn{m_s} are shared and carry the question,
#' while the per-assay levels \eqn{(A_0, C_1)} and the spike ratios are
#' nuisance: they differ between antigen targets and nothing is being asked
#' about them.
#'
#' Those nuisance coordinates must be **integrated out**, not competed against.
#' If they are left as separate bank atoms, a participant can be matched by any
#' candidate that happens to fit, whatever its \eqn{g}, and the class
#' coordinate stops being identified — a flat trajectory then scores equally
#' well at \eqn{g>0} and \eqn{g<0}. Averaging instead asks the right question:
#' *given this shared point, how well does a typical draw of the nuisance
#' coordinates explain this participant?* A flat participant scores badly at
#' \eqn{g>0} because most of the pool predicts a rise there.
#'
#' The result is the inner integral of the two-level construction, and is what
#' [fit_bank()] should be given: the bank's support is the shared grid, not the
#' candidate list.
#'
#' @param loglik A participant-by-candidate log-likelihood matrix, as returned
#'   by [trajectory_loglik()].
#' @param shared A vector of length `ncol(loglik)` giving the shared point each
#'   candidate belongs to. A factor fixes the output column order; anything
#'   else is converted with [factor()].
#' @param weights Optional pool weights, length `ncol(loglik)`, normalised
#'   within each shared point. The default is an equal-weight average, which is
#'   the Monte Carlo estimate of the integral when the pool is drawn from the
#'   nuisance prior.
#'
#' @return A matrix with one row per participant and one column per shared
#'   point, holding \eqn{\log \overline{p}(Y_i \mid s)}. Column names are the
#'   shared point labels.
#'
#' @seealso [jitter_marginal()] for integrating a calibration offset,
#'   [fit_bank()] for the outer problem.
#' @export
#' @examples
#' ll <- matrix(c(-1, -2, -3, -4, -5, -6), nrow = 2)
#' pool_marginal(ll, shared = c("a", "a", "b"))
pool_marginal <- function(loglik, shared, weights = NULL) {
  stopifnot(is.matrix(loglik), ncol(loglik) == length(shared))
  if (!is.factor(shared)) shared <- factor(shared)
  if (!is.null(weights)) {
    stopifnot(length(weights) == length(shared), all(weights >= 0))
  }
  idx <- split(seq_along(shared), shared)
  out <- matrix(NA_real_, nrow(loglik), length(idx),
                dimnames = list(rownames(loglik), names(idx)))
  for (s in seq_along(idx)) {
    j <- idx[[s]]
    if (is.null(weights)) {
      out[, s] <- row_lse(loglik[, j, drop = FALSE]) - log(length(j))
    } else {
      w <- weights[j]
      if (sum(w) <= 0) stop("Shared point ", names(idx)[s], " has no pool weight.")
      out[, s] <- row_lse(add_cols(loglik[, j, drop = FALSE],
                                   log(pmax(w / sum(w), 1e-300))))
    }
  }
  out
}

#' Integrate a calibration offset along a grid coordinate
#'
#' Several assays can measure the same underlying gain on scales that do not
#' quite agree. Rather than give each assay its own class coordinate, which
#' would stop them speaking to the same question, the shared coordinate is read
#' with an assay-specific offset: assay \eqn{t} sees \eqn{g + w\delta_t}, and
#' the offset is integrated out with equal weight.
#'
#' This requires the likelihood on a *wider* grid than the one kept, so that
#' every shifted read stays inside it. That is why the design has an outer grid
#' and an interior shared sub-range: `keep` indexes the latter within the
#' former.
#'
#' @param E A participant-by-column matrix whose columns are the grid
#'   coordinate (running fastest) crossed with any number of other
#'   coordinates.
#' @param n_grid The length of the grid coordinate, which must divide
#'   `ncol(E)`.
#' @param keep Integer indices into `seq_len(n_grid)` to retain: the interior
#'   sub-range on which the result is reported.
#' @param offsets Integer offsets, in grid steps, to average over. Every
#'   `keep + offsets` must lie inside the grid, which is checked.
#' @param weights Optional weights on `offsets`; equal by default.
#'
#' @return A matrix with `length(keep)` times the other coordinates' columns,
#'   the grid coordinate still running fastest.
#'
#' @seealso [pool_marginal()].
#' @export
#' @examples
#' E <- matrix(log(c(1, 2, 3, 4, 5)), nrow = 1)
#' jitter_marginal(E, n_grid = 5, keep = 2:4, offsets = c(-1, 0, 1))
jitter_marginal <- function(E, n_grid, keep, offsets, weights = NULL) {
  stopifnot(is.matrix(E), ncol(E) %% n_grid == 0)
  keep <- as.integer(keep); offsets <- as.integer(offsets)
  if (min(keep) + min(offsets) < 1L || max(keep) + max(offsets) > n_grid)
    stop("`keep + offsets` leaves the grid; widen the grid or shrink `keep`.")
  nd <- length(offsets)
  w  <- if (is.null(weights)) rep(1 / nd, nd) else weights / sum(weights)
  K  <- ncol(E) / n_grid
  out <- matrix(NA_real_, nrow(E), length(keep) * K)
  for (k in seq_len(K)) {
    base <- (k - 1L) * n_grid
    A <- array(NA_real_, c(nrow(E), length(keep), nd))
    for (j in seq_len(nd))
      A[, , j] <- E[, base + keep + offsets[j], drop = FALSE] + log(w[j])
    mx  <- apply(A, c(1, 2), max)
    acc <- mx + log(apply(exp(A - as.vector(mx)), c(1, 2), sum))
    out[, (k - 1L) * length(keep) + seq_along(keep)] <- acc
  }
  out
}
