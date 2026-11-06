#' Screen candidates to the ones participants actually support
#'
#' A candidate family large enough to cover the observed data is far larger
#' than the set any one participant is compatible with. Retaining the `nu`
#' best candidates per participant and pooling them keeps the support the data
#' speak to and makes the weight estimation tractable.
#'
#' The screen uses the measurements, so the pooled support is a function of
#' the donor data. That is a data-dependent step: it must be run on the donors
#' only, never on a participant who will later be scored against the result,
#' or the honesty condition fails. It also means the support is random, which
#' the fixed-support consistency theory does not cover.
#'
#' @param loglik A participant-by-candidate log-likelihood matrix.
#' @param nu Number of candidates retained per participant. Larger `nu` keeps
#'   more of the family at the cost of a bigger weight problem; the retained
#'   count is not usually binding, because the maximiser on a fixed support
#'   has at most one positive weight per participant.
#'
#' @return A list with the pooled candidate indices `keep`, the screened
#'   matrix `loglik`, and the number of candidates `before` and `after`.
#' @seealso [fit_bank()], [screen_positive()]
#' @export
#' @examples
#' set.seed(1)
#' LL <- matrix(rnorm(30 * 200), 30, 200)
#' s <- screen_candidates(LL, nu = 20)
#' c(before = s$before, after = s$after)
screen_candidates <- function(loglik, nu = 200L) {
  stopifnot(is.matrix(loglik), nu >= 1L)
  nu <- min(as.integer(nu), ncol(loglik))
  keep <- sort(unique(unlist(lapply(seq_len(nrow(loglik)), function(i) {
    order(loglik[i, ], decreasing = TRUE)[seq_len(nu)]
  }))))
  list(keep = keep, loglik = loglik[, keep, drop = FALSE],
       before = ncol(loglik), after = length(keep), nu = nu)
}

#' Reject candidates whose trajectory is not positive
#'
#' A sum of exponentials with mixed signs is not constrained to be positive,
#' so a candidate can dip below zero between observations. The likelihood
#' sees only the observation times, so such an excursion costs nothing there
#' and the candidate survives into the fitted population, where it appears in
#' posterior summaries as a dive to the floor and back.
#'
#' The grid must be dense enough to catch an excursion *between* visits;
#' evaluating on the observation times defeats the purpose.
#'
#' @param candidates A data frame of candidate parameter values.
#' @param trajectory A function `(time, theta)` as in [trajectory_loglik()].
#'   It must be called in a way that returns the raw solution, not one already
#'   floored at a positive value.
#' @param grid Dense time grid on which positivity is checked.
#' @param floor Smallest value treated as positive.
#'
#' @return A list with the logical vector `keep`, the screened `candidates`,
#'   and the `rejected` count.
#' @seealso [screen_candidates()]
#' @export
#' @examples
#' cand <- data.frame(a = c(1, 1), b = c(0.5, 5))
#' f <- function(t, th) th$a - th$b * t          # the second goes negative
#' screen_positive(cand, f, grid = seq(0, 1, by = 0.05))$keep
screen_positive <- function(candidates, trajectory,
                            grid = seq(0, 60, by = 1), floor = 0.01) {
  keep <- vapply(seq_len(nrow(candidates)), function(m) {
    v <- trajectory(grid, candidates[m, , drop = FALSE])
    all(is.finite(v)) && min(v) > floor
  }, logical(1))
  list(keep = keep, candidates = candidates[keep, , drop = FALSE],
       rejected = sum(!keep),
       fraction_rejected = mean(!keep))
}
