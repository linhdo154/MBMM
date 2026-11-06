#' Profile the residual scale
#'
#' The residual scale is not a nuisance that can be set aside: it enters every
#' class likelihood, and a smaller value sharpens them all, making the Bayes
#' factors more extreme and the responsibilities more confident. The gate is
#' then fitted against something closer to a hard label, which attenuates its
#' coefficient. Estimating the scale with the weights, rather than fixing it,
#' is what keeps that from being an analyst's choice.
#'
#' For each value on the grid the bank is refitted and its attained objective
#' recorded, giving a profile whose maximiser is the jointly estimated scale.
#'
#' @param obs,candidates,trajectory As in [trajectory_loglik()].
#' @param sigma_grid Values of the residual scale to profile over.
#' @param lambda Entropy penalty passed to [fit_bank()].
#' @param lags,lag_prior Optional lag integration, as in
#'   [trajectory_loglik()].
#' @param ... Further arguments to [fit_bank()].
#'
#' @return A list with the `profile` data frame (`sigma` and `objective`), the
#'   maximising `sigma`, and the `bank` fitted there.
#' @seealso [fit_bank()]
#' @export
#' @examples
#' set.seed(1)
#' obs <- data.frame(id = rep(1:8, each = 3), time = rep(c(1, 4, 9), 8),
#'                   value = pmax(rlnorm(24, log(200), 1), 1))
#' cand <- data.frame(peak = c(100, 300, 900), rate = c(.1, .2, .3))
#' decay <- function(t, th) th$peak * exp(-th$rate * t)
#' profile_sigma(obs, cand, decay, sigma_grid = c(0.3, 0.5, 0.8))$sigma
profile_sigma <- function(obs, candidates, trajectory, sigma_grid,
                          lambda = 1, lags = NULL, lag_prior = NULL, ...) {
  stopifnot(length(sigma_grid) >= 1, all(sigma_grid > 0))
  fits <- lapply(sigma_grid, function(s) {
    LL <- trajectory_loglik(obs, candidates, trajectory, sigma = s,
                            lags = lags, lag_prior = lag_prior)
    fit_bank(LL, lambda = lambda, ...)
  })
  obj <- vapply(fits, function(f) f$objective, numeric(1))
  k <- which.max(obj)
  list(profile = data.frame(sigma = sigma_grid, objective = obj),
       sigma = sigma_grid[k], bank = fits[[k]])
}
