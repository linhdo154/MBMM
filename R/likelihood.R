#' Log-likelihood of one measurement under a trajectory
#'
#' Observations are modelled on the \eqn{\log_{10}} scale. A detected
#' measurement contributes a normal density about the trajectory; a measurement
#' below the assay limit contributes the corresponding tail probability, so
#' left-censoring is handled exactly rather than by imputation.
#'
#' @param y Observed value on the natural (not log) scale.
#' @param mu Trajectory value at the observation time, natural scale.
#' @param sigma Residual standard deviation on the \eqn{\log_{10}} scale.
#' @param lod Assay limit of detection, natural scale. Used only when
#'   `censored` is `TRUE`.
#' @param censored Logical; `TRUE` if the measurement is below `lod`.
#' @param eps Small constant guarding `log10()` against zero.
#'
#' @return A numeric vector of log-likelihood contributions.
#' @export
#' @examples
#' obs_loglik(y = 1000, mu = 1200, sigma = 0.5)
#' obs_loglik(y = NA, mu = 50, sigma = 0.5, lod = 100, censored = TRUE)
obs_loglik <- function(y, mu, sigma, lod = NULL, censored = FALSE,
                       eps = 1e-12) {
  stopifnot(all(sigma > 0))
  n <- max(length(y), length(mu), length(censored),
           if (is.null(lod)) 0L else length(lod))
  y        <- rep_len(y, n)
  mu       <- rep_len(mu, n)
  sigma    <- rep_len(sigma, n)
  censored <- rep_len(censored, n)
  if (!is.null(lod)) lod <- rep_len(lod, n)
  out <- numeric(n)
  det <- !censored
  if (any(det)) {
    out[det] <- stats::dnorm(log10(y[det] + eps),
                             mean = log10(mu[det] + eps),
                             sd = sigma[det], log = TRUE)
  }
  if (any(censored)) {
    if (is.null(lod)) stop("`lod` is required when any observation is censored.")
    out[censored] <- stats::pnorm(log10(lod[censored]),
                                  mean = log10(mu[censored] + eps),
                                  sd = sigma[censored], log.p = TRUE)
  }
  out
}

#' Participant-by-candidate log-likelihood matrix
#'
#' Scores every participant against every candidate trajectory, giving the
#' matrix \eqn{\log p_k(Y_i \mid \theta_m)} that [fit_bank()] and
#' [class_loglik()] consume. This is the only place the mechanistic model
#' enters: `trajectory` can be any function of time and a parameter row.
#'
#' @param obs A data frame with one row per measurement and columns `id`,
#'   `time`, `value`, and optionally `censored` and `lod`.
#' @param candidates A data frame or matrix of candidate parameter values, one
#'   row per candidate.
#' @param trajectory A function `(time, theta)` returning the trajectory value
#'   at `time` for the single candidate `theta` (one row of `candidates`).
#' @param sigma Residual standard deviation on the \eqn{\log_{10}} scale.
#' @param lags Optional numeric vector of candidate time shifts, in the units
#'   of `obs$time`. When the clock on which measurements are recorded is not
#'   the biological clock -- days since a reported symptom onset rather than
#'   days since infection -- the offset is unknown and must be integrated out
#'   rather than fixed. Each observation time is evaluated at `time + lag`.
#' @param lag_prior Prior weights over `lags`, normalised internally. Defaults
#'   to uniform.
#'
#' @return A numeric matrix with one row per participant (named by `id`) and
#'   one column per candidate, with any lag already marginalised out.
#' @seealso [fit_bank()], [class_loglik()], [background_loglik()]
#' @export
#' @examples
#' obs <- data.frame(id = rep(1:2, each = 3), time = rep(c(1, 3, 6), 2),
#'                   value = c(100, 60, 20, 500, 300, 90))
#' cand <- data.frame(rate = c(0.1, 0.3), peak = c(200, 600))
#' decay <- function(t, th) th$peak * exp(-th$rate * t)
#' trajectory_loglik(obs, cand, decay, sigma = 0.3)
#'
#' ## with an unknown infection time, integrated over a one-day grid
#' trajectory_loglik(obs, cand, decay, sigma = 0.3, lags = -2:2)
trajectory_loglik <- function(obs, candidates, trajectory, sigma,
                              lags = NULL, lag_prior = NULL) {
  stopifnot(all(c("id", "time", "value") %in% names(obs)))
  if (is.null(obs$censored)) obs$censored <- FALSE
  if (is.null(obs$lod))      obs$lod      <- NA_real_
  if (is.null(lags)) lags <- 0
  lp <- if (is.null(lag_prior)) rep(1 / length(lags), length(lags)) else {
    stopifnot(length(lag_prior) == length(lags), all(lag_prior >= 0))
    lag_prior / sum(lag_prior)
  }
  ids <- unique(obs$id)
  M   <- nrow(candidates)
  LL  <- matrix(0, length(ids), M, dimnames = list(as.character(ids), NULL))
  idx <- lapply(ids, function(i) which(obs$id == i))
  for (m in seq_len(M)) {
    theta <- candidates[m, , drop = FALSE]
    bylag <- vapply(seq_along(lags), function(k) {
      mu <- trajectory(obs$time + lags[k], theta)
      ll <- obs_loglik(obs$value, mu, sigma, obs$lod, obs$censored)
      vapply(idx, function(ii) sum(ll[ii]), numeric(1)) + log(lp[k])
    }, numeric(length(ids)))
    LL[, m] <- if (length(lags) == 1L) as.vector(bylag) else
      row_lse(matrix(bylag, nrow = length(ids)))
  }
  LL
}

#' Background (class-0) log-likelihood at a fraction of the assay limit
#'
#' The absence of the modelled process is not itself a mechanistic
#' trajectory. In the viral application, class 0 is specified rather than
#' estimated: measurements are taken to fluctuate about a floor set by the
#' assay, by default half the limit of detection.
#'
#' Because the floor is a fixed fraction of the limit, every censored
#' measurement contributes the same constant
#' \eqn{\log \Phi(-\log_{10}(\mathrm{floor})/\sigma)} whatever limit was in
#' force. A record consisting only of censored measurements therefore has a
#' class-0 likelihood that depends on how many visits there were and not on
#' which visits they were, so its Bayes factor is driven entirely by the
#' class-1 term.
#'
#' @param obs As in [trajectory_loglik()]; needs `lod` for censored rows.
#' @param sigma Residual standard deviation on the \eqn{\log_{10}} scale.
#'   Use the same value as the class-1 likelihood: the two classes share one
#'   scale, so sensitivity to it acts on both.
#' @param floor Fraction of the limit of detection taken as the background
#'   level. `0.5` reproduces the analysis in the accompanying manuscript.
#'
#' @return A named numeric vector of log class-0 likelihoods.
#' @export
#' @examples
#' obs <- data.frame(id = c(1, 1, 2, 2), time = c(1, 3, 1, 3),
#'                   value = c(NA, 500, NA, NA), lod = 100,
#'                   censored = c(TRUE, FALSE, TRUE, TRUE))
#' background_loglik(obs, sigma = 1.19)
background_loglik <- function(obs, sigma, floor = 0.5) {
  stopifnot(all(c("id", "value") %in% names(obs)), floor > 0)
  if (is.null(obs$censored)) obs$censored <- FALSE
  if (is.null(obs$lod)) stop("`lod` is required to place the background floor.")
  mu  <- floor * obs$lod
  ll  <- obs_loglik(obs$value, mu, sigma, obs$lod, obs$censored)
  ids <- unique(obs$id)
  out <- vapply(ids, function(i) sum(ll[obs$id == i]), numeric(1))
  stats::setNames(out, as.character(ids))
}

#' Class likelihood from a bank
#'
#' Integrates the trajectory likelihood over an estimated kinetic population,
#' \eqn{\widehat L_{ik} = \sum_m \widehat p_{km} \, p_k(Y_i \mid \theta_{km})}.
#' Returned on the log scale.
#'
#' @param loglik A participant-by-candidate log-likelihood matrix, as returned
#'   by [trajectory_loglik()].
#' @param weights Bank weights summing to one, one per candidate.
#'
#' @return A named numeric vector of log class likelihoods, one per participant.
#' @export
class_loglik <- function(loglik, weights) {
  stopifnot(ncol(loglik) == length(weights),
            abs(sum(weights) - 1) < 1e-6)
  row_lse(add_cols(loglik, log(pmax(weights, 1e-300))))
}

## row-wise log-sum-exp, and add a vector across columns
row_lse <- function(x) {
  m <- x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]
  m + log(rowSums(exp(x - m)))
}
add_cols <- function(x, v) x + rep(v, each = nrow(x))
