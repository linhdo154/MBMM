#' Estimate bank weights by entropy-penalised maximum likelihood
#'
#' Given candidate trajectories and the log-likelihood of every participant
#' under each, estimates how often each candidate occurs in the population.
#' The objective is the marginal log-likelihood of the donor measurements plus
#' an entropy penalty,
#' \deqn{J_\lambda(p) = \sum_{j} \log \sum_m p_m \, p(Y_j \mid \theta_m)
#'       + \lambda H(p), \qquad H(p) = -\sum_m p_m \log p_m .}
#'
#' With `lambda = 0` this is the nonparametric maximum likelihood estimator on
#' a fixed support, whose solution is extremely sparse and unstable when many
#' near-identical trajectories compete for the same donors. The entropy term
#' spreads weight across such near-duplicates. For `lambda > 0`, the entropy
#' term is strictly concave and the maximum is unique in the relative interior
#' of the simplex. At `lambda = 0`, a maximum need not be unique.
#'
#' The objective is a **sum** over donors, not a mean. Raising `lambda` with a
#' fixed donor set therefore strengthens the penalty relative to a fixed
#' amount of data, but holding `lambda` fixed as donors accumulate lets the
#' penalty become negligible.
#'
#' @param loglik A donor-by-candidate log-likelihood matrix, as returned by
#'   [trajectory_loglik()]. Donors must not include any participant who will
#'   later be scored against the resulting bank (see Details).
#' @param lambda Entropy penalty, non-negative.
#' @param start Optional starting weights; defaults to uniform.
#' @param maxit,tol Optimiser controls.
#'
#' @details
#' **Honesty.** The bank used to score a participant must be built without
#' that participant's measurements. Pass a disjoint donor set, or use
#' [bank_loo()] to refit leaving each scored participant out in turn.
#'
#' @return An object of class `mbmm_bank`: a list with `weights`, the
#'   donor-by-candidate posterior `responsibility`, the attained `objective`,
#'   `lambda`, and convergence information.
#' @seealso [bank_loo()], [restrict_bank()], [class_loglik()]
#' @export
#' @examples
#' set.seed(1)
#' LL <- matrix(rnorm(20 * 6), 20, 6)
#' b <- fit_bank(LL, lambda = 1)
#' round(b$weights, 3)
fit_bank <- function(loglik, lambda = 1, start = NULL,
                     maxit = 5000L, tol = 1e-8) {
  stopifnot(is.matrix(loglik), lambda >= 0)
  N <- nrow(loglik); M <- ncol(loglik)
  if (M < 2L) stop("A bank needs at least two candidates.")
  w0 <- if (is.null(start)) rep(1 / M, M) else {
    stopifnot(length(start) == M, all(start > 0)); start / sum(start)
  }
  from_logits <- function(x) {
    z <- c(x, 0); z <- z - max(z); e <- exp(z); e / sum(e)
  }
  calc <- function(x) {
    w  <- from_logits(x)
    lw <- log(pmax(w, .Machine$double.xmin))
    joint <- add_cols(loglik, lw)
    mx  <- joint[cbind(seq_len(N), max.col(joint, ties.method = "first"))]
    sc  <- exp(joint - mx); tot <- rowSums(sc); post <- sc / tot
    obj <- sum(mx + log(tot)) + lambda * (-sum(w * lw))
    gw  <- colSums(post) / pmax(w, .Machine$double.xmin) - lambda * (lw + 1)
    list(objective = obj, gradient = (w * (gw - sum(w * gw)))[-M],
         weights = w, posterior = post)
  }
  cache <- new.env(parent = emptyenv())
  ev <- function(x) {
    if (is.null(cache$x) || !identical(cache$x, x)) {
      cache$x <- x; cache$v <- calc(x)
    }
    cache$v
  }
  opt <- stats::optim(
    log(w0[-M]) - log(w0[M]),
    function(x) -ev(x)$objective,
    function(x) -ev(x)$gradient,
    method = "L-BFGS-B",
    control = list(maxit = as.integer(maxit),
                   factr = tol / .Machine$double.eps, pgtol = 0, lmm = 20L))
  f <- calc(opt$par)
  structure(list(weights = f$weights, responsibility = f$posterior,
                 objective = f$objective, lambda = lambda,
                 converged = opt$convergence == 0L,
                 iterations = unname(opt$counts[["function"]])),
            class = "mbmm_bank")
}

#' @export
print.mbmm_bank <- function(x, ...) {
  eff <- exp(-sum(x$weights * log(pmax(x$weights, 1e-300))))
  cat("<mbmm_bank>\n")
  cat(sprintf("  candidates      : %d\n", length(x$weights)))
  cat(sprintf("  effective support: %.1f  (exp of the entropy)\n", eff))
  cat(sprintf("  lambda          : %g\n", x$lambda))
  cat(sprintf("  objective       : %.3f%s\n", x$objective,
              if (isTRUE(x$converged)) "" else "  [NOT CONVERGED]"))
  invisible(x)
}

#' Leave-one-out banks
#'
#' Refits the bank once per scored participant, each time excluding that
#' participant's own measurements. This restores the honesty condition when no
#' external donor sample exists, at the cost of one fit per participant.
#'
#' @param loglik A participant-by-candidate log-likelihood matrix.
#' @param lambda Entropy penalty passed to [fit_bank()].
#' @param nu Optional per-participant candidate retention count. When supplied,
#'   [screen_candidates()] is rerun inside every leave-one-out fold using only
#'   that fold's donors. This is required when support screening is
#'   data-dependent; screening once on the full cohort before calling
#'   `bank_loo()` would leak the held-out participant into its bank.
#' @param ... Further arguments to [fit_bank()].
#'
#' @return A list with `weights`, a participant-by-candidate matrix whose
#'   row `i` holds the bank fitted without participant `i`, and the
#'   corresponding `logLik`, a vector of each participant's log class
#'   likelihood under their own leave-one-out bank, and `support_size`, the
#'   number of candidates retained in each fold.
#'
#' @details
#' The columns of `loglik` must come from a candidate family fixed independently
#' of the participant being scored. If the family is reduced using participant
#' measurements, pass `nu` here so screening is repeated within each fold.
#' Positivity or other screens that depend only on the mechanistic design, and
#' not on participant measurements, may be applied before this function.
#' @seealso [fit_bank()]
#' @export
bank_loo <- function(loglik, lambda = 1, nu = NULL, ...) {
  n <- nrow(loglik)
  if (n < 3L) stop("Leave-one-out needs at least three participants.")
  if (!is.null(nu)) stopifnot(length(nu) == 1L, is.finite(nu), nu >= 1)
  W <- matrix(NA_real_, n, ncol(loglik),
              dimnames = list(rownames(loglik), NULL))
  ll <- numeric(n)
  support_size <- integer(n)
  for (i in seq_len(n)) {
    train <- loglik[-i, , drop = FALSE]
    keep <- if (is.null(nu)) seq_len(ncol(loglik)) else
      screen_candidates(train, nu = nu)$keep
    wi <- numeric(ncol(loglik))
    if (length(keep) == 1L) {
      wi[keep] <- 1
    } else {
      b <- fit_bank(train[, keep, drop = FALSE], lambda = lambda, ...)
      wi[keep] <- b$weights
    }
    W[i, ] <- wi
    ll[i]  <- class_loglik(loglik[i, , drop = FALSE], wi)
    support_size[i] <- length(keep)
  }
  names(ll) <- rownames(loglik)
  names(support_size) <- rownames(loglik)
  list(weights = W, logLik = ll, support_size = support_size, nu = nu)
}
