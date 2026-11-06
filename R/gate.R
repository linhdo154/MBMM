#' Fit the covariate gate by EM
#'
#' Combines the prior state probability \eqn{\pi_i(\gamma) = \mathrm{expit}
#' (x_i^\top \gamma)} with the two class likelihoods to form the posterior
#' responsibility
#' \deqn{w_i(\gamma) = \frac{\pi_i L_{i1}}{\pi_i L_{i1} + (1-\pi_i) L_{i0}},}
#' then refits \eqn{\gamma} by logistic regression using \eqn{w_i} as a
#' fractional outcome, iterating to convergence. A participant whose
#' trajectory does not discriminate the states keeps a responsibility close to
#' the gate probability and contributes at that weight rather than being
#' forced into a class.
#'
#' @param l1,l0 Log class likelihoods for state 1 and state 0, one per
#'   participant. See [class_loglik()].
#' @param X Gate design matrix, one row per participant, including an
#'   intercept column.
#' @param firth Use the Jeffreys-penalised score of Firth (1993). Recommended
#'   when any fitted probability approaches zero or one, which is routine with
#'   few participants per covariate level.
#' @param maxit,tol Iteration controls.
#'
#' @return An object of class `mbmm_gate`: `gamma`, the observed-data
#'   `loglik`, the responsibilities `w`, fitted probabilities `pi`, and, when
#'   `firth = TRUE`, the `penalised` objective
#'   \eqn{\ell + \tfrac12 \log |X^\top W X|}.
#' @references Firth D (1993). Bias reduction of maximum likelihood estimates.
#'   *Biometrika* 80, 27-38.
#' @seealso [gate_test()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 60; x <- rnorm(n); X <- cbind(1, x)
#' l1 <- rnorm(n, 1); l0 <- rnorm(n)
#' fit_gate(l1, l0, X)$gamma
fit_gate <- function(l1, l0, X, firth = FALSE, maxit = 2000L, tol = 1e-11) {
  stopifnot(length(l1) == length(l0), length(l1) == nrow(X))
  if (firth) return(gate_firth(l1, l0, X, maxit = maxit, tol = tol))
  g <- rep(0, ncol(X)); prev <- -Inf; n <- length(l1)
  for (it in seq_len(maxit)) {
    pp <- stats::plogis(as.vector(X %*% g))
    a  <- log(pmax(pp, 1e-300)) + l1
    b  <- log(pmax(1 - pp, 1e-300)) + l0
    m  <- pmax(a, b)
    ll <- sum(m + log(exp(a - m) + exp(b - m)))
    w  <- exp(a - m) / (exp(a - m) + exp(b - m))
    fitm <- suppressWarnings(
      stats::glm.fit(rbind(X, X), c(rep(1, n), rep(0, n)),
                     weights = c(w, 1 - w), family = stats::binomial()))
    g <- ifelse(is.na(fitm$coefficients), 0, fitm$coefficients)
    if (is.finite(prev) && abs(ll - prev) < tol) break
    prev <- ll
  }
  structure(list(gamma = g, loglik = ll, w = w,
                 pi = stats::plogis(as.vector(X %*% g)),
                 firth = FALSE, iter = it), class = "mbmm_gate")
}

gate_firth <- function(l1, l0, X, maxit = 500L, tol = 1e-10) {
  n <- length(l1); g <- rep(0, ncol(X)); prev <- -Inf
  w <- rep(0.5, n); obj <- -Inf
  pen_of <- function(gam) {
    pi <- stats::plogis(as.vector(X %*% gam))
    Wd <- pmax(pi * (1 - pi), 1e-12)
    0.5 * as.numeric(determinant(crossprod(X * sqrt(Wd)),
                                 logarithm = TRUE)$modulus)
  }
  for (it in seq_len(maxit)) {
    pp <- stats::plogis(as.vector(X %*% g))
    a  <- log(pmax(pp, 1e-300)) + l1
    b  <- log(pmax(1 - pp, 1e-300)) + l0
    m  <- pmax(a, b)
    ll <- sum(m + log(exp(a - m) + exp(b - m)))
    w  <- exp(a - m) / (exp(a - m) + exp(b - m))
    for (inner in seq_len(100L)) {
      pi  <- stats::plogis(as.vector(X %*% g))
      Wd  <- pmax(pi * (1 - pi), 1e-12)
      XtWX <- crossprod(X * sqrt(Wd))
      Inv <- tryCatch(solve(XtWX),
                      error = function(e) solve(XtWX + diag(1e-8, ncol(X))))
      Hd  <- rowSums((X %*% Inv) * X) * Wd
      step <- as.vector(Inv %*% crossprod(X, w - pi + Hd * (0.5 - pi)))
      step <- step * min(1, 5 / max(abs(step)))
      g <- g + step
      if (max(abs(step)) < 1e-11) break
    }
    obj <- ll + pen_of(g)
    if (is.finite(prev) && abs(obj - prev) < tol) break
    prev <- obj
  }
  structure(list(gamma = g, loglik = ll, penalised = obj, w = w,
                 pi = stats::plogis(as.vector(X %*% g)),
                 firth = TRUE, iter = it), class = "mbmm_gate")
}

#' @export
print.mbmm_gate <- function(x, ...) {
  cat("<mbmm_gate>", if (x$firth) " Firth-penalised" else "", "\n", sep = "")
  cat("  gamma  :", paste(sprintf("%.3f", x$gamma), collapse = "  "), "\n")
  cat(sprintf("  loglik : %.3f\n", x$loglik))
  if (!is.null(x$penalised))
    cat(sprintf("  penalised objective: %.3f\n", x$penalised))
  unc <- sum(x$w > 0.1 & x$w < 0.9)
  cat(sprintf("  responsibilities in (0.1, 0.9): %d of %d\n", unc, length(x$w)))
  invisible(x)
}
