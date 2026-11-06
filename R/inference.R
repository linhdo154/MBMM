#' Test a covariate effect in the gate
#'
#' Compares a full gate against a reduced gate in which the effect is removed.
#' With `firth = FALSE` the statistic is the usual likelihood ratio
#' \deqn{T = 2\{\hat\ell(\hat\gamma_{\mathrm{full}})
#'            - \hat\ell(\hat\gamma_{\mathrm{red}})\}.}
#'
#' With `firth = TRUE` it is instead twice the difference of **penalised
#' objectives**,
#' \deqn{T_{\mathrm{PO}} = 2\{\tilde\ell_X(\hat\gamma_{\mathrm{full}})
#'            - \tilde\ell_{X_0}(\hat\gamma_{\mathrm{red}})\},
#'       \qquad \tilde\ell_X(\gamma) = \ell(\gamma)
#'            + \tfrac12 \log |X^\top W(\gamma) X| .}
#' Because the two penalty determinants have different dimensions,
#' \eqn{T_{\mathrm{PO}}} is a penalised-objective difference and **not** a
#' likelihood ratio. No chi-squared reference is reported for it; use the
#' permutation null, which is valid for any fixed statistic.
#'
#' @param l1,l0 Log class likelihoods, one per participant.
#' @param X Full gate design matrix, including an intercept.
#' @param X0 Reduced design matrix. Defaults to the intercept alone.
#' @param firth Use the Firth-penalised gate and \eqn{T_{\mathrm{PO}}}.
#' @param B Number of permutations. `0` skips the permutation null. The class
#'   likelihoods are held fixed while the covariate rows are shuffled, which
#'   is valid because they do not use the covariate.
#' @param seed Optional seed for the permutation draw.
#'
#' @return An object of class `mbmm_test` with the observed `statistic`, its
#'   `df`, the permutation `p_perm` and null distribution when `B > 0`, the
#'   chi-squared `p_chisq` when `firth = FALSE`, and both fitted gates.
#' @seealso [fit_gate()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 80; x <- rnorm(n); X <- cbind(1, x)
#' l1 <- rnorm(n, 0.6 * x); l0 <- rnorm(n)
#' gate_test(l1, l0, X, B = 99, seed = 1)
gate_test <- function(l1, l0, X, X0 = NULL, firth = FALSE, B = 0L,
                      seed = NULL) {
  if (is.null(X0)) X0 <- matrix(1, nrow(X), 1)
  stat_of <- function(Xf, X0f) {
    ff <- fit_gate(l1, l0, Xf, firth = firth)
    nn <- fit_gate(l1, l0, X0f, firth = firth)
    if (firth) 2 * (ff$penalised - nn$penalised) else 2 * (ff$loglik - nn$loglik)
  }
  full <- fit_gate(l1, l0, X,  firth = firth)
  red  <- fit_gate(l1, l0, X0, firth = firth)
  stat <- if (firth) 2 * (full$penalised - red$penalised) else
                     2 * (full$loglik - red$loglik)
  df <- ncol(X) - ncol(X0)

  p_chisq <- if (!firth) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_
  null <- NULL; p_perm <- NA_real_
  if (B > 0L) {
    if (!is.null(seed)) set.seed(seed)
    keep <- seq_len(ncol(X))[!apply(X, 2, function(c) all(c == c[1]))]
    null <- vapply(seq_len(B), function(b) {
      idx <- sample.int(nrow(X))
      Xp <- X; Xp[, keep] <- X[idx, keep, drop = FALSE]
      stat_of(Xp, X0)
    }, numeric(1))
    p_perm <- (1 + sum(null >= stat)) / (B + 1)
  }
  structure(list(statistic = stat, df = df, firth = firth,
                 p_chisq = p_chisq, p_perm = p_perm, B = B, null = null,
                 full = full, reduced = red), class = "mbmm_test")
}

#' @export
print.mbmm_test <- function(x, ...) {
  nm <- if (x$firth) "T_PO (penalised-objective difference)" else
                     "T (likelihood ratio)"
  cat("<mbmm_test>\n")
  cat(sprintf("  %s = %.3f on %d df\n", nm, x$statistic, x$df))
  if (x$firth)
    cat("  not a likelihood ratio: the penalty determinants differ in dimension,\n",
        "  so no chi-squared reference is reported\n", sep = "")
  else
    cat(sprintf("  chi-squared p = %.4g\n", x$p_chisq))
  if (x$B > 0)
    cat(sprintf("  permutation p = %.5f  (B = %d, smallest attainable %.5f)\n",
                x$p_perm, x$B, 1 / (x$B + 1)))
  invisible(x)
}
