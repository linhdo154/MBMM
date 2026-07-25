#' Numerically stable log-sum-exp
#'
#' Computes \code{log(sum(exp(v)))} without numerical overflow or underflow.
#' Negative infinity values are retained as valid zero-likelihood
#' contributions. Missing values are not silently removed.
#'
#' @param v Numeric vector of log-scale values.
#' @return A scalar equal to \code{log(sum(exp(v)))}.
#' @export
mbmm_logsumexp <- function(v) {
     v <- as.numeric(v)
     
     if (length(v) == 0L) {
          return(-Inf)
     }
     
     if (anyNA(v) || any(is.nan(v))) {
          stop(
               "`v` contains NA or NaN values in `mbmm_logsumexp()`.",
               call. = FALSE
          )
     }
     
     # Any +Inf term makes the sum infinite.
     if (any(v == Inf)) {
          return(Inf)
     }
     
     # All terms equal -Inf means that every likelihood contribution is zero.
     if (all(v == -Inf)) {
          return(-Inf)
     }
     
     m <- max(v)
     
     m + log(sum(exp(v - m)))
}


#' Numerically stable log-mean-exp
#'
#' Used to Monte Carlo marginalize a class-conditional likelihood over a
#' reference bank of kinetic-parameter particles or bundles:
#'
#' \deqn{
#' \log\left\{
#' \frac{1}{M}
#' \sum_{m=1}^{M}\exp(v_m)
#' \right\}.
#' }
#'
#' Importantly, \code{-Inf} values remain part of the denominator \eqn{M};
#' they represent valid particles with zero likelihood and must not be
#' discarded.
#'
#' @param v Numeric vector of per-particle log-likelihoods.
#' @return Scalar log-mean-exp.
#' @export
mbmm_logmeanexp <- function(v) {
     v <- as.numeric(v)
     
     if (length(v) == 0L) {
          return(-Inf)
     }
     
     if (anyNA(v) || any(is.nan(v))) {
          stop(
               "`v` contains NA or NaN values in `mbmm_logmeanexp()`.",
               call. = FALSE
          )
     }
     
     mbmm_logsumexp(v) - log(length(v))
}


#' Inverse logit
#'
#' Uses \code{\link[stats]{plogis}} for numerical stability.
#'
#' @param x Numeric vector.
#' @return Numeric vector in the interval [0, 1].
#' @export
mbmm_inv_logit <- function(x) {
     stats::plogis(x)
}


#' Logit
#'
#' @param p Numeric vector of probabilities in [0, 1].
#' @return Numeric vector on the log-odds scale. Values 0 and 1 map to
#'   \code{-Inf} and \code{Inf}, respectively.
#' @export
mbmm_logit <- function(p) {
     if (!is.numeric(p)) {
          stop("`p` must be numeric.", call. = FALSE)
     }
     
     if (anyNA(p)) {
          stop("`p` contains missing values.", call. = FALSE)
     }
     
     if (any(p < 0 | p > 1)) {
          stop("All values of `p` must lie between 0 and 1.", call. = FALSE)
     }
     
     stats::qlogis(p)
}


#' Safe log10 with a positive floor
#'
#' Values below \code{floor}, including zero and negative values, are replaced
#' by \code{floor} before transformation. Missing values remain missing.
#'
#' @param x Numeric vector.
#' @param floor Smallest positive value passed to \code{log10}.
#' @return Numeric vector on the log10 scale.
#' @export
mbmm_safe_log10 <- function(x, floor = 0.01) {
     if (!is.numeric(floor) ||
         length(floor) != 1L ||
         !is.finite(floor) ||
         floor <= 0) {
          stop("`floor` must be one finite, strictly positive number.",
               call. = FALSE)
     }
     
     x <- as.numeric(x)
     
     log10(pmax(x, floor))
}


#' Gaussian log-density with optional left censoring at an LOD
#'
#' Detected observations contribute a Gaussian log-density on the measurement
#' scale. Left-censored observations contribute the Gaussian probability of
#' falling at or below the LOD.
#'
#' All arguments must already be on the same likelihood scale, normally the
#' log10 scale used by MBMM.
#'
#' @param y Observed value on the same scale as \code{mean}. Values of
#'   \code{y} are ignored where \code{censored} is TRUE.
#' @param mean Predicted mean on the likelihood scale.
#' @param sd Strictly positive measurement-error standard deviation.
#' @param censored Logical or 0/1 vector indicating left-censored observations.
#' @param lod LOD on the same scale as \code{y} and \code{mean}. Required for
#'   every censored observation.
#' @return Numeric vector of pointwise log-likelihood contributions.
#' @export
mbmm_dnorm_censored <- function(
          y,
          mean,
          sd,
          censored = FALSE,
          lod = NULL
) {
     n <- max(
          length(y),
          length(mean),
          length(sd),
          length(censored),
          if (is.null(lod)) 0L else length(lod)
     )
     
     if (n == 0L) {
          return(numeric(0))
     }
     
     y <- rep_len(as.numeric(y), n)
     mean <- rep_len(as.numeric(mean), n)
     sd <- rep_len(as.numeric(sd), n)
     censored <- rep_len(as.logical(censored), n)
     
     if (anyNA(censored)) {
          stop("`censored` contains missing values.", call. = FALSE)
     }
     
     if (any(!is.finite(mean))) {
          stop("`mean` must contain only finite values.", call. = FALSE)
     }
     
     if (any(!is.finite(sd)) || any(sd <= 0)) {
          stop("`sd` must contain finite, strictly positive values.",
               call. = FALSE)
     }
     
     detected <- !censored
     
     if (any(detected) && any(!is.finite(y[detected]))) {
          stop(
               "Detected observations require finite `y` values.",
               call. = FALSE
          )
     }
     
     if (any(censored)) {
          if (is.null(lod)) {
               stop(
                    "`lod` is required when any observation is censored.",
                    call. = FALSE
               )
          }
          
          lod <- rep_len(as.numeric(lod), n)
          
          if (any(!is.finite(lod[censored]))) {
               stop(
                    "Censored observations require finite `lod` values.",
                    call. = FALSE
               )
          }
     } else if (is.null(lod)) {
          lod <- rep(NA_real_, n)
     } else {
          lod <- rep_len(as.numeric(lod), n)
     }
     
     out <- numeric(n)
     
     if (any(detected)) {
          out[detected] <- stats::dnorm(
               x = y[detected],
               mean = mean[detected],
               sd = sd[detected],
               log = TRUE
          )
     }
     
     if (any(censored)) {
          out[censored] <- stats::pnorm(
               q = lod[censored],
               mean = mean[censored],
               sd = sd[censored],
               log.p = TRUE
          )
     }
     
     out
}
#' Fixed-step 4th-order Runge-Kutta ODE integrator (base R, no dependencies)
#'
#' A small, dependency-free ODE solver intended for building the
#' \code{simulate_fn} ingredient MBMM needs (see \code{\link{mbmm_abc_fit}}).
#' For stiff systems or very fine control, use \pkg{deSolve} instead -- this
#' function exists so the package (and its examples) do not require any
#' package beyond base R.
#'
#' @param f function(t, y, parms) returning the derivative dy/dt as a
#'   numeric vector the same length as y
#' @param y0 numeric vector of initial states
#' @param times numeric vector of output times (need not be evenly spaced);
#'   the first element is treated as t0 with state y0
#' @param parms a list or named vector passed through to \code{f}
#' @param n_substeps number of RK4 substeps to take between each pair of
#'   consecutive \code{times} (higher = more accurate, slower)
#' @return a matrix with one row per element of \code{times} and one column
#'   per state variable (plus a leading \code{time} column)
#' @export
mbmm_rk4_solve <- function(f, y0, times, parms, n_substeps = 10) {
  times <- sort(unique(times))
  nT <- length(times)
  nY <- length(y0)
  out <- matrix(NA_real_, nrow = nT, ncol = nY)
  colnames(out) <- names(y0)
  y <- y0
  out[1, ] <- y
  if (nT > 1) {
    for (i in 2:nT) {
      t0 <- times[i - 1]; t1 <- times[i]
      h <- (t1 - t0) / n_substeps
      t <- t0
      if (h > 0) {
        for (s in seq_len(n_substeps)) {
          k1 <- f(t, y, parms)
          k2 <- f(t + h / 2, y + h / 2 * k1, parms)
          k3 <- f(t + h / 2, y + h / 2 * k2, parms)
          k4 <- f(t + h, y + h * k3, parms)
          y <- y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
          t <- t + h
        }
      }
      out[i, ] <- y
    }
  }
  cbind(time = times, out)
}
