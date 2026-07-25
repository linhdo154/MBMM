#' Internal: observed-data marginal log-likelihood for a binary MBMM gate
#'
#' Computes
#'
#' \deqn{
#' \ell(\beta)
#' =
#' \sum_i
#' \log\left[
#' \pi_i(\beta)L_{i1}
#' +
#' \{1-\pi_i(\beta)\}L_{i0}
#' \right].
#' }
#'
#' This function returns the unpenalized observed-data likelihood.
#'
#' @keywords internal
.mbmm_obs_ll_binary <- function(
          beta,
          X,
          logL_pos,
          logL_neg
) {
     eta <- as.numeric(X %*% beta)
     
     # Numerically stable log gate probabilities.
     log_pi_pos <- stats::plogis(
          eta,
          log.p = TRUE
     )
     
     log_pi_neg <- stats::plogis(
          eta,
          lower.tail = FALSE,
          log.p = TRUE
     )
     
     terms <- vapply(
          seq_along(logL_pos),
          function(i) {
               mbmm_logsumexp(
                    c(
                         logL_pos[i] + log_pi_pos[i],
                         logL_neg[i] + log_pi_neg[i]
                    )
               )
          },
          numeric(1)
     )
     
     if (anyNA(terms)) {
          stop(
               "Binary gate likelihood produced NA values.",
               call. = FALSE
          )
     }
     
     sum(terms)
}


#' Internal: ridge penalty for a binary gate
#'
#' Applies an L2 penalty to non-intercept coefficients only.
#'
#' @keywords internal
.mbmm_ridge_penalty_binary <- function(
          beta,
          ridge
) {
     if (!is.numeric(ridge) ||
         length(ridge) != 1L ||
         !is.finite(ridge) ||
         ridge < 0) {
          stop(
               "`ridge` must be one finite, nonnegative number.",
               call. = FALSE
          )
     }
     
     if (ridge == 0 || length(beta) <= 1L) {
          return(0)
     }
     
     0.5 * ridge * sum(beta[-1L]^2)
}


#' Internal: observed-data marginal log-likelihood for a multinomial gate
#'
#' The final latent class is the reference class, with its coefficient vector
#' fixed at zero.
#'
#' @param beta_vec Numeric vector containing coefficients for the first
#'   \code{K - 1} classes.
#' @param X Gate design matrix.
#' @param logL_mat Subject-by-class matrix of class log-likelihoods.
#' @param K Number of latent classes.
#'
#' @keywords internal
.mbmm_obs_ll_multinom <- function(
          beta_vec,
          X,
          logL_mat,
          K
) {
     p <- ncol(X)
     
     beta_mat <- matrix(
          beta_vec,
          nrow = K - 1L,
          ncol = p,
          byrow = TRUE
     )
     
     eta <- X %*% t(beta_mat)
     
     # Final class is the multinomial reference class.
     eta_full <- cbind(eta, 0)
     
     row_max <- apply(
          eta_full,
          1L,
          max
     )
     
     eta_shift <- eta_full - row_max
     
     log_denom <- log(
          rowSums(exp(eta_shift))
     ) + row_max
     
     log_pi <- eta_full - log_denom
     
     terms <- vapply(
          seq_len(nrow(X)),
          function(i) {
               mbmm_logsumexp(
                    logL_mat[i, ] + log_pi[i, ]
               )
          },
          numeric(1)
     )
     
     if (anyNA(terms)) {
          stop(
               "Multinomial gate likelihood produced NA values.",
               call. = FALSE
          )
     }
     
     sum(terms)
}


#' Internal: ridge penalty for a multinomial gate
#'
#' Applies an L2 penalty to all non-intercept coefficients for the first
#' \code{K - 1} classes.
#'
#' @keywords internal
.mbmm_ridge_penalty_multinom <- function(
          beta_vec,
          p,
          K,
          ridge
) {
     if (!is.numeric(ridge) ||
         length(ridge) != 1L ||
         !is.finite(ridge) ||
         ridge < 0) {
          stop(
               "`ridge` must be one finite, nonnegative number.",
               call. = FALSE
          )
     }
     
     if (ridge == 0 || p <= 1L) {
          return(0)
     }
     
     beta_mat <- matrix(
          beta_vec,
          nrow = K - 1L,
          ncol = p,
          byrow = TRUE
     )
     
     0.5 * ridge * sum(
          beta_mat[, -1L, drop = FALSE]^2
     )
}


#' Internal: fit a binary MBMM gate by direct marginal optimization
#'
#' @keywords internal
.mbmm_fit_binary_direct <- function(
          X,
          logL_mat,
          penalty,
          ridge,
          start = NULL,
          maxit = 2000
) {
     p <- ncol(X)
     
     if (is.null(start)) {
          start <- rep(0, p)
     }
     
     if (!is.numeric(start) || length(start) != p) {
          stop(
               "Binary gate starting vector has the wrong length.",
               call. = FALSE
          )
     }
     
     ridge_use <- if (identical(penalty, "ridge")) {
          ridge
     } else {
          0
     }
     
     unpenalized_ll <- function(beta) {
          .mbmm_obs_ll_binary(
               beta = beta,
               X = X,
               logL_pos = logL_mat[, 1L],
               logL_neg = logL_mat[, 2L]
          )
     }
     
     penalized_objective <- function(beta) {
          unpenalized_ll(beta) -
               .mbmm_ridge_penalty_binary(
                    beta = beta,
                    ridge = ridge_use
               )
     }
     
     opt <- stats::optim(
          par = start,
          fn = function(beta) {
               -penalized_objective(beta)
          },
          method = "BFGS",
          control = list(
               maxit = maxit,
               reltol = 1e-10
          ),
          hessian = TRUE
     )
     
     beta_hat <- opt$par
     
     list(
          beta = beta_hat,
          X = X,
          logLik_unpenalized = unpenalized_ll(beta_hat),
          objective_penalized = penalized_objective(beta_hat),
          convergence = opt$convergence,
          message = opt$message,
          counts = opt$counts,
          hessian = opt$hessian
     )
}


#' Internal: fit a multinomial MBMM gate by direct marginal optimization
#'
#' @keywords internal
.mbmm_fit_multinom_direct <- function(
          X,
          logL_mat,
          K,
          penalty,
          ridge,
          start = NULL,
          maxit = 4000
) {
     p <- ncol(X)
     npar <- (K - 1L) * p
     
     if (is.null(start)) {
          start <- rep(0, npar)
     }
     
     if (!is.numeric(start) || length(start) != npar) {
          stop(
               "Multinomial gate starting vector has the wrong length.",
               call. = FALSE
          )
     }
     
     ridge_use <- if (identical(penalty, "ridge")) {
          ridge
     } else {
          0
     }
     
     unpenalized_ll <- function(beta_vec) {
          .mbmm_obs_ll_multinom(
               beta_vec = beta_vec,
               X = X,
               logL_mat = logL_mat,
               K = K
          )
     }
     
     penalized_objective <- function(beta_vec) {
          unpenalized_ll(beta_vec) -
               .mbmm_ridge_penalty_multinom(
                    beta_vec = beta_vec,
                    p = p,
                    K = K,
                    ridge = ridge_use
               )
     }
     
     opt <- stats::optim(
          par = start,
          fn = function(beta_vec) {
               -penalized_objective(beta_vec)
          },
          method = "BFGS",
          control = list(
               maxit = maxit,
               reltol = 1e-10
          ),
          hessian = TRUE
     )
     
     beta_hat <- matrix(
          opt$par,
          nrow = K - 1L,
          ncol = p,
          byrow = TRUE
     )
     
     list(
          beta = beta_hat,
          beta_vec = opt$par,
          X = X,
          logLik_unpenalized = unpenalized_ll(opt$par),
          objective_penalized = penalized_objective(opt$par),
          convergence = opt$convergence,
          message = opt$message,
          counts = opt$counts,
          hessian = opt$hessian
     )
}


#' Fit the MBMM gate model
#'
#' Fits the covariate-linked gate model using fixed class-conditional
#' mechanistic log-likelihoods.
#'
#' For a binary latent state, the observed-data marginal likelihood is
#'
#' \deqn{
#' \ell(\beta)
#' =
#' \sum_i
#' \log\left[
#' \pi_i(\beta)L_{i1}
#' +
#' \{1-\pi_i(\beta)\}L_{i0}
#' \right].
#' }
#'
#' The parameters are estimated by direct optimization of this marginal
#' likelihood. An optional ridge penalty can be applied to non-intercept
#' coefficients for numerical stabilization. Ridge is not part of the core
#' MBMM model.
#'
#' The reported likelihood-ratio statistic is evaluated using the
#' unpenalized observed-data marginal likelihood at the fitted full- and
#' null-model coefficient estimates:
#'
#' \deqn{
#' 2\{\ell(\widehat{\beta}_{full})
#' -
#' \ell(\widehat{\beta}_{null})\}.
#' }
#'
#' @param logL Data.frame or matrix of class-conditional log-likelihoods.
#'   Columns are normally named \code{logL_<class>}.
#' @param formula One-sided formula defining the gate covariates, such as
#'   \code{~ baseline_antibody} or \code{~ arm}.
#' @param data Data.frame containing gate covariates.
#' @param classes Character vector defining the class order. For two classes,
#'   \code{classes[1]} is the positive/event class and \code{classes[2]} is
#'   the reference class. For more than two classes, the final class is the
#'   multinomial reference class.
#' @param ridge Nonnegative L2 penalty applied to non-intercept coefficients
#'   when \code{penalty = "ridge"}.
#' @param penalty Either \code{"ridge"} or \code{"none"}.
#' @param n_perm Number of permutations for an optional permutation p-value.
#'   The implemented procedure permutes the complete full-model design matrix
#'   and is appropriate only when the tested non-intercept covariates are
#'   exchangeable as a block.
#' @param seed Optional random seed for permutation testing.
#' @param id_col Optional subject-ID column in \code{data}, copied to the
#'   posterior table.
#' @param maxit Maximum number of optimizer iterations.
#'
#' @return An object of class \code{"mbmm_gate"}.
#'
#' @examples
#' set.seed(1)
#'
#' n <- 150
#' x <- rnorm(n)
#'
#' z <- rbinom(
#'   n,
#'   size = 1,
#'   prob = mbmm_inv_logit(-0.5 + 1.5 * x)
#' )
#'
#' logL <- data.frame(
#'   logL_pos = ifelse(
#'     z == 1,
#'     log(0.98),
#'     log(0.02)
#'   ),
#'   logL_neg = ifelse(
#'     z == 1,
#'     log(0.02),
#'     log(0.98)
#'   )
#' )
#'
#' fit <- mbmm_fit_gate(
#'   logL = logL,
#'   formula = ~x,
#'   data = data.frame(x = x),
#'   classes = c("pos", "neg"),
#'   penalty = "ridge",
#'   ridge = 1e-4
#' )
#'
#' fit
#'
#' @export
mbmm_fit_gate <- function(
          logL,
          formula,
          data,
          classes = NULL,
          ridge = 1e-4,
          penalty = "ridge",
          n_perm = 0,
          seed = NULL,
          id_col = NULL,
          maxit = 4000
) {
     allowed_penalties <- c("ridge", "none")
     
     if (
          !is.character(penalty) ||
          length(penalty) != 1L ||
          is.na(penalty) ||
          !penalty %in% allowed_penalties
     ) {
          stop(
               "`penalty` must be one of \"ridge\" or \"none\".",
               call. = FALSE
          )
     }
     
     if (!inherits(formula, "formula")) {
          stop(
               "`formula` must be a formula.",
               call. = FALSE
          )
     }
     
     if (!is.data.frame(data)) {
          stop(
               "`data` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(ridge) ||
         length(ridge) != 1L ||
         !is.finite(ridge) ||
         ridge < 0) {
          stop(
               "`ridge` must be one finite, nonnegative number.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(n_perm) ||
         length(n_perm) != 1L ||
         !is.finite(n_perm) ||
         n_perm < 0 ||
         n_perm != as.integer(n_perm)) {
          stop(
               "`n_perm` must be a nonnegative integer.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(maxit) ||
         length(maxit) != 1L ||
         !is.finite(maxit) ||
         maxit < 1 ||
         maxit != as.integer(maxit)) {
          stop(
               "`maxit` must be a positive integer.",
               call. = FALSE
          )
     }
     
     logL <- as.data.frame(logL)
     
     if (is.null(classes)) {
          logL_columns <- grep(
               "^logL_",
               names(logL),
               value = TRUE
          )
          
          classes <- sub(
               "^logL_",
               "",
               logL_columns
          )
     }
     
     if (length(classes) < 2L) {
          stop(
               "At least two latent classes are required.",
               call. = FALSE
          )
     }
     
     if (any(classes == "") || anyDuplicated(classes)) {
          stop(
               "`classes` must contain unique, non-empty labels.",
               call. = FALSE
          )
     }
     
     logL_names <- paste0(
          "logL_",
          classes
     )
     
     missing_logL <- setdiff(
          logL_names,
          names(logL)
     )
     
     if (length(missing_logL) > 0L) {
          stop(
               "`logL` is missing class column(s): ",
               paste(missing_logL, collapse = ", "),
               call. = FALSE
          )
     }
     
     logL_mat <- as.matrix(
          logL[, logL_names, drop = FALSE]
     )
     
     storage.mode(logL_mat) <- "double"
     
     if (nrow(logL_mat) != nrow(data)) {
          stop(
               "`logL` and `data` must contain the same number of rows.",
               call. = FALSE
          )
     }
     
     K <- length(classes)
     
     # ---------------------------------------------------------------------------
     # Construct the full gate design matrix while preserving row alignment.
     # ---------------------------------------------------------------------------
     
     model_frame <- stats::model.frame(
          formula = formula,
          data = data,
          na.action = stats::na.pass
     )
     
     X_full_all <- stats::model.matrix(
          object = formula,
          data = model_frame
     )
     
     finite_logL <- apply(
          logL_mat,
          1L,
          function(x) {
               all(is.finite(x))
          }
     )
     
     complete_X <- apply(
          X_full_all,
          1L,
          function(x) {
               all(is.finite(x))
          }
     )
     
     keep <- finite_logL & complete_X
     
     if (!all(keep)) {
          message(
               sprintf(
                    paste0(
                         "Dropping %d subject(s): ",
                         "%d with non-finite class likelihoods and ",
                         "%d with incomplete gate covariates."
                    ),
                    sum(!keep),
                    sum(!finite_logL),
                    sum(!complete_X)
               )
          )
     }
     
     if (!any(keep)) {
          stop(
               "No subjects remain after likelihood and covariate validation.",
               call. = FALSE
          )
     }
     
     logL_mat <- logL_mat[
          keep,
          ,
          drop = FALSE
     ]
     
     data_used <- data[
          keep,
          ,
          drop = FALSE
     ]
     
     X_full <- X_full_all[
          keep,
          ,
          drop = FALSE
     ]
     
     X_null <- matrix(
          1,
          nrow = nrow(X_full),
          ncol = 1L,
          dimnames = list(
               NULL,
               "(Intercept)"
          )
     )
     
     if (nrow(X_full) <= ncol(X_full)) {
          warning(
               paste0(
                    "The full gate has at least as many coefficients as ",
                    "retained subjects."
               ),
               call. = FALSE
          )
     }
     
     # ---------------------------------------------------------------------------
     # Fit full and null models.
     # ---------------------------------------------------------------------------
     
     if (K == 2L) {
          full <- .mbmm_fit_binary_direct(
               X = X_full,
               logL_mat = logL_mat,
               penalty = penalty,
               ridge = ridge,
               maxit = maxit
          )
          
          null <- .mbmm_fit_binary_direct(
               X = X_null,
               logL_mat = logL_mat,
               penalty = penalty,
               ridge = ridge,
               maxit = maxit
          )
          
          df_lrt <- ncol(X_full) -
               ncol(X_null)
     } else {
          full <- .mbmm_fit_multinom_direct(
               X = X_full,
               logL_mat = logL_mat,
               K = K,
               penalty = penalty,
               ridge = ridge,
               maxit = maxit
          )
          
          null <- .mbmm_fit_multinom_direct(
               X = X_null,
               logL_mat = logL_mat,
               K = K,
               penalty = penalty,
               ridge = ridge,
               maxit = maxit
          )
          
          df_lrt <- (K - 1L) * (
               ncol(X_full) -
                    ncol(X_null)
          )
     }
     
     if (full$convergence != 0L) {
          warning(
               paste0(
                    "The full gate optimizer did not report convergence. Code: ",
                    full$convergence,
                    if (!is.null(full$message)) {
                         paste0("; ", full$message)
                    } else {
                         ""
                    }
               ),
               call. = FALSE
          )
     }
     
     if (null$convergence != 0L) {
          warning(
               paste0(
                    "The null gate optimizer did not report convergence. Code: ",
                    null$convergence,
                    if (!is.null(null$message)) {
                         paste0("; ", null$message)
                    } else {
                         ""
                    }
               ),
               call. = FALSE
          )
     }
     
     # ---------------------------------------------------------------------------
     # Likelihood-ratio statistic based on the unpenalized MBMM likelihood.
     # ---------------------------------------------------------------------------
     
     LRT <- 2 * (
          full$logLik_unpenalized -
               null$logLik_unpenalized
     )
     
     if (is.finite(LRT) &&
         LRT < 0 &&
         abs(LRT) < 1e-8) {
          LRT <- 0
     }
     
     if (!is.finite(LRT) || LRT < 0) {
          warning(
               paste0(
                    "The likelihood-ratio statistic is non-finite or negative. ",
                    "Inspect gate convergence and model specification."
               ),
               call. = FALSE
          )
     }
     
     p_value <- if (
          df_lrt > 0 &&
          is.finite(LRT) &&
          LRT >= 0
     ) {
          stats::pchisq(
               LRT,
               df = df_lrt,
               lower.tail = FALSE
          )
     } else {
          NA_real_
     }
     
     # ---------------------------------------------------------------------------
     # Posterior latent-class probabilities under the full model.
     # ---------------------------------------------------------------------------
     
     if (K == 2L) {
          eta <- as.numeric(
               X_full %*% full$beta
          )
          
          log_pi_pos <- stats::plogis(
               eta,
               log.p = TRUE
          )
          
          log_pi_neg <- stats::plogis(
               eta,
               lower.tail = FALSE,
               log.p = TRUE
          )
          
          log_post_pos <- logL_mat[, 1L] +
               log_pi_pos
          
          log_post_neg <- logL_mat[, 2L] +
               log_pi_neg
          
          log_denom <- vapply(
               seq_len(nrow(logL_mat)),
               function(i) {
                    mbmm_logsumexp(
                         c(
                              log_post_pos[i],
                              log_post_neg[i]
                         )
                    )
               },
               numeric(1)
          )
          
          prior_pos <- exp(log_pi_pos)
          
          posterior_pos <- exp(
               log_post_pos - log_denom
          )
          
          resp <- data.frame(
               gate_prior = prior_pos,
               posterior = posterior_pos,
               stringsAsFactors = FALSE
          )
          
          names(resp) <- c(
               paste0(
                    "gate_prior_",
                    classes[1L]
               ),
               paste0(
                    "posterior_",
                    classes[1L]
               )
          )
          
          resp$map_class <- ifelse(
               posterior_pos >= 0.5,
               classes[1L],
               classes[2L]
          )
     } else {
          eta <- X_full %*% t(full$beta)
          
          eta_full <- cbind(
               eta,
               0
          )
          
          row_max <- apply(
               eta_full,
               1L,
               max
          )
          
          log_denom_pi <- log(
               rowSums(
                    exp(eta_full - row_max)
               )
          ) + row_max
          
          log_pi <- eta_full -
               log_denom_pi
          
          log_post <- logL_mat +
               log_pi
          
          log_denom_post <- apply(
               log_post,
               1L,
               mbmm_logsumexp
          )
          
          post <- exp(
               log_post - log_denom_post
          )
          
          prior <- exp(log_pi)
          
          colnames(post) <- paste0(
               "posterior_",
               classes
          )
          
          colnames(prior) <- paste0(
               "gate_prior_",
               classes
          )
          
          resp <- as.data.frame(
               cbind(
                    prior,
                    post
               ),
               stringsAsFactors = FALSE
          )
          
          resp$map_class <- classes[
               apply(
                    post,
                    1L,
                    which.max
               )
          ]
     }
     
     if (!is.null(id_col)) {
          if (!id_col %in% names(data_used)) {
               stop(
                    "`id_col` was supplied but is not present in `data`.",
                    call. = FALSE
               )
          }
          
          resp <- cbind(
               id = data_used[[id_col]],
               resp
          )
     }
     
     # ---------------------------------------------------------------------------
     # Optional permutation test.
     #
     # This permutes the complete full-model design matrix. It should be used
     # only when all tested non-intercept covariates are exchangeable as a block.
     # ---------------------------------------------------------------------------
     
     p_perm <- NA_real_
     perm_LRTs <- NULL
     
     if (n_perm > 0L) {
          if (!is.null(seed)) {
               set.seed(seed)
          }
          
          n <- nrow(X_full)
          
          perm_LRTs <- rep(
               NA_real_,
               n_perm
          )
          
          for (b in seq_len(n_perm)) {
               perm_index <- sample.int(n)
               
               X_perm <- X_full[
                    perm_index,
                    ,
                    drop = FALSE
               ]
               
               if (K == 2L) {
                    perm_full <- .mbmm_fit_binary_direct(
                         X = X_perm,
                         logL_mat = logL_mat,
                         penalty = penalty,
                         ridge = ridge,
                         maxit = maxit
                    )
               } else {
                    perm_full <- .mbmm_fit_multinom_direct(
                         X = X_perm,
                         logL_mat = logL_mat,
                         K = K,
                         penalty = penalty,
                         ridge = ridge,
                         maxit = maxit
                    )
               }
               
               if (perm_full$convergence == 0L) {
                    perm_LRTs[b] <- 2 * (
                         perm_full$logLik_unpenalized -
                              null$logLik_unpenalized
                    )
               }
          }
          
          valid_perm <- is.finite(
               perm_LRTs
          )
          
          if (!any(valid_perm)) {
               warning(
                    "No permutation gate fits converged successfully.",
                    call. = FALSE
               )
          } else {
               p_perm <- (
                    1 +
                         sum(
                              perm_LRTs[valid_perm] >= LRT
                         )
               ) / (
                    1 +
                         sum(valid_perm)
               )
          }
     }
     
     # ---------------------------------------------------------------------------
     # Coefficient table.
     # ---------------------------------------------------------------------------
     
     coef_df <- if (K == 2L) {
          data.frame(
               term = colnames(X_full),
               estimate = as.numeric(
                    full$beta
               ),
               odds_ratio = exp(
                    as.numeric(full$beta)
               ),
               stringsAsFactors = FALSE
          )
     } else {
          do.call(
               rbind,
               lapply(
                    seq_len(K - 1L),
                    function(k) {
                         data.frame(
                              class = classes[k],
                              term = colnames(X_full),
                              estimate = full$beta[k, ],
                              odds_ratio = exp(
                                   full$beta[k, ]
                              ),
                              stringsAsFactors = FALSE
                         )
                    }
               )
          )
     }
     
     structure(
          list(
               K = K,
               classes = classes,
               method = "direct-marginal-optimization",
               formula = formula,
               
               penalty = penalty,
               ridge = if (
                    identical(penalty, "ridge")
               ) {
                    ridge
               } else {
                    0
               },
               
               coef = coef_df,
               beta = full$beta,
               
               logLik_full =
                    full$logLik_unpenalized,
               
               logLik_null =
                    null$logLik_unpenalized,
               
               objective_full_penalized =
                    full$objective_penalized,
               
               objective_null_penalized =
                    null$objective_penalized,
               
               # Backward-compatible aliases.
               ll_full =
                    full$logLik_unpenalized,
               
               ll_null =
                    null$logLik_unpenalized,
               
               LRT = LRT,
               df = df_lrt,
               p_value = p_value,
               
               n_perm = n_perm,
               p_perm = p_perm,
               perm_LRT = perm_LRTs,
               
               posterior = resp,
               
               X_full = X_full,
               X_null = X_null,
               logL = logL_mat,
               
               convergence_full =
                    full$convergence,
               
               convergence_null =
                    null$convergence,
               
               convergence_message_full =
                    full$message,
               
               convergence_message_null =
                    null$message,
               
               optimizer_counts_full =
                    full$counts,
               
               optimizer_counts_null =
                    null$counts,
               
               rows_used = unname(which(keep)),
               rows_dropped = unname(which(!keep)),
               
               n = nrow(X_full)
          ),
          class = "mbmm_gate"
     )
}


#' Print an mbmm_gate object
#'
#' @param x An \code{"mbmm_gate"} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.mbmm_gate <- function(
          x,
          ...
) {
     cat("MBMM gate model\n")
     
     cat(
          "  classes:",
          paste(
               x$classes,
               collapse = " vs "
          ),
          "\n"
     )
     
     cat(
          "  method:",
          x$method,
          "\n"
     )
     
     cat(
          "  penalty:",
          x$penalty,
          "\n"
     )
     
     if (identical(
          x$penalty,
          "ridge"
     )) {
          cat(
               "  ridge:",
               format(x$ridge),
               "\n"
          )
     }
     
     cat(
          "  n subjects:",
          x$n,
          "\n"
     )
     
     cat(
          "  optimizer convergence: full =",
          x$convergence_full,
          ", null =",
          x$convergence_null,
          "\n"
     )
     
     cat(
          "  unpenalized logLik: full =",
          signif(
               x$logLik_full,
               6
          ),
          ", null =",
          signif(
               x$logLik_null,
               6
          ),
          "\n"
     )
     
     cat(
          "  LRT =",
          signif(
               x$LRT,
               4
          ),
          "on df =",
          x$df,
          " -> p =",
          signif(
               x$p_value,
               4
          ),
          "\n"
     )
     
     if (!is.na(x$p_perm)) {
          cat(
               "  permutation p (valid n = ",
               sum(
                    is.finite(
                         x$perm_LRT
                    )
               ),
               ") = ",
               signif(
                    x$p_perm,
                    4
               ),
               "\n",
               sep = ""
          )
     }
     
     cat(
          "\nGate coefficients:\n"
     )
     
     print(
          x$coef,
          row.names = FALSE
     )
     
     invisible(x)
}


#' Summarize an mbmm_gate object
#'
#' @param object An \code{"mbmm_gate"} object.
#' @param ... Unused.
#' @return \code{object}, invisibly.
#' @export
summary.mbmm_gate <- function(
          object,
          ...
) {
     print(object)
     
     cat(
          "\nPosterior class-membership summary:\n"
     )
     
     print(
          table(
               object$posterior$map_class
          )
     )
     
     invisible(object)
}