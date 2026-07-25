#' Predict fitted MBMM gate probabilities
#'
#' Predicts the gate probabilities implied by the fitted covariate model before
#' incorporating subject-specific longitudinal mechanistic likelihoods. These
#' are therefore prior latent-class probabilities under the MBMM gate.
#'
#' @param fit An object of class \code{"mbmm_gate"} returned by
#'   \code{\link{mbmm_fit_gate}}.
#' @param newdata A data.frame containing the variables used in
#'   \code{fit$formula}.
#'
#' @return For a binary gate, a numeric vector giving the probability of
#'   \code{fit$classes[1]}. For a multinomial gate, a matrix with one column
#'   per latent class.
#'
#' @export
mbmm_predict_gate <- function(
          fit,
          newdata
) {
     if (!inherits(fit, "mbmm_gate")) {
          stop(
               "`fit` must be an object of class \"mbmm_gate\".",
               call. = FALSE
          )
     }
     
     if (!is.data.frame(newdata)) {
          stop(
               "`newdata` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (nrow(newdata) == 0L) {
          stop(
               "`newdata` must contain at least one row.",
               call. = FALSE
          )
     }
     
     required_variables <- all.vars(
          fit$formula
     )
     
     missing_variables <- setdiff(
          required_variables,
          names(newdata)
     )
     
     if (length(missing_variables) > 0L) {
          stop(
               "`newdata` is missing required gate variable(s): ",
               paste(missing_variables, collapse = ", "),
               call. = FALSE
          )
     }
     
     model_frame <- stats::model.frame(
          formula = fit$formula,
          data = newdata,
          na.action = stats::na.pass
     )
     
     X <- stats::model.matrix(
          object = fit$formula,
          data = model_frame
     )
     
     fitted_columns <- colnames(
          fit$X_full
     )
     
     new_columns <- colnames(X)
     
     missing_design_columns <- setdiff(
          fitted_columns,
          new_columns
     )
     
     extra_design_columns <- setdiff(
          new_columns,
          fitted_columns
     )
     
     if (length(missing_design_columns) > 0L ||
         length(extra_design_columns) > 0L) {
          stop(
               paste0(
                    "The design matrix created from `newdata` does not match the fitted ",
                    "gate model. This may be caused by different factor levels, contrasts, ",
                    "or transformed variables."
               ),
               call. = FALSE
          )
     }
     
     # Reorder columns to exactly match the fitted design matrix.
     X <- X[
          ,
          fitted_columns,
          drop = FALSE
     ]
     
     if (any(!is.finite(X))) {
          stop(
               "`newdata` produces missing or non-finite gate covariates.",
               call. = FALSE
          )
     }
     
     if (fit$K == 2L) {
          beta <- as.numeric(
               fit$beta
          )
          
          if (length(beta) != ncol(X)) {
               stop(
                    "The fitted binary gate coefficients do not match the design matrix.",
                    call. = FALSE
               )
          }
          
          eta <- as.numeric(
               X %*% beta
          )
          
          probability <- stats::plogis(
               eta
          )
          
          names(probability) <- rownames(newdata)
          
          return(probability)
     }
     
     beta <- fit$beta
     
     if (!is.matrix(beta) ||
         ncol(beta) != ncol(X) ||
         nrow(beta) != fit$K - 1L) {
          stop(
               "The fitted multinomial gate coefficients do not match the design matrix.",
               call. = FALSE
          )
     }
     
     eta <- X %*% t(beta)
     
     # Final class is the reference class.
     eta_full <- cbind(
          eta,
          0
     )
     
     row_max <- apply(
          eta_full,
          1L,
          max
     )
     
     exp_shifted <- exp(
          eta_full - row_max
     )
     
     probabilities <- exp_shifted /
          rowSums(exp_shifted)
     
     colnames(probabilities) <- fit$classes
     rownames(probabilities) <- rownames(newdata)
     
     probabilities
}


#' Fitted MBMM gate threshold for a continuous covariate
#'
#' For a binary gate, calculates the value of a focal continuous covariate at
#' which the fitted prior probability of the positive class equals
#' \code{level}.
#'
#' For the common intercept-plus-one-covariate model,
#'
#' \deqn{
#' \mathrm{logit}\{\pi(x)\}
#' =
#' \beta_0 + \beta_1 x,
#' }
#'
#' the threshold is
#'
#' \deqn{
#' x^*
#' =
#' \frac{\mathrm{logit}(\mathrm{level})-\beta_0}{\beta_1}.
#' }
#'
#' When additional covariates are included, their values must be supplied
#' through \code{hold_at}. The function constructs the fitted design matrix at
#' two focal values to identify the linear-predictor slope associated with the
#' focal variable. This also correctly incorporates fixed contributions from
#' other additive covariates.
#'
#' This function returns a threshold on the scale used in the gate model. For
#' standardized biomarkers, back-transformation and bootstrap uncertainty
#' should be handled separately.
#'
#' @param fit A binary \code{"mbmm_gate"} object.
#' @param focal Name of the continuous focal variable.
#' @param hold_at Named list providing fixed values for every other variable
#'   in the gate formula.
#' @param level Probability level at which to calculate the crossing.
#' @param slope_tolerance Small positive value used to determine whether the
#'   fitted gate is effectively flat in the focal covariate.
#'
#' @return A numeric threshold on the fitted gate-covariate scale. Returns
#'   \code{NA_real_} with a message when the focal variable has no estimable
#'   linear effect.
#'
#' @examples
#' set.seed(1)
#'
#' n <- 150
#' x <- rnorm(n)
#'
#' z <- rbinom(
#'   n,
#'   1,
#'   mbmm_inv_logit(-0.5 + 1.5 * x)
#' )
#'
#' logL <- data.frame(
#'   logL_pos = ifelse(z == 1, log(0.98), log(0.02)),
#'   logL_neg = ifelse(z == 1, log(0.02), log(0.98))
#' )
#'
#' fit <- mbmm_fit_gate(
#'   logL = logL,
#'   formula = ~x,
#'   data = data.frame(x = x),
#'   classes = c("pos", "neg")
#' )
#'
#' mbmm_gate_threshold(
#'   fit,
#'   focal = "x"
#' )
#'
#' @export
mbmm_gate_threshold <- function(
          fit,
          focal,
          hold_at = list(),
          level = 0.5,
          slope_tolerance = 1e-10
) {
     if (!inherits(fit, "mbmm_gate")) {
          stop(
               "`fit` must be an object of class \"mbmm_gate\".",
               call. = FALSE
          )
     }
     
     if (fit$K != 2L) {
          stop(
               "`mbmm_gate_threshold()` is currently defined only for binary gates.",
               call. = FALSE
          )
     }
     
     if (!is.character(focal) ||
         length(focal) != 1L ||
         is.na(focal) ||
         focal == "") {
          stop(
               "`focal` must be one non-empty variable name.",
               call. = FALSE
          )
     }
     
     if (!is.list(hold_at)) {
          stop(
               "`hold_at` must be a named list.",
               call. = FALSE
          )
     }
     
     if (length(hold_at) > 0L &&
         (is.null(names(hold_at)) ||
          any(names(hold_at) == "") ||
          anyDuplicated(names(hold_at)))) {
          stop(
               "`hold_at` must have unique, non-empty names.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(level) ||
         length(level) != 1L ||
         !is.finite(level) ||
         level <= 0 ||
         level >= 1) {
          stop(
               "`level` must be one finite probability strictly between 0 and 1.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(slope_tolerance) ||
         length(slope_tolerance) != 1L ||
         !is.finite(slope_tolerance) ||
         slope_tolerance <= 0) {
          stop(
               "`slope_tolerance` must be one finite, strictly positive number.",
               call. = FALSE
          )
     }
     
     terms_used <- all.vars(
          fit$formula
     )
     
     if (!focal %in% terms_used) {
          stop(
               "`focal` is not a variable in the fitted gate formula: ",
               paste(deparse(fit$formula), collapse = " "),
               call. = FALSE
          )
     }
     
     other_variables <- setdiff(
          terms_used,
          focal
     )
     
     missing_hold_values <- setdiff(
          other_variables,
          names(hold_at)
     )
     
     if (length(missing_hold_values) > 0L) {
          stop(
               "No fixed value was supplied in `hold_at` for covariate(s): ",
               paste(missing_hold_values, collapse = ", "),
               call. = FALSE
          )
     }
     
     extra_hold_values <- setdiff(
          names(hold_at),
          other_variables
     )
     
     if (length(extra_hold_values) > 0L) {
          stop(
               "`hold_at` contains variable(s) not needed by the gate formula: ",
               paste(extra_hold_values, collapse = ", "),
               call. = FALSE
          )
     }
     
     # Construct two rows differing only by one unit in the focal covariate.
     newdata_zero <- as.list(
          hold_at
     )
     
     newdata_one <- as.list(
          hold_at
     )
     
     newdata_zero[[focal]] <- 0
     newdata_one[[focal]] <- 1
     
     newdata_pair <- rbind(
          as.data.frame(
               newdata_zero,
               stringsAsFactors = FALSE
          ),
          as.data.frame(
               newdata_one,
               stringsAsFactors = FALSE
          )
     )
     
     model_frame <- stats::model.frame(
          formula = fit$formula,
          data = newdata_pair,
          na.action = stats::na.pass
     )
     
     X_pair <- stats::model.matrix(
          object = fit$formula,
          data = model_frame
     )
     
     fitted_columns <- colnames(
          fit$X_full
     )
     
     if (!setequal(colnames(X_pair), fitted_columns)) {
          stop(
               paste0(
                    "The design matrix constructed for the threshold does not match the ",
                    "fitted gate. Check factor levels, contrasts, transformations, and ",
                    "`hold_at` values."
               ),
               call. = FALSE
          )
     }
     
     X_pair <- X_pair[
          ,
          fitted_columns,
          drop = FALSE
     ]
     
     if (any(!is.finite(X_pair))) {
          stop(
               "The supplied threshold covariate values produce a non-finite design matrix.",
               call. = FALSE
          )
     }
     
     beta <- as.numeric(
          fit$beta
     )
     
     eta_at_zero <- as.numeric(
          X_pair[1L, , drop = FALSE] %*% beta
     )
     
     eta_at_one <- as.numeric(
          X_pair[2L, , drop = FALSE] %*% beta
     )
     
     focal_slope <- eta_at_one -
          eta_at_zero
     
     if (!is.finite(focal_slope) ||
         abs(focal_slope) <= slope_tolerance) {
          message(
               paste0(
                    "The fitted gate has no appreciable linear change in the focal ",
                    "covariate; returning NA."
               )
          )
          
          return(NA_real_)
     }
     
     target_eta <- stats::qlogis(
          level
     )
     
     threshold <- (
          target_eta -
               eta_at_zero
     ) / focal_slope
     
     as.numeric(threshold)
}