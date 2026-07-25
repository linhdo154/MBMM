#' Back-transform an MBMM gate threshold
#'
#' Converts a threshold from the covariate scale used in the gate model back
#' to the original biomarker scale.
#'
#' If the fitted covariate is a standardized variable
#'
#' \deqn{
#' z = \frac{x - \mu}{s},
#' }
#'
#' the model-scale threshold is first converted to
#'
#' \deqn{
#' x = \mu + sz.
#' }
#'
#' An optional inverse transformation is then applied. For example, when the
#' standardized covariate was constructed from \code{log10(biomarker)}, use
#' \code{transform = "log10"}.
#'
#' @param threshold Numeric threshold on the fitted gate-covariate scale.
#' @param center Center used when standardizing the covariate.
#' @param scale Scale used when standardizing the covariate.
#' @param transform Transformation applied before standardization. One of
#'   \code{"identity"}, \code{"log10"}, or \code{"log"}.
#'
#' @return Numeric threshold on the original biomarker scale.
#'
#' @examples
#' mbmm_backtransform_threshold(
#'   threshold = 1.2,
#'   center = 2,
#'   scale = 0.5,
#'   transform = "log10"
#' )
#'
#' @export
mbmm_backtransform_threshold <- function(
          threshold,
          center = 0,
          scale = 1,
          transform = c("identity", "log10", "log")
) {
     transform <- match.arg(transform)
     
     if (!is.numeric(threshold)) {
          stop(
               "`threshold` must be numeric.",
               call. = FALSE
          )
     }
     
     if (
          !is.numeric(center) ||
          length(center) != 1L ||
          !is.finite(center)
     ) {
          stop(
               "`center` must be one finite numeric value.",
               call. = FALSE
          )
     }
     
     if (
          !is.numeric(scale) ||
          length(scale) != 1L ||
          !is.finite(scale) ||
          scale <= 0
     ) {
          stop(
               "`scale` must be one finite, strictly positive numeric value.",
               call. = FALSE
          )
     }
     
     transformed_scale <- center + scale * threshold
     
     switch(
          transform,
          identity = transformed_scale,
          log10 = 10^transformed_scale,
          log = exp(transformed_scale)
     )
}


#' Bootstrap an MBMM gate threshold
#'
#' Estimates a fitted gate threshold and its uncertainty by resampling subjects,
#' refitting the MBMM gate, and recalculating the threshold in each bootstrap
#' sample.
#'
#' The function reports thresholds both on the covariate scale used in the gate
#' model and, optionally, on the original biomarker scale.
#'
#' A bootstrap estimate is classified as extrapolated when its model-scale
#' threshold lies outside the observed range of the focal gate covariate.
#'
#' @param logL Data.frame containing the class-conditional log-likelihood
#'   columns used by \code{\link{mbmm_fit_gate}}.
#' @param formula One-sided gate-model formula.
#' @param data Subject-level covariate data.frame.
#' @param focal Name of the continuous covariate whose threshold is requested.
#' @param classes Optional latent-class order passed to
#'   \code{\link{mbmm_fit_gate}}.
#' @param hold_at Named list fixing all other covariates in the gate formula.
#' @param level Gate probability at which the threshold is calculated.
#' @param center Center used to standardize the focal covariate.
#' @param scale Scale used to standardize the focal covariate.
#' @param transform Transformation applied before standardization. One of
#'   \code{"identity"}, \code{"log10"}, or \code{"log"}.
#' @param observed_range Optional two-element numeric vector giving the
#'   observed range of the focal covariate on the fitted model scale. By
#'   default, the range is calculated from \code{data[[focal]]}.
#' @param n_boot Number of subject-level bootstrap samples.
#' @param conf Confidence level for the percentile interval.
#' @param ridge Ridge parameter passed to \code{\link{mbmm_fit_gate}}.
#' @param penalty Either \code{"ridge"} or \code{"none"}.
#' @param id_col Subject-ID column in \code{data}.
#' @param seed Optional bootstrap seed.
#' @param maxit Maximum gate optimizer iterations.
#' @param verbose Print bootstrap progress.
#'
#' @return An object of class \code{"mbmm_threshold_summary"} containing the
#'   point estimate, percentile confidence interval, bootstrap estimates,
#'   bootstrap success rate, and extrapolation percentage.
#'
#' @examples
#' \dontrun{
#' result <- mbmm_bootstrap_gate_threshold(
#'   logL = logL,
#'   formula = ~ biomarker_z,
#'   data = biomarker_data,
#'   focal = "biomarker_z",
#'   classes = c("productive", "nonestablishment"),
#'   center = mean(log10(biomarker_data$biomarker)),
#'   scale = sd(log10(biomarker_data$biomarker)),
#'   transform = "log10",
#'   n_boot = 1000,
#'   seed = 1
#' )
#'
#' result
#' }
#'
#' @export
mbmm_bootstrap_gate_threshold <- function(
          logL,
          formula,
          data,
          focal,
          classes = NULL,
          hold_at = list(),
          level = 0.5,
          center = 0,
          scale = 1,
          transform = c("identity", "log10", "log"),
          observed_range = NULL,
          n_boot = 1000,
          conf = 0.95,
          ridge = 1e-4,
          penalty = "ridge",
          id_col = "id",
          seed = NULL,
          maxit = 4000,
          verbose = TRUE
) {
     transform <- match.arg(transform)
     
     if (!is.data.frame(logL)) {
          stop(
               "`logL` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (!is.data.frame(data)) {
          stop(
               "`data` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (nrow(logL) != nrow(data)) {
          stop(
               "`logL` and `data` must contain the same number of subjects.",
               call. = FALSE
          )
     }
     
     if (nrow(data) < 2L) {
          stop(
               "At least two subjects are required.",
               call. = FALSE
          )
     }
     
     if (
          !is.character(focal) ||
          length(focal) != 1L ||
          is.na(focal) ||
          focal == ""
     ) {
          stop(
               "`focal` must be one non-empty variable name.",
               call. = FALSE
          )
     }
     
     if (!focal %in% names(data)) {
          stop(
               "`focal` is not a column in `data`: ",
               focal,
               call. = FALSE
          )
     }
     
     if (!is.numeric(data[[focal]])) {
          stop(
               "The focal threshold covariate must be numeric.",
               call. = FALSE
          )
     }
     
     if (!id_col %in% names(data)) {
          stop(
               "`data` does not contain subject-ID column `",
               id_col,
               "`.",
               call. = FALSE
          )
     }
     
     if (!"id" %in% names(logL)) {
          stop(
               "`logL` must contain an `id` column.",
               call. = FALSE
          )
     }
     
     if (
          !is.numeric(n_boot) ||
          length(n_boot) != 1L ||
          !is.finite(n_boot) ||
          n_boot < 1 ||
          n_boot != floor(n_boot)
     ) {
          stop(
               "`n_boot` must be one positive integer.",
               call. = FALSE
          )
     }
     
     n_boot <- as.integer(n_boot)
     
     if (
          !is.numeric(conf) ||
          length(conf) != 1L ||
          !is.finite(conf) ||
          conf <= 0 ||
          conf >= 1
     ) {
          stop(
               "`conf` must be strictly between 0 and 1.",
               call. = FALSE
          )
     }
     
     allowed_penalties <- c("ridge", "none")
     
     if (
          !is.character(penalty) ||
          length(penalty) != 1L ||
          !penalty %in% allowed_penalties
     ) {
          stop(
               "`penalty` must be one of \"ridge\" or \"none\".",
               call. = FALSE
          )
     }
     
     if (is.null(observed_range)) {
          focal_values <- data[[focal]]
          
          if (any(!is.finite(focal_values))) {
               stop(
                    "`data[[focal]]` contains missing or non-finite values.",
                    call. = FALSE
               )
          }
          
          observed_range <- range(focal_values)
     } else {
          if (
               !is.numeric(observed_range) ||
               length(observed_range) != 2L ||
               any(!is.finite(observed_range)) ||
               observed_range[1L] > observed_range[2L]
          ) {
               stop(
                    "`observed_range` must be a finite numeric vector of length two.",
                    call. = FALSE
               )
          }
     }
     
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     # ---------------------------------------------------------------------------
     # Fit the original sample.
     # ---------------------------------------------------------------------------
     
     original_fit <- mbmm_fit_gate(
          logL = logL,
          formula = formula,
          data = data,
          classes = classes,
          ridge = ridge,
          penalty = penalty,
          n_perm = 0,
          id_col = id_col,
          maxit = maxit
     )
     
     classes_used <- original_fit$classes
     
     point_model <- mbmm_gate_threshold(
          fit = original_fit,
          focal = focal,
          hold_at = hold_at,
          level = level
     )
     
     point_original <- mbmm_backtransform_threshold(
          threshold = point_model,
          center = center,
          scale = scale,
          transform = transform
     )
     
     # ---------------------------------------------------------------------------
     # Subject-level nonparametric bootstrap.
     # ---------------------------------------------------------------------------
     
     n_subjects <- nrow(data)
     
     boot_model <- rep(
          NA_real_,
          n_boot
     )
     
     fit_success <- rep(
          FALSE,
          n_boot
     )
     
     optimizer_success <- rep(
          FALSE,
          n_boot
     )
     
     for (b in seq_len(n_boot)) {
          sampled_rows <- sample.int(
               n = n_subjects,
               size = n_subjects,
               replace = TRUE
          )
          
          logL_b <- logL[
               sampled_rows,
               ,
               drop = FALSE
          ]
          
          data_b <- data[
               sampled_rows,
               ,
               drop = FALSE
          ]
          
          # Bootstrap samples contain repeated subjects. Assign synthetic IDs so
          # each sampled row is treated as a distinct bootstrap observation.
          bootstrap_ids <- paste0(
               "bootstrap_",
               b,
               "_",
               seq_len(n_subjects)
          )
          
          logL_b$id <- bootstrap_ids
          data_b[[id_col]] <- bootstrap_ids
          
          fit_b <- tryCatch(
               mbmm_fit_gate(
                    logL = logL_b,
                    formula = formula,
                    data = data_b,
                    classes = classes_used,
                    ridge = ridge,
                    penalty = penalty,
                    n_perm = 0,
                    id_col = id_col,
                    maxit = maxit
               ),
               error = function(e) {
                    NULL
               }
          )
          
          if (!is.null(fit_b)) {
               fit_success[b] <- TRUE
               
               convergence_full <- fit_b$convergence_full
               convergence_null <- fit_b$convergence_null
               
               optimizer_success[b] <-
                    isTRUE(convergence_full == 0) &&
                    isTRUE(convergence_null == 0)
               
               threshold_b <- tryCatch(
                    mbmm_gate_threshold(
                         fit = fit_b,
                         focal = focal,
                         hold_at = hold_at,
                         level = level
                    ),
                    error = function(e) {
                         NA_real_
                    }
               )
               
               if (
                    length(threshold_b) == 1L &&
                    is.finite(threshold_b)
               ) {
                    boot_model[b] <- threshold_b
               }
          }
          
          if (
               isTRUE(verbose) &&
               (b %% 50L == 0L || b == n_boot)
          ) {
               message(
                    sprintf(
                         "Completed %d of %d bootstrap samples.",
                         b,
                         n_boot
                    )
               )
          }
     }
     
     successful <- is.finite(
          boot_model
     )
     
     n_successful <- sum(
          successful
     )
     
     if (n_successful == 0L) {
          stop(
               "No bootstrap sample produced a finite threshold.",
               call. = FALSE
          )
     }
     
     boot_original <- rep(
          NA_real_,
          n_boot
     )
     
     boot_original[successful] <-
          mbmm_backtransform_threshold(
               threshold = boot_model[successful],
               center = center,
               scale = scale,
               transform = transform
          )
     
     alpha <- 1 - conf
     
     ci_probabilities <- c(
          alpha / 2,
          1 - alpha / 2
     )
     
     ci_model <- as.numeric(
          stats::quantile(
               boot_model[successful],
               probs = ci_probabilities,
               names = FALSE,
               na.rm = TRUE,
               type = 7
          )
     )
     
     ci_original <- as.numeric(
          stats::quantile(
               boot_original[successful],
               probs = ci_probabilities,
               names = FALSE,
               na.rm = TRUE,
               type = 7
          )
     )
     
     extrapolated <- successful &
          (
               boot_model < observed_range[1L] |
                    boot_model > observed_range[2L]
          )
     
     pct_extrapolated <- 100 *
          sum(extrapolated) /
          n_successful
     
     bootstrap_table <- data.frame(
          replicate = seq_len(n_boot),
          threshold_model = boot_model,
          threshold_original = boot_original,
          fit_success = fit_success,
          optimizer_success = optimizer_success,
          extrapolated = ifelse(
               successful,
               extrapolated,
               NA
          )
     )
     
     structure(
          list(
               point = list(
                    threshold_model = point_model,
                    threshold_original = point_original
               ),
               confidence_interval = list(
                    level = conf,
                    model_scale = ci_model,
                    original_scale = ci_original
               ),
               observed_range = observed_range,
               pct_extrapolated = pct_extrapolated,
               bootstrap = bootstrap_table,
               diagnostics = list(
                    n_requested = n_boot,
                    n_fit_success = sum(fit_success),
                    n_optimizer_success = sum(optimizer_success),
                    n_threshold_success = n_successful,
                    threshold_success_rate = n_successful / n_boot
               ),
               settings = list(
                    focal = focal,
                    hold_at = hold_at,
                    probability_level = level,
                    center = center,
                    scale = scale,
                    transform = transform,
                    classes = classes_used,
                    penalty = penalty,
                    ridge = ridge
               ),
               gate = original_fit
          ),
          class = "mbmm_threshold_summary"
     )
}


#' Print an MBMM threshold summary
#'
#' @param x An object of class \code{"mbmm_threshold_summary"}.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.mbmm_threshold_summary <- function(
          x,
          digits = 4,
          ...
) {
     cat("MBMM gate-threshold summary\n")
     
     cat(
          "  focal covariate:",
          x$settings$focal,
          "\n"
     )
     
     cat(
          "  probability level:",
          format(x$settings$probability_level, digits = digits),
          "\n"
     )
     
     cat(
          "  threshold, model scale:",
          format(x$point$threshold_model, digits = digits),
          "\n"
     )
     
     cat(
          "  threshold, original scale:",
          format(x$point$threshold_original, digits = digits),
          "\n"
     )
     
     cat(
          sprintf(
               "  %.1f%% CI, original scale: [%s, %s]\n",
               100 * x$confidence_interval$level,
               format(
                    x$confidence_interval$original_scale[1L],
                    digits = digits
               ),
               format(
                    x$confidence_interval$original_scale[2L],
                    digits = digits
               )
          )
     )
     
     cat(
          "  bootstrap thresholds retained:",
          x$diagnostics$n_threshold_success,
          "of",
          x$diagnostics$n_requested,
          "\n"
     )
     
     cat(
          "  extrapolated bootstrap thresholds:",
          format(
               x$pct_extrapolated,
               digits = digits
          ),
          "%\n"
     )
     
     invisible(x)
}


#' Summarize an MBMM threshold summary
#'
#' @param object An object of class \code{"mbmm_threshold_summary"}.
#' @param ... Unused.
#'
#' @return \code{object}, invisibly.
#'
#' @export
summary.mbmm_threshold_summary <- function(
          object,
          ...
) {
     print(object)
     
     invisible(object)
}