#' Reference-bank class-conditional log-likelihood for one subject
#'
#' Computes the class-conditional mechanistic likelihood for one subject by
#' averaging trajectory likelihoods over a class-specific reference bank.
#'
#' For equal donor weighting, the function first averages over the retained
#' particles/bundles contributed by each donor and then averages across donors:
#'
#' \deqn{
#' L_{ik}
#' =
#' \frac{1}{D_k}
#' \sum_{d=1}^{D_k}
#' \left[
#' \frac{1}{M_{kd}}
#' \sum_{m=1}^{M_{kd}}
#' p(Y_i \mid Z_i = k, \theta_{kdm})
#' \right].
#' }
#'
#' This two-stage averaging ensures that a donor with more retained particles
#' does not receive greater weight in the reference-bank likelihood.
#'
#' Detected observations contribute a Gaussian density. Left-censored
#' observations contribute the Gaussian probability of falling below the LOD.
#'
#' @param obs This subject's observations. A data.frame containing
#'   \code{time} and \code{y}, and optionally \code{target},
#'   \code{censored}, and \code{lod}.
#' @param bank_pred Long-format bank predictions for one class or phenotype.
#'   Must contain the columns named by \code{donor_col} and
#'   \code{bundle_col}, together with \code{time}, \code{pred}, and optionally
#'   \code{target}.
#' @param sigma Measurement-error SD on the Gaussian likelihood scale. Either
#'   one positive numeric value or, for multitarget data, a named numeric vector
#'   with one value for every target.
#' @param exclude_donor Optional donor ID or vector of donor IDs to remove
#'   before scoring. Use this for leave-one-subject/donor-out scoring.
#' @param donor_col Name of the donor column in \code{bank_pred}.
#' @param bundle_col Name of the within-donor bundle/particle column in
#'   \code{bank_pred}.
#' @param log_scale If \code{TRUE}, \code{y}, \code{pred}, and \code{lod} are
#'   assumed to already be on the log10 scale. If \code{FALSE}, all three are
#'   log10-transformed internally.
#' @param donor_weighting Either \code{"equal"} or \code{"particle"}.
#'   \code{"equal"} first averages particles within donor and then donors,
#'   matching the manuscript reference-bank construction.
#'   \code{"particle"} averages all donor-particle bundles together.
#'
#' @return A scalar class-conditional log-likelihood. Returns \code{-Inf} when
#'   no eligible bank draw overlaps the subject's observations.
#'
#' @examples
#' obs <- data.frame(
#'   time = c(0, 1, 2),
#'   y = c(1, 2, 3)
#' )
#'
#' close_bank <- data.frame(
#'   id = "donorA",
#'   bundle_id = 1,
#'   time = c(0, 1, 2),
#'   pred = c(1.05, 1.95, 3.1)
#' )
#'
#' far_bank <- data.frame(
#'   id = "donorB",
#'   bundle_id = 1,
#'   time = c(0, 1, 2),
#'   pred = c(5, 5, 5)
#' )
#'
#' mbmm_class_loglik(obs, close_bank, sigma = 0.2)
#' mbmm_class_loglik(obs, far_bank, sigma = 0.2)
#'
#' @export
mbmm_class_loglik <- function(
          obs,
          bank_pred,
          sigma,
          exclude_donor = NULL,
          donor_col = "id",
          bundle_col = "bundle_id",
          log_scale = TRUE,
          donor_weighting = c("equal", "particle")
) {
     donor_weighting <- match.arg(donor_weighting)
     
     # ---------------------------------------------------------------------------
     # 1. Validate inputs
     # ---------------------------------------------------------------------------
     
     if (!is.data.frame(obs)) {
          stop("`obs` must be a data.frame.", call. = FALSE)
     }
     
     if (!is.data.frame(bank_pred)) {
          stop("`bank_pred` must be a data.frame.", call. = FALSE)
     }
     
     obs_required <- c("time", "y")
     bank_required <- c(donor_col, bundle_col, "time", "pred")
     
     missing_obs <- setdiff(obs_required, names(obs))
     missing_bank <- setdiff(bank_required, names(bank_pred))
     
     if (length(missing_obs) > 0L) {
          stop(
               "`obs` is missing required column(s): ",
               paste(missing_obs, collapse = ", "),
               call. = FALSE
          )
     }
     
     if (length(missing_bank) > 0L) {
          stop(
               "`bank_pred` is missing required column(s): ",
               paste(missing_bank, collapse = ", "),
               call. = FALSE
          )
     }
     
     obs_has_target <- "target" %in% names(obs)
     bank_has_target <- "target" %in% names(bank_pred)
     
     if (xor(obs_has_target, bank_has_target)) {
          stop(
               "`target` must appear in both `obs` and `bank_pred`, or in neither.",
               call. = FALSE
          )
     }
     
     has_target <- obs_has_target && bank_has_target
     
     if (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale)) {
          stop("`log_scale` must be TRUE or FALSE.", call. = FALSE)
     }
     
     if (!is.numeric(sigma) || length(sigma) < 1L) {
          stop("`sigma` must be numeric and non-empty.", call. = FALSE)
     }
     
     if (any(!is.finite(as.numeric(sigma))) || any(as.numeric(sigma) <= 0)) {
          stop("All values of `sigma` must be finite and strictly positive.",
               call. = FALSE)
     }
     
     # ---------------------------------------------------------------------------
     # 2. Apply donor exclusion
     # ---------------------------------------------------------------------------
     
     if (!is.null(exclude_donor)) {
          keep <- !(as.character(bank_pred[[donor_col]]) %in%
                         as.character(exclude_donor))
          
          bank_pred <- bank_pred[keep, , drop = FALSE]
     }
     
     if (nrow(bank_pred) == 0L) {
          return(-Inf)
     }
     
     # ---------------------------------------------------------------------------
     # 3. Validate observation and bank keys
     # ---------------------------------------------------------------------------
     
     obs_key_cols <- if (has_target) c("target", "time") else "time"
     
     bank_key_cols <- c(
          donor_col,
          bundle_col,
          if (has_target) "target",
          "time"
     )
     
     if (anyDuplicated(obs[obs_key_cols])) {
          stop(
               "Duplicate observation keys detected. Each subject must have at most ",
               "one observation for each ",
               if (has_target) "`target`-`time`" else "`time`",
               " combination.",
               call. = FALSE
          )
     }
     
     if (anyDuplicated(bank_pred[bank_key_cols])) {
          stop(
               "Duplicate bank prediction keys detected. Each donor-bundle",
               if (has_target) "-target" else "",
               "-time combination must be unique.",
               call. = FALSE
          )
     }
     
     # ---------------------------------------------------------------------------
     # 4. Prepare likelihood-scale values
     # ---------------------------------------------------------------------------
     
     obs2 <- obs
     bank2 <- bank_pred
     
     if (log_scale) {
          obs2$y_use <- obs2$y
          bank2$pred_use <- bank2$pred
          
          if ("lod" %in% names(obs2)) {
               obs2$lod_use <- obs2$lod
          }
     } else {
          obs2$y_use <- mbmm_safe_log10(obs2$y)
          bank2$pred_use <- mbmm_safe_log10(bank2$pred)
          
          if ("lod" %in% names(obs2)) {
               obs2$lod_use <- mbmm_safe_log10(obs2$lod)
          }
     }
     
     if (any(!is.finite(bank2$pred_use))) {
          stop(
               "`bank_pred$pred` contains non-finite values on the likelihood scale.",
               call. = FALSE
          )
     }
     
     # Detected observations must have a usable y value.
     censored_obs <- if ("censored" %in% names(obs2)) {
          as.logical(obs2$censored)
     } else {
          rep(FALSE, nrow(obs2))
     }
     
     censored_obs[is.na(censored_obs)] <- FALSE
     
     if (any(!censored_obs & !is.finite(obs2$y_use))) {
          stop(
               "Detected observations contain non-finite `y` values on the ",
               "likelihood scale.",
               call. = FALSE
          )
     }
     
     if (any(censored_obs)) {
          if (!"lod_use" %in% names(obs2)) {
               stop(
                    "Censored observations require an `lod` column.",
                    call. = FALSE
               )
          }
          
          if (any(!is.finite(obs2$lod_use[censored_obs]))) {
               stop(
                    "All censored observations must have a finite LOD.",
                    call. = FALSE
               )
          }
     }
     
     # ---------------------------------------------------------------------------
     # 5. Merge observations with every eligible bank trajectory
     # ---------------------------------------------------------------------------
     
     join_cols <- if (has_target) c("target", "time") else "time"
     
     merged <- merge(
          obs2,
          bank2,
          by = join_cols,
          suffixes = c("", ".bank"),
          all = FALSE,
          sort = FALSE
     )
     
     if (nrow(merged) == 0L) {
          return(-Inf)
     }
     
     # ---------------------------------------------------------------------------
     # 6. Assign measurement-error SD
     # ---------------------------------------------------------------------------
     
     if (has_target) {
          target_values <- as.character(merged$target)
          
          if (length(sigma) == 1L) {
               sig_vec <- rep(as.numeric(sigma), nrow(merged))
          } else {
               if (is.null(names(sigma)) || any(names(sigma) == "")) {
                    stop(
                         "For multitarget data, a multi-value `sigma` must be named by target.",
                         call. = FALSE
                    )
               }
               
               missing_sigma_targets <- setdiff(unique(target_values), names(sigma))
               
               if (length(missing_sigma_targets) > 0L) {
                    stop(
                         "Missing `sigma` value(s) for target(s): ",
                         paste(missing_sigma_targets, collapse = ", "),
                         call. = FALSE
                    )
               }
               
               sig_vec <- as.numeric(sigma[target_values])
          }
     } else {
          if (length(sigma) != 1L) {
               stop(
                    "For single-target data, `sigma` must contain exactly one value.",
                    call. = FALSE
               )
          }
          
          sig_vec <- rep(as.numeric(sigma), nrow(merged))
     }
     
     if (any(!is.finite(sig_vec)) || any(sig_vec <= 0)) {
          stop(
               "Resolved measurement-error SD values must be finite and positive.",
               call. = FALSE
          )
     }
     
     # ---------------------------------------------------------------------------
     # 7. Evaluate detected and censored pointwise likelihood contributions
     # ---------------------------------------------------------------------------
     
     censored <- if ("censored" %in% names(merged)) {
          as.logical(merged$censored)
     } else {
          rep(FALSE, nrow(merged))
     }
     
     censored[is.na(censored)] <- FALSE
     
     lod_use <- if ("lod_use" %in% names(merged)) {
          merged$lod_use
     } else {
          NULL
     }
     
     merged$ll_point <- mbmm_dnorm_censored(
          y = merged$y_use,
          mean = merged$pred_use,
          sd = sig_vec,
          censored = censored,
          lod = lod_use
     )
     
     if (anyNA(merged$ll_point)) {
          stop(
               "The pointwise mechanistic likelihood produced NA values.",
               call. = FALSE
          )
     }
     
     # ---------------------------------------------------------------------------
     # 8. Define a unique mechanistic draw using donor AND bundle
     # ---------------------------------------------------------------------------
     
     merged$.donor_id <- as.character(merged[[donor_col]])
     
     merged$.draw_id <- interaction(
          merged[[donor_col]],
          merged[[bundle_col]],
          drop = TRUE,
          lex.order = TRUE
     )
     
     # Sum the longitudinal/target-specific contributions within each complete
     # mechanistic draw.
     draw_ll <- tapply(
          merged$ll_point,
          INDEX = merged$.draw_id,
          FUN = sum
     )
     
     if (length(draw_ll) == 0L) {
          return(-Inf)
     }
     
     # Map every globally unique draw back to its donor.
     draw_donor <- tapply(
          merged$.donor_id,
          INDEX = merged$.draw_id,
          FUN = function(x) x[[1L]]
     )
     
     # ---------------------------------------------------------------------------
     # 9. Monte Carlo averaging
     # ---------------------------------------------------------------------------
     
     if (donor_weighting == "particle") {
          # Ordinary Monte Carlo mean over all available donor-particle draws.
          return(mbmm_logmeanexp(as.numeric(draw_ll)))
     }
     
     # Manuscript-oriented equal donor weighting:
     #   1. average particles within each donor;
     #   2. average donor-specific likelihoods across donors.
     draw_names <- names(draw_ll)
     donor_for_draw <- as.character(draw_donor[draw_names])
     
     draw_ll_by_donor <- split(
          as.numeric(draw_ll),
          donor_for_draw
     )
     
     donor_loglik <- vapply(
          draw_ll_by_donor,
          FUN = mbmm_logmeanexp,
          FUN.VALUE = numeric(1)
     )
     
     mbmm_logmeanexp(donor_loglik)
}


#' Analytic non-establishment class likelihood for COVID5001-style analyses
#'
#' Computes the class-conditional likelihood for a non-establishment or
#' nonproductive class whose expected viral load is a fixed fraction of the
#' observation-specific LOD. In the COVID5001 manuscript analysis, the default
#' mean was LOD / 2 on the natural viral-load scale.
#'
#' This class is analytic and does not require a reference trajectory bank.
#'
#' @param obs Subject observations with columns \code{y}, \code{lod}, and
#'   optionally \code{target} and \code{censored}.
#' @param sigma Measurement-error SD on the log10 likelihood scale. Either a
#'   scalar or a target-named numeric vector.
#' @param lod_fraction Expected non-establishment mean as a fraction of the LOD.
#'   The manuscript default is \code{0.5}.
#' @param log_scale Whether \code{obs$y} and \code{obs$lod} are already on the
#'   log10 scale. When \code{FALSE}, the non-establishment mean is constructed
#'   as \code{lod_fraction * lod} on the natural scale and then log10-transformed.
#'
#' @return A scalar log-likelihood under the analytic non-establishment model.
#'
#' @export
mbmm_noise_floor_loglik <- function(
          obs,
          sigma,
          lod_fraction = 0.5,
          log_scale = FALSE
) {
     if (!is.data.frame(obs)) {
          stop("`obs` must be a data.frame.", call. = FALSE)
     }
     
     required <- c("y", "lod")
     missing_cols <- setdiff(required, names(obs))
     
     if (length(missing_cols) > 0L) {
          stop(
               "`obs` is missing required column(s): ",
               paste(missing_cols, collapse = ", "),
               call. = FALSE
          )
     }
     
     if (!is.numeric(lod_fraction) ||
         length(lod_fraction) != 1L ||
         !is.finite(lod_fraction) ||
         lod_fraction <= 0) {
          stop(
               "`lod_fraction` must be one finite positive number.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(sigma) ||
         length(sigma) < 1L ||
         any(!is.finite(as.numeric(sigma))) ||
         any(as.numeric(sigma) <= 0)) {
          stop(
               "`sigma` must contain finite, strictly positive numeric value(s).",
               call. = FALSE
          )
     }
     
     has_target <- "target" %in% names(obs)
     
     if (has_target) {
          target_values <- as.character(obs$target)
          
          if (length(sigma) == 1L) {
               sig_vec <- rep(as.numeric(sigma), nrow(obs))
          } else {
               if (is.null(names(sigma)) || any(names(sigma) == "")) {
                    stop(
                         "For multitarget data, a multi-value `sigma` must be named by target.",
                         call. = FALSE
                    )
               }
               
               missing_sigma_targets <- setdiff(unique(target_values), names(sigma))
               
               if (length(missing_sigma_targets) > 0L) {
                    stop(
                         "Missing `sigma` value(s) for target(s): ",
                         paste(missing_sigma_targets, collapse = ", "),
                         call. = FALSE
                    )
               }
               
               sig_vec <- as.numeric(sigma[target_values])
          }
     } else {
          if (length(sigma) != 1L) {
               stop(
                    "For single-target data, `sigma` must contain exactly one value.",
                    call. = FALSE
               )
          }
          
          sig_vec <- rep(as.numeric(sigma), nrow(obs))
     }
     
     censored <- if ("censored" %in% names(obs)) {
          as.logical(obs$censored)
     } else {
          rep(FALSE, nrow(obs))
     }
     
     censored[is.na(censored)] <- FALSE
     
     if (log_scale) {
          y_use <- obs$y
          lod_use <- obs$lod
          
          # If LOD is already log10 transformed, multiplication by lod_fraction must
          # be represented as addition on the log10 scale.
          mean_use <- lod_use + log10(lod_fraction)
     } else {
          y_use <- mbmm_safe_log10(obs$y)
          lod_use <- mbmm_safe_log10(obs$lod)
          mean_use <- mbmm_safe_log10(lod_fraction * obs$lod)
     }
     
     if (any(!censored & !is.finite(y_use))) {
          stop(
               "Detected observations contain non-finite `y` values.",
               call. = FALSE
          )
     }
     
     if (any(!is.finite(lod_use))) {
          stop(
               "All observations must have a finite LOD for the noise-floor model.",
               call. = FALSE
          )
     }
     
     ll_point <- mbmm_dnorm_censored(
          y = y_use,
          mean = mean_use,
          sd = sig_vec,
          censored = censored,
          lod = lod_use
     )
     
     if (anyNA(ll_point)) {
          stop(
               "The analytic non-establishment likelihood produced NA values.",
               call. = FALSE
          )
     }
     
     sum(ll_point)
}


#' Score every subject against class-specific reference banks
#'
#' Builds the subject-by-class matrix of class-conditional log-likelihoods used
#' by the MBMM gate model.
#'
#' Leave-one-subject/donor-out scoring is controlled by \code{loo}; it is not
#' applied automatically unless requested. This allows the same function to
#' reproduce application-specific choices:
#'
#' \itemize{
#'   \item SDY2354: typically \code{loo = TRUE};
#'   \item COVID5001 productive bank: use the setting specified by the final
#'   manuscript analysis.
#' }
#'
#' @param obs_list Named list of subject observation data.frames.
#' @param bank_pred_list Named list of reference-bank prediction tables, one
#'   per bank-based class.
#' @param sigma Measurement-error SD passed to
#'   \code{\link{mbmm_class_loglik}}.
#' @param donor_col Name of the donor column in each bank table.
#' @param bundle_col Name of the bundle column in each bank table.
#' @param log_scale Whether observations, predictions, and LODs are already on
#'   the log10 likelihood scale.
#' @param donor_weighting Either \code{"equal"} or \code{"particle"}.
#' @param loo Logical scalar or named logical vector. If a scalar, the same LOSO
#'   rule is applied to every class. If named, names must correspond to
#'   \code{bank_pred_list}.
#' @param verbose Print progress every ten subjects.
#'
#' @return A data.frame with one row per subject and one
#'   \code{logL_<class>} column per bank.
#'
#' @export
mbmm_score_subjects <- function(
          obs_list,
          bank_pred_list,
          sigma,
          donor_col = "id",
          bundle_col = "bundle_id",
          log_scale = TRUE,
          donor_weighting = c("equal", "particle"),
          loo = TRUE,
          verbose = TRUE
) {
     donor_weighting <- match.arg(donor_weighting)
     
     if (!is.list(obs_list) || length(obs_list) == 0L) {
          stop("`obs_list` must be a non-empty named list.", call. = FALSE)
     }
     
     if (!is.list(bank_pred_list) || length(bank_pred_list) == 0L) {
          stop("`bank_pred_list` must be a non-empty named list.", call. = FALSE)
     }
     
     ids <- names(obs_list)
     classes <- names(bank_pred_list)
     
     if (is.null(ids) || any(ids == "") || anyDuplicated(ids)) {
          stop(
               "`obs_list` must have unique, non-empty subject-ID names.",
               call. = FALSE
          )
     }
     
     if (is.null(classes) || any(classes == "") || anyDuplicated(classes)) {
          stop(
               "`bank_pred_list` must have unique, non-empty class names.",
               call. = FALSE
          )
     }
     
     # Resolve class-specific LOSO behavior.
     if (length(loo) == 1L) {
          loo_by_class <- stats::setNames(
               rep(as.logical(loo), length(classes)),
               classes
          )
     } else {
          if (is.null(names(loo))) {
               stop(
                    "When `loo` has more than one value, it must be named by class.",
                    call. = FALSE
               )
          }
          
          missing_loo_classes <- setdiff(classes, names(loo))
          
          if (length(missing_loo_classes) > 0L) {
               stop(
                    "Missing `loo` setting for class(es): ",
                    paste(missing_loo_classes, collapse = ", "),
                    call. = FALSE
               )
          }
          
          loo_by_class <- as.logical(loo[classes])
          names(loo_by_class) <- classes
     }
     
     if (anyNA(loo_by_class)) {
          stop("All `loo` values must be TRUE or FALSE.", call. = FALSE)
     }
     
     logL <- matrix(
          NA_real_,
          nrow = length(ids),
          ncol = length(classes),
          dimnames = list(ids, classes)
     )
     
     for (i in seq_along(ids)) {
          sid <- ids[[i]]
          obs_i <- obs_list[[sid]]
          
          for (k in classes) {
               bank_k <- bank_pred_list[[k]]
               
               if (!is.data.frame(bank_k)) {
                    stop(
                         "Bank for class `", k, "` is not a data.frame.",
                         call. = FALSE
                    )
               }
               
               if (!donor_col %in% names(bank_k)) {
                    stop(
                         "Bank for class `", k, "` does not contain donor column `",
                         donor_col, "`.",
                         call. = FALSE
                    )
               }
               
               exclude <- NULL
               
               if (
                    isTRUE(loo_by_class[[k]]) &&
                    sid %in% as.character(unique(bank_k[[donor_col]]))
               ) {
                    exclude <- sid
               }
               
               logL[sid, k] <- mbmm_class_loglik(
                    obs = obs_i,
                    bank_pred = bank_k,
                    sigma = sigma,
                    exclude_donor = exclude,
                    donor_col = donor_col,
                    bundle_col = bundle_col,
                    log_scale = log_scale,
                    donor_weighting = donor_weighting
               )
          }
          
          if (isTRUE(verbose) && i %% 10L == 0L) {
               message(
                    sprintf("Scored %d of %d subjects.", i, length(ids))
               )
          }
     }
     
     out <- as.data.frame(logL, stringsAsFactors = FALSE)
     names(out) <- paste0("logL_", classes)
     out$id <- ids
     
     out <- out[, c("id", paste0("logL_", classes)), drop = FALSE]
     rownames(out) <- NULL
     
     out
}