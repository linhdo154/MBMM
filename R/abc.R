#' Default ABC discrepancy: censoring-aware weighted RMSE on the log10 scale
#'
#' Detected observations contribute squared error on the log10 scale.
#' Left-censored observations are penalized only when the predicted trajectory
#' exceeds the observation-specific LOD.
#'
#' @param y_obs Observed values on the natural scale.
#' @param pred Predicted values on the natural scale.
#' @param censored Logical or 0/1 vector indicating left-censored observations.
#' @param lod Observation-specific LOD on the natural scale.
#' @param w_detect Weight for detected-observation squared errors.
#' @param w_censored Weight for censored-observation squared errors.
#' @param floor Smallest positive value passed to log10.
#'
#' @return A nonnegative scalar discrepancy.
#' @export
mbmm_default_distance <- function(
          y_obs,
          pred,
          censored = NULL,
          lod = NULL,
          w_detect = 5,
          w_censored = 1,
          floor = 1e-2
) {
     y_obs <- as.numeric(y_obs)
     pred <- as.numeric(pred)
     
     n <- length(y_obs)
     
     if (length(pred) != n) {
          stop(
               "`pred` must have the same length as `y_obs`.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(w_detect) ||
         length(w_detect) != 1L ||
         !is.finite(w_detect) ||
         w_detect < 0) {
          stop(
               "`w_detect` must be one finite, nonnegative number.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(w_censored) ||
         length(w_censored) != 1L ||
         !is.finite(w_censored) ||
         w_censored < 0) {
          stop(
               "`w_censored` must be one finite, nonnegative number.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(floor) ||
         length(floor) != 1L ||
         !is.finite(floor) ||
         floor <= 0) {
          stop(
               "`floor` must be one finite, strictly positive number.",
               call. = FALSE
          )
     }
     
     if (is.null(censored)) {
          censored <- rep(FALSE, n)
     }
     
     censored <- as.logical(censored)
     
     if (length(censored) != n || anyNA(censored)) {
          stop(
               "`censored` must have the same length as `y_obs` and contain no missing values.",
               call. = FALSE
          )
     }
     
     detected <- !censored
     
     if (any(detected) && any(!is.finite(y_obs[detected]))) {
          stop(
               "Detected observations must have finite `y_obs` values.",
               call. = FALSE
          )
     }
     
     if (any(!is.finite(pred))) {
          stop(
               "`pred` must contain only finite values.",
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
          
          lod <- as.numeric(lod)
          
          if (length(lod) != n) {
               stop(
                    "`lod` must have the same length as `y_obs`.",
                    call. = FALSE
               )
          }
          
          if (any(!is.finite(lod[censored]))) {
               stop(
                    "Censored observations must have finite LOD values.",
                    call. = FALSE
               )
          }
     } else if (is.null(lod)) {
          lod <- rep(NA_real_, n)
     } else {
          lod <- rep_len(as.numeric(lod), n)
     }
     
     pred_use <- pmax(pred, floor)
     
     err_detected <- if (any(detected)) {
          (
               mbmm_safe_log10(y_obs[detected], floor = floor) -
                    mbmm_safe_log10(pred_use[detected], floor = floor)
          )^2
     } else {
          numeric(0)
     }
     
     err_censored <- if (any(censored)) {
          excess <- pmax(
               mbmm_safe_log10(pred_use[censored], floor = floor) -
                    mbmm_safe_log10(lod[censored], floor = floor),
               0
          )
          
          excess^2
     } else {
          numeric(0)
     }
     
     numerator <-
          w_detect * sum(err_detected) +
          w_censored * sum(err_censored)
     
     denominator <-
          w_detect * length(err_detected) +
          w_censored * length(err_censored)
     
     if (denominator <= 0) {
          stop(
               "The discrepancy has zero total weight.",
               call. = FALSE
          )
     }
     
     sqrt(numerator / denominator)
}


#' Fit an ABC rejection-sampling posterior for one subject
#'
#' Draws kinetic parameter vectors from the prior, simulates the mechanistic
#' trajectory, evaluates a censoring-aware discrepancy, and retains the
#' lowest-discrepancy candidates.
#'
#' The number retained can be specified directly using \code{n_keep}, or by
#' an acceptance fraction using \code{accept_frac}. The manuscript analyses
#' generally used a 2.5 percent acceptance fraction before donor balancing.
#'
#' @param id Subject identifier.
#' @param obs Subject observation data.frame containing \code{time} and
#'   \code{y}, and optionally \code{censored} and \code{lod}.
#' @param simulate_fn Function with signature
#'   \code{simulate_fn(theta_row, times)}.
#' @param prior_sample_fn Function with signature
#'   \code{prior_sample_fn(n)} returning a data.frame.
#' @param n_draws Number of prior draws.
#' @param n_keep Optional number of particles to retain.
#' @param accept_frac Optional fraction of valid draws to retain. Used when
#'   \code{n_keep = NULL}.
#' @param distance_fn Optional custom function with signature
#'   \code{distance_fn(theta_row, obs, pred)}.
#' @param seed Optional random seed.
#' @param failure Either \code{"discard"} or \code{"error"}. With
#'   \code{"discard"}, failed simulations receive infinite discrepancy.
#'
#' @return A data.frame of accepted particles with \code{distance} and
#'   \code{id} columns. ABC diagnostic information is stored as attributes.
#'
#' @export
mbmm_abc_fit <- function(
          id,
          obs,
          simulate_fn,
          prior_sample_fn,
          n_draws = 10000,
          n_keep = NULL,
          accept_frac = 0.025,
          distance_fn = NULL,
          seed = NULL,
          failure = c("discard", "error")
) {
     failure <- match.arg(failure)
     
     if (!is.data.frame(obs)) {
          stop(
               "`obs` must be a data.frame.",
               call. = FALSE
          )
     }
     
     missing_obs <- setdiff(
          c("time", "y"),
          names(obs)
     )
     
     if (length(missing_obs) > 0L) {
          stop(
               "`obs` is missing required column(s): ",
               paste(missing_obs, collapse = ", "),
               call. = FALSE
          )
     }
     
     if (nrow(obs) == 0L) {
          stop(
               "`obs` must contain at least one observation.",
               call. = FALSE
          )
     }
     
     if (any(!is.finite(obs$time))) {
          stop(
               "`obs$time` must contain only finite values.",
               call. = FALSE
          )
     }
     
     if (!is.function(simulate_fn)) {
          stop(
               "`simulate_fn` must be a function.",
               call. = FALSE
          )
     }
     
     if (!is.function(prior_sample_fn)) {
          stop(
               "`prior_sample_fn` must be a function.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(n_draws) ||
         length(n_draws) != 1L ||
         !is.finite(n_draws) ||
         n_draws < 1 ||
         n_draws != as.integer(n_draws)) {
          stop(
               "`n_draws` must be one positive integer.",
               call. = FALSE
          )
     }
     
     if (!is.null(n_keep)) {
          if (!is.numeric(n_keep) ||
              length(n_keep) != 1L ||
              !is.finite(n_keep) ||
              n_keep < 1 ||
              n_keep != as.integer(n_keep)) {
               stop(
                    "`n_keep` must be NULL or one positive integer.",
                    call. = FALSE
               )
          }
          
          if (n_keep > n_draws) {
               stop(
                    "`n_keep` cannot exceed `n_draws`.",
                    call. = FALSE
               )
          }
     }
     
     if (!is.numeric(accept_frac) ||
         length(accept_frac) != 1L ||
         !is.finite(accept_frac) ||
         accept_frac <= 0 ||
         accept_frac > 1) {
          stop(
               "`accept_frac` must lie in (0, 1].",
               call. = FALSE
          )
     }
     
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     prior_draws <- prior_sample_fn(n_draws)
     
     if (!is.data.frame(prior_draws)) {
          stop(
               "`prior_sample_fn` must return a data.frame.",
               call. = FALSE
          )
     }
     
     if (nrow(prior_draws) != n_draws) {
          stop(
               "`prior_sample_fn` returned ",
               nrow(prior_draws),
               " rows; expected ",
               n_draws,
               ".",
               call. = FALSE
          )
     }
     
     if (ncol(prior_draws) == 0L) {
          stop(
               "`prior_sample_fn` returned no parameter columns.",
               call. = FALSE
          )
     }
     
     if (is.null(distance_fn)) {
          distance_fn <- function(theta_row, obs, pred) {
               mbmm_default_distance(
                    y_obs = obs$y,
                    pred = pred,
                    censored = if ("censored" %in% names(obs)) {
                         obs$censored
                    } else {
                         NULL
                    },
                    lod = if ("lod" %in% names(obs)) {
                         obs$lod
                    } else {
                         NULL
                    },
                    w_detect = 5,
                    w_censored = 1
               )
          }
     }
     
     if (!is.function(distance_fn)) {
          stop(
               "`distance_fn` must be NULL or a function.",
               call. = FALSE
          )
     }
     
     distances <- rep(Inf, n_draws)
     failure_messages <- character(n_draws)
     
     for (j in seq_len(n_draws)) {
          theta_row <- prior_draws[
               j,
               ,
               drop = FALSE
          ]
          
          result <- tryCatch(
               {
                    pred <- simulate_fn(
                         theta_row,
                         obs$time
                    )
                    
                    if (!is.numeric(pred)) {
                         stop(
                              "Simulator output is not numeric."
                         )
                    }
                    
                    if (length(pred) != nrow(obs)) {
                         stop(
                              "Simulator returned ",
                              length(pred),
                              " predictions; expected ",
                              nrow(obs),
                              "."
                         )
                    }
                    
                    if (any(!is.finite(pred))) {
                         stop(
                              "Simulator returned non-finite predictions."
                         )
                    }
                    
                    distance <- distance_fn(
                         theta_row,
                         obs,
                         pred
                    )
                    
                    if (!is.numeric(distance) ||
                        length(distance) != 1L ||
                        !is.finite(distance) ||
                        distance < 0) {
                         stop(
                              "Distance function must return one finite nonnegative number."
                         )
                    }
                    
                    list(
                         distance = distance,
                         error = NULL
                    )
               },
               error = function(e) {
                    list(
                         distance = Inf,
                         error = conditionMessage(e)
                    )
               }
          )
          
          distances[j] <- result$distance
          
          if (!is.null(result$error)) {
               failure_messages[j] <- result$error
               
               if (identical(failure, "error")) {
                    stop(
                         "ABC simulation failed for subject `",
                         id,
                         "` at prior draw ",
                         j,
                         ": ",
                         result$error,
                         call. = FALSE
                    )
               }
          }
     }
     
     valid <- is.finite(distances)
     n_valid <- sum(valid)
     n_failed <- n_draws - n_valid
     
     if (n_valid == 0L) {
          stop(
               "All ABC simulations failed for subject `",
               id,
               "`.",
               call. = FALSE
          )
     }
     
     if (is.null(n_keep)) {
          n_keep_use <- ceiling(
               accept_frac * n_valid
          )
     } else {
          n_keep_use <- n_keep
     }
     
     n_keep_use <- min(
          n_keep_use,
          n_valid
     )
     
     valid_indices <- which(valid)
     
     ordered_valid <- valid_indices[
          order(
               distances[valid_indices],
               na.last = NA
          )
     ]
     
     accepted_indices <- ordered_valid[
          seq_len(n_keep_use)
     ]
     
     accepted <- prior_draws[
          accepted_indices,
          ,
          drop = FALSE
     ]
     
     accepted$distance <- distances[
          accepted_indices
     ]
     
     accepted$id <- id
     
     rownames(accepted) <- NULL
     
     attr(accepted, "n_draws") <- n_draws
     attr(accepted, "n_valid") <- n_valid
     attr(accepted, "n_failed") <- n_failed
     attr(accepted, "n_keep") <- n_keep_use
     attr(accepted, "acceptance_fraction_valid") <-
          n_keep_use / n_valid
     attr(accepted, "acceptance_fraction_total") <-
          n_keep_use / n_draws
     attr(accepted, "distance_cutoff") <-
          max(accepted$distance)
     attr(accepted, "failed_draws") <-
          which(!valid)
     
     accepted
}


#' Fit ABC posteriors for multiple subjects or donors
#'
#' Applies \code{\link{mbmm_abc_fit}} to each subject in a named observation
#' list and stacks the accepted particles.
#'
#' @param obs_list Named list of subject observation data.frames.
#' @param simulate_fn See \code{\link{mbmm_abc_fit}}.
#' @param prior_sample_fn See \code{\link{mbmm_abc_fit}}.
#' @param n_draws See \code{\link{mbmm_abc_fit}}.
#' @param n_keep See \code{\link{mbmm_abc_fit}}.
#' @param accept_frac See \code{\link{mbmm_abc_fit}}.
#' @param distance_fn See \code{\link{mbmm_abc_fit}}.
#' @param seed Base random seed. Subject \code{i} uses \code{seed + i}.
#' @param failure See \code{\link{mbmm_abc_fit}}.
#' @param verbose Print progress.
#'
#' @return A stacked particle data.frame with donor IDs. Subject-level ABC
#'   diagnostics are stored in the \code{"abc_diagnostics"} attribute.
#'
#' @export
mbmm_abc_bank <- function(
          obs_list,
          simulate_fn,
          prior_sample_fn,
          n_draws = 10000,
          n_keep = NULL,
          accept_frac = 0.025,
          distance_fn = NULL,
          seed = 1,
          failure = c("discard", "error"),
          verbose = TRUE
) {
     failure <- match.arg(failure)
     
     if (!is.list(obs_list) ||
         length(obs_list) == 0L) {
          stop(
               "`obs_list` must be a non-empty named list.",
               call. = FALSE
          )
     }
     
     ids <- names(obs_list)
     
     if (is.null(ids) ||
         any(ids == "") ||
         anyDuplicated(ids)) {
          stop(
               "`obs_list` must have unique, non-empty subject-ID names.",
               call. = FALSE
          )
     }
     
     output <- vector(
          "list",
          length(ids)
     )
     
     diagnostics <- vector(
          "list",
          length(ids)
     )
     
     names(diagnostics) <- ids
     
     for (i in seq_along(ids)) {
          subject_id <- ids[[i]]
          
          if (isTRUE(verbose)) {
               message(
                    sprintf(
                         "[%d/%d] ABC fitting subject %s.",
                         i,
                         length(ids),
                         subject_id
                    )
               )
          }
          
          fit_i <- mbmm_abc_fit(
               id = subject_id,
               obs = obs_list[[subject_id]],
               simulate_fn = simulate_fn,
               prior_sample_fn = prior_sample_fn,
               n_draws = n_draws,
               n_keep = n_keep,
               accept_frac = accept_frac,
               distance_fn = distance_fn,
               seed = if (is.null(seed)) {
                    NULL
               } else {
                    seed + i
               },
               failure = failure
          )
          
          output[[i]] <- fit_i
          
          diagnostics[[subject_id]] <- data.frame(
               id = subject_id,
               n_draws = attr(fit_i, "n_draws"),
               n_valid = attr(fit_i, "n_valid"),
               n_failed = attr(fit_i, "n_failed"),
               n_keep = attr(fit_i, "n_keep"),
               acceptance_fraction_valid =
                    attr(fit_i, "acceptance_fraction_valid"),
               acceptance_fraction_total =
                    attr(fit_i, "acceptance_fraction_total"),
               distance_cutoff =
                    attr(fit_i, "distance_cutoff"),
               stringsAsFactors = FALSE
          )
     }
     
     result <- do.call(
          rbind,
          output
     )
     
     rownames(result) <- NULL
     
     attr(result, "abc_diagnostics") <- do.call(
          rbind,
          diagnostics
     )
     
     result
}