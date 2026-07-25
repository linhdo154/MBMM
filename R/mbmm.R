#' Fit a bank-based MBMM analysis
#'
#' Convenience wrapper that first scores every subject against class-specific
#' mechanistic reference banks and then fits the covariate-linked MBMM gate.
#'
#' This wrapper is appropriate when every latent class is represented by a
#' mechanistic trajectory bank, such as a boost-responsive versus weak/stable
#' phenotype analysis.
#'
#' For productive-infection analyses in which the non-establishment class uses
#' an analytic likelihood rather than a second trajectory bank, use
#' \code{\link{mbmm_productive_infection}}.
#'
#' For binary models, gate coefficients describe the log-odds of the first
#' class in \code{classes} relative to the second class.
#'
#' @param obs_list Named list of subject observation data.frames. Names must
#'   contain unique subject IDs.
#' @param bank_pred_list Named list of class-specific bank-prediction tables.
#'   Names define the latent-class labels.
#' @param sigma Measurement-error SD passed to
#'   \code{\link{mbmm_class_loglik}}.
#' @param formula One-sided gate-model formula, such as
#'   \code{~ baseline_antibody} or \code{~ arm}.
#' @param data Covariate data.frame containing one row per subject.
#' @param classes Optional character vector fixing the latent-class order.
#'   Defaults to \code{names(bank_pred_list)}.
#' @param ridge Nonnegative ridge penalty used when
#'   \code{penalty = "ridge"}.
#' @param penalty Either \code{"ridge"} or \code{"none"}.
#' @param n_perm Number of optional gate permutations.
#' @param seed Optional permutation-test seed.
#' @param id_col Name of the subject-ID column in \code{data}.
#' @param donor_col Name of the donor-ID column in the bank tables.
#' @param bundle_col Name of the bundle-ID column in the bank tables.
#' @param log_scale Whether observations, predictions, and LOD values are
#'   already on the log10 likelihood scale.
#' @param donor_weighting Either \code{"equal"} or \code{"particle"}.
#'   Equal donor weighting is the manuscript-oriented default.
#' @param loo Logical scalar or class-named logical vector controlling
#'   leave-one-subject/donor-out scoring.
#' @param verbose Print scoring progress.
#' @param maxit Maximum gate optimizer iterations.
#'
#' @return An object of class \code{"mbmm_fit"} containing the scored
#'   class-conditional likelihoods and fitted gate.
#'
#' @examples
#' obs_list <- list(
#'   S1 = data.frame(time = 0:1, y = c(1.0, 1.1)),
#'   S2 = data.frame(time = 0:1, y = c(1.2, 0.9)),
#'   S3 = data.frame(time = 0:1, y = c(1.1, 1.0)),
#'   S4 = data.frame(time = 0:1, y = c(4.8, 5.1)),
#'   S5 = data.frame(time = 0:1, y = c(5.2, 4.9)),
#'   S6 = data.frame(time = 0:1, y = c(5.0, 5.3))
#' )
#'
#' bank_low <- data.frame(
#'   id = "donorA",
#'   bundle_id = 1,
#'   time = 0:1,
#'   pred = c(1, 1)
#' )
#'
#' bank_high <- data.frame(
#'   id = "donorB",
#'   bundle_id = 1,
#'   time = 0:1,
#'   pred = c(5, 5)
#' )
#'
#' covariates <- data.frame(
#'   id = paste0("S", 1:6),
#'   x = c(-1.2, -0.8, -0.4, 0.4, 0.8, 1.2)
#' )
#'
#' fit <- mbmm(
#'   obs_list = obs_list,
#'   bank_pred_list = list(
#'     low = bank_low,
#'     high = bank_high
#'   ),
#'   sigma = 0.3,
#'   formula = ~x,
#'   data = covariates,
#'   classes = c("low", "high"),
#'   id_col = "id",
#'   loo = FALSE,
#'   verbose = FALSE
#' )
#'
#' fit
#'
#' @export
mbmm <- function(
          obs_list,
          bank_pred_list,
          sigma,
          formula,
          data,
          classes = NULL,
          ridge = 1e-4,
          penalty = "ridge",
          n_perm = 0,
          seed = NULL,
          id_col = "id",
          donor_col = "id",
          bundle_col = "bundle_id",
          log_scale = TRUE,
          donor_weighting = c("equal", "particle"),
          loo = TRUE,
          verbose = TRUE,
          maxit = 4000
) {
     donor_weighting <- match.arg(donor_weighting)
     
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
     
     if (!is.list(obs_list) || length(obs_list) == 0L) {
          stop(
               "`obs_list` must be a non-empty named list.",
               call. = FALSE
          )
     }
     
     subject_ids <- names(obs_list)
     
     if (
          is.null(subject_ids) ||
          any(subject_ids == "") ||
          anyDuplicated(subject_ids)
     ) {
          stop(
               "`obs_list` must have unique, non-empty subject-ID names.",
               call. = FALSE
          )
     }
     
     if (!is.list(bank_pred_list) || length(bank_pred_list) == 0L) {
          stop(
               "`bank_pred_list` must be a non-empty named list.",
               call. = FALSE
          )
     }
     
     bank_classes <- names(bank_pred_list)
     
     if (
          is.null(bank_classes) ||
          any(bank_classes == "") ||
          anyDuplicated(bank_classes)
     ) {
          stop(
               "`bank_pred_list` must have unique, non-empty class names.",
               call. = FALSE
          )
     }
     
     if (is.null(classes)) {
          classes <- bank_classes
     }
     
     if (
          any(classes == "") ||
          anyDuplicated(classes) ||
          !setequal(classes, bank_classes)
     ) {
          stop(
               "`classes` must contain exactly the names of `bank_pred_list`.",
               call. = FALSE
          )
     }
     
     bank_pred_list <- bank_pred_list[classes]
     
     data_used <- .mbmm_align_subject_data(
          data = data,
          subject_ids = subject_ids,
          id_col = id_col
     )
     
     logL <- mbmm_score_subjects(
          obs_list = obs_list,
          bank_pred_list = bank_pred_list,
          sigma = sigma,
          donor_col = donor_col,
          bundle_col = bundle_col,
          log_scale = log_scale,
          donor_weighting = donor_weighting,
          loo = loo,
          verbose = verbose
     )
     
     logL_match <- match(
          subject_ids,
          as.character(logL$id)
     )
     
     if (anyNA(logL_match)) {
          stop(
               "Internal error: scored likelihoods are missing one or more subjects.",
               call. = FALSE
          )
     }
     
     logL <- logL[
          logL_match,
          ,
          drop = FALSE
     ]
     
     rownames(logL) <- NULL
     
     if (!identical(
          as.character(logL$id),
          as.character(data_used[[id_col]])
     )) {
          stop(
               "Internal subject alignment failed between `logL` and `data`.",
               call. = FALSE
          )
     }
     
     gate <- mbmm_fit_gate(
          logL = logL,
          formula = formula,
          data = data_used,
          classes = classes,
          ridge = ridge,
          penalty = penalty,
          n_perm = n_perm,
          seed = seed,
          id_col = id_col,
          maxit = maxit
     )
     
     structure(
          list(
               logL = logL,
               gate = gate,
               obs_list = obs_list,
               bank_pred_list = bank_pred_list,
               data = data_used,
               classes = classes,
               scoring = list(
                    type = "bank-versus-bank",
                    donor_col = donor_col,
                    bundle_col = bundle_col,
                    log_scale = log_scale,
                    donor_weighting = donor_weighting,
                    loo = loo
               )
          ),
          class = "mbmm_fit"
     )
}


#' Fit a productive-infection MBMM analysis
#'
#' Fits an MBMM analysis in which the productive-infection class is represented
#' by a mechanistic trajectory bank and the non-establishment class is
#' represented by an analytic observation model below the assay LOD.
#'
#' This structure is appropriate for infection-establishment applications in
#' which productive infection produces structured longitudinal kinetics,
#' whereas non-establishment or transient/abortive detection is represented
#' by a low-level analytic likelihood rather than a second ODE trajectory bank.
#'
#' The expected measurement under non-establishment is defined as a specified
#' fraction of the observation-specific LOD. A value of \code{0.5}
#' corresponds to an expected natural-scale measurement of LOD divided by two.
#'
#' For binary models, gate coefficients describe the log-odds of the
#' \code{productive_class} relative to the
#' \code{nonestablishment_class}.
#'
#' @param obs_list Named list of subject longitudinal observations.
#' @param productive_bank Predicted productive-infection reference bank.
#' @param sigma Measurement-error SD on the Gaussian likelihood scale.
#' @param formula One-sided gate-model formula.
#' @param data Subject-level gate covariate data.frame.
#' @param lod_fraction Expected non-establishment mean as a fraction of LOD.
#' @param productive_class Label for the productive-infection class.
#' @param nonestablishment_class Label for the non-establishment class.
#' @param ridge Nonnegative ridge stabilization parameter.
#' @param penalty Either \code{"ridge"} or \code{"none"}.
#' @param n_perm Number of optional gate permutations.
#' @param seed Optional permutation-test seed.
#' @param id_col Subject-ID column in \code{data}.
#' @param donor_col Donor-ID column in \code{productive_bank}.
#' @param bundle_col Bundle-ID column in \code{productive_bank}.
#' @param log_scale Whether observations, predictions, and LOD values are
#'   already on the log10 likelihood scale.
#' @param donor_weighting Either \code{"equal"} or \code{"particle"}.
#' @param loo_productive Whether to exclude the scored subject from the
#'   productive bank when that subject is also a bank donor.
#' @param verbose Print scoring progress.
#' @param maxit Maximum gate optimizer iterations.
#'
#' @return An object of class \code{"mbmm_fit"}.
#'
#' @examples
#' \dontrun{
#' # See vignette("mbmm-intro") for a complete productive-infection example.
#' }
#'
#' @export
mbmm_productive_infection <- function(
          obs_list,
          productive_bank,
          sigma,
          formula,
          data,
          lod_fraction = 0.5,
          productive_class = "productive",
          nonestablishment_class = "nonestablishment",
          ridge = 1e-4,
          penalty = "none",
          n_perm = 0,
          seed = NULL,
          id_col = "id",
          donor_col = "id",
          bundle_col = "bundle_id",
          log_scale = FALSE,
          donor_weighting = c("equal", "particle"),
          loo_productive = FALSE,
          verbose = TRUE,
          maxit = 4000
) {
     donor_weighting <- match.arg(donor_weighting)
     
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
     
     if (!is.list(obs_list) || length(obs_list) == 0L) {
          stop(
               "`obs_list` must be a non-empty named list.",
               call. = FALSE
          )
     }
     
     subject_ids <- names(obs_list)
     
     if (
          is.null(subject_ids) ||
          any(subject_ids == "") ||
          anyDuplicated(subject_ids)
     ) {
          stop(
               "`obs_list` must have unique, non-empty subject-ID names.",
               call. = FALSE
          )
     }
     
     if (!is.data.frame(productive_bank)) {
          stop(
               "`productive_bank` must be a data.frame.",
               call. = FALSE
          )
     }
     
     required_bank_columns <- c(
          donor_col,
          bundle_col,
          "time",
          "pred"
     )
     
     missing_bank_columns <- setdiff(
          required_bank_columns,
          names(productive_bank)
     )
     
     if (length(missing_bank_columns) > 0L) {
          stop(
               "`productive_bank` is missing required column(s): ",
               paste(missing_bank_columns, collapse = ", "),
               call. = FALSE
          )
     }
     
     class_labels <- c(
          productive_class,
          nonestablishment_class
     )
     
     if (
          any(class_labels == "") ||
          anyDuplicated(class_labels)
     ) {
          stop(
               "Productive and non-establishment class labels must be unique and non-empty.",
               call. = FALSE
          )
     }
     
     data_used <- .mbmm_align_subject_data(
          data = data,
          subject_ids = subject_ids,
          id_col = id_col
     )
     
     logL_productive <- numeric(
          length(subject_ids)
     )
     
     logL_nonestablishment <- numeric(
          length(subject_ids)
     )
     
     for (i in seq_along(subject_ids)) {
          subject_id <- subject_ids[[i]]
          obs_i <- obs_list[[subject_id]]
          
          exclude_donor <- NULL
          
          if (
               isTRUE(loo_productive) &&
               subject_id %in%
               as.character(unique(productive_bank[[donor_col]]))
          ) {
               exclude_donor <- subject_id
          }
          
          logL_productive[i] <- mbmm_class_loglik(
               obs = obs_i,
               bank_pred = productive_bank,
               sigma = sigma,
               exclude_donor = exclude_donor,
               donor_col = donor_col,
               bundle_col = bundle_col,
               log_scale = log_scale,
               donor_weighting = donor_weighting
          )
          
          logL_nonestablishment[i] <- mbmm_noise_floor_loglik(
               obs = obs_i,
               sigma = sigma,
               lod_fraction = lod_fraction,
               log_scale = log_scale
          )
          
          if (isTRUE(verbose) && i %% 10L == 0L) {
               message(
                    sprintf(
                         "Scored %d of %d subjects.",
                         i,
                         length(subject_ids)
                    )
               )
          }
     }
     
     logL <- data.frame(
          id = subject_ids,
          stringsAsFactors = FALSE
     )
     
     logL[[paste0("logL_", productive_class)]] <-
          logL_productive
     
     logL[[paste0("logL_", nonestablishment_class)]] <-
          logL_nonestablishment
     
     gate <- mbmm_fit_gate(
          logL = logL,
          formula = formula,
          data = data_used,
          classes = class_labels,
          ridge = ridge,
          penalty = penalty,
          n_perm = n_perm,
          seed = seed,
          id_col = id_col,
          maxit = maxit
     )
     
     structure(
          list(
               logL = logL,
               gate = gate,
               obs_list = obs_list,
               bank_pred_list = list(
                    productive = productive_bank
               ),
               data = data_used,
               classes = class_labels,
               scoring = list(
                    type = "productive-bank-versus-analytic-nonestablishment",
                    donor_col = donor_col,
                    bundle_col = bundle_col,
                    log_scale = log_scale,
                    donor_weighting = donor_weighting,
                    loo_productive = loo_productive,
                    lod_fraction = lod_fraction
               )
          ),
          class = "mbmm_fit"
     )
}


#' Internal: align subject-level covariate data with observation-list IDs
#'
#' @keywords internal
.mbmm_align_subject_data <- function(
          data,
          subject_ids,
          id_col
) {
     if (!is.data.frame(data)) {
          stop(
               "`data` must be a data.frame.",
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
     
     data_ids <- as.character(
          data[[id_col]]
     )
     
     if (anyNA(data_ids) || any(data_ids == "")) {
          stop(
               "Subject IDs in `data` must be non-missing and non-empty.",
               call. = FALSE
          )
     }
     
     if (anyDuplicated(data_ids)) {
          duplicated_ids <- unique(
               data_ids[duplicated(data_ids)]
          )
          
          stop(
               "`data` contains duplicate subject ID(s): ",
               paste(duplicated_ids, collapse = ", "),
               call. = FALSE
          )
     }
     
     missing_data_ids <- setdiff(
          subject_ids,
          data_ids
     )
     
     if (length(missing_data_ids) > 0L) {
          stop(
               "The following subject ID(s) from `obs_list` are missing from `data`: ",
               paste(missing_data_ids, collapse = ", "),
               call. = FALSE
          )
     }
     
     extra_data_ids <- setdiff(
          data_ids,
          subject_ids
     )
     
     if (length(extra_data_ids) > 0L) {
          message(
               "Ignoring ",
               length(extra_data_ids),
               " row(s) in `data` that are not present in `obs_list`."
          )
     }
     
     matched_rows <- match(
          subject_ids,
          data_ids
     )
     
     data_used <- data[
          matched_rows,
          ,
          drop = FALSE
     ]
     
     rownames(data_used) <- NULL
     
     data_used
}


#' Print an mbmm_fit object
#'
#' @param x An \code{"mbmm_fit"} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.mbmm_fit <- function(
          x,
          ...
) {
     cat("MBMM fit\n")
     
     cat(
          "  classes:",
          paste(x$classes, collapse = " vs "),
          "\n"
     )
     
     if (!is.null(x$scoring$type)) {
          cat(
               "  scoring:",
               x$scoring$type,
               "\n"
          )
     }
     
     cat(
          "  n subjects scored:",
          nrow(x$logL),
          "\n"
     )
     
     print(x$gate)
     
     invisible(x)
}


#' Summarize an mbmm_fit object
#'
#' @param object An \code{"mbmm_fit"} object.
#' @param ... Unused.
#' @return \code{object}, invisibly.
#' @export
summary.mbmm_fit <- function(
          object,
          ...
) {
     summary(object$gate)
     
     invisible(object)
}