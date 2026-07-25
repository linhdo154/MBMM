#' Pool per-donor ABC particles into a mechanistic reference bank
#'
#' Constructs a balanced mechanistic reference bank by retaining the same
#' number of accepted ABC particles from every eligible donor.
#'
#' When \code{n_per_donor = NULL}, the smallest available donor-specific
#' particle count is used, so that every donor contributes equally.
#'
#' If \code{n_per_donor} is supplied and some donor has fewer particles than
#' requested, the default behavior is to stop rather than silently construct
#' an unbalanced bank. Set \code{allow_fewer = TRUE} only when unequal
#' donor-specific particle counts are intentional.
#'
#' @param particles A data.frame containing accepted ABC particles, with one
#'   row per particle.
#' @param n_per_donor Number of particles to retain per donor. If NULL, uses
#'   the minimum available count across eligible donors.
#' @param donor_col Name of the donor ID column.
#' @param eligible_donors Optional vector of donor IDs to include. This can be
#'   used for application-specific eligibility rules, such as requiring a
#'   minimum number of detected measurements.
#' @param allow_fewer If FALSE, stop when a donor has fewer particles than
#'   \code{n_per_donor}. If TRUE, retain all available particles for that donor.
#' @param seed Optional random seed.
#'
#' @return A pooled particle bank with a balanced donor contribution when
#'   \code{allow_fewer = FALSE}.
#'
#' @export
mbmm_pool_bank <- function(
          particles,
          n_per_donor = NULL,
          donor_col = "id",
          eligible_donors = NULL,
          allow_fewer = FALSE,
          seed = NULL
) {
     if (!is.data.frame(particles)) {
          stop(
               "`particles` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (!donor_col %in% names(particles)) {
          stop(
               "`particles` does not contain donor column `",
               donor_col,
               "`.",
               call. = FALSE
          )
     }
     
     if (nrow(particles) == 0L) {
          stop(
               "`particles` is empty.",
               call. = FALSE
          )
     }
     
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     donor_values <- as.character(
          particles[[donor_col]]
     )
     
     if (anyNA(donor_values) || any(donor_values == "")) {
          stop(
               "Donor IDs must be non-missing and non-empty.",
               call. = FALSE
          )
     }
     
     if (!is.null(eligible_donors)) {
          eligible_donors <- as.character(
               eligible_donors
          )
          
          particles <- particles[
               donor_values %in% eligible_donors,
               ,
               drop = FALSE
          ]
          
          donor_values <- as.character(
               particles[[donor_col]]
          )
          
          if (nrow(particles) == 0L) {
               stop(
                    "No particles remain after applying `eligible_donors`.",
                    call. = FALSE
               )
          }
     }
     
     donors <- unique(donor_values)
     
     counts <- vapply(
          donors,
          function(donor_id) {
               sum(donor_values == donor_id)
          },
          integer(1)
     )
     
     names(counts) <- donors
     
     if (is.null(n_per_donor)) {
          n_per_donor <- min(counts)
     }
     
     if (!is.numeric(n_per_donor) ||
         length(n_per_donor) != 1L ||
         !is.finite(n_per_donor) ||
         n_per_donor < 1 ||
         n_per_donor != as.integer(n_per_donor)) {
          stop(
               "`n_per_donor` must be NULL or one positive integer.",
               call. = FALSE
          )
     }
     
     insufficient <- names(counts)[
          counts < n_per_donor
     ]
     
     if (
          length(insufficient) > 0L &&
          !isTRUE(allow_fewer)
     ) {
          stop(
               "The following donor(s) have fewer than ",
               n_per_donor,
               " accepted particles: ",
               paste(
                    paste0(
                         insufficient,
                         " (",
                         counts[insufficient],
                         ")"
                    ),
                    collapse = ", "
               ),
               ". Use a smaller `n_per_donor`, set `n_per_donor = NULL`, or set ",
               "`allow_fewer = TRUE` if an unbalanced bank is intentional.",
               call. = FALSE
          )
     }
     
     selected <- lapply(
          donors,
          function(donor_id) {
               donor_rows <- particles[
                    as.character(particles[[donor_col]]) == donor_id,
                    ,
                    drop = FALSE
               ]
               
               n_keep <- if (isTRUE(allow_fewer)) {
                    min(n_per_donor, nrow(donor_rows))
               } else {
                    n_per_donor
               }
               
               sampled_rows <- sample.int(
                    n = nrow(donor_rows),
                    size = n_keep,
                    replace = FALSE
               )
               
               donor_rows[
                    sampled_rows,
                    ,
                    drop = FALSE
               ]
          }
     )
     
     result <- do.call(
          rbind,
          selected
     )
     
     rownames(result) <- NULL
     
     attr(result, "particles_per_donor") <- table(
          result[[donor_col]]
     )
     
     result
}


#' Join target-specific particle banks into complete multitarget bundles
#'
#' Constructs joint mechanistic bundles for multitarget applications. For each
#' donor and each bundle, one accepted particle is sampled independently from
#' every target-specific donor pool.
#'
#' Only donors represented in every target are eligible. Every returned
#' donor-bundle combination is guaranteed to contain exactly one row for every
#' target.
#'
#' @param particles_by_target Named list of particle data.frames, one per
#'   target.
#' @param n_bundles Number of joint bundles to construct per donor.
#' @param donor_col Name of the donor ID column.
#' @param replace Whether to sample particles with replacement within each
#'   target-specific donor pool. The SDY2354 construction typically uses
#'   replacement because the requested number of joint bundles may exceed the
#'   number of retained particles for an individual target.
#' @param seed Optional random seed.
#'
#' @return A long data.frame with one row per donor, bundle, and target.
#'
#' @export
mbmm_join_targets <- function(
          particles_by_target,
          n_bundles = 80,
          donor_col = "id",
          replace = TRUE,
          seed = NULL
) {
     if (!is.list(particles_by_target) ||
         length(particles_by_target) == 0L) {
          stop(
               "`particles_by_target` must be a non-empty named list.",
               call. = FALSE
          )
     }
     
     targets <- names(particles_by_target)
     
     if (is.null(targets) ||
         any(targets == "") ||
         anyDuplicated(targets)) {
          stop(
               "`particles_by_target` must have unique, non-empty target names.",
               call. = FALSE
          )
     }
     
     if (!is.numeric(n_bundles) ||
         length(n_bundles) != 1L ||
         !is.finite(n_bundles) ||
         n_bundles < 1 ||
         n_bundles != as.integer(n_bundles)) {
          stop(
               "`n_bundles` must be one positive integer.",
               call. = FALSE
          )
     }
     
     if (!is.logical(replace) ||
         length(replace) != 1L ||
         is.na(replace)) {
          stop(
               "`replace` must be TRUE or FALSE.",
               call. = FALSE
          )
     }
     
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     for (target_name in targets) {
          target_particles <- particles_by_target[[target_name]]
          
          if (!is.data.frame(target_particles)) {
               stop(
                    "Particle object for target `",
                    target_name,
                    "` is not a data.frame.",
                    call. = FALSE
               )
          }
          
          if (!donor_col %in% names(target_particles)) {
               stop(
                    "Particle object for target `",
                    target_name,
                    "` does not contain donor column `",
                    donor_col,
                    "`.",
                    call. = FALSE
               )
          }
          
          donor_ids <- as.character(
               target_particles[[donor_col]]
          )
          
          if (anyNA(donor_ids) || any(donor_ids == "")) {
               stop(
                    "Target `",
                    target_name,
                    "` contains missing or empty donor IDs.",
                    call. = FALSE
               )
          }
     }
     
     donors_by_target <- lapply(
          particles_by_target,
          function(target_particles) {
               unique(
                    as.character(
                         target_particles[[donor_col]]
                    )
               )
          }
     )
     
     common_donors <- Reduce(
          intersect,
          donors_by_target
     )
     
     if (length(common_donors) == 0L) {
          stop(
               "No donors are represented across every requested target.",
               call. = FALSE
          )
     }
     
     if (!replace) {
          minimum_counts <- vapply(
               common_donors,
               function(donor_id) {
                    min(
                         vapply(
                              particles_by_target,
                              function(target_particles) {
                                   sum(
                                        as.character(
                                             target_particles[[donor_col]]
                                        ) == donor_id
                                   )
                              },
                              integer(1)
                         )
                    )
               },
               integer(1)
          )
          
          insufficient <- names(minimum_counts)[
               minimum_counts < n_bundles
          ]
          
          if (length(insufficient) > 0L) {
               stop(
                    "Cannot construct ",
                    n_bundles,
                    " bundles without replacement for donor(s): ",
                    paste(insufficient, collapse = ", "),
                    ". Use fewer bundles or set `replace = TRUE`.",
                    call. = FALSE
               )
          }
     }
     
     output <- vector(
          "list",
          length(common_donors) *
               n_bundles *
               length(targets)
     )
     
     output_index <- 1L
     
     for (donor_id in common_donors) {
          donor_target_pools <- lapply(
               particles_by_target,
               function(target_particles) {
                    target_particles[
                         as.character(
                              target_particles[[donor_col]]
                         ) == donor_id,
                         ,
                         drop = FALSE
                    ]
               }
          )
          
          for (bundle_number in seq_len(n_bundles)) {
               for (target_name in targets) {
                    donor_target_pool <- donor_target_pools[[target_name]]
                    
                    if (nrow(donor_target_pool) == 0L) {
                         stop(
                              "Unexpected empty donor-target pool for donor `",
                              donor_id,
                              "` and target `",
                              target_name,
                              "`.",
                              call. = FALSE
                         )
                    }
                    
                    sampled_index <- sample.int(
                         n = nrow(donor_target_pool),
                         size = 1L,
                         replace = replace
                    )
                    
                    selected_row <- donor_target_pool[
                         sampled_index,
                         ,
                         drop = FALSE
                    ]
                    
                    selected_row[[donor_col]] <- donor_id
                    selected_row$target <- target_name
                    
                    # bundle_id is unique within donor. The donor + bundle pair is the
                    # complete mechanistic draw identity.
                    selected_row$bundle_id <- bundle_number
                    
                    output[[output_index]] <- selected_row
                    output_index <- output_index + 1L
               }
          }
     }
     
     result <- do.call(
          rbind,
          output
     )
     
     rownames(result) <- NULL
     
     # Confirm every donor-bundle contains exactly one row for every target.
     bundle_key <- interaction(
          result[[donor_col]],
          result$bundle_id,
          drop = TRUE
     )
     
     target_count <- tapply(
          as.character(result$target),
          bundle_key,
          function(x) {
               length(unique(x))
          }
     )
     
     if (any(target_count != length(targets))) {
          stop(
               "Incomplete multitarget bundles were created.",
               call. = FALSE
          )
     }
     
     duplicate_key <- data.frame(
          donor = result[[donor_col]],
          bundle_id = result$bundle_id,
          target = result$target,
          stringsAsFactors = FALSE
     )
     
     if (anyDuplicated(duplicate_key)) {
          stop(
               "Duplicate donor-bundle-target rows were created.",
               call. = FALSE
          )
     }
     
     result
}


#' Predict mechanistic bank trajectories at required observation times
#'
#' Applies a user-supplied mechanistic simulator to every particle or
#' multitarget bundle row in a reference bank.
#'
#' @param bank Reference-bank data.frame.
#' @param simulate_fn Function with signature
#'   \code{simulate_fn(theta_row, times)} returning one numeric prediction per
#'   requested time.
#' @param times Numeric vector of evaluation times for single-target banks.
#' @param times_by_target Named list of time vectors for multitarget banks.
#' @param donor_col Name of the donor ID column.
#' @param bundle_col Name of the bundle ID column.
#' @param failure Either \code{"error"} or \code{"drop"}. With
#'   \code{"error"}, any failed or invalid simulation stops. With
#'   \code{"drop"}, failed bank rows are omitted and reported in attributes.
#'
#' @return Long-format bank predictions containing donor ID, bundle ID,
#'   optional target, time, and prediction.
#'
#' @export
mbmm_predict_bank <- function(
          bank,
          simulate_fn,
          times = NULL,
          times_by_target = NULL,
          donor_col = "id",
          bundle_col = "bundle_id",
          failure = c("error", "drop")
) {
     failure <- match.arg(failure)
     
     if (!is.data.frame(bank)) {
          stop(
               "`bank` must be a data.frame.",
               call. = FALSE
          )
     }
     
     if (nrow(bank) == 0L) {
          stop(
               "`bank` is empty.",
               call. = FALSE
          )
     }
     
     if (!is.function(simulate_fn)) {
          stop(
               "`simulate_fn` must be a function.",
               call. = FALSE
          )
     }
     
     if (!donor_col %in% names(bank)) {
          stop(
               "`bank` does not contain donor column `",
               donor_col,
               "`.",
               call. = FALSE
          )
     }
     
     if (!(bundle_col %in% names(bank))) {
          bank[[bundle_col]] <- seq_len(
               nrow(bank)
          )
     }
     
     donor_values <- as.character(
          bank[[donor_col]]
     )
     
     if (anyNA(donor_values) || any(donor_values == "")) {
          stop(
               "Bank donor IDs must be non-missing and non-empty.",
               call. = FALSE
          )
     }
     
     if (anyNA(bank[[bundle_col]])) {
          stop(
               "Bank bundle IDs must be non-missing.",
               call. = FALSE
          )
     }
     
     has_target <- "target" %in% names(bank)
     
     if (has_target) {
          if (is.null(times_by_target) ||
              !is.list(times_by_target) ||
              is.null(names(times_by_target))) {
               stop(
                    "A named `times_by_target` list is required for a multitarget bank.",
                    call. = FALSE
               )
          }
          
          bank_targets <- unique(
               as.character(bank$target)
          )
          
          missing_times <- setdiff(
               bank_targets,
               names(times_by_target)
          )
          
          if (length(missing_times) > 0L) {
               stop(
                    "No evaluation times were supplied for target(s): ",
                    paste(missing_times, collapse = ", "),
                    call. = FALSE
               )
          }
     } else {
          if (is.null(times)) {
               stop(
                    "`times` must be supplied for a single-target bank.",
                    call. = FALSE
               )
          }
     }
     
     rows <- vector(
          "list",
          nrow(bank)
     )
     
     failed_rows <- integer(0)
     
     for (i in seq_len(nrow(bank))) {
          theta_row <- bank[
               i,
               ,
               drop = FALSE
          ]
          
          if (has_target) {
               target_name <- as.character(
                    theta_row$target
               )
               
               evaluation_times <- times_by_target[[target_name]]
          } else {
               target_name <- NULL
               evaluation_times <- times
          }
          
          if (!is.numeric(evaluation_times) ||
              length(evaluation_times) == 0L ||
              any(!is.finite(evaluation_times))) {
               stop(
                    "Evaluation times must be a non-empty finite numeric vector.",
                    call. = FALSE
               )
          }
          
          prediction <- tryCatch(
               simulate_fn(
                    theta_row,
                    evaluation_times
               ),
               error = function(e) {
                    structure(
                         NA_real_,
                         simulation_error = conditionMessage(e)
                    )
               }
          )
          
          valid_prediction <- is.numeric(prediction) &&
               length(prediction) == length(evaluation_times) &&
               all(is.finite(prediction))
          
          if (!valid_prediction) {
               failed_rows <- c(
                    failed_rows,
                    i
               )
               
               if (identical(failure, "error")) {
                    error_message <- attr(
                         prediction,
                         "simulation_error"
                    )
                    
                    stop(
                         "Simulation failed for bank row ",
                         i,
                         ", donor `",
                         theta_row[[donor_col]],
                         "`, bundle `",
                         theta_row[[bundle_col]],
                         "`",
                         if (has_target) {
                              paste0(
                                   ", target `",
                                   target_name,
                                   "`"
                              )
                         } else {
                              ""
                         },
                         if (!is.null(error_message)) {
                              paste0(
                                   ": ",
                                   error_message
                              )
                         } else {
                              paste0(
                                   ". Expected ",
                                   length(evaluation_times),
                                   " finite predictions."
                              )
                         },
                         call. = FALSE
                    )
               }
               
               rows[[i]] <- NULL
               next
          }
          
          prediction_row <- data.frame(
               donor = as.character(
                    theta_row[[donor_col]]
               ),
               bundle_id = theta_row[[bundle_col]],
               time = as.numeric(
                    evaluation_times
               ),
               pred = as.numeric(
                    prediction
               ),
               stringsAsFactors = FALSE
          )
          
          if (has_target) {
               prediction_row$target <- target_name
          }
          
          rows[[i]] <- prediction_row
     }
     
     rows <- rows[
          !vapply(
               rows,
               is.null,
               logical(1)
          )
     ]
     
     if (length(rows) == 0L) {
          stop(
               "No valid bank predictions were produced.",
               call. = FALSE
          )
     }
     
     result <- do.call(
          rbind,
          rows
     )
     
     names(result)[
          names(result) == "donor"
     ] <- donor_col
     
     rownames(result) <- NULL
     
     prediction_key <- data.frame(
          donor = result[[donor_col]],
          bundle = result$bundle_id,
          target = if ("target" %in% names(result)) {
               result$target
          } else {
               ""
          },
          time = result$time,
          stringsAsFactors = FALSE
     )
     
     if (anyDuplicated(prediction_key)) {
          stop(
               "Duplicate donor-bundle-target-time predictions were produced.",
               call. = FALSE
          )
     }
     
     attr(result, "failed_bank_rows") <- failed_rows
     attr(result, "n_failed_bank_rows") <- length(
          failed_rows
     )
     
     result
}