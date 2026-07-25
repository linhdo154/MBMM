#' Reduced manuscript-style simulation helpers
#'
#' These helpers provide a small, fully synthetic simulation that preserves the
#' main structure of the MBMM manuscript simulation:
#'
#' 1. a continuous baseline biomarker;
#' 2. a biomarker-dependent latent productive-infection state;
#' 3. productive trajectories generated from the package's eclipse-phase ODE;
#' 4. protected subjects with no establishment or transient/abortive detections;
#' 5. irregular visit schedules, an assay LOD, and censoring.
#'
#' No participant-level COVID5001 data are included.

# ------------------------------------------------------------------------------
# Internal utilities
# ------------------------------------------------------------------------------
#' @noRd
.mbmm_demo_clip <- function(x, lower, upper) {
     pmin(pmax(x, lower), upper)
}

#' @noRd
.mbmm_demo_calibrate_alpha <- function(z, beta, target_prevalence) {
     if (!is.numeric(z) || length(z) < 2L || any(!is.finite(z))) {
          stop("`z` must contain at least two finite numeric values.", call. = FALSE)
     }
     
     if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta)) {
          stop("`beta` must be one finite numeric value.", call. = FALSE)
     }
     
     if (
          !is.numeric(target_prevalence) ||
          length(target_prevalence) != 1L ||
          !is.finite(target_prevalence) ||
          target_prevalence <= 0 ||
          target_prevalence >= 1
     ) {
          stop(
               "`target_prevalence` must lie strictly between 0 and 1.",
               call. = FALSE
          )
     }
     
     objective <- function(alpha) {
          mean(stats::plogis(alpha + beta * z)) - target_prevalence
     }
     
     stats::uniroot(
          objective,
          interval = c(-30, 30)
     )$root
}


# ------------------------------------------------------------------------------
# Public synthetic visit schedules
# ------------------------------------------------------------------------------

#' Public visit schedules for the reduced MBMM simulation
#'
#' These schedules are fully synthetic. They resemble a longitudinal infection
#' study with occasional missed visits but are not copied from participant-level
#' data.
#'
#' @return A named list of numeric visit-time vectors.
#' @export
mbmm_demo_visit_schedules <- function() {
     list(
          complete = c(0, 2, 4, 7, 10, 14, 21),
          missed_day_4 = c(0, 2, 7, 10, 14, 21),
          missed_day_10 = c(0, 2, 4, 7, 14, 21),
          sparse_early = c(0, 2, 4, 7, 14),
          sparse_late = c(0, 4, 7, 10, 21)
     )
}


# ------------------------------------------------------------------------------
# Fully synthetic biomarker
# ------------------------------------------------------------------------------

#' Simulate a fully synthetic baseline antibody measurement
#'
#' The public demonstration uses a bounded two-component Gaussian mixture on the
#' natural-log scale. This produces a broad, right-skewed antibody distribution
#' on the natural scale without using empirical COVID5001 measurements.
#'
#' @param n Number of subjects.
#' @param high_fraction Probability of belonging to the higher-response component.
#' @param low_mean,high_mean Component means on the natural-log scale.
#' @param low_sd,high_sd Component standard deviations on the natural-log scale.
#' @param clip Lower and upper bounds applied on the natural-log scale.
#' @param calibration_mean Optional fixed mean used to standardize the biomarker.
#' @param calibration_sd Optional fixed standard deviation used to standardize it.
#' @param seed Optional random seed.
#'
#' @return A data frame containing natural-scale, log-scale, and standardized
#'   antibody values.
#' @export
mbmm_simulate_demo_antibody <- function(
          n,
          high_fraction = 0.25,
          low_mean = 1.25,
          high_mean = 2.35,
          low_sd = 0.55,
          high_sd = 0.40,
          clip = c(0.0, 3.3),
          calibration_mean = NULL,
          calibration_sd = NULL,
          seed = NULL
) {
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 2) {
          stop("`n` must be an integer of at least 2.", call. = FALSE)
     }
     n <- as.integer(n)
     
     if (
          !is.numeric(high_fraction) ||
          length(high_fraction) != 1L ||
          !is.finite(high_fraction) ||
          high_fraction < 0 ||
          high_fraction > 1
     ) {
          stop("`high_fraction` must lie between 0 and 1.", call. = FALSE)
     }
     
     numeric_args <- c(
          low_mean = low_mean,
          high_mean = high_mean,
          low_sd = low_sd,
          high_sd = high_sd
     )
     if (any(!is.finite(numeric_args))) {
          stop("Biomarker mixture parameters must be finite.", call. = FALSE)
     }
     if (low_sd <= 0 || high_sd <= 0) {
          stop("Biomarker component standard deviations must be positive.", call. = FALSE)
     }
     
     if (
          !is.numeric(clip) ||
          length(clip) != 2L ||
          any(!is.finite(clip)) ||
          clip[1] >= clip[2]
     ) {
          stop("`clip` must contain increasing finite bounds.", call. = FALSE)
     }
     
     high_component <- stats::rbinom(
          n,
          size = 1L,
          prob = high_fraction
     )
     
     log_antibody <- stats::rnorm(
          n,
          mean = ifelse(high_component == 1L, high_mean, low_mean),
          sd = ifelse(high_component == 1L, high_sd, low_sd)
     )
     
     log_antibody <- .mbmm_demo_clip(
          log_antibody,
          lower = clip[1],
          upper = clip[2]
     )
     
     if (is.null(calibration_mean)) {
          calibration_mean <- mean(log_antibody)
     }
     
     if (is.null(calibration_sd)) {
          calibration_sd <- stats::sd(log_antibody)
     }
     
     if (
          !is.numeric(calibration_mean) ||
          length(calibration_mean) != 1L ||
          !is.finite(calibration_mean)
     ) {
          stop("`calibration_mean` must be one finite number.", call. = FALSE)
     }
     
     if (
          !is.numeric(calibration_sd) ||
          length(calibration_sd) != 1L ||
          !is.finite(calibration_sd) ||
          calibration_sd <= 0
     ) {
          stop("`calibration_sd` must be one positive finite number.", call. = FALSE)
     }
     
     z_antibody <- (log_antibody - calibration_mean) / calibration_sd
     
     data.frame(
          antibody = exp(log_antibody),
          log_antibody = log_antibody,
          z_antibody = z_antibody,
          high_component = high_component,
          calibration_mean = calibration_mean,
          calibration_sd = calibration_sd,
          stringsAsFactors = FALSE
     )
}


# ------------------------------------------------------------------------------
# Transient / abortive patterns
# ------------------------------------------------------------------------------

#' Simulate one transient or abortive detection pattern
#'
#' The four patterns reflect the reduced manuscript simulation:
#' \itemize{
#'   \item \code{noise_like}: no intended signal above the LOD;
#'   \item \code{single_blip}: one intended detectable visit;
#'   \item \code{two_point}: two intended near-LOD detections;
#'   \item \code{fast_rise_drop}: one higher isolated peak with neighboring
#'   visits remaining near background.
#' }
#'
#' @param times Observation times.
#' @param lod Natural-scale assay limit of detection.
#' @param sigma_obs Observation-noise SD on the log10 scale.
#' @param subtype Optional fixed subtype. When \code{NULL}, one is sampled.
#' @param subtype_prob Named subtype probabilities.
#' @param max_log10_above_lod Maximum signal above log10(LOD).
#' @param target_detect_prob Probability used to protect intended detections from
#'   being lost after observation noise.
#' @param seed Optional random seed.
#'
#' @return A data frame containing natural-scale observations, censoring,
#'   latent log10 signal, and subtype.
#' @export
mbmm_simulate_abortive_pattern <- function(
          times,
          lod = 100,
          sigma_obs = 0.5,
          subtype = NULL,
          subtype_prob = c(
               noise_like = 0.20,
               single_blip = 0.35,
               two_point = 0.30,
               fast_rise_drop = 0.15
          ),
          max_log10_above_lod = 1.5,
          target_detect_prob = 0.90,
          seed = NULL
) {
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     times <- sort(unique(as.numeric(times)))
     
     if (length(times) == 0L || any(!is.finite(times))) {
          stop("`times` must contain finite numeric values.", call. = FALSE)
     }
     
     if (
          !is.numeric(lod) ||
          length(lod) != 1L ||
          !is.finite(lod) ||
          lod <= 0
     ) {
          stop("`lod` must be one positive finite number.", call. = FALSE)
     }
     
     if (
          !is.numeric(sigma_obs) ||
          length(sigma_obs) != 1L ||
          !is.finite(sigma_obs) ||
          sigma_obs < 0
     ) {
          stop("`sigma_obs` must be one nonnegative finite number.", call. = FALSE)
     }
     
     if (
          !is.numeric(max_log10_above_lod) ||
          length(max_log10_above_lod) != 1L ||
          !is.finite(max_log10_above_lod) ||
          max_log10_above_lod <= 0
     ) {
          stop("`max_log10_above_lod` must be positive.", call. = FALSE)
     }
     
     if (
          !is.numeric(target_detect_prob) ||
          length(target_detect_prob) != 1L ||
          !is.finite(target_detect_prob) ||
          target_detect_prob <= 0 ||
          target_detect_prob >= 1
     ) {
          stop("`target_detect_prob` must lie strictly between 0 and 1.", call. = FALSE)
     }
     
     allowed <- c(
          "noise_like",
          "single_blip",
          "two_point",
          "fast_rise_drop"
     )
     
     if (is.null(subtype)) {
          if (
               !is.numeric(subtype_prob) ||
               is.null(names(subtype_prob)) ||
               !all(allowed %in% names(subtype_prob)) ||
               any(!is.finite(subtype_prob)) ||
               any(subtype_prob < 0) ||
               sum(subtype_prob) <= 0
          ) {
               stop(
                    "`subtype_prob` must be a valid named probability vector.",
                    call. = FALSE
               )
          }
          
          subtype_prob <- subtype_prob[allowed]
          subtype_prob <- subtype_prob / sum(subtype_prob)
          
          subtype <- sample(
               allowed,
               size = 1L,
               prob = subtype_prob
          )
     } else {
          subtype <- match.arg(subtype, allowed)
     }
     
     lod_log10 <- log10(lod)
     n <- length(times)
     
     background <- stats::rnorm(
          n,
          mean = lod_log10 - 1.1,
          sd = 0.25
     )
     
     log10_true <- background
     
     margin_standard <- stats::qnorm(target_detect_prob) * sigma_obs
     margin_high <- stats::qnorm(0.97) * sigma_obs
     maximum <- lod_log10 + max_log10_above_lod
     
     if (subtype == "single_blip") {
          peak <- sample(seq_len(n), size = 1L)
          
          log10_true[peak] <- stats::runif(
               1L,
               min = lod_log10 + margin_standard,
               max = min(
                    lod_log10 + margin_standard + 0.5,
                    maximum
               )
          )
     }
     
     if (subtype == "two_point") {
          if (n == 1L) {
               log10_true[1L] <- stats::runif(
                    1L,
                    min = lod_log10 + margin_standard,
                    max = min(
                         lod_log10 + margin_standard + 0.5,
                         maximum
                    )
               )
          } else {
               margin <- if (n <= 3L) margin_high else margin_standard
               
               detected <- sample(
                    seq_len(n),
                    size = 2L,
                    replace = FALSE
               )
               
               log10_true[detected[1L]] <- stats::runif(
                    1L,
                    min = lod_log10 + margin + 0.1,
                    max = min(
                         lod_log10 + margin + 0.6,
                         maximum
                    )
               )
               
               log10_true[detected[2L]] <- stats::runif(
                    1L,
                    min = lod_log10 + margin,
                    max = min(
                         lod_log10 + margin + 0.3,
                         maximum
                    )
               )
          }
     }
     
     if (subtype == "fast_rise_drop") {
          peak <- sample(seq_len(n), size = 1L)
          
          log10_true[peak] <- stats::runif(
               1L,
               min = lod_log10 + 0.8,
               max = maximum
          )
     }
     
     log10_observed <- log10_true +
          stats::rnorm(n, mean = 0, sd = sigma_obs)
     
     censored <- log10_observed < lod_log10
     
     observed <- ifelse(
          censored,
          lod,
          10^log10_observed
     )
     
     data.frame(
          time = times,
          y = observed,
          censored = censored,
          lod = lod,
          log10_true = log10_true,
          log10_observed = log10_observed,
          subtype = subtype,
          stringsAsFactors = FALSE
     )
}


# ------------------------------------------------------------------------------
# Complete public manuscript-style simulation
# ------------------------------------------------------------------------------

#' Simulate a reduced manuscript-style MBMM dataset
#'
#' This public simulation preserves the main manuscript structure while replacing
#' private empirical biomarker and visit-pattern inputs with fully synthetic
#' analogues.
#'
#' Productive subjects are generated from
#' \code{\link{mbmm_example_eclipse_simulate}} using parameter draws from
#' \code{\link{mbmm_example_eclipse_prior}} unless a custom productive parameter
#' bank is supplied.
#'
#' @param n_subjects Number of subjects.
#' @param prevalence Target marginal productive-infection prevalence for the
#'   realized synthetic biomarker distribution.
#' @param beta_true True standardized biomarker effect in the logistic gate.
#'   Negative values represent protection.
#' @param p_abortive Probability that a nonproductive subject receives an
#'   abortive/transient pattern rather than a purely near-background trajectory.
#' @param n_bank_donors Number of productive subjects assigned the complete visit
#'   schedule and flagged as bank eligible.
#' @param productive_parameter_bank Optional data frame of productive ODE
#'   parameters. When omitted, a synthetic public bank is generated from
#'   \code{\link{mbmm_example_eclipse_prior}}.
#' @param visit_schedules Named list of public visit schedules.
#' @param sigma_obs Observation-noise SD on the log10 scale.
#' @param lod Natural-scale assay limit of detection.
#' @param productive_min_detect Minimum number of detections requested for a
#'   productive trajectory.
#' @param productive_max_tries Maximum attempts used to obtain that minimum.
#' @param seed Optional random seed.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{subject_data}}{Subject-level covariates and truth.}
#'   \item{\code{truth}}{Alias of \code{subject_data}.}
#'   \item{\code{observations}}{Long-format observations.}
#'   \item{\code{obs_list}}{Named per-subject observation data frames.}
#'   \item{\code{covariates}}{Data frame containing \code{id} and
#'   \code{baseline_antibody}.}
#'   \item{\code{productive_parameter_bank}}{Parameter bank used for productive
#'   trajectory generation.}
#'   \item{\code{settings}}{Generating settings.}
#' }
#' @export
mbmm_simulate_manuscript_demo <- function(
          n_subjects = 50,
          prevalence = 0.30,
          beta_true = -1.5,
          p_abortive = 0.40,
          n_bank_donors = 8L,
          productive_parameter_bank = NULL,
          visit_schedules = mbmm_demo_visit_schedules(),
          sigma_obs = 0.5,
          lod = 100,
          productive_min_detect = 1L,
          productive_max_tries = 50L,
          seed = NULL
) {
     if (!is.null(seed)) {
          set.seed(seed)
     }
     
     if (
          !is.numeric(n_subjects) ||
          length(n_subjects) != 1L ||
          !is.finite(n_subjects) ||
          n_subjects < 10
     ) {
          stop("`n_subjects` must be an integer of at least 10.", call. = FALSE)
     }
     n_subjects <- as.integer(n_subjects)
     
     if (
          !is.numeric(prevalence) ||
          length(prevalence) != 1L ||
          !is.finite(prevalence) ||
          prevalence <= 0 ||
          prevalence >= 1
     ) {
          stop("`prevalence` must lie strictly between 0 and 1.", call. = FALSE)
     }
     
     if (
          !is.numeric(beta_true) ||
          length(beta_true) != 1L ||
          !is.finite(beta_true)
     ) {
          stop("`beta_true` must be one finite number.", call. = FALSE)
     }
     
     if (
          !is.numeric(p_abortive) ||
          length(p_abortive) != 1L ||
          !is.finite(p_abortive) ||
          p_abortive < 0 ||
          p_abortive > 1
     ) {
          stop("`p_abortive` must lie between 0 and 1.", call. = FALSE)
     }
     
     if (
          !is.numeric(n_bank_donors) ||
          length(n_bank_donors) != 1L ||
          !is.finite(n_bank_donors) ||
          n_bank_donors < 0
     ) {
          stop("`n_bank_donors` must be a nonnegative integer.", call. = FALSE)
     }
     n_bank_donors <- as.integer(n_bank_donors)
     
     if (
          !is.numeric(sigma_obs) ||
          length(sigma_obs) != 1L ||
          !is.finite(sigma_obs) ||
          sigma_obs < 0
     ) {
          stop("`sigma_obs` must be one nonnegative finite number.", call. = FALSE)
     }
     
     if (
          !is.numeric(lod) ||
          length(lod) != 1L ||
          !is.finite(lod) ||
          lod <= 0
     ) {
          stop("`lod` must be one positive finite number.", call. = FALSE)
     }
     
     productive_min_detect <- as.integer(productive_min_detect)
     productive_max_tries <- as.integer(productive_max_tries)
     
     if (productive_min_detect < 0L) {
          stop("`productive_min_detect` must be nonnegative.", call. = FALSE)
     }
     
     if (productive_max_tries < 1L) {
          stop("`productive_max_tries` must be at least 1.", call. = FALSE)
     }
     
     if (!is.list(visit_schedules) || length(visit_schedules) == 0L) {
          stop("`visit_schedules` must be a nonempty list.", call. = FALSE)
     }
     
     if (is.null(names(visit_schedules))) {
          names(visit_schedules) <- paste0(
               "schedule_",
               seq_along(visit_schedules)
          )
     }
     
     if (!"complete" %in% names(visit_schedules)) {
          stop(
               "`visit_schedules` must contain a schedule named `complete`.",
               call. = FALSE
          )
     }
     
     visit_schedules <- lapply(
          visit_schedules,
          function(x) sort(unique(as.numeric(x)))
     )
     
     valid_schedule <- vapply(
          visit_schedules,
          function(x) length(x) > 0L && all(is.finite(x)) && all(x >= 0),
          logical(1)
     )
     
     if (!all(valid_schedule)) {
          stop(
               "Every visit schedule must contain finite nonnegative times.",
               call. = FALSE
          )
     }
     
     # Step 1: fully synthetic baseline biomarker.
     antibody <- mbmm_simulate_demo_antibody(
          n = n_subjects
     )
     
     z <- antibody$z_antibody
     
     # Step 2: calibrate and generate the latent productive state.
     alpha_true <- .mbmm_demo_calibrate_alpha(
          z = z,
          beta = beta_true,
          target_prevalence = prevalence
     )
     
     productive_probability <- stats::plogis(
          alpha_true + beta_true * z
     )
     
     productive_true <- stats::rbinom(
          n_subjects,
          size = 1L,
          prob = productive_probability
     )
     
     # Step 3: choose productive bank donors.
     productive_indices <- which(productive_true == 1L)
     
     n_bank_available <- min(
          n_bank_donors,
          length(productive_indices)
     )
     
     bank_indices <- if (n_bank_available > 0L) {
          sample(
               productive_indices,
               size = n_bank_available,
               replace = FALSE
          )
     } else {
          integer(0)
     }
     
     bank_eligible <- seq_len(n_subjects) %in% bank_indices
     
     # Step 4: public productive parameter bank.
     if (is.null(productive_parameter_bank)) {
          productive_parameter_bank <- mbmm_example_eclipse_prior(
               n = 1000L
          )
     }
     
     if (!is.data.frame(productive_parameter_bank)) {
          productive_parameter_bank <- as.data.frame(
               productive_parameter_bank,
               stringsAsFactors = FALSE
          )
     }
     
     required_parameters <- c(
          "beta",
          "delta",
          "p",
          "V0",
          "offset"
     )
     
     missing_parameters <- setdiff(
          required_parameters,
          names(productive_parameter_bank)
     )
     
     if (length(missing_parameters) > 0L) {
          stop(
               "`productive_parameter_bank` is missing: ",
               paste(missing_parameters, collapse = ", "),
               call. = FALSE
          )
     }
     
     if (nrow(productive_parameter_bank) < 1L) {
          stop(
               "`productive_parameter_bank` must contain at least one row.",
               call. = FALSE
          )
     }
     
     parameter_matrix <- as.matrix(
          productive_parameter_bank[, required_parameters, drop = FALSE]
     )
     
     if (any(!is.finite(parameter_matrix))) {
          stop(
               "`productive_parameter_bank` contains non-finite required parameters.",
               call. = FALSE
          )
     }
     
     # Step 5: simulate observations.
     ids <- sprintf(
          "SUBJ%03d",
          seq_len(n_subjects)
     )
     
     subject_rows <- vector(
          "list",
          n_subjects
     )
     
     obs_list <- stats::setNames(
          vector("list", n_subjects),
          ids
     )
     
     schedule_names <- names(visit_schedules)
     
     for (i in seq_len(n_subjects)) {
          id <- ids[i]
          
          if (bank_eligible[i]) {
               schedule_name <- "complete"
          } else {
               schedule_name <- sample(
                    schedule_names,
                    size = 1L
               )
          }
          
          times <- visit_schedules[[schedule_name]]
          
          subtype <- NA_character_
          parameter_row <- NA_integer_
          simulation_attempts <- NA_integer_
          
          if (productive_true[i] == 1L) {
               best <- NULL
               best_n_detect <- -1L
               
               for (attempt in seq_len(productive_max_tries)) {
                    parameter_row_candidate <- sample(
                         seq_len(nrow(productive_parameter_bank)),
                         size = 1L
                    )
                    
                    theta <- productive_parameter_bank[
                         parameter_row_candidate,
                         ,
                         drop = FALSE
                    ]
                    
                    trajectory <- tryCatch(
                         mbmm_example_eclipse_simulate(
                              theta_row = theta,
                              times = times
                         ),
                         error = function(e) NULL
                    )
                    
                    if (
                         is.null(trajectory) ||
                         length(trajectory) != length(times) ||
                         any(!is.finite(trajectory))
                    ) {
                         next
                    }
                    
                    log10_observed <- mbmm_safe_log10(trajectory) +
                         stats::rnorm(
                              length(times),
                              mean = 0,
                              sd = sigma_obs
                         )
                    
                    censored <- log10_observed < log10(lod)
                    
                    candidate <- data.frame(
                         id = id,
                         time = times,
                         y = ifelse(
                              censored,
                              lod,
                              10^log10_observed
                         ),
                         censored = censored,
                         lod = lod,
                         stringsAsFactors = FALSE
                    )
                    
                    n_detect_candidate <- sum(!candidate$censored)
                    
                    if (n_detect_candidate > best_n_detect) {
                         best <- candidate
                         best_n_detect <- n_detect_candidate
                         parameter_row <- parameter_row_candidate
                         simulation_attempts <- attempt
                    }
                    
                    if (n_detect_candidate >= productive_min_detect) {
                         break
                    }
               }
               
               if (is.null(best)) {
                    stop(
                         "Unable to simulate a valid productive trajectory for ",
                         id,
                         ".",
                         call. = FALSE
                    )
               }
               
               obs_i <- best
               true_class <- "productive"
          } else {
               has_abortive_pattern <- stats::rbinom(
                    1L,
                    size = 1L,
                    prob = p_abortive
               ) == 1L
               
               if (has_abortive_pattern) {
                    abortive <- mbmm_simulate_abortive_pattern(
                         times = times,
                         lod = lod,
                         sigma_obs = sigma_obs
                    )
                    
                    subtype <- unique(abortive$subtype)
                    
                    obs_i <- data.frame(
                         id = id,
                         time = abortive$time,
                         y = abortive$y,
                         censored = abortive$censored,
                         lod = abortive$lod,
                         stringsAsFactors = FALSE
                    )
               } else {
                    subtype <- "noise_like"
                    
                    log10_observed <- stats::rnorm(
                         length(times),
                         mean = log10(lod) - 1,
                         sd = sigma_obs
                    )
                    
                    censored <- log10_observed < log10(lod)
                    
                    obs_i <- data.frame(
                         id = id,
                         time = times,
                         y = ifelse(
                              censored,
                              lod,
                              10^log10_observed
                         ),
                         censored = censored,
                         lod = lod,
                         stringsAsFactors = FALSE
                    )
               }
               
               true_class <- "nonestablishment"
          }
          
          obs_list[[id]] <- obs_i
          
          n_detect <- sum(!obs_i$censored)
          
          subject_rows[[i]] <- data.frame(
               id = id,
               antibody = antibody$antibody[i],
               log_antibody = antibody$log_antibody[i],
               baseline_antibody = z[i],
               alpha_true = alpha_true,
               beta_true = beta_true,
               productive_probability = productive_probability[i],
               productive_true = productive_true[i],
               true_class = true_class,
               bank_eligible = bank_eligible[i],
               visit_schedule = schedule_name,
               abortive_subtype = subtype,
               n_detect = n_detect,
               any_detect = n_detect >= 1L,
               two_detect = n_detect >= 2L,
               productive_parameter_row = parameter_row,
               productive_simulation_attempts = simulation_attempts,
               stringsAsFactors = FALSE
          )
     }
     
     subject_data <- do.call(
          rbind,
          subject_rows
     )
     
     observations <- do.call(
          rbind,
          obs_list
     )
     rownames(observations) <- NULL
     
     covariates <- subject_data[
          ,
          c("id", "baseline_antibody"),
          drop = FALSE
     ]
     
     list(
          subject_data = subject_data,
          truth = subject_data,
          observations = observations,
          obs_list = obs_list,
          covariates = covariates,
          productive_parameter_bank = productive_parameter_bank,
          settings = list(
               n_subjects = n_subjects,
               prevalence = prevalence,
               realized_mean_probability = mean(productive_probability),
               realized_productive_fraction = mean(productive_true),
               alpha_true = alpha_true,
               beta_true = beta_true,
               p_abortive = p_abortive,
               n_bank_donors_requested = n_bank_donors,
               n_bank_donors_realized = sum(bank_eligible),
               sigma_obs = sigma_obs,
               lod = lod,
               productive_min_detect = productive_min_detect,
               productive_max_tries = productive_max_tries
          )
     )
}


# ------------------------------------------------------------------------------
# Reduced demonstration grid
# ------------------------------------------------------------------------------

#' Create the default reduced manuscript-style scenario grid
#'
#' This grid is deliberately tiny and is intended for a vignette demonstration,
#' not for estimating statistical power or type-I error.
#'
#' @param replicates Number of demonstration replicates per setting.
#' @param n_subjects Number of subjects per generated dataset.
#' @param prevalence Target marginal productive-infection prevalence.
#' @param beta_grid True biomarker effects.
#' @param p_abortive_grid Abortive-pattern probabilities.
#' @param base_seed First seed in the grid.
#'
#' @return A data frame.
#' @export
mbmm_demo_scenario_grid <- function(
          replicates = 2L,
          n_subjects = 50L,
          prevalence = 0.30,
          beta_grid = c(0, -1.5),
          p_abortive_grid = c(0, 0.40),
          base_seed = 20260724L
) {
     replicates <- as.integer(replicates)
     n_subjects <- as.integer(n_subjects)
     base_seed <- as.integer(base_seed)
     
     if (replicates < 1L) {
          stop("`replicates` must be at least 1.", call. = FALSE)
     }
     
     if (n_subjects < 10L) {
          stop("`n_subjects` must be at least 10.", call. = FALSE)
     }
     
     if (
          !is.numeric(prevalence) ||
          length(prevalence) != 1L ||
          !is.finite(prevalence) ||
          prevalence <= 0 ||
          prevalence >= 1
     ) {
          stop("`prevalence` must lie strictly between 0 and 1.", call. = FALSE)
     }
     
     if (
          !is.numeric(beta_grid) ||
          length(beta_grid) < 1L ||
          any(!is.finite(beta_grid))
     ) {
          stop("`beta_grid` must contain finite numeric values.", call. = FALSE)
     }
     
     if (
          !is.numeric(p_abortive_grid) ||
          length(p_abortive_grid) < 1L ||
          any(!is.finite(p_abortive_grid)) ||
          any(p_abortive_grid < 0) ||
          any(p_abortive_grid > 1)
     ) {
          stop(
               "`p_abortive_grid` must contain values between 0 and 1.",
               call. = FALSE
          )
     }
     
     grid <- expand.grid(
          n_subjects = n_subjects,
          prevalence = prevalence,
          beta_true = beta_grid,
          p_abortive = p_abortive_grid,
          replicate = seq_len(replicates),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
     )
     
     grid$scenario <- paste0(
          "beta=",
          grid$beta_true,
          "; abortive=",
          grid$p_abortive
     )
     
     grid$seed <- base_seed +
          seq_len(nrow(grid)) -
          1L
     
     grid
}