#' Example mechanistic trajectory model: eclipse-phase viral kinetics
#'
#' A target-cell-limited within-host viral kinetic model with an eclipse
#' phase (Eq. 11-14 of Do, McCarthy & Rogers / Baccam et al. 2006), provided
#' as a ready-to-use example of the \code{simulate_fn} ingredient MBMM needs
#' (see \code{\link{mbmm_abc_fit}}). Integration uses \pkg{deSolve}
#' (\code{lsoda}), because this system is numerically stiff enough that a
#' naive fixed-step integrator (like \code{\link{mbmm_rk4_solve}}) can
#' become unstable at the eclipse/clearance rates used here; \pkg{deSolve}
#' is a package dependency for exactly this reason. Fixed constants
#' (\code{k}, \code{c}) and prior locations below are illustrative, not
#' pathogen-specific -- re-derive them from the literature for a real
#' analysis (as Sec. 4.2/App. A.1 of the manuscript do from UKC19 human
#' challenge data).
#'
#' Free parameters (columns expected in \code{theta_row}): \code{beta}
#' (infection rate), \code{delta} (clearance rate of infected cells),
#' \code{p} (viral production rate), \code{V0} (initial viral burden), and
#' \code{offset} (days between the model's internal time origin and the
#' cohort's observed day-0, allowing unobserved infection timing to vary
#' across subjects -- see Sec. 4.1 of the manuscript).
#'
#' @param theta_row a one-row data.frame/list with elements \code{beta},
#'   \code{delta}, \code{p}, \code{V0}, \code{offset}
#' @param times observed times (e.g. "day post enrollment"), any
#'   non-negative numeric vector
#' @param T0 fixed initial target-cell count (default 1e7)
#' @param k fixed eclipse-to-productive transition rate, 1/day (default 4)
#' @param c fixed free-virus clearance rate, 1/day (default 10)
#' @return numeric vector of predicted viral load (natural scale, e.g.
#'   copies/mL) at each element of \code{times}
#' @export
mbmm_example_eclipse_simulate <- function(theta_row, times, T0 = 1e7,
                                           k = 4, c = 10) {
  beta <- as.numeric(theta_row$beta)
  delta <- as.numeric(theta_row$delta)
  p <- as.numeric(theta_row$p)
  V0 <- as.numeric(theta_row$V0)
  offset <- as.numeric(theta_row$offset)

  internal_t <- pmax(times - offset, 0)
  grid <- sort(unique(c(0, internal_t)))

  y0 <- c(T = T0, I1 = 0, I2 = 0, V = V0)
  deriv <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
      dT <- -beta * T * V
      dI1 <- beta * T * V - k * I1
      dI2 <- k * I1 - delta * I2
      dV <- p * I2 - c * V
      list(c(dT, dI1, dI2, dV))
    })
  }

  sol <- deSolve::ode(y = y0, times = grid, func = deriv,
                       parms = c(beta = beta, k = k, delta = delta, p = p, c = c),
                       method = "lsoda")
  idx <- match(internal_t, sol[, "time"])
  pmax(sol[idx, "V"], 1e-3)
}

#' Example ABC prior sampler for the eclipse-phase model
#'
#' Draws candidate parameters for \code{\link{mbmm_example_eclipse_simulate}}
#' from independent priors calibrated to give a growth/peak/decline
#' trajectory over roughly 3 weeks (basic reproduction number
#' \eqn{R_0 = \beta T_0 p / (\delta c)} centered in the 5-15 range at
#' \code{T0 = 1e7}, \code{c = 10}). Intended as a template: real analyses
#' should replace these with priors informed by the literature or an
#' external challenge-study posterior, as in Sec. 4.2/App. A.1 of the
#' manuscript.
#'
#' @param n number of draws
#' @return a data.frame with columns \code{beta}, \code{delta}, \code{p},
#'   \code{V0}, \code{offset}
#' @export
mbmm_example_eclipse_prior <- function(n) {
  data.frame(
    beta = 10^stats::rnorm(n, mean = -6.0, sd = 0.3),
    delta = pmax(stats::rnorm(n, mean = 1.0, sd = 0.25), 0.2),
    p = 10^stats::rnorm(n, mean = 1.1, sd = 0.3),
    V0 = 10^stats::runif(n, 0, 1.5),
    offset = stats::runif(n, 0, 5)
  )
}

#' Simulate a synthetic correlate-of-protection cohort
#'
#' Generates a small synthetic dataset in the spirit of the simulation
#' study in Sec. 3.2 of Do, McCarthy & Rogers: each subject has a baseline
#' (standardized) antibody value; a latent productive-infection state
#' \code{Z} is drawn from an antibody-dependent logistic gate; subjects with
#' \code{Z = 1} get a viral-load trajectory from the eclipse-phase model
#' (\code{\link{mbmm_example_eclipse_simulate}}) under subject-specific true
#' kinetic parameters, observed with log-normal noise and a limit of
#' detection; subjects with \code{Z = 0} get a flat near-LOD trajectory with
#' an occasional single-point transient/abortive "blip" (mimicking a
#' non-productive PCR-detectable event). A subset of the \code{Z = 1}
#' subjects is given a denser visit schedule and flagged as bank-eligible,
#' mirroring the donor-richness threshold used for real reference banks
#' (Sec. 4.3).
#'
#' This function only *generates data*; it does not call any MBMM fitting
#' functions. See \code{vignette("mbmm-intro")} or
#' \code{system.file("examples/01_simulated_correlate_of_protection.R", package = "mbmm")}
#' for the full analysis pipeline built on top of it.
#'
#' @param n_subjects total number of subjects to simulate
#' @param prevalence marginal probability that a "typical" (zAb = 0) subject
#'   is truly infected (Z = 1); combined with \code{beta_true} to set the
#'   gate intercept
#' @param beta_true true log-odds effect of standardized baseline antibody
#'   on infection probability (negative = protective)
#' @param n_bank_donors number of Z = 1 subjects given a dense visit
#'   schedule and flagged \code{bank_eligible = TRUE}
#' @param sigma_obs log10-scale observation noise SD
#' @param lod limit of detection (natural scale, e.g. copies/mL)
#' @param seed RNG seed
#' @return a list with elements \code{obs_list} (named list of per-subject
#'   data.frames with columns \code{time}, \code{y}, \code{censored},
#'   \code{lod}), \code{covariates} (data.frame with \code{id} and
#'   \code{baseline_antibody}), and \code{truth} (data.frame with the true
#'   generating values, for checking recovery -- never used by the fitting
#'   functions themselves)
#' @export
mbmm_simulate_cop_demo <- function(n_subjects = 90, prevalence = 0.35,
                                    beta_true = -1.2, n_bank_donors = 18,
                                    sigma_obs = 0.35, lod = 3, seed = 1) {
  set.seed(seed)
  ids <- sprintf("SUBJ%03d", seq_len(n_subjects))
  zAb <- stats::rnorm(n_subjects)
  alpha0 <- mbmm_logit(prevalence)
  pZ <- mbmm_inv_logit(alpha0 + beta_true * zAb)
  Z <- stats::rbinom(n_subjects, 1, pZ)

  bank_idx <- sample(which(Z == 1), min(n_bank_donors, sum(Z == 1)))
  bank_eligible <- seq_len(n_subjects) %in% bank_idx

  # All subjects are scheduled from the SAME fixed protocol visit grid (as in
  # real trials, and as required by the exact-time join `mbmm_class_loglik`
  # uses between observed days and bank-predicted days); non-bank subjects
  # simply miss some scheduled visits rather than visiting on arbitrary days.
  visit_grid <- c(0, 2, 4, 6, 9, 12, 16, 21)

  obs_list <- stats::setNames(vector("list", n_subjects), ids)
  truth_rows <- vector("list", n_subjects)

  for (i in seq_len(n_subjects)) {
    if (bank_eligible[i]) {
      visit_days <- visit_grid
    } else if (Z[i] == 1) {
      n_visits <- sample(3:5, 1)
      visit_days <- sort(sample(visit_grid, n_visits))
    } else {
      n_visits <- sample(3:6, 1)
      visit_days <- sort(sample(visit_grid, n_visits))
    }

    if (Z[i] == 1) {
      theta_true <- list(
        beta = 10^stats::rnorm(1, -6.0, 0.2),
        delta = pmax(stats::rnorm(1, 1.0, 0.15), 0.3),
        p = 10^stats::rnorm(1, 1.1, 0.2),
        V0 = 10^stats::runif(1, 0.3, 1.2),
        offset = stats::runif(1, 1, 4)
      )
      mu <- mbmm_example_eclipse_simulate(theta_true, visit_days)
      y_true <- 10^(mbmm_safe_log10(mu) + stats::rnorm(length(mu), 0, sigma_obs))
    } else {
      theta_true <- list(beta = NA, delta = NA, p = NA, V0 = NA, offset = NA)
      y_true <- rep(lod * 0.5, length(visit_days))
      if (stats::runif(1) < 0.25) {
        blip_idx <- sample(seq_along(visit_days), 1)
        y_true[blip_idx] <- lod * stats::runif(1, 1.5, 6)
      }
      y_true <- pmax(y_true * 10^stats::rnorm(length(y_true), 0, 0.2), 0.1)
    }

    censored <- y_true < lod
    y_obs <- ifelse(censored, lod * stats::runif(length(y_true), 0.3, 0.99), y_true)

    obs_list[[ids[i]]] <- data.frame(
      time = visit_days, y = y_obs, censored = censored, lod = lod
    )
    truth_rows[[i]] <- data.frame(
      id = ids[i], Z_true = Z[i], bank_eligible = bank_eligible[i],
      beta_true = theta_true$beta, delta_true = theta_true$delta,
      p_true = theta_true$p, V0_true = theta_true$V0, offset_true = theta_true$offset
    )
  }

  list(
    obs_list = obs_list,
    covariates = data.frame(id = ids, baseline_antibody = zAb),
    truth = do.call(rbind, truth_rows),
    gate_truth = list(alpha = alpha0, beta = beta_true)
  )
}
