#' Simulate a cohort with a latent kinetic state
#'
#' Generates the setting the method is built for: a baseline marker shifts the
#' probability of a latent state, participants in that state follow a
#' mechanistic trajectory, and some participants *not* in it still produce a
#' brief detectable event. A rule that counts any detection as the state
#' therefore mislabels part of the cohort, which is what attenuates a
#' covariate effect and costs power.
#'
#' Because the state is generated rather than inferred, the simulation is the
#' one setting in which recovery can be checked directly. It is also the only
#' worked example the package can ship: the analyses in the accompanying
#' manuscript use access-controlled cohorts.
#'
#' @param n Number of participants.
#' @param beta_ab Covariate effect on the log-odds of the state. Negative
#'   means a higher marker lowers the probability.
#' @param prevalence Target marginal prevalence of the state. The gate
#'   intercept is solved to hit it against a large reference draw, so the
#'   realised prevalence does not drift with `n`.
#' @param p_abortive Fraction of participants *not* in the state who still
#'   produce a transient detectable event. This is the contamination the
#'   method is meant to absorb; at zero, a detection rule is already faithful
#'   and there is little for a latent-state model to recover.
#' @param times Visit schedule.
#' @param lod Assay limit of detection.
#' @param sigma_obs Observation noise, \eqn{\log_{10}} scale.
#' @param n_pool Size of the pre-solved trajectory population drawn from
#'   [eclipse_design()].
#' @param seed Optional seed.
#'
#' @return A list with `obs` (long-format measurements, with `censored` and
#'   `lod`), `truth` (one row per participant: the latent `z`, the covariate
#'   `x`, and the generating parameters), and the `alpha` that was solved for.
#' @seealso [eclipse_design()], [trajectory_loglik()]
#' @export
#' @examples
#' if (requireNamespace("deSolve", quietly = TRUE)) {
#'   sim <- simulate_mbmm(n = 12, beta_ab = -1, n_pool = 8, seed = 1)
#'   table(sim$truth$z)
#' }
simulate_mbmm <- function(n, beta_ab = -1, prevalence = 0.2,
                          p_abortive = 0.4,
                          times = c(0, 2, 4, 7, 10, 14, 21),
                          lod = 100, sigma_obs = 0.5, n_pool = 200L,
                          seed = NULL) {
  if (!requireNamespace("deSolve", quietly = TRUE))
    stop("Package 'deSolve' is required for simulate_mbmm().")
  if (!is.null(seed)) set.seed(seed)

  ## gate: solve the intercept against a large reference draw so the realised
  ## prevalence does not drift with n
  zref  <- stats::rnorm(2e5)
  alpha <- stats::uniroot(
    function(a) mean(stats::plogis(a + beta_ab * zref)) - prevalence,
    c(-50, 50))$root

  x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, stats::plogis(alpha + beta_ab * x))

  ## one pre-solved population of establishing trajectories
  pool <- eclipse_design(n_pool, seed = NULL)
  grid <- seq(0, max(times) + 1, by = 0.25)
  curves <- vapply(seq_len(n_pool), function(m) {
    eclipse_trajectory(grid, pool$beta[m], pool$delta[m], pool$p[m],
                       pool$V0[m], kappa = pool$kappa[m],
                       c_clear = pool$c_clear[m], T0 = pool$T0[m],
                       grid = grid)
  }, numeric(length(grid)))

  pick <- integer(n); pick[z == 1] <- sample.int(n_pool, sum(z == 1), TRUE)
  abortive <- rep(FALSE, n)
  abortive[z == 0] <- stats::runif(sum(z == 0)) < p_abortive

  rows <- lapply(seq_len(n), function(i) {
    if (z[i] == 1) {
      mu <- stats::approx(grid, curves[, pick[i]], xout = times, rule = 2)$y
    } else if (abortive[i]) {
      ## a brief low blip early, then nothing: the records a detection rule
      ## must call infection and a trajectory likelihood need not
      peak <- 10^stats::runif(1, log10(lod), log10(lod) + 1.2)
      mu <- peak * exp(-1.2 * pmax(times - times[1], 0))
    } else {
      mu <- rep(lod / 10, length(times))
    }
    v <- 10^(log10(pmax(mu, 1e-6)) + stats::rnorm(length(times), 0, sigma_obs))
    cen <- v < lod
    data.frame(id = i, time = times, value = ifelse(cen, NA_real_, v),
               censored = cen, lod = lod)
  })

  list(obs = do.call(rbind, rows),
       truth = data.frame(id = seq_len(n), z = z, x = x,
                          abortive = abortive, pool_index = pick),
       alpha = alpha)
}
