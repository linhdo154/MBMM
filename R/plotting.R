utils::globalVariables(c(
  ".mbmm_x", ".mbmm_y", ".mbmm_weight", "state", "responsibility",
  "fitted", "probability"
))

#' Plot a weighted mechanistic trajectory bank
#'
#' Shows where an estimated or known kinetic population places its mass on two
#' candidate coordinates. Point area and opacity both encode the candidate
#' weight. The function returns a ggplot object, so labels, themes and scales
#' can be changed with ordinary ggplot2 layers.
#'
#' @param candidates A data frame with one row per candidate.
#' @param weights Non-negative candidate weights. They are normalised to sum
#'   to one for display.
#' @param x,y Names of the candidate columns to plot.
#' @param log10_x,log10_y Transform the corresponding coordinate to base-10
#'   logarithms before plotting. All transformed values must be positive.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' candidates <- expand.grid(peak = 10^(2:4), rate = c(0.1, 0.3, 0.5))
#' weights <- seq_len(nrow(candidates))
#' plot_bank(candidates, weights, "peak", "rate", log10_x = TRUE)
plot_bank <- function(candidates, weights, x, y,
                      log10_x = FALSE, log10_y = FALSE, title = NULL) {
  if (!is.data.frame(candidates)) stop("`candidates` must be a data frame.")
  if (!all(c(x, y) %in% names(candidates)))
    stop("`x` and `y` must name columns in `candidates`.")
  if (length(weights) != nrow(candidates))
    stop("`weights` must have one value per candidate.")
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0)
    stop("`weights` must be finite, non-negative, and have positive sum.")

  xx <- candidates[[x]]
  yy <- candidates[[y]]
  if (!is.numeric(xx) || !is.numeric(yy))
    stop("The plotted candidate columns must be numeric.")
  if (log10_x) {
    if (any(xx <= 0)) stop("`x` values must be positive for log10 display.")
    xx <- log10(xx)
  }
  if (log10_y) {
    if (any(yy <= 0)) stop("`y` values must be positive for log10 display.")
    yy <- log10(yy)
  }

  dat <- data.frame(.mbmm_x = xx, .mbmm_y = yy,
                    .mbmm_weight = weights / sum(weights))
  xlab <- if (log10_x) paste0("log10(", x, ")") else x
  ylab <- if (log10_y) paste0("log10(", y, ")") else y

  ggplot2::ggplot(dat, ggplot2::aes(x = .mbmm_x, y = .mbmm_y)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .mbmm_weight, alpha = .mbmm_weight),
      shape = 21, fill = "#2C7FB8", colour = "#174A69", stroke = 0.35
    ) +
    ggplot2::scale_size_continuous(range = c(1, 10), name = "Weight") +
    ggplot2::scale_alpha_continuous(range = c(0.2, 0.95), guide = "none") +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot posterior responsibilities against simulated truth
#'
#' Intended for simulations, where the true latent state is known. In real
#' data there is no truth column to supply: the responsibility is the inferred
#' probability of the state rather than a calibrated diagnostic label.
#'
#' @param fit An object returned by [fit_gate()].
#' @param truth Binary vector of known simulated states (0 or 1).
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' n <- 30; x <- rnorm(n); X <- cbind(1, x)
#' truth <- rbinom(n, 1, plogis(-0.5 - x))
#' fit <- fit_gate(rnorm(n, truth), rnorm(n, 1 - truth), X)
#' plot_responsibility(fit, truth)
plot_responsibility <- function(fit, truth, title = NULL) {
  if (is.null(fit$w)) stop("`fit` must contain posterior responsibilities.")
  if (length(truth) != length(fit$w) || any(!truth %in% c(0, 1)))
    stop("`truth` must be a binary vector matching the fitted cohort.")

  dat <- data.frame(
    state = factor(truth, levels = 0:1, labels = c("State 0", "State 1")),
    responsibility = fit$w
  )
  ggplot2::ggplot(dat, ggplot2::aes(x = state, y = responsibility,
                                    fill = state)) +
    ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 1.4) +
    ggplot2::scale_fill_manual(values = c("#BDBDBD", "#FDB863"),
                               guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Posterior responsibility", title = title) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Plot a fitted covariate gate
#'
#' Displays the fitted prior state probability against one covariate. For a
#' simulation, optional true states and generating probabilities make recovery
#' visible. For real data, omit both and the fitted gate alone is shown.
#'
#' @param fit An object returned by [fit_gate()].
#' @param x Numeric covariate used in the gate.
#' @param truth Optional binary vector of known simulated latent states.
#' @param true_probability Optional vector of generating state probabilities.
#' @param xlab Horizontal-axis label.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' set.seed(1)
#' n <- 40; x <- rnorm(n); prob <- plogis(-0.5 - x)
#' truth <- rbinom(n, 1, prob); X <- cbind(1, x)
#' fit <- fit_gate(rnorm(n, truth), rnorm(n, 1 - truth), X)
#' plot_gate(fit, x, truth, prob)
plot_gate <- function(fit, x, truth = NULL, true_probability = NULL,
                      xlab = "Covariate", title = NULL) {
  if (is.null(fit$pi)) stop("`fit` must contain fitted gate probabilities.")
  n <- length(fit$pi)
  if (!is.numeric(x) || length(x) != n || any(!is.finite(x)))
    stop("`x` must be a finite numeric vector matching the fitted cohort.")
  if (!is.null(truth) &&
      (length(truth) != n || any(!truth %in% c(0, 1))))
    stop("`truth` must be binary and match the fitted cohort.")
  if (!is.null(true_probability) &&
      (length(true_probability) != n || any(!is.finite(true_probability)) ||
       any(true_probability < 0 | true_probability > 1)))
    stop("`true_probability` must lie in [0, 1] and match the cohort.")

  dat <- data.frame(x = x, fitted = fit$pi)
  dat <- dat[order(dat$x), , drop = FALSE]
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x))

  if (!is.null(truth)) {
    points <- data.frame(x = x, truth = truth)
    p <- p + ggplot2::geom_jitter(
      data = points, ggplot2::aes(x = x, y = truth),
      width = 0, height = 0.035, alpha = 0.25, size = 1.5,
      inherit.aes = FALSE
    )
  }
  if (!is.null(true_probability)) {
    truth_line <- data.frame(x = x, probability = true_probability)
    truth_line <- truth_line[order(truth_line$x), , drop = FALSE]
    p <- p + ggplot2::geom_line(
      data = truth_line,
      ggplot2::aes(x = x, y = probability, colour = "Truth"),
      linewidth = 1, inherit.aes = FALSE
    )
  }

  p +
    ggplot2::geom_line(ggplot2::aes(y = fitted, colour = "MBMM estimate"),
                       linewidth = 1, linetype = 2) +
    ggplot2::scale_colour_manual(
      values = c("Truth" = "#1B9E77", "MBMM estimate" = "#D95F02"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = xlab, y = "State probability", title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}
