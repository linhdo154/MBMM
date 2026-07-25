#' Plot a fitted MBMM gate curve against one continuous covariate
#'
#' For a binary (K = 2) gate fit with a single continuous covariate (plus
#' optional intercept), draws the fitted gate probability curve
#' \eqn{P(Z = \code{classes[1]} \mid x)} over the observed range of that
#' covariate, with subjects' posterior class probabilities overlaid as
#' points (jittered on the y-axis at their MAP class for readability). This
#' is the plot behind Fig. 4/Fig. 6-style panels in Do, McCarthy & Rogers.
#'
#' @param x an \code{"mbmm_gate"} object from \code{\link{mbmm_fit_gate}}
#' @param focal name of the continuous covariate to plot on the x-axis. If
#'   NULL (default), the first non-intercept term in the gate formula is
#'   used.
#' @param hold_at named list of fixed values for any other covariates in
#'   the gate formula (only needed if there is more than one covariate)
#' @param n_grid number of points used to draw the fitted curve
#' @param ... additional arguments passed to the underlying \code{plot()}
#'   call (e.g. \code{main}, \code{col})
#' @return the fitted gate object, invisibly
#' @export
plot.mbmm_gate <- function(x, focal = NULL, hold_at = list(), n_grid = 200, ...) {
  stopifnot(inherits(x, "mbmm_gate"), x$K == 2)
  terms_used <- all.vars(x$formula)
  if (is.null(focal)) {
    if (length(terms_used) == 0) {
      stop("Gate formula has no covariates to plot against; supply `focal`.")
    }
    focal <- terms_used[1]
  }

  obs_x <- x$X_full[, colnames(x$X_full) != "(Intercept)", drop = FALSE]
  xr <- if (focal %in% colnames(obs_x)) range(obs_x[, focal]) else c(-3, 3)
  grid_x <- seq(xr[1], xr[2], length.out = n_grid)

  nd <- as.list(hold_at)
  nd[[focal]] <- grid_x
  for (v in terms_used) if (is.null(nd[[v]])) nd[[v]] <- 0
  nd <- as.data.frame(nd)
  fitted_p <- mbmm_predict_gate(x, nd)

  graphics::plot(grid_x, fitted_p, type = "l", lwd = 2,
                  xlab = focal, ylab = paste0("P(", x$classes[1], ")"),
                  ylim = c(0, 1), ...)
  graphics::abline(h = 0.5, lty = 3, col = "grey50")

  if (focal %in% colnames(obs_x)) {
    post_col <- paste0("posterior_", x$classes[1])
    graphics::points(obs_x[, focal], x$posterior[[post_col]], pch = 16,
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.5))
  }
  invisible(x)
}

#' Plot the bank log-likelihood-ratio diagnostic for an MBMM fit
#'
#' For a K = 2 MBMM fit, plots each subject's
#' \eqn{\log L_{i,\code{classes[1]}} - \log L_{i,\code{classes[2]}}}
#' (positive values favor the first class's reference bank), grouped and
#' jittered by MAP posterior class -- the diagnostic used throughout Do,
#' McCarthy & Rogers to sanity-check that the reference banks separate
#' subjects before any gate/covariate information is used.
#'
#' @param x an \code{"mbmm_fit"} object from \code{\link{mbmm}}
#' @param ... additional arguments passed to the underlying \code{plot()}
#'   call
#' @return \code{x}, invisibly
#' @export
plot.mbmm_fit <- function(x, ...) {
  gate <- x$gate
  stopifnot(gate$K == 2)
  logL <- x$logL
  lc <- paste0("logL_", gate$classes)
  logLR <- logL[[lc[1]]] - logL[[lc[2]]]
  grp <- factor(gate$posterior$map_class, levels = gate$classes)
  jitter_x <- as.numeric(grp) + stats::runif(length(grp), -0.12, 0.12)

  graphics::plot(jitter_x, logLR, xaxt = "n", pch = 16,
                  col = grDevices::adjustcolor("steelblue", alpha.f = 0.6),
                  xlab = "", xlim = c(0.5, 2.5),
                  ylab = paste0("log L(", gate$classes[1], ") - log L(",
                                 gate$classes[2], ")"), ...)
  graphics::axis(1, at = c(1, 2), labels = gate$classes)
  graphics::abline(h = 0, lty = 2, col = "grey40")
  invisible(x)
}
