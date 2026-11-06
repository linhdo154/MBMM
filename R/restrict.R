#' Restrict a bank to a region of parameter space (Route II)
#'
#' Defines a latent class as a prespecified region of kinetic-parameter space
#' and recovers its kinetic distribution by restricting an estimated marginal
#' population to that region and renormalising. Because the class prevalence
#' cancels, the population used to build the bank need not consist only of
#' class members, nor share the analysis population's class prevalence.
#'
#' The region must be fixed before any outcome is observed, and must carry
#' population mass bounded away from zero: restriction divides by that mass,
#' so a thin region amplifies any error in the marginal.
#'
#' @param bank An `mbmm_bank` from [fit_bank()], or a numeric weight vector.
#' @param in_region Logical vector, one entry per candidate, `TRUE` for
#'   candidates inside the region.
#'
#' @return A list with the renormalised `weights`, the `mass` the region
#'   carried in the marginal, and `in_region`.
#' @references Result D of the accompanying manuscript.
#' @export
#' @examples
#' set.seed(1)
#' b <- fit_bank(matrix(rnorm(20 * 6), 20, 6))
#' restrict_bank(b, in_region = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))$mass
restrict_bank <- function(bank, in_region) {
  w <- if (inherits(bank, "mbmm_bank")) bank$weights else bank
  stopifnot(length(w) == length(in_region), is.logical(in_region))
  mass <- sum(w[in_region])
  if (mass <= 0) stop("The region carries no population mass.")
  if (mass < 0.01)
    warning(sprintf("Region mass is %.3f; restriction amplifies error in the marginal.", mass))
  out <- numeric(length(w)); out[in_region] <- w[in_region] / mass
  list(weights = out, mass = mass, in_region = in_region)
}
