test_that("logsumexp / logmeanexp are numerically correct and stable", {
     v <- c(1, 2, 3)
     
     expect_equal(
          mbmm_logsumexp(v),
          log(sum(exp(v))),
          tolerance = 1e-12
     )
     
     expect_equal(
          mbmm_logmeanexp(v),
          log(mean(exp(v))),
          tolerance = 1e-12
     )
     
     # Large values should remain finite rather than overflow.
     big <- c(1000, 1000.5, 999)
     
     expect_true(is.finite(mbmm_logsumexp(big)))
     expect_true(is.finite(mbmm_logmeanexp(big)))
     
     # All -Inf values represent zero total likelihood.
     expect_equal(
          mbmm_logsumexp(c(-Inf, -Inf)),
          -Inf
     )
     
     expect_equal(
          mbmm_logmeanexp(c(-Inf, -Inf)),
          -Inf
     )
     
     # Empty input returns -Inf.
     expect_equal(
          mbmm_logsumexp(numeric(0)),
          -Inf
     )
     
     expect_equal(
          mbmm_logmeanexp(numeric(0)),
          -Inf
     )
     
     # A +Inf term dominates the sum and mean.
     expect_equal(
          mbmm_logsumexp(c(0, Inf)),
          Inf
     )
     
     expect_equal(
          mbmm_logmeanexp(c(0, Inf)),
          Inf
     )
})


test_that("logmeanexp retains -Inf particles in the Monte Carlo denominator", {
     # The likelihoods are 1, 0, and 0, so their mean is 1/3.
     v <- c(0, -Inf, -Inf)
     
     expect_equal(
          mbmm_logmeanexp(v),
          log(1 / 3),
          tolerance = 1e-12
     )
     
     # This guards against accidentally dropping -Inf values before dividing
     # by the number of bank particles.
     expect_false(
          isTRUE(all.equal(
               mbmm_logmeanexp(v),
               0
          ))
     )
})


test_that("logsumexp and logmeanexp do not silently remove missing values", {
     expect_error(
          mbmm_logsumexp(c(0, NA_real_)),
          "NA or NaN"
     )
     
     expect_error(
          mbmm_logmeanexp(c(0, NA_real_)),
          "NA or NaN"
     )
     
     expect_error(
          mbmm_logsumexp(c(0, NaN)),
          "NA or NaN"
     )
     
     expect_error(
          mbmm_logmeanexp(c(0, NaN)),
          "NA or NaN"
     )
})


test_that("inv_logit and logit are inverses and correctly bounded", {
     p <- c(0.001, 0.25, 0.5, 0.75, 0.999)
     
     expect_equal(
          mbmm_inv_logit(mbmm_logit(p)),
          p,
          tolerance = 1e-8
     )
     
     inv_values <- mbmm_inv_logit(c(-50, 50))
     
     expect_true(all(inv_values >= 0 & inv_values <= 1))
     expect_equal(mbmm_inv_logit(0), 0.5)
     
     # Boundary probabilities have the mathematically correct infinite logits.
     expect_equal(mbmm_logit(0), -Inf)
     expect_equal(mbmm_logit(1), Inf)
     
     # Invalid probabilities should not be accepted.
     expect_error(
          mbmm_logit(c(-0.1, 0.5)),
          "between 0 and 1"
     )
     
     expect_error(
          mbmm_logit(c(0.5, 1.1)),
          "between 0 and 1"
     )
})


test_that("mbmm_safe_log10 floors values before taking log10", {
     expect_equal(
          mbmm_safe_log10(0, floor = 0.01),
          log10(0.01)
     )
     
     expect_equal(
          mbmm_safe_log10(-5, floor = 0.01),
          log10(0.01)
     )
     
     expect_equal(
          mbmm_safe_log10(100),
          log10(100)
     )
     
     # Missing values remain missing rather than being replaced by the floor.
     expect_true(
          is.na(mbmm_safe_log10(NA_real_))
     )
     
     # The floor must be finite and strictly positive.
     expect_error(
          mbmm_safe_log10(1, floor = 0),
          "strictly positive"
     )
     
     expect_error(
          mbmm_safe_log10(1, floor = -0.1),
          "strictly positive"
     )
})


test_that("mbmm_dnorm_censored matches dnorm and pnorm", {
     y <- c(2.0, 1.0)
     mean <- c(2.1, 1.5)
     sd <- c(0.3, 0.3)
     censored <- c(FALSE, TRUE)
     lod <- c(NA_real_, 1.2)
     
     out <- mbmm_dnorm_censored(
          y = y,
          mean = mean,
          sd = sd,
          censored = censored,
          lod = lod
     )
     
     expect_equal(
          out[1],
          dnorm(
               y[1],
               mean = mean[1],
               sd = sd[1],
               log = TRUE
          )
     )
     
     expect_equal(
          out[2],
          pnorm(
               lod[2],
               mean = mean[2],
               sd = sd[2],
               log.p = TRUE
          )
     )
})


test_that("censored observations do not require a finite y value", {
     out <- mbmm_dnorm_censored(
          y = NA_real_,
          mean = 1.5,
          sd = 0.3,
          censored = TRUE,
          lod = 1.2
     )
     
     expected <- pnorm(
          1.2,
          mean = 1.5,
          sd = 0.3,
          log.p = TRUE
     )
     
     expect_equal(out, expected)
})


test_that("detected observations require finite y values", {
     expect_error(
          mbmm_dnorm_censored(
               y = NA_real_,
               mean = 1.5,
               sd = 0.3,
               censored = FALSE
          ),
          "Detected observations require finite"
     )
})


test_that("censored observations require a finite LOD", {
     expect_error(
          mbmm_dnorm_censored(
               y = NA_real_,
               mean = 1.5,
               sd = 0.3,
               censored = TRUE,
               lod = NULL
          ),
          "lod.*required"
     )
     
     expect_error(
          mbmm_dnorm_censored(
               y = NA_real_,
               mean = 1.5,
               sd = 0.3,
               censored = TRUE,
               lod = NA_real_
          ),
          "finite.*lod"
     )
})


test_that("mbmm_dnorm_censored rejects invalid SD values", {
     expect_error(
          mbmm_dnorm_censored(
               y = 1,
               mean = 1,
               sd = 0,
               censored = FALSE
          ),
          "strictly positive"
     )
     
     expect_error(
          mbmm_dnorm_censored(
               y = 1,
               mean = 1,
               sd = -0.1,
               censored = FALSE
          ),
          "strictly positive"
     )
     
     expect_error(
          mbmm_dnorm_censored(
               y = 1,
               mean = 1,
               sd = NA_real_,
               censored = FALSE
          ),
          "finite"
     )
})


test_that("mbmm_dnorm_censored recycles scalar arguments correctly", {
     y <- c(1, 2, 3)
     
     out <- mbmm_dnorm_censored(
          y = y,
          mean = 2,
          sd = 0.5,
          censored = FALSE
     )
     
     expected <- dnorm(
          y,
          mean = 2,
          sd = 0.5,
          log = TRUE
     )
     
     expect_equal(out, expected)
})


test_that("mbmm_rk4_solve recovers the analytic solution of dy/dt = -k y", {
     f <- function(t, y, parms) {
          -parms$k * y
     }
     
     times <- seq(0, 5, by = 0.5)
     
     sol <- mbmm_rk4_solve(
          f = f,
          y0 = c(y = 10),
          times = times,
          parms = list(k = 0.7)
     )
     
     analytic <- 10 * exp(-0.7 * times)
     
     expect_equal(
          as.numeric(sol[, "y"]),
          analytic,
          tolerance = 1e-4
     )
})