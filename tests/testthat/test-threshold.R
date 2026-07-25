make_binary_gate_data <- function(n = 150, seed = 1) {
     set.seed(seed)
     
     x <- stats::rnorm(n)
     
     z <- stats::rbinom(
          n,
          size = 1,
          prob = mbmm_inv_logit(-0.5 + 1.5 * x)
     )
     
     data.frame(
          x = x,
          logL_pos = ifelse(
               z == 1,
               log(0.98),
               log(0.02)
          ),
          logL_neg = ifelse(
               z == 1,
               log(0.02),
               log(0.98)
          )
     )
}

test_that("mbmm_predict_gate validates missing variables", {
     d <- make_binary_gate_data()
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg")
     )
     
     expect_error(
          mbmm_predict_gate(
               fit,
               data.frame(y = 1)
          ),
          "missing required gate variable"
     )
})


test_that("mbmm_gate_threshold matches the analytic binary threshold", {
     d <- make_binary_gate_data()
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg")
     )
     
     beta_0 <- fit$coef$estimate[
          fit$coef$term == "(Intercept)"
     ]
     
     beta_1 <- fit$coef$estimate[
          fit$coef$term == "x"
     ]
     
     expected <- -beta_0 / beta_1
     
     observed <- mbmm_gate_threshold(
          fit,
          focal = "x"
     )
     
     expect_equal(
          observed,
          expected,
          tolerance = 1e-10
     )
})


test_that("mbmm_gate_threshold supports adjusted additive gates", {
     set.seed(21)
     
     n <- 200
     x <- rnorm(n)
     adjustment <- rnorm(n)
     
     latent_probability <- mbmm_inv_logit(
          -0.4 +
               1.2 * x +
               0.5 * adjustment
     )
     
     z <- rbinom(
          n,
          1,
          latent_probability
     )
     
     logL <- data.frame(
          logL_pos = ifelse(z == 1, log(0.95), log(0.05)),
          logL_neg = ifelse(z == 1, log(0.05), log(0.95))
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x + adjustment,
          data = data.frame(
               x = x,
               adjustment = adjustment
          ),
          classes = c("pos", "neg")
     )
     
     beta_0 <- fit$coef$estimate[
          fit$coef$term == "(Intercept)"
     ]
     
     beta_x <- fit$coef$estimate[
          fit$coef$term == "x"
     ]
     
     beta_adjustment <- fit$coef$estimate[
          fit$coef$term == "adjustment"
     ]
     
     adjustment_value <- 0.25
     
     expected <- -(
          beta_0 +
               beta_adjustment * adjustment_value
     ) / beta_x
     
     observed <- mbmm_gate_threshold(
          fit,
          focal = "x",
          hold_at = list(
               adjustment = adjustment_value
          )
     )
     
     expect_equal(
          observed,
          expected,
          tolerance = 1e-8
     )
})


test_that("mbmm_gate_threshold rejects invalid probability levels", {
     d <- make_binary_gate_data()
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg")
     )
     
     expect_error(
          mbmm_gate_threshold(
               fit,
               focal = "x",
               level = 1
          ),
          "strictly between 0 and 1"
     )
})