make_binary_gate_data <- function(
          n = 300,
          alpha = -0.5,
          beta = 1.5,
          seed = 1
) {
     set.seed(seed)
     
     x <- rnorm(n)
     pi <- mbmm_inv_logit(alpha + beta * x)
     z <- rbinom(n, 1, pi)
     
     # Near-perfectly informative class likelihoods, so the gate should
     # recover something close to the true covariate effect.
     logL_pos <- ifelse(
          z == 1,
          log(0.98),
          log(0.02)
     )
     
     logL_neg <- ifelse(
          z == 1,
          log(0.02),
          log(0.98)
     )
     
     data.frame(
          id = seq_len(n),
          x = x,
          logL_pos = logL_pos,
          logL_neg = logL_neg
     )
}


test_that("mbmm_fit_gate recovers the sign and rough magnitude of a true effect", {
     d <- make_binary_gate_data()
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          id_col = "id",
          penalty = "ridge",
          ridge = 1e-4
     )
     
     expect_s3_class(fit, "mbmm_gate")
     expect_equal(fit$K, 2)
     expect_equal(fit$penalty, "ridge")
     
     beta_hat <- fit$coef$estimate[
          fit$coef$term == "x"
     ]
     
     expect_gt(beta_hat, 0)
     expect_equal(
          beta_hat,
          1.5,
          tolerance = 0.5
     )
     
     expect_true(fit$p_value < 0.01)
     expect_equal(fit$convergence_full, 0)
     expect_equal(fit$convergence_null, 0)
})


test_that("mbmm_fit_gate LRT p-value is large under a null effect", {
     set.seed(3)
     
     n <- 200
     x <- rnorm(n)
     z <- rbinom(n, 1, 0.4)
     
     logL <- data.frame(
          logL_pos = ifelse(
               z == 1,
               log(0.9),
               log(0.1)
          ),
          logL_neg = ifelse(
               z == 1,
               log(0.1),
               log(0.9)
          )
     )
     
     d <- data.frame(x = x)
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "none"
     )
     
     expect_true(fit$p_value > 0.05)
})


test_that("reported LRT uses the unpenalized marginal likelihood", {
     d <- make_binary_gate_data(
          n = 200,
          seed = 11
     )
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          ridge = 0.1
     )
     
     expect_equal(
          fit$LRT,
          2 * (
               fit$logLik_full -
                    fit$logLik_null
          ),
          tolerance = 1e-10
     )
     
     expect_equal(
          fit$ll_full,
          fit$logLik_full
     )
     
     expect_equal(
          fit$ll_null,
          fit$logLik_null
     )
     
     expect_lte(
          fit$objective_full_penalized,
          fit$logLik_full + 1e-12
     )
     
     expect_lte(
          fit$objective_null_penalized,
          fit$logLik_null + 1e-12
     )
})


test_that("ridge and unpenalized fits are nearly identical when ridge is tiny", {
     d <- make_binary_gate_data(
          n = 250,
          seed = 12
     )
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit_none <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "none"
     )
     
     fit_ridge <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          ridge = 1e-8
     )
     
     expect_equal(
          fit_ridge$coef$estimate,
          fit_none$coef$estimate,
          tolerance = 1e-3
     )
     
     expect_equal(
          fit_ridge$logLik_full,
          fit_none$logLik_full,
          tolerance = 1e-5
     )
})


test_that("permutation p-value is between zero and one and detects a strong effect", {
     d <- make_binary_gate_data(
          n = 150
     )
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          n_perm = 100,
          seed = 42
     )
     
     expect_true(
          fit$p_perm >= 0 &&
               fit$p_perm <= 1
     )
     
     expect_true(fit$p_perm < 0.05)
     expect_equal(length(fit$perm_LRT), 100)
     expect_true(any(is.finite(fit$perm_LRT)))
})


test_that("mbmm_predict_gate and mbmm_gate_threshold agree with fitted coefficients", {
     d <- make_binary_gate_data()
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge"
     )
     
     b0 <- fit$coef$estimate[
          fit$coef$term == "(Intercept)"
     ]
     
     b1 <- fit$coef$estimate[
          fit$coef$term == "x"
     ]
     
     x_star <- -b0 / b1
     
     expect_equal(
          mbmm_gate_threshold(
               fit,
               "x"
          ),
          x_star,
          tolerance = 1e-4
     )
     
     p_at_star <- mbmm_predict_gate(
          fit,
          data.frame(x = x_star)
     )
     
     expect_equal(
          as.numeric(p_at_star),
          0.5,
          tolerance = 1e-4
     )
})


test_that("mbmm_fit_gate handles more than two classes via the multinomial path", {
     set.seed(5)
     
     n <- 240
     x <- rnorm(n)
     
     logL <- data.frame(
          logL_A = c(
               rep(log(0.9), n / 3),
               rep(log(0.05), 2 * n / 3)
          ),
          logL_B = c(
               rep(log(0.05), n / 3),
               rep(log(0.9), n / 3),
               rep(log(0.05), n / 3)
          ),
          logL_C = c(
               rep(log(0.05), 2 * n / 3),
               rep(log(0.9), n / 3)
          )
     )
     
     d <- data.frame(x = x)
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~1,
          data = d,
          classes = c("A", "B", "C"),
          penalty = "ridge"
     )
     
     expect_equal(fit$K, 3)
     expect_equal(nrow(fit$posterior), n)
     
     expect_true(
          all(
               fit$posterior$map_class %in%
                    c("A", "B", "C")
          )
     )
     
     expect_equal(fit$convergence_full, 0)
     expect_equal(fit$convergence_null, 0)
})


test_that("ridge keeps coefficients finite under near-complete separation", {
     set.seed(7)
     
     n <- 60
     x <- rnorm(n)
     z <- as.integer(x > 0)
     
     flip_ids <- sample(
          seq_len(n),
          2
     )
     
     z[flip_ids] <- 1L - z[flip_ids]
     
     logL <- data.frame(
          logL_pos = ifelse(
               z == 1,
               log(0.99),
               log(0.01)
          ),
          logL_neg = ifelse(
               z == 1,
               log(0.01),
               log(0.99)
          )
     )
     
     d <- data.frame(x = x)
     
     fit_ridge <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          ridge = 1e-4
     )
     
     expect_equal(
          fit_ridge$penalty,
          "ridge"
     )
     
     expect_true(
          all(
               is.finite(
                    fit_ridge$coef$estimate
               )
          )
     )
     
     beta_x <- fit_ridge$coef$estimate[
          fit_ridge$coef$term == "x"
     ]
     
     expect_gt(beta_x, 0)
     expect_equal(fit_ridge$convergence_full, 0)
})


test_that("larger ridge penalty produces stronger coefficient shrinkage", {
     set.seed(17)
     
     n <- 80
     x <- rnorm(n)
     z <- as.integer(x > 0)
     
     logL <- data.frame(
          logL_pos = ifelse(
               z == 1,
               log(0.995),
               log(0.005)
          ),
          logL_neg = ifelse(
               z == 1,
               log(0.005),
               log(0.995)
          )
     )
     
     d <- data.frame(x = x)
     
     fit_small <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          ridge = 1e-6
     )
     
     fit_large <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "ridge",
          ridge = 0.1
     )
     
     beta_small <- fit_small$coef$estimate[
          fit_small$coef$term == "x"
     ]
     
     beta_large <- fit_large$coef$estimate[
          fit_large$coef$term == "x"
     ]
     
     expect_true(
          abs(beta_large) <
               abs(beta_small)
     )
     
     expect_gt(beta_large, 0)
})


test_that("unsupported Firth penalty is rejected", {
     logL <- data.frame(
          logL_pos = c(-1, -2, -1),
          logL_neg = c(-2, -1, -2)
     )
     
     d <- data.frame(
          x = c(-1, 0, 1)
     )
     
     expect_error(
          mbmm_fit_gate(
               logL = logL,
               formula = ~x,
               data = d,
               classes = c("pos", "neg"),
               penalty = "firth"
          ),
          "`penalty` must be one of"
     )
})


test_that("mbmm_fit_gate drops subjects with non-finite class likelihoods", {
     d <- make_binary_gate_data(
          n = 50
     )
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     logL$logL_pos[1] <- -Inf
     
     expect_message(
          fit <- mbmm_fit_gate(
               logL = logL,
               formula = ~x,
               data = d,
               classes = c("pos", "neg")
          ),
          "Dropping 1 subject"
     )
     
     expect_equal(fit$n, 49)
     expect_equal(fit$rows_dropped, 1)
})


test_that("mbmm_fit_gate drops missing covariates without misaligning subjects", {
     d <- make_binary_gate_data(
          n = 50
     )
     
     d$x[c(3, 8)] <- NA_real_
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     expect_message(
          fit <- mbmm_fit_gate(
               logL = logL,
               formula = ~x,
               data = d,
               classes = c("pos", "neg"),
               id_col = "id"
          ),
          "incomplete gate covariates"
     )
     
     expect_equal(fit$n, 48)
     expect_equal(fit$rows_dropped, c(3, 8))
     
     expect_equal(
          fit$posterior$id,
          d$id[-c(3, 8)]
     )
})


test_that("penalty none stores zero ridge and unpenalized objective", {
     d <- make_binary_gate_data(
          n = 120,
          seed = 31
     )
     
     logL <- data.frame(
          logL_pos = d$logL_pos,
          logL_neg = d$logL_neg
     )
     
     fit <- mbmm_fit_gate(
          logL = logL,
          formula = ~x,
          data = d,
          classes = c("pos", "neg"),
          penalty = "none"
     )
     
     expect_equal(fit$penalty, "none")
     expect_equal(fit$ridge, 0)
     
     expect_equal(
          fit$objective_full_penalized,
          fit$logLik_full,
          tolerance = 1e-12
     )
     
     expect_equal(
          fit$objective_null_penalized,
          fit$logLik_null,
          tolerance = 1e-12
     )
})