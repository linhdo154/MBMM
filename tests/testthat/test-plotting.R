test_that("plot_bank returns a ggplot and validates weights", {
  candidates <- expand.grid(peak = 10^(2:3), rate = c(0.1, 0.3))
  p <- plot_bank(candidates, 1:4, "peak", "rate", log10_x = TRUE)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_error(plot_bank(candidates, c(1, -1, 1, 1), "peak", "rate"),
               "non-negative")
  expect_error(plot_bank(candidates, 1:3, "peak", "rate"),
               "one value per candidate")
})

test_that("recovery plots return ggplots and validate simulation truth", {
  fit <- list(w = c(0.1, 0.8, 0.3, 0.9),
              pi = c(0.2, 0.7, 0.4, 0.8))
  truth <- c(0, 1, 0, 1)

  p_resp <- plot_responsibility(fit, truth)
  p_gate <- plot_gate(fit, 1:4, truth, c(0.1, 0.4, 0.7, 0.9))
  expect_s3_class(p_resp, "ggplot")
  expect_s3_class(p_gate, "ggplot")
  expect_silent(ggplot2::ggplot_build(p_resp))
  expect_silent(ggplot2::ggplot_build(p_gate))
  expect_error(plot_responsibility(fit, c(0, 1)), "matching")
  expect_error(plot_gate(fit, 1:3), "matching")
})
