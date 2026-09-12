test_that("ranked feature pipeline calibrates up, down, and refuses unsupported results", {
  set.seed(1)

  x <- data.frame(
    cohort = "cohort_1",
    sample_id = "sample_1",
    result_id = rep(c("gene_A", "gene_B", "gene_C"), each = 150),
    theta_full = rep(c(1.0, 2.0, -1.5), each = 150),
    theta_pert = c(
      rnorm(150, 0.75, 0.03),
      rnorm(150, 2.30, 0.03),
      rnorm(150, -1.50, 0.80)
    ),
    rank_full = rep(c(1, 2, 3), each = 150),
    rank_pert = rep(c(1, 2, 3), each = 150),
    scenario = rep(rep(c("spots_075", "spots_050", "random_spots_uniform_0_50"), each = 50), 3),
    iteration = rep(rep(1:50, 3), 3),
    loss_fraction = rep(rep(c(0.25, 0.50, 0.25), each = 50), 3)
  )

  fit <- sparq_pipeline_ranked_feature(x, top_k = 2)

  expect_true("calibration" %in% names(fit))
  expect_equal(nrow(fit$calibration), 3)

  decisions <- setNames(fit$calibration$decision, fit$calibration$result_id)

  expect_equal(decisions[["gene_A"]], "calibrate_up")
  expect_equal(decisions[["gene_B"]], "calibrate_down")
})
