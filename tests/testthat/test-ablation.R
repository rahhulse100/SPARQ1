test_that("ablation scores remove only the requested support component", {
  rows <- data.frame(
    support_quality = c(0.4, 0.3),
    iteration_support = c(0.8, 0.6),
    bias_support = c(0.5, 0.5),
    uncertainty_support = c(0.5, 0.5),
    reproducibility_support = c(1, 1)
  )

  expect_equal(sparq_ablation_score(rows, "full"), c(0.4, 0.3))
  expect_equal(sparq_ablation_score(rows, "no_uncertainty"), c(0.4, 0.3))
  expect_equal(sparq_ablation_score(rows, "no_reproducibility"), c(0.2, 0.15))
})

test_that("continuous ablation metrics expose threshold-independent metrics", {
  metrics <- sparq_ablation_continuous_metrics(
    truth = c(1, 1, 0, 0),
    score = c(0.9, 0.8, 0.2, 0.1)
  )
  expect_equal(metrics$auroc, 1)
  expect_equal(metrics$auprc, 1)
  expect_true(is.na(metrics$sensitivity))
})
