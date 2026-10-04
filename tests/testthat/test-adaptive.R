test_that("adaptive assessment stops after a conclusive support interval", {
  set.seed(2)
  data <- data.frame(value = rnorm(60), x = runif(60), y = runif(60))

  fit <- sparq_assess_adaptive(
    data = data,
    analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar,
    retention = 0.75,
    support_threshold = 0.80,
    batch_size = 5,
    min_iterations = 10,
    max_iterations = 20,
    support_bootstrap_B = 20,
    instability_bootstrap_B = 20,
    seed = 3
  )

  expect_true(nrow(fit$adaptive_history) >= 1)
  expect_true(fit$adaptive_decision$n_attempted <= 20)
  expect_equal(nrow(fit$comparisons), fit$adaptive_decision$n_attempted)
  expect_true(fit$adaptive_decision$decision %in%
    c("above_threshold", "below_threshold", "inconclusive"))
})

test_that("adaptive assessment requires a user-specified threshold", {
  data <- data.frame(value = 1:10)
  expect_error(
    sparq_assess_adaptive(
      data = data,
      analysis_function = function(x) mean(x$value),
      comparator = sparq_compare_scalar
    ),
    "support_threshold"
  )
})
