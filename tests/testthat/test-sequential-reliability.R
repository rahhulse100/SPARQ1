test_that("confidence sequence records every look and conservatively counts failures", {
  comparisons <- data.frame(
    status = c(rep("ok", 99), "failed"),
    instability = c(rep(0.01, 99), NA_real_)
  )
  sequence <- sparq_reliability_confidence_sequence(
    comparisons,
    tolerance = 0.05,
    required_reliability = 0.80,
    confidence = 0.95
  )

  expect_equal(nrow(sequence), 100)
  expect_equal(sequence$n_attempted[100], 100)
  expect_equal(sequence$n_failed[100], 1)
  expect_equal(sequence$n_preserved[100], 99)
  expect_true(all(sequence$lower_confidence_bound >= 0))
  expect_true(all(sequence$upper_confidence_bound <= 1))
  expect_equal(sequence$method[100], "beta_binomial_mixture_confidence_sequence")
})

test_that("adaptive reliability uses an anytime-valid confidence sequence", {
  set.seed(102)
  spots <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
  fit <- sparq_assess_adaptive(
    data = spots,
    analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar,
    reliability_tolerance = 1,
    required_reliability = 0.5,
    reliability_confidence = 0.8,
    preset = "quick",
    max_iterations = 6,
    min_iterations = 3,
    batch_size = 3,
    instability_bootstrap_B = 10,
    verbose = FALSE,
    seed = 12
  )

  expect_true(is.data.frame(fit$reliability_confidence_sequence))
  expect_equal(
    fit$reliability_decision$method,
    "beta_binomial_mixture_confidence_sequence"
  )
  expect_equal(fit$adaptive_decision$decision_mode, "reliability_certificate")
})
