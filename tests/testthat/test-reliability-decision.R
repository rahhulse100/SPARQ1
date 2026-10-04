test_that("reliability decision uses conservative exact intervals", {
  comparisons <- data.frame(
    instability = c(rep(0.02, 19), 0.20),
    status = rep("ok", 20),
    iteration = seq_len(20)
  )
  decision <- sparq_reliability_decision(
    comparisons, tolerance = 0.05,
    required_reliability = 0.70, confidence = 0.95
  )
  expect_equal(decision$n_preserved, 19)
  expect_equal(decision$preservation_rate, 0.95)
  expect_true(decision$lower_confidence_bound < decision$preservation_rate)
  expect_equal(decision$decision, "supported")

  comparisons$status[20] <- "failed"
  decision_with_failure <- sparq_reliability_decision(
    comparisons, tolerance = 0.05,
    required_reliability = 0.70, confidence = 0.95
  )
  expect_equal(decision_with_failure$n_preserved, 19)
  expect_equal(decision_with_failure$n_failed, 1)
})

test_that("workflow assessment attaches a reliability decision without changing support calibration", {
  set.seed(71)
  spots <- data.frame(value = rnorm(30), x = runif(30), y = runif(30))
  fit <- sparq_assess_workflow(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    preset = "quick",
    n_iterations = 8,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    reliability_tolerance = 1,
    required_reliability = 0.50,
    verbose = FALSE,
    seed = 7
  )
  expect_true(is.data.frame(fit$reliability_decision))
  expect_true(all(c("support_quality", "reliability_decision") %in% names(fit$summary)))
  expect_equal(fit$settings$reliability_tolerance, 1)
})

test_that("adaptive reliability mode does not require a calibrated support threshold", {
  set.seed(72)
  spots <- data.frame(value = rnorm(30), x = runif(30), y = runif(30))
  fit <- sparq_assess_adaptive(
    data = spots,
    analysis_function = function(x) mean(x$value),
    comparator = sparq_compare_scalar,
    preset = "quick",
    max_iterations = 6,
    min_iterations = 3,
    batch_size = 3,
    instability_bootstrap_B = 20,
    reliability_tolerance = 1,
    required_reliability = 0.50,
    reliability_confidence = 0.80,
    verbose = FALSE,
    seed = 8
  )
  expect_equal(fit$adaptive_decision$decision_mode, "reliability_certificate")
  expect_true(is.data.frame(fit$reliability_decision))
})
