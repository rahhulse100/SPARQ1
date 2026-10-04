test_that("frontier records each evaluated retention without interpolation", {
  set.seed(8)
  spots <- data.frame(value = rnorm(40), x = runif(40), y = runif(40))
  frontier <- sparq_assess_frontier(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    retentions = c(0.5, 0.75, 1),
    reliability_tolerance = 1,
    required_reliability = 0.5,
    confidence = 0.8,
    preset = "quick",
    n_iterations = 6,
    min_iterations = 3,
    instability_bootstrap_B = 10,
    verbose = FALSE,
    seed = 9
  )

  expect_s3_class(frontier, "sparq_reliability_frontier")
  expect_equal(frontier$frontier$retention, c(0.5, 0.75, 1))
  expect_equal(frontier$per_retention_confidence, 1 - (1 - 0.8) / 3)
  expect_true(all(frontier$frontier$per_retention_confidence ==
    frontier$per_retention_confidence))
  expect_true(is.na(frontier$critical_retention) ||
    frontier$critical_retention %in% frontier$frontier$retention)
})

test_that("cohort reliability pools section certificates without treating spots as sections", {
  assessments <- list(
    section_a = list(reliability_decision = data.frame(
      n_attempted = 20L, n_preserved = 19L, n_completed = 20L, n_failed = 0L,
      preservation_rate = 0.95, tolerance = 0.1, required_reliability = 0.8,
      confidence = 0.95, reference_scale = 1, decision = "supported"
    )),
    section_b = list(reliability_decision = data.frame(
      n_attempted = 20L, n_preserved = 12L, n_completed = 20L, n_failed = 0L,
      preservation_rate = 0.60, tolerance = 0.1, required_reliability = 0.8,
      confidence = 0.95, reference_scale = 1, decision = "inconclusive"
    ))
  )
  fit <- sparq_infer_cohort_reliability(
    assessments,
    confidence = 0.9,
    section_metadata = data.frame(sample_id = c("section_a", "section_b"),
      quality = c(1000, 500))
  )

  expect_s3_class(fit, "sparq_cohort_reliability")
  expect_equal(nrow(fit$section_reliability), 2)
  expect_equal(fit$cohort_summary$n_sections, 2)
  expect_true("quality" %in% names(fit$section_reliability))
})
