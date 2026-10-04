test_that("direct, serial-cohort, and parallel-cohort assessments return one result per sample", {
  set.seed(4)
  samples <- list(
    section_a = data.frame(value = rnorm(30)),
    section_b = data.frame(value = rnorm(30))
  )
  analysis <- function(data) mean(data$value)

  direct <- sparq_assess(
    full_result = analysis(samples$section_a),
    perturbed_results = list(0.1, 0.2, 0.3),
    comparator = sparq_compare_scalar,
    min_iterations = 3,
    instability_bootstrap_B = 20
  )
  expect_equal(nrow(direct$comparisons), 3)

  serial <- sparq_assess_cohort(
    sample_data = samples,
    analysis_function = analysis,
    comparator = sparq_compare_scalar,
    retention = 0.75,
    n_iterations = 4,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    seed = 5
  )
  expect_equal(nrow(serial), 2)
  expect_equal(sort(serial$result_id), sort(names(samples)))

  parallel <- sparq_assess_cohort_parallel(
    sample_data = samples,
    analysis_function = analysis,
    comparator = sparq_compare_scalar,
    retention = 0.75,
    n_iterations = 4,
    workers = 1,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    seed = 5
  )
  expect_equal(nrow(parallel$results), 2)
  expect_equal(nrow(parallel$errors), 0)

  summary <- sparq_run_with_stress_model(
    data = samples$section_a,
    analysis_function = analysis,
    comparator = sparq_compare_scalar,
    retention = 0.75,
    n_iterations = 4,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    seed = 5
  )
  expect_true("support_quality" %in% names(summary))
})
