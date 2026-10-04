test_that("presets resolve known protocols and preserve explicit overrides", {
  quick <- sparq_preset("quick")
  expect_equal(quick$n_iterations, 25L)
  expect_equal(quick$retention, 0.75)

  resolved <- sparq_resolve_preset(
    "quick", n_iterations = 7, min_iterations = 3,
    instability_bootstrap_B = 20, verbose = FALSE
  )
  expect_equal(resolved$n_iterations, 7)
  expect_equal(resolved$min_iterations, 3)
  expect_equal(resolved$instability_bootstrap_B, 20)
  expect_false(resolved$verbose)
})

test_that("standard methods expose assessment summaries", {
  set.seed(51)
  spots <- data.frame(value = rnorm(20), x = runif(20), y = runif(20))
  fit <- sparq_assess_workflow(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    preset = "quick",
    min_iterations = 3,
    instability_bootstrap_B = 20,
    verbose = FALSE,
    seed = 8
  )
  expect_s3_class(summary(fit), "summary.sparq_assessment")
  expect_true(is.data.frame(as.data.frame(fit)))
  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(fit)
  })
})

test_that("workflow cohort writes a complete audit manifest", {
  set.seed(52)
  sections <- list(
    section_a = data.frame(value = rnorm(20), x = runif(20), y = runif(20)),
    section_b = data.frame(value = rnorm(20), x = runif(20), y = runif(20))
  )
  output_dir <- tempfile("sparq_cohort_")
  fit <- sparq_assess_cohort_workflow(
    sample_data = sections,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    preset = "quick",
    n_iterations = 4,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    output_dir = output_dir,
    verbose = FALSE,
    seed = 9
  )
  expect_s3_class(fit, "sparq_cohort_assessment")
  expect_equal(nrow(fit$manifest), 2)
  expect_true(all(fit$manifest$status == "completed"))
  expect_true(file.exists(fit$output_files$manifest))
  expect_true(file.exists(fit$output_files$settings))
})
