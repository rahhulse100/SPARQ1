test_that("quality-weighted loss has fixed retained count and reproducible audit", {
  spots <- data.frame(
    quality = seq_len(40),
    x = seq_len(40),
    y = rev(seq_len(40))
  )
  model <- sparq_define_loss_model(
    type = "quality_weighted",
    quality_col = "quality",
    quality_direction = "higher_quality_retained"
  )
  first <- sparq_apply_loss_model(spots, model, retention = 0.75, seed = 4)
  second <- sparq_apply_loss_model(spots, model, retention = 0.75, seed = 4)

  expect_equal(nrow(first$data), 30)
  expect_equal(first$retained_index, second$retained_index)
  expect_equal(sum(first$audit$retained), 30)
  expect_true(mean(first$data$quality) > mean(spots$quality))
})

test_that("edge-weighted loss validates coordinates and produces a spot audit", {
  spots <- data.frame(x = runif(30), y = runif(30))
  model <- sparq_define_loss_model("edge_weighted", edge_strength = 2)
  application <- sparq_apply_loss_model(
    spots, model, retention = 0.7, x_col = "x", y_col = "y", seed = 2
  )

  expect_equal(nrow(application$data), 21)
  expect_true(all(c("retention_weight", "relative_loss_weight", "retained") %in%
    names(application$audit)))
  expect_error(
    sparq_apply_loss_model(spots, model, retention = 0.7),
    "x_col and y_col"
  )
})

test_that("workflow records a supplied spatial loss model", {
  spots <- data.frame(value = rnorm(30), quality = seq_len(30))
  model <- sparq_define_loss_model("quality_weighted", quality_col = "quality")
  fit <- sparq_assess_workflow(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    loss_model = model,
    preset = "quick",
    n_iterations = 6,
    min_iterations = 3,
    instability_bootstrap_B = 10,
    verbose = FALSE,
    seed = 2
  )
  expect_equal(fit$settings$loss_model$type, "quality_weighted")
  expect_true(all(fit$comparisons$loss_model == "quality_weighted"))
  expect_true(is.list(fit$loss_audits))
  expect_true(all(vapply(fit$loss_audits, is.data.frame, logical(1))))
})

test_that("reports export loss and reliability audit tables", {
  spots <- data.frame(value = rnorm(30), quality = seq_len(30))
  model <- sparq_define_loss_model("quality_weighted", quality_col = "quality")
  assessment <- sparq_assess_workflow(
    data = spots,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    loss_model = model,
    n_iterations = 6,
    min_iterations = 3,
    instability_bootstrap_B = 10,
    reliability_tolerance = 1,
    required_reliability = 0.5,
    verbose = FALSE,
    seed = 20
  )
  output_dir <- tempfile("sparq-loss-report-")
  report <- sparq_report(assessment, output_dir = output_dir,
    include_plots = FALSE, verbose = FALSE)

  expect_true(file.exists(report$output_files$spot_loss_audit))
  expect_true(file.exists(report$output_files$reliability_decision))
})
