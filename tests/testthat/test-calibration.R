test_that("scalar calibration honors correction policy", {
  perturbed <- seq(0.89, 0.91, length.out = 60)
  raw <- sparq_calibrate_scalar(
    observed_value = 1,
    perturbed_values = perturbed,
    reference_type = "pathology",
    apply_correction = TRUE,
    B = 100
  )
  expect_equal(raw$calibration_status, "validated")
  expect_true(raw$apply_correction)
  expect_true(is.finite(raw$calibrated_value))
})

test_that("threshold calibration never uses held-out labels to choose its threshold", {
  train_score <- c(0.1, 0.2, 0.8, 0.9)
  train_truth <- c(0, 0, 1, 1)
  first <- sparq_calibrate_threshold(train_score, train_truth, c(0.2, 0.8), c(0, 1), B = 20)
  second <- sparq_calibrate_threshold(train_score, train_truth, c(0.2, 0.8), c(1, 0), B = 20)
  expect_equal(first$threshold, second$threshold)
})
