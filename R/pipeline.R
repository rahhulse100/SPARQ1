sparq_pipeline <- function(data,
                           min_iterations = 50,
                           min_direction_consistency = 0.70,
                           max_sign_flip_rate = 0.30,
                           max_correction_fraction = 0.50,
                           estimator_aggregation = c("median", "mean")) {
  estimator_aggregation <- match.arg(estimator_aggregation)
  sparq_validate_results(data)

  reproducibility <- sparq_reproducibility_feature(data)
  bias <- sparq_bias(data)

  support <- sparq_support(
    bias_table = bias,
    reproducibility_table = reproducibility,
    min_iterations = min_iterations,
    min_direction_consistency = min_direction_consistency,
    max_sign_flip_rate = max_sign_flip_rate
  )

  calibration <- sparq_calibrate(
    data = data,
    support_table = support,
    estimator_aggregation = estimator_aggregation,
    max_correction_fraction = max_correction_fraction
  )

  list(
    reproducibility = reproducibility,
    bias = bias,
    support = support,
    calibration = calibration
  )
}
