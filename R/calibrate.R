sparq_calibration_estimators <- function(theta_full, theta_pert) {
  delta <- theta_pert - theta_full
  delta <- delta[is.finite(delta)]

  if (length(delta) < 3) {
    return(c(
      median_bias = NA_real_,
      trimmed_bias = NA_real_,
      severity_regression = NA_real_,
      empirical_bayes = NA_real_,
      signed_quantile = NA_real_
    ))
  }

  median_bias <- stats::median(delta, na.rm = TRUE)
  trimmed_bias <- mean(delta, trim = 0.20, na.rm = TRUE)

  severity <- abs(delta)
  severity_regression <- if (stats::sd(severity, na.rm = TRUE) > 0) {
    fit <- stats::lm(delta ~ severity)
    as.numeric(stats::predict(fit, newdata = data.frame(severity = stats::median(severity, na.rm = TRUE))))
  } else {
    median_bias
  }

  between_variance <- stats::var(delta, na.rm = TRUE)
  within_variance <- stats::mad(delta, constant = 1.4826, na.rm = TRUE)^2
  eb_weight <- if (is.finite(between_variance + within_variance) && between_variance + within_variance > 0) {
    between_variance / (between_variance + within_variance)
  } else {
    0
  }
  empirical_bayes <- eb_weight * median_bias

  signed_quantile <- if (median_bias < 0) {
    stats::quantile(delta, probs = 0.25, na.rm = TRUE, names = FALSE)
  } else {
    stats::quantile(delta, probs = 0.75, na.rm = TRUE, names = FALSE)
  }

  c(
    median_bias = median_bias,
    trimmed_bias = trimmed_bias,
    severity_regression = severity_regression,
    empirical_bayes = empirical_bayes,
    signed_quantile = signed_quantile
  )
}

sparq_calibrate <- function(data,
                            support_table,
                            result_id_col = "result_id",
                            theta_full_col = "theta_full",
                            theta_pert_col = "theta_pert",
                            estimator_aggregation = c("median", "mean"),
                            max_correction_fraction = 0.50,
                            min_estimator_direction_agreement = 0.60) {
  estimator_aggregation <- match.arg(estimator_aggregation)
  sparq_validate_results(data, required_cols = c(result_id_col, theta_full_col, theta_pert_col))

  needed <- c("result_id", "support_quality", "insufficient_support")
  missing <- setdiff(needed, colnames(support_table))
  if (length(missing)) {
    stop("Missing required support columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  ids <- unique(as.character(data[[result_id_col]]))

  out <- lapply(ids, function(id) {
    d <- data[as.character(data[[result_id_col]]) == id, , drop = FALSE]
    s <- support_table[as.character(support_table$result_id) == id, , drop = FALSE]

    theta_full <- stats::median(d[[theta_full_col]], na.rm = TRUE)
    est <- sparq_calibration_estimators(
      theta_full = d[[theta_full_col]],
      theta_pert = d[[theta_pert_col]]
    )

    est_finite <- est[is.finite(est)]

    if (!nrow(s) || !length(est_finite)) {
      return(data.frame(
        result_id = id,
        theta_full = theta_full,
        theta_calibrated = theta_full,
        correction = 0,
        calibration_bias = NA_real_,
        estimator_direction_agreement = NA_real_,
        support_quality = NA_real_,
        decision = "insufficient_support",
        t(est),
        check.names = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    main_sign <- sign(stats::median(est_finite, na.rm = TRUE))
    direction_agreement <- mean(sign(est_finite) == main_sign)

    calibration_bias <- if (estimator_aggregation == "median") {
      stats::median(est_finite, na.rm = TRUE)
    } else {
      mean(est_finite, na.rm = TRUE)
    }

    correction <- -sparq_bound01(s$support_quality[1]) * calibration_bias

    bound <- max_correction_fraction * max(abs(theta_full), .Machine$double.eps)
    correction <- max(-bound, min(bound, correction))

    unsafe <- isTRUE(s$insufficient_support[1]) ||
      is.na(direction_agreement) ||
      direction_agreement < min_estimator_direction_agreement

    if (unsafe) {
      decision <- "insufficient_support"
      theta_calibrated <- theta_full
      correction <- 0
    } else {
      theta_calibrated <- theta_full + correction
      if (correction > 0) {
        decision <- "calibrate_up"
      } else if (correction < 0) {
        decision <- "calibrate_down"
      } else {
        decision <- "hold"
      }
    }

    data.frame(
      result_id = id,
      theta_full = theta_full,
      theta_calibrated = theta_calibrated,
      correction = correction,
      calibration_bias = calibration_bias,
      estimator_direction_agreement = direction_agreement,
      support_quality = s$support_quality[1],
      decision = decision,
      t(est),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

sparq_decide <- function(calibration_table) {
  needed <- c("result_id", "theta_full", "theta_calibrated", "correction", "decision")
  missing <- setdiff(needed, colnames(calibration_table))
  if (length(missing)) {
    stop("Missing required calibration columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  calibration_table
}
