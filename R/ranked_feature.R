sparq_prepare_ranked_feature_input <- function(
  data,
  result_id_col = "result_id",
  theta_full_col = "theta_full",
  theta_pert_col = "theta_pert",
  rank_full_col = "rank_full",
  rank_pert_col = "rank_pert",
  scenario_col = "scenario",
  iteration_col = "iteration",
  loss_fraction_col = "loss_fraction",
  cohort_col = "cohort",
  sample_id_col = "sample_id"
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("SPARQ requires data.table.", call. = FALSE)
  }

  x <- data.table::as.data.table(data)

  required_source <- c(
    result_id_col,
    theta_full_col,
    theta_pert_col,
    rank_full_col,
    rank_pert_col,
    scenario_col,
    iteration_col
  )

  missing <- setdiff(required_source, names(x))
  if (length(missing) > 0) {
    stop(
      "SPARQ input missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- data.table::data.table(
    result_id = as.character(x[[result_id_col]]),
    theta_full = as.numeric(x[[theta_full_col]]),
    theta_pert = as.numeric(x[[theta_pert_col]]),
    rank_full = as.numeric(x[[rank_full_col]]),
    rank_pert = as.numeric(x[[rank_pert_col]]),
    scenario = as.character(x[[scenario_col]]),
    iteration = as.integer(x[[iteration_col]])
  )

  if (!is.null(loss_fraction_col) && loss_fraction_col %in% names(x)) {
    out[, loss_fraction := as.numeric(x[[loss_fraction_col]])]
  } else {
    out[, loss_fraction := NA_real_]
  }

  if (!is.null(cohort_col) && cohort_col %in% names(x)) {
    out[, cohort := as.character(x[[cohort_col]])]
  } else {
    out[, cohort := "cohort_1"]
  }

  if (!is.null(sample_id_col) && sample_id_col %in% names(x)) {
    out[, sample_id := as.character(x[[sample_id_col]])]
  } else {
    out[, sample_id := "sample_1"]
  }

  out <- out[
    !is.na(result_id) &
      result_id != "" &
      is.finite(theta_full) &
      is.finite(theta_pert) &
      is.finite(rank_full) &
      is.finite(rank_pert) &
      !is.na(scenario) &
      scenario != "" &
      !is.na(iteration)
  ]

  if (nrow(out) == 0) {
    stop("SPARQ input has zero usable rows after validation.", call. = FALSE)
  }

  data.table::setcolorder(
    out,
    c(
      "cohort", "sample_id", "result_id",
      "theta_full", "theta_pert",
      "rank_full", "rank_pert",
      "scenario", "iteration", "loss_fraction"
    )
  )

  out[]
}

sparq_validate_ranked_feature <- function(data, ...) {
  invisible(sparq_prepare_ranked_feature_input(data, ...))
}

sparq_trimmed_median_internal <- function(x, trim = 0.2) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n == 0) return(NA_real_)

  lo <- floor(n * trim) + 1L
  hi <- ceiling(n * (1 - trim))

  if (lo > hi) return(stats::median(x))
  stats::median(x[lo:hi])
}

sparq_signed_quantile_internal <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)

  med <- stats::median(x)
  if (med < 0) return(as.numeric(stats::quantile(x, 0.25, names = FALSE)))
  if (med > 0) return(as.numeric(stats::quantile(x, 0.75, names = FALSE)))
  0
}

sparq_empirical_bayes_internal <- function(x, prior_strength = 20) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x) * (length(x) / (length(x) + prior_strength))
}

sparq_severity_regression_internal <- function(loss_fraction, delta) {
  ok <- is.finite(loss_fraction) & is.finite(delta)

  if (sum(ok) < 10 || length(unique(loss_fraction[ok])) < 2) {
    return(stats::median(delta[is.finite(delta)]))
  }

  fit <- try(stats::lm(delta[ok] ~ loss_fraction[ok]), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(stats::median(delta[is.finite(delta)]))
  }

  pred_loss <- stats::median(loss_fraction[ok])
  as.numeric(stats::coef(fit)[1] + stats::coef(fit)[2] * pred_loss)
}

sparq_safe_cor_internal <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (length(unique(x[ok])) < 2) return(NA_real_)
  if (length(unique(y[ok])) < 2) return(NA_real_)

  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

sparq_rank_turnover <- function(
  data,
  top_k = 100,
  ...
) {
  dt <- sparq_prepare_ranked_feature_input(data, ...)

  dt[, {
    full_ids <- result_id[order(rank_full, result_id)]
    pert_ids <- result_id[order(rank_pert, result_id)]

    full_top <- unique(full_ids[seq_len(min(top_k, length(full_ids)))])
    pert_top <- unique(pert_ids[seq_len(min(top_k, length(pert_ids)))])

    top_intersect <- length(intersect(full_top, pert_top))
    top_union <- length(union(full_top, pert_top))
    top_lost <- length(setdiff(full_top, pert_top))
    top_gained <- length(setdiff(pert_top, full_top))

    list(
      top_k = top_k,
      n_features_compared = .N,
      spearman_rank_correlation = sparq_safe_cor_internal(rank_full, rank_pert, "spearman"),
      topk_jaccard = if (top_union == 0) NA_real_ else top_intersect / top_union,
      topk_lost = top_lost,
      topk_gained = top_gained,
      topk_replaced_fraction = if (length(full_top) == 0) NA_real_ else top_lost / length(full_top)
    )
  }, by = .(cohort, sample_id, scenario, iteration)]
}

sparq_bias_ranked_feature <- function(
  data,
  ...
) {
  dt <- sparq_prepare_ranked_feature_input(data, ...)
  dt[, delta := theta_pert - theta_full]
  dt[, abs_delta := abs(delta)]

  dt[, {
    mb <- stats::median(delta, na.rm = TRUE)
    tb <- sparq_trimmed_median_internal(delta)
    sr <- sparq_severity_regression_internal(loss_fraction, delta)
    eb <- sparq_empirical_bayes_internal(delta)
    sq <- sparq_signed_quantile_internal(delta)

    estimators <- c(
      median_bias = mb,
      trimmed_bias = tb,
      severity_regression = sr,
      empirical_bayes = eb,
      signed_quantile = sq
    )
    estimators <- estimators[is.finite(estimators)]

    calibration_bias <- if (length(estimators) > 0) {
      stats::median(estimators)
    } else {
      NA_real_
    }

    estimator_signs <- sign(estimators[abs(estimators) > 1e-12])

    estimator_direction_agreement <- if (length(estimator_signs) > 0) {
      mean(estimator_signs == sign(stats::median(estimators)))
    } else {
      0
    }

    direction_consistency <- if (is.finite(calibration_bias) && abs(calibration_bias) > 1e-12) {
      mean(sign(delta) == sign(calibration_bias), na.rm = TRUE)
    } else {
      0
    }

    sign_flip_rate <- {
      base_sign <- sign(theta_full[1])
      if (!is.finite(base_sign) || base_sign == 0) {
        NA_real_
      } else {
        mean(sign(theta_pert) != base_sign, na.rm = TRUE)
      }
    }

    uncertainty_mad <- stats::median(
      abs(delta - stats::median(delta, na.rm = TRUE)),
      na.rm = TRUE
    )

    signal_to_uncertainty <- abs(calibration_bias) / (uncertainty_mad + 1e-8)
    relative_bias <- abs(calibration_bias) / (abs(theta_full[1]) + 1e-8)

    n_runs <- data.table::uniqueN(paste(scenario, iteration, sep = "__"))

    list(
      theta_full = theta_full[1],
      rank_full = rank_full[1],
      n_runs = n_runs,
      calibration_bias = calibration_bias,
      median_bias = mb,
      trimmed_bias = tb,
      severity_regression = sr,
      empirical_bayes = eb,
      signed_quantile = sq,
      direction_consistency = direction_consistency,
      estimator_direction_agreement = estimator_direction_agreement,
      sign_flip_rate = sign_flip_rate,
      uncertainty_mad = uncertainty_mad,
      signal_to_uncertainty = signal_to_uncertainty,
      relative_bias = relative_bias
    )
  }, by = .(cohort, sample_id, result_id)]
}

sparq_calibrate_ranked_feature <- function(
  bias_table,
  expected_runs = 150,
  max_correction_fraction = 0.50,
  min_runs = 30,
  min_support_quality = 0.35,
  min_direction_consistency = 0.65,
  max_sign_flip_rate = 0.25
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("SPARQ requires data.table.", call. = FALSE)
  }

  x <- data.table::as.data.table(bias_table)

  required <- c(
    "cohort", "sample_id", "result_id",
    "theta_full", "rank_full", "n_runs",
    "calibration_bias", "direction_consistency",
    "estimator_direction_agreement",
    "sign_flip_rate", "signal_to_uncertainty"
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "SPARQ bias table missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  x[, run_support := pmin(n_runs / expected_runs, 1)]
  x[, direction_support := pmax(0, pmin(direction_consistency, 1))]
  x[, estimator_support := pmax(0, pmin(estimator_direction_agreement, 1))]
  x[, uncertainty_support := pmax(0, pmin(signal_to_uncertainty / 2, 1))]
  x[, flip_support := 1 - pmax(0, pmin(ifelse(is.na(sign_flip_rate), 1, sign_flip_rate), 1))]

  x[, support_quality := rowMeans(
    cbind(
      run_support,
      direction_support,
      estimator_support,
      uncertainty_support,
      flip_support
    ),
    na.rm = TRUE
  )]

  x[, lambda := support_quality * direction_support * estimator_support]
  x[, correction_raw := -lambda * calibration_bias]
  x[, correction_bound := max_correction_fraction * abs(theta_full)]

  x[, correction := pmax(-correction_bound, pmin(correction_bound, correction_raw))]
  x[, theta_calibrated := theta_full + correction]

  x[, decision := data.table::fcase(
    n_runs < min_runs, "insufficient_support",
    support_quality < min_support_quality, "insufficient_support",
    direction_consistency < min_direction_consistency, "insufficient_support",
    !is.na(sign_flip_rate) & sign_flip_rate > max_sign_flip_rate, "insufficient_support",
    abs(correction) < 1e-12, "hold",
    correction > 0, "calibrate_up",
    correction < 0, "calibrate_down",
    default = "hold"
  )]

  x[decision == "insufficient_support", correction := 0]
  x[decision == "insufficient_support", theta_calibrated := theta_full]

  x[, rank_calibrated := data.table::frank(-theta_calibrated, ties.method = "min"), by = .(cohort, sample_id)]
  x[, rank_shift_after_calibration := rank_calibrated - rank_full]

  keep <- c(
    "cohort", "sample_id", "result_id",
    "theta_full", "theta_calibrated", "correction",
    "calibration_bias", "rank_full", "rank_calibrated",
    "rank_shift_after_calibration", "decision",
    "support_quality", "n_runs",
    "direction_consistency", "estimator_direction_agreement",
    "sign_flip_rate", "uncertainty_mad",
    "signal_to_uncertainty", "relative_bias",
    "median_bias", "trimmed_bias", "severity_regression",
    "empirical_bayes", "signed_quantile"
  )

  as.data.frame(x[, ..keep])
}

sparq_sample_fragility_ranked_feature <- function(turnover_table) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("SPARQ requires data.table.", call. = FALSE)
  }

  x <- data.table::as.data.table(turnover_table)

  required <- c(
    "cohort", "sample_id",
    "spearman_rank_correlation",
    "topk_jaccard",
    "topk_replaced_fraction"
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "SPARQ turnover table missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  out <- x[, .(
    n_runs = .N,
    median_spearman_rank_correlation = stats::median(spearman_rank_correlation, na.rm = TRUE),
    q25_spearman_rank_correlation = as.numeric(stats::quantile(spearman_rank_correlation, 0.25, na.rm = TRUE)),
    median_topk_jaccard = stats::median(topk_jaccard, na.rm = TRUE),
    median_topk_replaced_fraction = stats::median(topk_replaced_fraction, na.rm = TRUE),
    fragility_score = (1 - stats::median(spearman_rank_correlation, na.rm = TRUE)) +
      stats::median(topk_replaced_fraction, na.rm = TRUE)
  ), by = .(cohort, sample_id)]

  as.data.frame(out)
}

sparq_pipeline_ranked_feature <- function(
  data,
  top_k = 100,
  expected_runs = 150,
  max_correction_fraction = 0.50,
  min_runs = 30,
  min_support_quality = 0.35,
  min_direction_consistency = 0.65,
  max_sign_flip_rate = 0.25,
  ...
) {
  input <- sparq_prepare_ranked_feature_input(data, ...)

  input_audit <- data.frame(
    n_rows = nrow(input),
    n_cohorts = length(unique(input$cohort)),
    n_samples = length(unique(paste(input$cohort, input$sample_id, sep = "__"))),
    n_results = length(unique(input$result_id)),
    n_scenarios = length(unique(input$scenario)),
    n_iterations = length(unique(input$iteration)),
    stringsAsFactors = FALSE
  )

  turnover <- sparq_rank_turnover(input, top_k = top_k)
  bias <- sparq_bias_ranked_feature(input)

  calibration <- sparq_calibrate_ranked_feature(
    bias,
    expected_runs = expected_runs,
    max_correction_fraction = max_correction_fraction,
    min_runs = min_runs,
    min_support_quality = min_support_quality,
    min_direction_consistency = min_direction_consistency,
    max_sign_flip_rate = max_sign_flip_rate
  )

  sample_summary <- sparq_sample_fragility_ranked_feature(turnover)

  list(
    input_audit = input_audit,
    turnover = turnover,
    bias = as.data.frame(bias),
    calibration = calibration,
    sample_summary = sample_summary
  )
}
