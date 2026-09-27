sparq_apply_builtin_calibration <- function(support_quality, calibration_model) {
    probability <- predict(calibration_model$model, newdata = data.frame(support_quality = support_quality), 
        type = "response")
    support_class <- ifelse(support_quality >= calibration_model$threshold, "supported", "insufficient")
    data.frame(support_quality = support_quality, support_probability = probability, support_threshold = calibration_model$threshold, 
        support_class = support_class, stringsAsFactors = FALSE)
}

sparq_bias <- function(data, result_id_col = "result_id", theta_full_col = "theta_full", theta_pert_col = "theta_pert") {
    sparq_validate_results(data, required_cols = c(result_id_col, theta_full_col, theta_pert_col))
    ids <- unique(as.character(data[[result_id_col]]))
    out <- lapply(ids, function(id) {
        d <- data[as.character(data[[result_id_col]]) == id, , drop = FALSE]
        theta_full <- stats::median(d[[theta_full_col]], na.rm = TRUE)
        delta <- d[[theta_pert_col]] - d[[theta_full_col]]
        delta <- delta[is.finite(delta)]
        bias <- stats::median(delta, na.rm = TRUE)
        uncertainty <- stats::mad(delta, constant = 1.4826, na.rm = TRUE)
        direction_consistency <- sparq_mode_sign_fraction(delta)
        data.frame(result_id = id, theta_full = theta_full, n_perturbations = length(delta), bias = bias, 
            uncertainty = uncertainty, direction_consistency = direction_consistency, sign_flip_rate = 1 - 
                direction_consistency, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

sparq_calibrate <- function(data, support_table, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", estimator_aggregation = c("median", "mean"), max_correction_fraction = 0.5, 
    min_estimator_direction_agreement = 0.6) {
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
        est <- sparq_calibration_estimators(theta_full = d[[theta_full_col]], theta_pert = d[[theta_pert_col]])
        est_finite <- est[is.finite(est)]
        if (!nrow(s) || !length(est_finite)) {
            return(data.frame(result_id = id, theta_full = theta_full, theta_calibrated = theta_full, 
                correction = 0, calibration_bias = NA_real_, estimator_direction_agreement = NA_real_, 
                support_quality = NA_real_, decision = "insufficient_support", t(est), check.names = FALSE, 
                stringsAsFactors = FALSE))
        }
        main_sign <- sign(stats::median(est_finite, na.rm = TRUE))
        direction_agreement <- mean(sign(est_finite) == main_sign)
        calibration_bias <- if (estimator_aggregation == "median") {
            stats::median(est_finite, na.rm = TRUE)
        }
        else {
            mean(est_finite, na.rm = TRUE)
        }
        correction <- -sparq_bound01(s$support_quality[1]) * calibration_bias
        bound <- max_correction_fraction * max(abs(theta_full), .Machine$double.eps)
        correction <- max(-bound, min(bound, correction))
        unsafe <- isTRUE(s$insufficient_support[1]) || is.na(direction_agreement) || direction_agreement < 
            min_estimator_direction_agreement
        if (unsafe) {
            decision <- "insufficient_support"
            theta_calibrated <- theta_full
            correction <- 0
        }
        else {
            theta_calibrated <- theta_full + correction
            if (correction > 0) {
                decision <- "calibrate_up"
            }
            else if (correction < 0) {
                decision <- "calibrate_down"
            }
            else {
                decision <- "hold"
            }
        }
        data.frame(result_id = id, theta_full = theta_full, theta_calibrated = theta_calibrated, correction = correction, 
            calibration_bias = calibration_bias, estimator_direction_agreement = direction_agreement, 
            support_quality = s$support_quality[1], decision = decision, t(est), check.names = FALSE, 
            stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

sparq_calibrate_scalar_raw <- function(observed_value, perturbed_values, reference_type = c("none", "simulation", 
    "replicate", "pathology"), apply_correction = FALSE, alpha = 0.05, B = 2000, min_iterations = 20, 
    min_direction_consistency = 0.8, min_reproducibility = 0.7, max_relative_uncertainty = 0.25) {
    reference_type <- match.arg(reference_type)
    observed_value <- as.numeric(observed_value)
    perturbed_values <- as.numeric(unlist(perturbed_values))
    perturbed_values <- perturbed_values[is.finite(perturbed_values)]
    n_iterations <- length(perturbed_values)
    base_output <- data.frame(calibration_status = "not_supported", calibration_reference = reference_type, 
        calibration_direction = "none", calibrated_value = observed_value, apply_correction = FALSE, 
        calibration_lower = NA_real_, calibration_upper = NA_real_, calibration_bias = NA_real_, calibration_bias_lower = NA_real_, 
        calibration_bias_upper = NA_real_, calibration_direction_consistency = NA_real_, calibration_reproducibility = NA_real_, 
        calibration_relative_uncertainty = NA_real_, calibration_reason = "not_assessable", stringsAsFactors = FALSE)
    if (length(observed_value) != 1 || !is.finite(observed_value)) {
        stop("observed_value must be one finite numeric value.", call. = FALSE)
    }
    if (n_iterations < min_iterations) {
        base_output$calibration_reason <- "few_iterations"
        return(base_output)
    }
    deltas <- perturbed_values - observed_value
    bias <- median(deltas, na.rm = TRUE)
    robust_uncertainty <- mad(deltas, constant = 1.4826, na.rm = TRUE)
    nonzero <- deltas[deltas != 0]
    if (!length(nonzero)) {
        base_output$calibration_status <- "not_supported"
        base_output$calibration_reason <- "no_directional_bias"
        base_output$calibrated_value <- observed_value
        base_output$calibration_lower <- observed_value
        base_output$calibration_upper <- observed_value
        return(base_output)
    }
    positive_fraction <- mean(nonzero > 0)
    negative_fraction <- mean(nonzero < 0)
    direction_consistency <- max(positive_fraction, negative_fraction)
    direction <- if (positive_fraction >= negative_fraction) {
        "down"
    }
    else {
        "up"
    }
    set.seed(20260911)
    bootstrap_bias <- replicate(B, median(sample(deltas, replace = TRUE), na.rm = TRUE))
    bias_ci <- quantile(bootstrap_bias, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    scale <- max(abs(observed_value), mad(perturbed_values, constant = 1.4826, na.rm = TRUE), 1e-06)
    relative_uncertainty <- robust_uncertainty/scale
    reproducibility_score <- sparq_bound01(1 - median(abs(deltas), na.rm = TRUE)/scale)
    supported <- (!(bias_ci[1] <= 0 && bias_ci[2] >= 0) && direction_consistency >= min_direction_consistency && 
        reproducibility_score >= min_reproducibility && relative_uncertainty <= max_relative_uncertainty)
    if (!supported) {
        reason_parts <- character()
        if (bias_ci[1] <= 0 && bias_ci[2] >= 0) {
            reason_parts <- c(reason_parts, "bias_interval_overlaps_zero")
        }
        if (direction_consistency < min_direction_consistency) {
            reason_parts <- c(reason_parts, "inconsistent_direction")
        }
        if (reproducibility_score < min_reproducibility) {
            reason_parts <- c(reason_parts, "low_reproducibility")
        }
        if (relative_uncertainty > max_relative_uncertainty) {
            reason_parts <- c(reason_parts, "high_uncertainty")
        }
        base_output$calibration_reference <- reference_type
        base_output$calibration_reason <- paste(reason_parts, collapse = ";")
        base_output$calibration_bias <- bias
        base_output$calibration_bias_lower <- bias_ci[1]
        base_output$calibration_bias_upper <- bias_ci[2]
        base_output$calibration_direction_consistency <- direction_consistency
        base_output$calibration_reproducibility <- reproducibility_score
        base_output$calibration_relative_uncertainty <- relative_uncertainty
        return(base_output)
    }
    policy <- sparq_calibration_policy(reference_type = reference_type, correction_supported = TRUE, 
        apply_correction = apply_correction)
    out <- data.frame(calibration_status = policy$calibration_status, calibration_reference = reference_type, 
        calibration_direction = direction, calibrated_value = observed_value, apply_correction = policy$apply_correction, 
        calibration_lower = observed_value - bias_ci[2], calibration_upper = observed_value - bias_ci[1], 
        calibration_bias = bias, calibration_bias_lower = bias_ci[1], calibration_bias_upper = bias_ci[2], 
        calibration_direction_consistency = direction_consistency, calibration_reproducibility = reproducibility_score, 
        calibration_relative_uncertainty = relative_uncertainty, calibration_reason = policy$calibration_reason, 
        stringsAsFactors = FALSE)
    if (policy$apply_correction) {
        out$calibrated_value <- observed_value - bias
    }
    else {
        out$calibrated_value <- observed_value
    }
    out
}

sparq_calibration_estimators <- function(theta_full, theta_pert) {
    delta <- theta_pert - theta_full
    delta <- delta[is.finite(delta)]
    if (length(delta) < 3) {
        return(c(median_bias = NA_real_, trimmed_bias = NA_real_, severity_regression = NA_real_, empirical_bayes = NA_real_, 
            signed_quantile = NA_real_))
    }
    median_bias <- stats::median(delta, na.rm = TRUE)
    trimmed_bias <- mean(delta, trim = 0.2, na.rm = TRUE)
    severity <- abs(delta)
    severity_regression <- if (stats::sd(severity, na.rm = TRUE) > 0) {
        fit <- stats::lm(delta ~ severity)
        as.numeric(stats::predict(fit, newdata = data.frame(severity = stats::median(severity, na.rm = TRUE))))
    }
    else {
        median_bias
    }
    between_variance <- stats::var(delta, na.rm = TRUE)
    within_variance <- stats::mad(delta, constant = 1.4826, na.rm = TRUE)^2
    eb_weight <- if (is.finite(between_variance + within_variance) && between_variance + within_variance > 
        0) {
        between_variance/(between_variance + within_variance)
    }
    else {
        0
    }
    empirical_bayes <- eb_weight * median_bias
    signed_quantile <- if (median_bias < 0) {
        stats::quantile(delta, probs = 0.25, na.rm = TRUE, names = FALSE)
    }
    else {
        stats::quantile(delta, probs = 0.75, na.rm = TRUE, names = FALSE)
    }
    c(median_bias = median_bias, trimmed_bias = trimmed_bias, severity_regression = severity_regression, 
        empirical_bayes = empirical_bayes, signed_quantile = signed_quantile)
}

sparq_decide <- function(calibration_table) {
    needed <- c("result_id", "theta_full", "theta_calibrated", "correction", "decision")
    missing <- setdiff(needed, colnames(calibration_table))
    if (length(missing)) {
        stop("Missing required calibration columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    calibration_table
}

sparq_fit_builtin_calibration <- function(benchmark, train_fraction = 0.7, seed = 1) {
    required <- c("support_quality", "truth")
    missing <- setdiff(required, names(benchmark))
    if (length(missing)) {
        stop("Benchmark is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    benchmark <- benchmark[is.finite(benchmark$support_quality) & !is.na(benchmark$truth), , drop = FALSE]
    benchmark$truth <- as.integer(benchmark$truth)
    set.seed(seed)
    train_index <- sample(seq_len(nrow(benchmark)), size = floor(train_fraction * nrow(benchmark)))
    train <- benchmark[train_index, , drop = FALSE]
    test <- benchmark[-train_index, , drop = FALSE]
    fit <- glm(truth ~ support_quality, data = train, family = binomial())
    threshold_grid <- seq(0.01, 0.99, length.out = 500)
    train_performance <- do.call(rbind, lapply(threshold_grid, function(threshold) {
        predicted <- train$support_quality >= threshold
        sensitivity <- mean(predicted[train$truth == 1])
        specificity <- mean(!predicted[train$truth == 0])
        data.frame(threshold = threshold, sensitivity = sensitivity, specificity = specificity, balanced_accuracy = (sensitivity + 
            specificity)/2)
    }))
    best_threshold <- train_performance[which.max(train_performance$balanced_accuracy), "threshold"]
    test_probability <- predict(fit, newdata = test, type = "response")
    test_prediction <- test_probability >= 0.5
    test_support_prediction <- test$support_quality >= best_threshold
    test_sensitivity <- mean(test_support_prediction[test$truth == 1])
    test_specificity <- mean(!test_support_prediction[test$truth == 0])
    test_balanced_accuracy <- mean(c(test_sensitivity, test_specificity))
    list(model = fit, threshold = best_threshold, train_performance = train_performance, heldout_sensitivity = test_sensitivity, 
        heldout_specificity = test_specificity, heldout_balanced_accuracy = test_balanced_accuracy, heldout_predictions = data.frame(support_quality = test$support_quality, 
            truth = test$truth, predicted_probability = test_probability, predicted_class = test_prediction, 
            predicted_support = test_support_prediction))
}

sparq_make_truth_benchmark <- function(n_stable = 200, n_unstable = 200, seed = 1) {
    set.seed(seed)
    stable_support <- rbeta(n_stable, shape1 = 18, shape2 = 2)
    unstable_support <- rbeta(n_unstable, shape1 = 5, shape2 = 5)
    data.frame(result_id = paste0("benchmark_", seq_len(n_stable + n_unstable)), support_quality = c(stable_support, 
        unstable_support), truth = c(rep(1, n_stable), rep(0, n_unstable)), stringsAsFactors = FALSE)
}

sparq_pipeline <- function(data, min_iterations = 50, min_direction_consistency = 0.7, max_sign_flip_rate = 0.3, 
    max_correction_fraction = 0.5, estimator_aggregation = c("median", "mean")) {
    estimator_aggregation <- match.arg(estimator_aggregation)
    sparq_validate_results(data)
    reproducibility <- sparq_reproducibility_feature(data)
    bias <- sparq_bias(data)
    support <- sparq_support(bias_table = bias, reproducibility_table = reproducibility, min_iterations = min_iterations, 
        min_direction_consistency = min_direction_consistency, max_sign_flip_rate = max_sign_flip_rate)
    calibration <- sparq_calibrate(data = data, support_table = support, estimator_aggregation = estimator_aggregation, 
        max_correction_fraction = max_correction_fraction)
    list(reproducibility = reproducibility, bias = bias, support = support, calibration = calibration)
}

sparq_reproducibility_feature <- function(data, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", scenario_col = "scenario", iteration_col = "iteration") {
    sparq_validate_results(data, required_cols = c(result_id_col, theta_full_col, theta_pert_col, scenario_col, 
        iteration_col))
    ids <- unique(as.character(data[[result_id_col]]))
    out <- lapply(ids, function(id) {
        d <- data[as.character(data[[result_id_col]]) == id, , drop = FALSE]
        theta_full <- d[[theta_full_col]]
        theta_pert <- d[[theta_pert_col]]
        delta <- theta_pert - theta_full
        data.frame(result_id = id, n_perturbations = nrow(d), n_scenarios = length(unique(as.character(d[[scenario_col]]))), 
            pearson = sparq_safe_cor(theta_full, theta_pert, method = "pearson"), spearman = sparq_safe_cor(theta_full, 
                theta_pert, method = "spearman"), rmse = sqrt(mean(delta^2, na.rm = TRUE)), mae = mean(abs(delta), 
                na.rm = TRUE), cv_perturbed = stats::sd(theta_pert, na.rm = TRUE)/max(abs(stats::median(theta_pert, 
                na.rm = TRUE)), .Machine$double.eps), sign_agreement = mean(sign(theta_full) == sign(theta_pert), 
                na.rm = TRUE), stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

