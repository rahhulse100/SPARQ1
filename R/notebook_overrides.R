sparq_apply_builtin_calibration <-
function (support_quality, calibration_model) 
{
    probability <- predict(calibration_model$model, newdata = data.frame(support_quality = support_quality), 
        type = "response")
    support_class <- ifelse(support_quality >= calibration_model$threshold, 
        "supported", "insufficient")
    data.frame(support_quality = support_quality, support_probability = probability, 
        support_threshold = calibration_model$threshold, support_class = support_class, 
        stringsAsFactors = FALSE)
}
sparq_assess <-
function (full_result, perturbed_results, comparator, result_id = "result", 
    reference_scale = 1, ...) 
{
    if (!length(perturbed_results)) {
        stop("perturbed_results is empty.")
    }
    comparison_list <- lapply(perturbed_results, function(x) {
        comparator(full_result, x)
    })
    comparison_table <- do.call(rbind, comparison_list)
    rownames(comparison_table) <- NULL
    summary <- sparq_support_from_comparisons(comparison_table = comparison_table, 
        n_iterations = nrow(comparison_table), reference_scale = reference_scale, 
        ...)
    summary$result_id <- result_id
    list(result_id = result_id, comparisons = comparison_table, 
        summary = summary)
}
sparq_assess_cohort <-
function (sample_data, analysis_function, comparator, stress_model = "uniform_random", 
    retention = 0.75, n_iterations = 20, x_col = NULL, y_col = NULL, 
    custom_function = NULL, seed = 1, ...) 
{
    if (is.null(names(sample_data))) {
        names(sample_data) <- paste0("sample_", seq_along(sample_data))
    }
    results <- lapply(names(sample_data), function(sample_id) {
        sparq_run_with_stress_model(data = sample_data[[sample_id]], 
            analysis_function = analysis_function, comparator = comparator, 
            stress_model = stress_model, retention = retention, 
            n_iterations = n_iterations, x_col = x_col, y_col = y_col, 
            custom_function = custom_function, result_id = sample_id, 
            seed = seed, ...)
    })
    output <- do.call(rbind, results)
    rownames(output) <- NULL
    output
}
sparq_assess_precomputed <-
function (full_results, perturbed_results, comparator, reference_scales = NULL, 
    ...) 
{
    if (is.null(names(full_results))) {
        names(full_results) <- paste0("sample_", seq_along(full_results))
    }
    if (is.null(names(perturbed_results))) {
        stop("perturbed_results must be a named list.")
    }
    sample_ids <- intersect(names(full_results), names(perturbed_results))
    if (!length(sample_ids)) {
        stop("No matching sample IDs between full_results and ", 
            "perturbed_results.")
    }
    summaries <- lapply(sample_ids, function(sample_id) {
        full_result <- full_results[[sample_id]]
        perturbed <- perturbed_results[[sample_id]]
        if (!length(perturbed)) {
            stop("No perturbed results for sample: ", sample_id)
        }
        reference_scale <- 1
        if (!is.null(reference_scales) && sample_id %in% names(reference_scales)) {
            reference_scale <- reference_scales[[sample_id]]
        }
        assessed <- sparq_assess(full_result = full_result, perturbed_results = perturbed, 
            comparator = comparator, result_id = sample_id, reference_scale = reference_scale, 
            ...)
        assessed$summary
    })
    output <- do.call(rbind, summaries)
    rownames(output) <- NULL
    output$assessment_mode <- "precomputed"
    output
}
sparq_bias <-
function (data, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert") 
{
    sparq_validate_results(data, required_cols = c(result_id_col, 
        theta_full_col, theta_pert_col))
    ids <- unique(as.character(data[[result_id_col]]))
    out <- lapply(ids, function(id) {
        d <- data[as.character(data[[result_id_col]]) == id, 
            , drop = FALSE]
        theta_full <- stats::median(d[[theta_full_col]], na.rm = TRUE)
        delta <- d[[theta_pert_col]] - d[[theta_full_col]]
        delta <- delta[is.finite(delta)]
        bias <- stats::median(delta, na.rm = TRUE)
        uncertainty <- stats::mad(delta, constant = 1.4825999999999999, 
            na.rm = TRUE)
        direction_consistency <- sparq_mode_sign_fraction(delta)
        data.frame(result_id = id, theta_full = theta_full, n_perturbations = length(delta), 
            bias = bias, uncertainty = uncertainty, direction_consistency = direction_consistency, 
            sign_flip_rate = 1 - direction_consistency, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}
sparq_bias_ranked_feature <-
function (data, ...) 
{
    dt <- sparq_prepare_ranked_feature_input(data, ...)
    dt[, `:=`(delta, theta_pert - theta_full)]
    dt[, `:=`(abs_delta, abs(delta))]
    dt[, {
        mb <- stats::median(delta, na.rm = TRUE)
        tb <- sparq_trimmed_median_internal(delta)
        sr <- sparq_severity_regression_internal(loss_fraction, 
            delta)
        eb <- sparq_empirical_bayes_internal(delta)
        sq <- sparq_signed_quantile_internal(delta)
        estimators <- c(median_bias = mb, trimmed_bias = tb, 
            severity_regression = sr, empirical_bayes = eb, signed_quantile = sq)
        estimators <- estimators[is.finite(estimators)]
        calibration_bias <- if (length(estimators) > 0) {
            stats::median(estimators)
        }
        else {
            NA_real_
        }
        estimator_signs <- sign(estimators[abs(estimators) > 
            9.9999999999999998e-13])
        estimator_direction_agreement <- if (length(estimator_signs) > 
            0) {
            mean(estimator_signs == sign(stats::median(estimators)))
        }
        else {
            0
        }
        direction_consistency <- if (is.finite(calibration_bias) && 
            abs(calibration_bias) > 9.9999999999999998e-13) {
            mean(sign(delta) == sign(calibration_bias), na.rm = TRUE)
        }
        else {
            0
        }
        sign_flip_rate <- {
            base_sign <- sign(theta_full[1])
            if (!is.finite(base_sign) || base_sign == 0) {
                NA_real_
            }
            else {
                mean(sign(theta_pert) != base_sign, na.rm = TRUE)
            }
        }
        uncertainty_mad <- stats::median(abs(delta - stats::median(delta, 
            na.rm = TRUE)), na.rm = TRUE)
        signal_to_uncertainty <- abs(calibration_bias)/(uncertainty_mad + 
            1e-08)
        relative_bias <- abs(calibration_bias)/(abs(theta_full[1]) + 
            1e-08)
        n_runs <- data.table::uniqueN(paste(scenario, iteration, 
            sep = "__"))
        list(theta_full = theta_full[1], rank_full = rank_full[1], 
            n_runs = n_runs, calibration_bias = calibration_bias, 
            median_bias = mb, trimmed_bias = tb, severity_regression = sr, 
            empirical_bayes = eb, signed_quantile = sq, direction_consistency = direction_consistency, 
            estimator_direction_agreement = estimator_direction_agreement, 
            sign_flip_rate = sign_flip_rate, uncertainty_mad = uncertainty_mad, 
            signal_to_uncertainty = signal_to_uncertainty, relative_bias = relative_bias)
    }, by = .(cohort, sample_id, result_id)]
}
sparq_bootstrap_ci <-
function (x, statistic = median, B = 2000, conf = 0.94999999999999996, 
    seed = 1) 
{
    x <- x[is.finite(x)]
    if (length(x) < 3) {
        return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    set.seed(seed)
    boot_values <- replicate(B, statistic(sample(x, replace = TRUE)))
    alpha <- (1 - conf)/2
    c(estimate = statistic(x), lower = as.numeric(quantile(boot_values, 
        alpha, na.rm = TRUE)), upper = as.numeric(quantile(boot_values, 
        1 - alpha, na.rm = TRUE)))
}
sparq_bound01 <-
function (x) 
{
    pmax(0, pmin(1, x))
}
sparq_calibrate <-
function (data, support_table, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", estimator_aggregation = c("median", 
        "mean"), max_correction_fraction = 0.5, min_estimator_direction_agreement = 0.59999999999999998) 
{
    estimator_aggregation <- match.arg(estimator_aggregation)
    sparq_validate_results(data, required_cols = c(result_id_col, 
        theta_full_col, theta_pert_col))
    needed <- c("result_id", "support_quality", "insufficient_support")
    missing <- setdiff(needed, colnames(support_table))
    if (length(missing)) {
        stop("Missing required support columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    ids <- unique(as.character(data[[result_id_col]]))
    out <- lapply(ids, function(id) {
        d <- data[as.character(data[[result_id_col]]) == id, 
            , drop = FALSE]
        s <- support_table[as.character(support_table$result_id) == 
            id, , drop = FALSE]
        theta_full <- stats::median(d[[theta_full_col]], na.rm = TRUE)
        est <- sparq_calibration_estimators(theta_full = d[[theta_full_col]], 
            theta_pert = d[[theta_pert_col]])
        est_finite <- est[is.finite(est)]
        if (!nrow(s) || !length(est_finite)) {
            return(data.frame(result_id = id, theta_full = theta_full, 
                theta_calibrated = theta_full, correction = 0, 
                calibration_bias = NA_real_, estimator_direction_agreement = NA_real_, 
                support_quality = NA_real_, decision = "insufficient_support", 
                t(est), check.names = FALSE, stringsAsFactors = FALSE))
        }
        main_sign <- sign(stats::median(est_finite, na.rm = TRUE))
        direction_agreement <- mean(sign(est_finite) == main_sign)
        calibration_bias <- if (estimator_aggregation == "median") {
            stats::median(est_finite, na.rm = TRUE)
        }
        else {
            mean(est_finite, na.rm = TRUE)
        }
        correction <- -sparq_bound01(s$support_quality[1]) * 
            calibration_bias
        bound <- max_correction_fraction * max(abs(theta_full), 
            .Machine$double.eps)
        correction <- max(-bound, min(bound, correction))
        unsafe <- isTRUE(s$insufficient_support[1]) || is.na(direction_agreement) || 
            direction_agreement < min_estimator_direction_agreement
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
        data.frame(result_id = id, theta_full = theta_full, theta_calibrated = theta_calibrated, 
            correction = correction, calibration_bias = calibration_bias, 
            estimator_direction_agreement = direction_agreement, 
            support_quality = s$support_quality[1], decision = decision, 
            t(est), check.names = FALSE, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}
sparq_calibrate_ranked_feature <-
function (bias_table, expected_runs = 150, max_correction_fraction = 0.5, 
    min_runs = 30, min_support_quality = 0.34999999999999998, 
    min_direction_consistency = 0.65000000000000002, max_sign_flip_rate = 0.25) 
{
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("SPARQ requires data.table.", call. = FALSE)
    }
    x <- data.table::as.data.table(bias_table)
    required <- c("cohort", "sample_id", "result_id", "theta_full", 
        "rank_full", "n_runs", "calibration_bias", "direction_consistency", 
        "estimator_direction_agreement", "sign_flip_rate", "signal_to_uncertainty")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
        stop("SPARQ bias table missing required columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    x[, `:=`(run_support, pmin(n_runs/expected_runs, 1))]
    x[, `:=`(direction_support, pmax(0, pmin(direction_consistency, 
        1)))]
    x[, `:=`(estimator_support, pmax(0, pmin(estimator_direction_agreement, 
        1)))]
    x[, `:=`(uncertainty_support, pmax(0, pmin(signal_to_uncertainty/2, 
        1)))]
    x[, `:=`(flip_support, 1 - pmax(0, pmin(ifelse(is.na(sign_flip_rate), 
        1, sign_flip_rate), 1)))]
    x[, `:=`(support_quality, rowMeans(cbind(run_support, direction_support, 
        estimator_support, uncertainty_support, flip_support), 
        na.rm = TRUE))]
    x[, `:=`(lambda, support_quality * direction_support * estimator_support)]
    x[, `:=`(correction_raw, -lambda * calibration_bias)]
    x[, `:=`(correction_bound, max_correction_fraction * abs(theta_full))]
    x[, `:=`(correction, pmax(-correction_bound, pmin(correction_bound, 
        correction_raw)))]
    x[, `:=`(theta_calibrated, theta_full + correction)]
    x[, `:=`(decision, data.table::fcase(n_runs < min_runs, "insufficient_support", 
        support_quality < min_support_quality, "insufficient_support", 
        direction_consistency < min_direction_consistency, "insufficient_support", 
        !is.na(sign_flip_rate) & sign_flip_rate > max_sign_flip_rate, 
        "insufficient_support", abs(correction) < 9.9999999999999998e-13, 
        "hold", correction > 0, "calibrate_up", correction < 
            0, "calibrate_down", default = "hold"))]
    x[decision == "insufficient_support", `:=`(correction, 0)]
    x[decision == "insufficient_support", `:=`(theta_calibrated, 
        theta_full)]
    x[, `:=`(rank_calibrated, data.table::frank(-theta_calibrated, 
        ties.method = "min")), by = .(cohort, sample_id)]
    x[, `:=`(rank_shift_after_calibration, rank_calibrated - 
        rank_full)]
    keep <- c("cohort", "sample_id", "result_id", "theta_full", 
        "theta_calibrated", "correction", "calibration_bias", 
        "rank_full", "rank_calibrated", "rank_shift_after_calibration", 
        "decision", "support_quality", "n_runs", "direction_consistency", 
        "estimator_direction_agreement", "sign_flip_rate", "uncertainty_mad", 
        "signal_to_uncertainty", "relative_bias", "median_bias", 
        "trimmed_bias", "severity_regression", "empirical_bayes", 
        "signed_quantile")
    as.data.frame(x[, ..keep])
}
sparq_calibrate_scalar <-
function (observed_value, perturbed_values, reference_type = c("none", 
    "simulation", "replicate", "pathology"), apply_correction = FALSE, 
    ...) 
{
    reference_type <- match.arg(reference_type)
    raw <- sparq_calibrate_scalar_raw(observed_value = observed_value, 
        perturbed_values = perturbed_values, ...)
    correction_supported <- raw$calibration_status == "supported"
    policy <- sparq_calibration_policy(reference_type = reference_type, 
        correction_supported = correction_supported, apply_correction = apply_correction)
    out <- raw
    out$calibration_reference <- reference_type
    out$apply_correction <- policy$apply_correction
    if (correction_supported && policy$apply_correction) {
        out$calibration_status <- "validated"
        out$calibration_direction <- raw$calibration_direction
        out$calibrated_value <- raw$calibrated_value
        out$calibration_reason <- policy$calibration_reason
    }
    else if (correction_supported && reference_type == "none") {
        out$calibration_status <- "conditional"
        out$calibration_direction <- raw$calibration_direction
        out$calibrated_value <- observed_value
        out$apply_correction <- FALSE
        out$calibration_reason <- policy$calibration_reason
    }
    else {
        out$calibration_status <- "not_supported"
        out$calibration_direction <- "none"
        out$calibrated_value <- observed_value
        out$apply_correction <- FALSE
        out$calibration_reason <- policy$calibration_reason
    }
    out
}
sparq_calibrate_scalar_raw <-
function (observed_value, perturbed_values, reference_type = c("none", 
    "simulation", "replicate", "pathology"), apply_correction = FALSE, 
    alpha = 0.050000000000000003, B = 2000, min_iterations = 20, 
    min_direction_consistency = 0.80000000000000004, min_reproducibility = 0.69999999999999996, 
    max_relative_uncertainty = 0.25) 
{
    reference_type <- match.arg(reference_type)
    observed_value <- as.numeric(observed_value)
    perturbed_values <- as.numeric(unlist(perturbed_values))
    perturbed_values <- perturbed_values[is.finite(perturbed_values)]
    n_iterations <- length(perturbed_values)
    base_output <- data.frame(calibration_status = "not_supported", 
        calibration_reference = reference_type, calibration_direction = "none", 
        calibrated_value = observed_value, apply_correction = FALSE, 
        calibration_lower = NA_real_, calibration_upper = NA_real_, 
        calibration_bias = NA_real_, calibration_bias_lower = NA_real_, 
        calibration_bias_upper = NA_real_, calibration_direction_consistency = NA_real_, 
        calibration_reproducibility = NA_real_, calibration_relative_uncertainty = NA_real_, 
        calibration_reason = "not_assessable", stringsAsFactors = FALSE)
    if (length(observed_value) != 1 || !is.finite(observed_value)) {
        stop("observed_value must be one finite numeric value.", 
            call. = FALSE)
    }
    if (n_iterations < min_iterations) {
        base_output$calibration_reason <- "few_iterations"
        return(base_output)
    }
    deltas <- perturbed_values - observed_value
    bias <- median(deltas, na.rm = TRUE)
    robust_uncertainty <- mad(deltas, constant = 1.4825999999999999, 
        na.rm = TRUE)
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
    bootstrap_bias <- replicate(B, median(sample(deltas, replace = TRUE), 
        na.rm = TRUE))
    bias_ci <- quantile(bootstrap_bias, probs = c(alpha/2, 1 - 
        alpha/2), na.rm = TRUE)
    scale <- max(abs(observed_value), mad(perturbed_values, constant = 1.4825999999999999, 
        na.rm = TRUE), 9.9999999999999995e-07)
    relative_uncertainty <- robust_uncertainty/scale
    reproducibility_score <- sparq_bound01(1 - median(abs(deltas), 
        na.rm = TRUE)/scale)
    supported <- (!(bias_ci[1] <= 0 && bias_ci[2] >= 0) && direction_consistency >= 
        min_direction_consistency && reproducibility_score >= 
        min_reproducibility && relative_uncertainty <= max_relative_uncertainty)
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
        base_output$calibration_reason <- paste(reason_parts, 
            collapse = ";")
        base_output$calibration_bias <- bias
        base_output$calibration_bias_lower <- bias_ci[1]
        base_output$calibration_bias_upper <- bias_ci[2]
        base_output$calibration_direction_consistency <- direction_consistency
        base_output$calibration_reproducibility <- reproducibility_score
        base_output$calibration_relative_uncertainty <- relative_uncertainty
        return(base_output)
    }
    policy <- sparq_calibration_policy(reference_type = reference_type, 
        correction_supported = TRUE, apply_correction = apply_correction)
    out <- data.frame(calibration_status = policy$calibration_status, 
        calibration_reference = reference_type, calibration_direction = direction, 
        calibrated_value = observed_value, apply_correction = policy$apply_correction, 
        calibration_lower = observed_value - bias_ci[2], calibration_upper = observed_value - 
            bias_ci[1], calibration_bias = bias, calibration_bias_lower = bias_ci[1], 
        calibration_bias_upper = bias_ci[2], calibration_direction_consistency = direction_consistency, 
        calibration_reproducibility = reproducibility_score, 
        calibration_relative_uncertainty = relative_uncertainty, 
        calibration_reason = policy$calibration_reason, stringsAsFactors = FALSE)
    if (policy$apply_correction) {
        out$calibrated_value <- observed_value - bias
    }
    else {
        out$calibrated_value <- observed_value
    }
    out
}
sparq_calibrate_threshold <-
function (support_quality, truth) 
{
    ok <- is.finite(support_quality) & !is.na(truth)
    support_quality <- support_quality[ok]
    truth <- as.integer(truth[ok])
    thresholds <- seq(min(support_quality), max(support_quality), 
        length.out = 200)
    performance <- lapply(thresholds, function(threshold) {
        predicted <- support_quality >= threshold
        sensitivity <- mean(predicted[truth == 1])
        specificity <- mean(!predicted[truth == 0])
        data.frame(threshold = threshold, sensitivity = sensitivity, 
            specificity = specificity, balanced_accuracy = (sensitivity + 
                specificity)/2)
    })
    performance <- do.call(rbind, performance)
    performance[which.max(performance$balanced_accuracy), , drop = FALSE]
}
sparq_calibration_estimators <-
function (theta_full, theta_pert) 
{
    delta <- theta_pert - theta_full
    delta <- delta[is.finite(delta)]
    if (length(delta) < 3) {
        return(c(median_bias = NA_real_, trimmed_bias = NA_real_, 
            severity_regression = NA_real_, empirical_bayes = NA_real_, 
            signed_quantile = NA_real_))
    }
    median_bias <- stats::median(delta, na.rm = TRUE)
    trimmed_bias <- mean(delta, trim = 0.20000000000000001, na.rm = TRUE)
    severity <- abs(delta)
    severity_regression <- if (stats::sd(severity, na.rm = TRUE) > 
        0) {
        fit <- stats::lm(delta ~ severity)
        as.numeric(stats::predict(fit, newdata = data.frame(severity = stats::median(severity, 
            na.rm = TRUE))))
    }
    else {
        median_bias
    }
    between_variance <- stats::var(delta, na.rm = TRUE)
    within_variance <- stats::mad(delta, constant = 1.4825999999999999, 
        na.rm = TRUE)^2
    eb_weight <- if (is.finite(between_variance + within_variance) && 
        between_variance + within_variance > 0) {
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
    c(median_bias = median_bias, trimmed_bias = trimmed_bias, 
        severity_regression = severity_regression, empirical_bayes = empirical_bayes, 
        signed_quantile = signed_quantile)
}
sparq_calibration_policy <-
function (reference_type = c("none", "simulation", "replicate", 
    "pathology"), correction_supported, apply_correction = FALSE) 
{
    reference_type <- match.arg(reference_type)
    if (!is.logical(correction_supported) || length(correction_supported) != 
        1) {
        stop("correction_supported must be TRUE or FALSE.", call. = FALSE)
    }
    if (reference_type == "none") {
        return(list(calibration_status = ifelse(correction_supported, 
            "conditional", "not_supported"), calibration_reference = "none", 
            apply_correction = FALSE, calibration_reason = ifelse(correction_supported, 
                "Perturbation evidence supports a conditional correction, but no external reference is available.", 
                "No validated calibration evidence is available.")))
    }
    if (!correction_supported) {
        return(list(calibration_status = "not_supported", calibration_reference = reference_type, 
            apply_correction = FALSE, calibration_reason = "Reference data do not support correction."))
    }
    list(calibration_status = "validated", calibration_reference = reference_type, 
        apply_correction = apply_correction, calibration_reason = paste("Correction supported using", 
            reference_type, "reference data."))
}
sparq_compare_graph <-
function (full_edges, perturbed_edges) 
{
    full_edges <- unique(as.character(full_edges))
    perturbed_edges <- unique(as.character(perturbed_edges))
    union_edges <- union(full_edges, perturbed_edges)
    jaccard <- if (!length(union_edges)) {
        NA_real_
    }
    else {
        length(intersect(full_edges, perturbed_edges))/length(union_edges)
    }
    data.frame(edge_jaccard = jaccard, instability = 1 - jaccard)
}
sparq_compare_partition <-
function (full_labels, perturbed_labels) 
{
    full_labels <- as.character(full_labels)
    perturbed_labels <- as.character(perturbed_labels)
    ok <- !is.na(full_labels) & !is.na(perturbed_labels)
    full_labels <- full_labels[ok]
    perturbed_labels <- perturbed_labels[ok]
    tab <- table(full_labels, perturbed_labels)
    n <- sum(tab)
    if (n < 2) {
        return(data.frame(ari = NA_real_, instability = NA_real_))
    }
    choose2 <- function(x) {
        x * (x - 1)/2
    }
    nij <- sum(choose2(tab))
    ai <- sum(choose2(rowSums(tab)))
    bj <- sum(choose2(colSums(tab)))
    total <- choose2(n)
    expected <- if (total == 0) {
        0
    }
    else {
        ai * bj/total
    }
    max_index <- 0.5 * (ai + bj)
    ari <- if ((max_index - expected) == 0) {
        1
    }
    else {
        (nij - expected)/(max_index - expected)
    }
    data.frame(ari = ari, instability = 1 - ari)
}
sparq_compare_ranked <-
function (full_rank, perturbed_rank, k = 100) 
{
    full_rank <- as.character(full_rank)
    perturbed_rank <- as.character(perturbed_rank)
    full_top <- unique(full_rank[seq_len(min(k, length(full_rank)))])
    pert_top <- unique(perturbed_rank[seq_len(min(k, length(perturbed_rank)))])
    union_ids <- union(full_top, pert_top)
    jaccard <- if (!length(union_ids)) {
        NA_real_
    }
    else {
        length(intersect(full_top, pert_top))/length(union_ids)
    }
    data.frame(spearman = suppressWarnings(cor(match(full_rank, 
        union_ids), match(perturbed_rank, union_ids), method = "spearman", 
        use = "complete.obs")), jaccard = jaccard, replacement = 1 - 
        jaccard, instability = 1 - jaccard, stringsAsFactors = FALSE)
}
sparq_compare_scalar <-
function (full_result, perturbed_result) 
{
    delta <- as.numeric(perturbed_result) - as.numeric(full_result)
    data.frame(delta = delta, instability = abs(delta), stringsAsFactors = FALSE)
}
sparq_decide <-
function (calibration_table) 
{
    needed <- c("result_id", "theta_full", "theta_calibrated", 
        "correction", "decision")
    missing <- setdiff(needed, colnames(calibration_table))
    if (length(missing)) {
        stop("Missing required calibration columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    calibration_table
}
sparq_effect_scale <-
function (reference, floor = 9.9999999999999995e-07) 
{
    reference <- reference[is.finite(reference)]
    if (!length(reference)) {
        return(floor)
    }
    max(stats::mad(reference, na.rm = TRUE), floor)
}
sparq_empirical_bayes_internal <-
function (x, prior_strength = 20) 
{
    x <- x[is.finite(x)]
    if (length(x) == 0) 
        return(NA_real_)
    stats::median(x) * (length(x)/(length(x) + prior_strength))
}
sparq_finalize_support <-
function (support_table, support_threshold) 
{
    support_table$support_interval_low <- sparq_bound01(1 - support_table$instability_ci_high)
    support_table$support_interval_high <- sparq_bound01(1 - 
        support_table$instability_ci_low)
    support_table$support_interval <- sprintf("%.3f–%.3f", 
        support_table$support_interval_low, support_table$support_interval_high)
    support_table$support_class <- ifelse(support_table$support_interval_low >= 
        support_threshold, "supported", ifelse(support_table$support_interval_high < 
        support_threshold, "insufficient", "borderline"))
    support_table$support_reason <- ifelse(support_table$support_class == 
        "supported", "interval_above_threshold", ifelse(support_table$support_class == 
        "insufficient", "interval_below_threshold", "interval_overlaps_threshold"))
    rownames(support_table) <- NULL
    support_table
}
sparq_fit_builtin_calibration <-
function (benchmark, train_fraction = 0.69999999999999996, seed = 1) 
{
    required <- c("support_quality", "truth")
    missing <- setdiff(required, names(benchmark))
    if (length(missing)) {
        stop("Benchmark is missing: ", paste(missing, collapse = ", "), 
            call. = FALSE)
    }
    benchmark <- benchmark[is.finite(benchmark$support_quality) & 
        !is.na(benchmark$truth), , drop = FALSE]
    benchmark$truth <- as.integer(benchmark$truth)
    set.seed(seed)
    train_index <- sample(seq_len(nrow(benchmark)), size = floor(train_fraction * 
        nrow(benchmark)))
    train <- benchmark[train_index, , drop = FALSE]
    test <- benchmark[-train_index, , drop = FALSE]
    fit <- glm(truth ~ support_quality, data = train, family = binomial())
    threshold_grid <- seq(0.01, 0.98999999999999999, length.out = 500)
    train_performance <- do.call(rbind, lapply(threshold_grid, 
        function(threshold) {
            predicted <- train$support_quality >= threshold
            sensitivity <- mean(predicted[train$truth == 1])
            specificity <- mean(!predicted[train$truth == 0])
            data.frame(threshold = threshold, sensitivity = sensitivity, 
                specificity = specificity, balanced_accuracy = (sensitivity + 
                  specificity)/2)
        }))
    best_threshold <- train_performance[which.max(train_performance$balanced_accuracy), 
        "threshold"]
    test_probability <- predict(fit, newdata = test, type = "response")
    test_prediction <- test_probability >= 0.5
    test_support_prediction <- test$support_quality >= best_threshold
    test_sensitivity <- mean(test_support_prediction[test$truth == 
        1])
    test_specificity <- mean(!test_support_prediction[test$truth == 
        0])
    test_balanced_accuracy <- mean(c(test_sensitivity, test_specificity))
    list(model = fit, threshold = best_threshold, train_performance = train_performance, 
        heldout_sensitivity = test_sensitivity, heldout_specificity = test_specificity, 
        heldout_balanced_accuracy = test_balanced_accuracy, heldout_predictions = data.frame(support_quality = test$support_quality, 
            truth = test$truth, predicted_probability = test_probability, 
            predicted_class = test_prediction, predicted_support = test_support_prediction))
}
sparq_make_truth_benchmark <-
function (n_stable = 200, n_unstable = 200, seed = 1) 
{
    set.seed(seed)
    stable_support <- rbeta(n_stable, shape1 = 18, shape2 = 2)
    unstable_support <- rbeta(n_unstable, shape1 = 5, shape2 = 5)
    data.frame(result_id = paste0("benchmark_", seq_len(n_stable + 
        n_unstable)), support_quality = c(stable_support, unstable_support), 
        truth = c(rep(1, n_stable), rep(0, n_unstable)), stringsAsFactors = FALSE)
}
sparq_mode_sign_fraction <-
function (x) 
{
    x <- x[is.finite(x)]
    x <- sign(x[x != 0])
    if (length(x) == 0) 
        return(0)
    max(mean(x > 0), mean(x < 0))
}
sparq_pipeline <-
function (data, min_iterations = 50, min_direction_consistency = 0.69999999999999996, 
    max_sign_flip_rate = 0.29999999999999999, max_correction_fraction = 0.5, 
    estimator_aggregation = c("median", "mean")) 
{
    estimator_aggregation <- match.arg(estimator_aggregation)
    sparq_validate_results(data)
    reproducibility <- sparq_reproducibility_feature(data)
    bias <- sparq_bias(data)
    support <- sparq_support(bias_table = bias, reproducibility_table = reproducibility, 
        min_iterations = min_iterations, min_direction_consistency = min_direction_consistency, 
        max_sign_flip_rate = max_sign_flip_rate)
    calibration <- sparq_calibrate(data = data, support_table = support, 
        estimator_aggregation = estimator_aggregation, max_correction_fraction = max_correction_fraction)
    list(reproducibility = reproducibility, bias = bias, support = support, 
        calibration = calibration)
}
sparq_pipeline_ranked_feature <-
function (data, top_k = 100, expected_runs = 150, max_correction_fraction = 0.5, 
    min_runs = 30, min_support_quality = 0.34999999999999998, 
    min_direction_consistency = 0.65000000000000002, max_sign_flip_rate = 0.25, 
    ...) 
{
    input <- sparq_prepare_ranked_feature_input(data, ...)
    input_audit <- data.frame(n_rows = nrow(input), n_cohorts = length(unique(input$cohort)), 
        n_samples = length(unique(paste(input$cohort, input$sample_id, 
            sep = "__"))), n_results = length(unique(input$result_id)), 
        n_scenarios = length(unique(input$scenario)), n_iterations = length(unique(input$iteration)), 
        stringsAsFactors = FALSE)
    turnover <- sparq_rank_turnover(input, top_k = top_k)
    bias <- sparq_bias_ranked_feature(input)
    calibration <- sparq_calibrate_ranked_feature(bias, expected_runs = expected_runs, 
        max_correction_fraction = max_correction_fraction, min_runs = min_runs, 
        min_support_quality = min_support_quality, min_direction_consistency = min_direction_consistency, 
        max_sign_flip_rate = max_sign_flip_rate)
    sample_summary <- sparq_sample_fragility_ranked_feature(turnover)
    list(input_audit = input_audit, turnover = turnover, bias = as.data.frame(bias), 
        calibration = calibration, sample_summary = sample_summary)
}
sparq_prepare_ranked_feature_input <-
function (data, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", rank_full_col = "rank_full", 
    rank_pert_col = "rank_pert", scenario_col = "scenario", iteration_col = "iteration", 
    loss_fraction_col = "loss_fraction", cohort_col = "cohort", 
    sample_id_col = "sample_id") 
{
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("SPARQ requires data.table.", call. = FALSE)
    }
    x <- data.table::as.data.table(data)
    required_source <- c(result_id_col, theta_full_col, theta_pert_col, 
        rank_full_col, rank_pert_col, scenario_col, iteration_col)
    missing <- setdiff(required_source, names(x))
    if (length(missing) > 0) {
        stop("SPARQ input missing required columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    out <- data.table::data.table(result_id = as.character(x[[result_id_col]]), 
        theta_full = as.numeric(x[[theta_full_col]]), theta_pert = as.numeric(x[[theta_pert_col]]), 
        rank_full = as.numeric(x[[rank_full_col]]), rank_pert = as.numeric(x[[rank_pert_col]]), 
        scenario = as.character(x[[scenario_col]]), iteration = as.integer(x[[iteration_col]]))
    if (!is.null(loss_fraction_col) && loss_fraction_col %in% 
        names(x)) {
        out[, `:=`(loss_fraction, as.numeric(x[[loss_fraction_col]]))]
    }
    else {
        out[, `:=`(loss_fraction, NA_real_)]
    }
    if (!is.null(cohort_col) && cohort_col %in% names(x)) {
        out[, `:=`(cohort, as.character(x[[cohort_col]]))]
    }
    else {
        out[, `:=`(cohort, "cohort_1")]
    }
    if (!is.null(sample_id_col) && sample_id_col %in% names(x)) {
        out[, `:=`(sample_id, as.character(x[[sample_id_col]]))]
    }
    else {
        out[, `:=`(sample_id, "sample_1")]
    }
    out <- out[!is.na(result_id) & result_id != "" & is.finite(theta_full) & 
        is.finite(theta_pert) & is.finite(rank_full) & is.finite(rank_pert) & 
        !is.na(scenario) & scenario != "" & !is.na(iteration)]
    if (nrow(out) == 0) {
        stop("SPARQ input has zero usable rows after validation.", 
            call. = FALSE)
    }
    data.table::setcolorder(out, c("cohort", "sample_id", "result_id", 
        "theta_full", "theta_pert", "rank_full", "rank_pert", 
        "scenario", "iteration", "loss_fraction"))
    out[]
}
sparq_rank_turnover <-
function (data, top_k = 100, ...) 
{
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
        list(top_k = top_k, n_features_compared = .N, spearman_rank_correlation = sparq_safe_cor_internal(rank_full, 
            rank_pert, "spearman"), topk_jaccard = if (top_union == 
            0) NA_real_ else top_intersect/top_union, topk_lost = top_lost, 
            topk_gained = top_gained, topk_replaced_fraction = if (length(full_top) == 
                0) NA_real_ else top_lost/length(full_top))
    }, by = .(cohort, sample_id, scenario, iteration)]
}
sparq_reliability_curve <-
function (support_quality, true_error, bins = 10) 
{
    ok <- is.finite(support_quality) & is.finite(true_error)
    df <- data.frame(support_quality = support_quality[ok], true_error = true_error[ok])
    breaks <- seq(0, 1, length.out = bins + 1)
    df$bin <- cut(df$support_quality, breaks = breaks, include.lowest = TRUE)
    df <- df[!is.na(df$bin), , drop = FALSE]
    split_df <- split(df, droplevels(df$bin))
    out <- lapply(split_df, function(x) {
        data.frame(support_bin = as.character(x$bin[1]), n = nrow(x), 
            median_support = median(x$support_quality), median_true_error = median(x$true_error), 
            q25_true_error = as.numeric(quantile(x$true_error, 
                0.25)), q75_true_error = as.numeric(quantile(x$true_error, 
                0.75)))
    })
    do.call(rbind, out)
}
sparq_reproducibility_feature <-
function (data, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", scenario_col = "scenario", 
    iteration_col = "iteration") 
{
    sparq_validate_results(data, required_cols = c(result_id_col, 
        theta_full_col, theta_pert_col, scenario_col, iteration_col))
    ids <- unique(as.character(data[[result_id_col]]))
    out <- lapply(ids, function(id) {
        d <- data[as.character(data[[result_id_col]]) == id, 
            , drop = FALSE]
        theta_full <- d[[theta_full_col]]
        theta_pert <- d[[theta_pert_col]]
        delta <- theta_pert - theta_full
        data.frame(result_id = id, n_perturbations = nrow(d), 
            n_scenarios = length(unique(as.character(d[[scenario_col]]))), 
            pearson = sparq_safe_cor(theta_full, theta_pert, 
                method = "pearson"), spearman = sparq_safe_cor(theta_full, 
                theta_pert, method = "spearman"), rmse = sqrt(mean(delta^2, 
                na.rm = TRUE)), mae = mean(abs(delta), na.rm = TRUE), 
            cv_perturbed = stats::sd(theta_pert, na.rm = TRUE)/max(abs(stats::median(theta_pert, 
                na.rm = TRUE)), .Machine$double.eps), sign_agreement = mean(sign(theta_full) == 
                sign(theta_pert), na.rm = TRUE), stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}
sparq_run <-
function (data, analysis_function, perturbation_function, comparator, 
    n_iterations = 50, result_id = "result", ...) 
{
    full_result <- analysis_function(data)
    perturbed_results <- vector(mode = "list", length = n_iterations)
    for (i in seq_len(n_iterations)) {
        perturbed_data <- perturbation_function(data, iteration = i)
        perturbed_results[[i]] <- analysis_function(perturbed_data)
    }
    sparq_assess(full_result = full_result, perturbed_results = perturbed_results, 
        comparator = comparator, result_id = result_id, ...)
}
sparq_run_with_stress_model <-
function (data, analysis_function, comparator, stress_model = "uniform_random", 
    retention = 0.75, n_iterations = 20, x_col = NULL, y_col = NULL, 
    custom_function = NULL, result_id = "result", reference_scale = 1, 
    seed = 1, ...) 
{
    if (stress_model == "none") {
        full_result <- analysis_function(data)
        return(data.frame(result_id = result_id, assessment_mode = "observed_baseline_only", 
            stress_model = "none", n_iterations = 0, support_quality = NA_real_, 
            support_class = "not_assessable", support_reason = "No perturbations were supplied", 
            stringsAsFactors = FALSE))
    }
    full_result <- analysis_function(data)
    perturbed_results <- lapply(seq_len(n_iterations), function(i) {
        perturbed_data <- sparq_stress_model(data = data, retention = retention, 
            model = stress_model, x_col = x_col, y_col = y_col, 
            custom_function = custom_function, iteration = i, 
            seed = seed)
        analysis_function(perturbed_data)
    })
    assessed <- sparq_assess(full_result = full_result, perturbed_results = perturbed_results, 
        comparator = comparator, result_id = result_id, reference_scale = reference_scale, 
        ...)
    output <- assessed$summary
    output$assessment_mode <- "automatic_rerun"
    output$stress_model <- stress_model
    output$retention <- retention
    output
}
sparq_safe_cor <-
function (x, y, method = "spearman") 
{
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) 
        return(NA_real_)
    if (length(unique(x[ok])) < 2) 
        return(NA_real_)
    if (length(unique(y[ok])) < 2) 
        return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}
sparq_safe_cor_internal <-
function (x, y, method = "spearman") 
{
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) 
        return(NA_real_)
    if (length(unique(x[ok])) < 2) 
        return(NA_real_)
    if (length(unique(y[ok])) < 2) 
        return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}
sparq_sample_fragility_ranked_feature <-
function (turnover_table) 
{
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("SPARQ requires data.table.", call. = FALSE)
    }
    x <- data.table::as.data.table(turnover_table)
    required <- c("cohort", "sample_id", "spearman_rank_correlation", 
        "topk_jaccard", "topk_replaced_fraction")
    missing <- setdiff(required, names(x))
    if (length(missing) > 0) {
        stop("SPARQ turnover table missing required columns: ", 
            paste(missing, collapse = ", "), call. = FALSE)
    }
    out <- x[, .(n_runs = .N, median_spearman_rank_correlation = stats::median(spearman_rank_correlation, 
        na.rm = TRUE), q25_spearman_rank_correlation = as.numeric(stats::quantile(spearman_rank_correlation, 
        0.25, na.rm = TRUE)), median_topk_jaccard = stats::median(topk_jaccard, 
        na.rm = TRUE), median_topk_replaced_fraction = stats::median(topk_replaced_fraction, 
        na.rm = TRUE), fragility_score = (1 - stats::median(spearman_rank_correlation, 
        na.rm = TRUE)) + stats::median(topk_replaced_fraction, 
        na.rm = TRUE)), by = .(cohort, sample_id)]
    as.data.frame(out)
}
sparq_severity_regression_internal <-
function (loss_fraction, delta) 
{
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
sparq_signed_quantile_internal <-
function (x) 
{
    x <- x[is.finite(x)]
    if (length(x) == 0) 
        return(NA_real_)
    med <- stats::median(x)
    if (med < 0) 
        return(as.numeric(stats::quantile(x, 0.25, names = FALSE)))
    if (med > 0) 
        return(as.numeric(stats::quantile(x, 0.75, names = FALSE)))
    0
}
sparq_stress_model <-
function (data, retention = 0.75, model = c("uniform_random", 
    "contiguous_hole", "none", "custom"), x_col = NULL, y_col = NULL, 
    custom_function = NULL, iteration = 1, seed = NULL) 
{
    model <- match.arg(model)
    if (!is.null(seed)) {
        set.seed(seed + iteration)
    }
    if (model == "none") {
        return(data)
    }
    if (model == "custom") {
        if (is.null(custom_function)) {
            stop("custom_function must be supplied when model = 'custom'.")
        }
        return(custom_function(data = data, retention = retention, 
            iteration = iteration))
    }
    if (retention <= 0 || retention > 1) {
        stop("retention must be greater than 0 and at most 1.")
    }
    n <- nrow(data)
    n_keep <- max(2, floor(n * retention))
    if (n_keep >= n) {
        return(data)
    }
    if (model == "uniform_random") {
        keep <- sample(seq_len(n), size = n_keep, replace = FALSE)
        return(data[keep, , drop = FALSE])
    }
    if (model == "contiguous_hole") {
        if (is.null(x_col) || is.null(y_col)) {
            stop("x_col and y_col are required for contiguous_hole.")
        }
        if (!all(c(x_col, y_col) %in% names(data))) {
            stop("Coordinate columns were not found in data.")
        }
        coordinates <- data[, c(x_col, y_col), drop = FALSE]
        coordinates[[x_col]] <- as.numeric(coordinates[[x_col]])
        coordinates[[y_col]] <- as.numeric(coordinates[[y_col]])
        valid <- complete.cases(coordinates)
        if (sum(valid) < n_keep) {
            stop("Not enough valid coordinates for contiguous-hole loss.")
        }
        center_index <- sample(which(valid), size = 1)
        dx <- coordinates[[x_col]] - coordinates[[x_col]][center_index]
        dy <- coordinates[[y_col]] - coordinates[[y_col]][center_index]
        distance <- sqrt(dx^2 + dy^2)
        distance[!valid] <- Inf
        keep <- order(distance, decreasing = TRUE)[seq_len(n_keep)]
        return(data[keep, , drop = FALSE])
    }
    stop("Unknown stress model.")
}
sparq_support <-
function (bias_table, reproducibility_table = NULL, min_iterations = 50, 
    min_direction_consistency = 0.69999999999999996, max_sign_flip_rate = 0.29999999999999999, 
    min_reproducibility = 0.5, max_relative_bias = 0.5, max_relative_uncertainty = 0.25) 
{
    needed <- c("result_id", "n_perturbations", "theta_full", 
        "bias", "uncertainty", "direction_consistency", "sign_flip_rate")
    missing <- setdiff(needed, colnames(bias_table))
    if (length(missing)) {
        stop("Missing required bias columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    out <- bias_table
    out$iteration_support <- sparq_bound01(out$n_perturbations/min_iterations)
    scale <- pmax(abs(out$theta_full), .Machine$double.eps)
    out$relative_bias <- abs(out$bias)/scale
    out$relative_uncertainty <- out$uncertainty/scale
    out$bias_support <- 1/(1 + out$relative_bias)
    out$uncertainty_support <- 1/(1 + out$relative_uncertainty)
    out$reproducibility_score <- 1
    if (!is.null(reproducibility_table)) {
        rep_cols <- intersect(c("spearman", "pearson", "sign_agreement"), 
            colnames(reproducibility_table))
        if (length(rep_cols)) {
            rep_score <- rowMeans(reproducibility_table[, rep_cols, 
                drop = FALSE], na.rm = TRUE)
            rep_df <- data.frame(result_id = reproducibility_table$result_id, 
                reproducibility_score = rep_score, stringsAsFactors = FALSE)
            out <- merge(out, rep_df, by = "result_id", all.x = TRUE, 
                suffixes = c("", ".from_rep"))
            out$reproducibility_score <- out$reproducibility_score.from_rep
            out$reproducibility_score.from_rep <- NULL
        }
    }
    out$reproducibility_support <- sparq_bound01(out$reproducibility_score/min_reproducibility)
    out$support_quality <- sparq_bound01(out$iteration_support * 
        out$bias_support * out$uncertainty_support * out$reproducibility_support)
    out$insufficient_support <- with(out, n_perturbations < min_iterations | 
        is.na(reproducibility_score) | reproducibility_score < 
        min_reproducibility | relative_bias > max_relative_bias | 
        relative_uncertainty > max_relative_uncertainty)
    out
}
sparq_support_from_comparisons <-
function (comparison_table, n_iterations, reference_scale = 1, 
    min_iterations = 50, min_reproducibility = 0.5, max_relative_bias = 0.5, 
    max_relative_uncertainty = 0.25) 
{
    if (!"instability" %in% names(comparison_table)) {
        stop("Comparator must return an 'instability' column.")
    }
    instability <- comparison_table$instability
    instability <- instability[is.finite(instability)]
    if (!length(instability)) {
        stop("No finite instability values were returned.")
    }
    bias <- if ("delta" %in% names(comparison_table)) {
        mean(comparison_table$delta, na.rm = TRUE)
    }
    else {
        median(instability, na.rm = TRUE)
    }
    uncertainty <- stats::sd(instability, na.rm = TRUE)
    reproducibility_score <- 1 - median(instability, na.rm = TRUE)
    reproducibility_score <- sparq_bound01(reproducibility_score)
    reference_scale <- max(abs(reference_scale), 9.9999999999999995e-07)
    relative_bias <- abs(bias)/reference_scale
    relative_uncertainty <- uncertainty/reference_scale
    iteration_support <- sparq_bound01(n_iterations/min_iterations)
    bias_support <- 1/(1 + relative_bias)
    uncertainty_support <- 1/(1 + relative_uncertainty)
    reproducibility_support <- sparq_bound01(reproducibility_score/min_reproducibility)
    support_quality <- sparq_bound01(iteration_support * bias_support * 
        uncertainty_support * reproducibility_support)
    insufficient_support <- n_iterations < min_iterations | reproducibility_score < 
        min_reproducibility | relative_bias > max_relative_bias | 
        relative_uncertainty > max_relative_uncertainty
    reason <- sparq_support_reason(n_iterations = n_iterations, 
        reproducibility_score = reproducibility_score, relative_bias = relative_bias, 
        relative_uncertainty = relative_uncertainty, min_iterations = min_iterations, 
        min_reproducibility = min_reproducibility, max_relative_bias = max_relative_bias, 
        max_relative_uncertainty = max_relative_uncertainty)
    ci <- sparq_bootstrap_ci(instability)
    data.frame(n_iterations = n_iterations, bias = bias, uncertainty = uncertainty, 
        relative_bias = relative_bias, relative_uncertainty = relative_uncertainty, 
        reproducibility_score = reproducibility_score, support_quality = support_quality, 
        instability_median = ci["estimate"], instability_ci_low = ci["lower"], 
        instability_ci_high = ci["upper"], insufficient_support = insufficient_support, 
        support_reason = reason, stringsAsFactors = FALSE)
}
sparq_support_reason <-
function (n_iterations, reproducibility_score, relative_bias, 
    relative_uncertainty, min_iterations = 50, min_reproducibility = 0.5, 
    max_relative_bias = 0.5, max_relative_uncertainty = 0.25) 
{
    reasons <- character()
    if (n_iterations < min_iterations) {
        reasons <- c(reasons, "few_iterations")
    }
    if (is.na(reproducibility_score)) {
        reasons <- c(reasons, "missing_reproducibility")
    }
    else if (reproducibility_score < min_reproducibility) {
        reasons <- c(reasons, "low_reproducibility")
    }
    if (is.finite(relative_bias) && relative_bias > max_relative_bias) {
        reasons <- c(reasons, "high_bias")
    }
    if (is.finite(relative_uncertainty) && relative_uncertainty > 
        max_relative_uncertainty) {
        reasons <- c(reasons, "high_uncertainty")
    }
    if (!length(reasons)) {
        "adequate_support"
    }
    else {
        paste(reasons, collapse = ";")
    }
}
sparq_topk_jaccard <-
function (full_ids, pert_ids, k = 100) 
{
    full_ids <- as.character(full_ids)
    pert_ids <- as.character(pert_ids)
    full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
    pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])
    if (length(full_top) == 0 && length(pert_top) == 0) 
        return(NA_real_)
    length(intersect(full_top, pert_top))/length(union(full_top, 
        pert_top))
}
sparq_topk_replacement <-
function (full_ids, pert_ids, k = 100) 
{
    full_ids <- as.character(full_ids)
    pert_ids <- as.character(pert_ids)
    full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
    pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])
    if (length(full_top) == 0) 
        return(NA_real_)
    1 - (length(intersect(full_top, pert_top))/length(full_top))
}
sparq_trimmed_median_internal <-
function (x, trim = 0.20000000000000001) 
{
    x <- sort(x[is.finite(x)])
    n <- length(x)
    if (n == 0) 
        return(NA_real_)
    lo <- floor(n * trim) + 1L
    hi <- ceiling(n * (1 - trim))
    if (lo > hi) 
        return(stats::median(x))
    stats::median(x[lo:hi])
}
sparq_validate_columns <-
function (data, required_cols) 
{
    missing <- setdiff(required_cols, names(data))
    if (length(missing) > 0) {
        stop("SPARQ input missing required columns: ", paste(missing, 
            collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}
sparq_validate_ranked_feature <-
function (data, result_id_col = "result_id", theta_full_col = "theta_full", 
    theta_pert_col = "theta_pert", rank_full_col = "rank_full", 
    rank_pert_col = "rank_pert", scenario_col = "scenario", iteration_col = "iteration") 
{
    required <- c(result_id_col, theta_full_col, theta_pert_col, 
        rank_full_col, rank_pert_col, scenario_col, iteration_col)
    sparq_validate_columns(data, required)
    if (!is.numeric(data[[theta_full_col]])) {
        stop(theta_full_col, " must be numeric.", call. = FALSE)
    }
    if (!is.numeric(data[[theta_pert_col]])) {
        stop(theta_pert_col, " must be numeric.", call. = FALSE)
    }
    if (!is.numeric(data[[rank_full_col]]) && !is.integer(data[[rank_full_col]])) {
        stop(rank_full_col, " must be numeric or integer.", call. = FALSE)
    }
    if (!is.numeric(data[[rank_pert_col]]) && !is.integer(data[[rank_pert_col]])) {
        stop(rank_pert_col, " must be numeric or integer.", call. = FALSE)
    }
    if (nrow(data) == 0) {
        stop("SPARQ input has zero rows.", call. = FALSE)
    }
    invisible(TRUE)
}
sparq_validate_results <-
function (data, required_cols = c("result_id", "theta_full", 
    "theta_pert", "scenario", "iteration")) 
{
    missing <- setdiff(required_cols, colnames(data))
    if (length(missing)) {
        stop("Missing required columns: ", paste(missing, collapse = ", "), 
            call. = FALSE)
    }
    if (!is.numeric(data$theta_full)) {
        stop("theta_full must be numeric.", call. = FALSE)
    }
    if (!is.numeric(data$theta_pert)) {
        stop("theta_pert must be numeric.", call. = FALSE)
    }
    invisible(TRUE)
}
sparq_assess_precomputed_files <-
function (full_file, perturbed_file, comparator = sparq_compare_scalar, 
    reader = read.delim, result_id_col = "result_id", full_value_col = "value", 
    perturbed_value_col = "value", ...) 
{
    if (!file.exists(full_file)) {
        stop("Full-result file does not exist: ", full_file, 
            call. = FALSE)
    }
    if (!file.exists(perturbed_file)) {
        stop("Perturbed-result file does not exist: ", perturbed_file, 
            call. = FALSE)
    }
    full_table <- reader(full_file, sep = "\t", stringsAsFactors = FALSE, 
        check.names = FALSE, ...)
    perturbed_table <- reader(perturbed_file, sep = "\t", stringsAsFactors = FALSE, 
        check.names = FALSE, ...)
    required_full <- c(result_id_col, full_value_col)
    required_perturbed <- c(result_id_col, perturbed_value_col)
    missing_full <- setdiff(required_full, names(full_table))
    missing_perturbed <- setdiff(required_perturbed, names(perturbed_table))
    if (length(missing_full)) {
        stop("Full-result file is missing: ", paste(missing_full, 
            collapse = ", "), call. = FALSE)
    }
    if (length(missing_perturbed)) {
        stop("Perturbed-result file is missing: ", paste(missing_perturbed, 
            collapse = ", "), call. = FALSE)
    }
    full_table[[result_id_col]] <- as.character(full_table[[result_id_col]])
    perturbed_table[[result_id_col]] <- as.character(perturbed_table[[result_id_col]])
    full_table[[full_value_col]] <- as.numeric(full_table[[full_value_col]])
    perturbed_table[[perturbed_value_col]] <- as.numeric(perturbed_table[[perturbed_value_col]])
    full_results <- lapply(split(full_table[[full_value_col]], 
        full_table[[result_id_col]]), function(x) median(x[is.finite(x)], 
        na.rm = TRUE))
    perturbed_results <- lapply(split(perturbed_table[[perturbed_value_col]], 
        perturbed_table[[result_id_col]]), function(x) x[is.finite(x)])
    sparq_assess_precomputed(full_results = full_results, perturbed_results = perturbed_results, 
        comparator = comparator)
}
