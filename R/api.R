sparq_run <-
function (data, analysis_function, comparator, stress_model = "uniform_random", 
    retention = 0.75, n_iterations = 50, result_id = "result", 
    x_col = NULL, y_col = NULL, custom_function = NULL, reference_scale = 1, 
    seed = 1, cache_dir = NULL, resume = TRUE, ...) 
{
    if (!is.data.frame(data)) {
        stop("data must be a data.frame.", call. = FALSE)
    }
    if (!is.function(analysis_function)) {
        stop("analysis_function must be a function.", call. = FALSE)
    }
    if (!is.function(comparator)) {
        stop("comparator must be a function.", call. = FALSE)
    }
    if (length(n_iterations) != 1 || !is.finite(n_iterations) || 
        n_iterations < 1) {
        stop("n_iterations must be a positive integer.", call. = FALSE)
    }
    n_iterations <- as.integer(n_iterations)
    if (!is.null(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    set.seed(seed)
    full_result <- analysis_function(data)
    comparison_list <- vector(mode = "list", length = n_iterations)
    for (i in seq_len(n_iterations)) {
        cache_file <- NULL
        if (!is.null(cache_dir)) {
            cache_file <- file.path(cache_dir, paste0(result_id, 
                "_iteration_", i, ".rds"))
        }
        if (resume && !is.null(cache_file) && file.exists(cache_file)) {
            comparison_list[[i]] <- readRDS(cache_file)
            next
        }
        perturbed_data <- sparq_stress_model(data = data, retention = retention, 
            model = stress_model, x_col = x_col, y_col = y_col, 
            custom_function = custom_function, iteration = i, 
            seed = seed)
        perturbed_result <- analysis_function(perturbed_data)
        comparison <- comparator(full_result, perturbed_result)
        if (!is.data.frame(comparison)) {
            comparison <- as.data.frame(comparison)
        }
        if (!"instability" %in% names(comparison)) {
            stop("Comparator must return an 'instability' column.", 
                call. = FALSE)
        }
        comparison$iteration <- i
        comparison_list[[i]] <- comparison
        if (!is.null(cache_file)) {
            saveRDS(comparison, cache_file)
        }
    }
    comparison_table <- do.call(rbind, comparison_list)
    rownames(comparison_table) <- NULL
    summary <- sparq_support_from_comparisons(comparison_table = comparison_table, 
        n_iterations = n_iterations, reference_scale = reference_scale, 
        ...)
    summary$result_id <- result_id
    summary$assessment_mode <- "automatic_rerun"
    summary$stress_model <- stress_model
    summary$retention <- retention
    list(result_id = result_id, full_result = full_result, comparisons = comparison_table, 
        summary = summary)
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
sparq_read_result_file <-
function (path, reader = NULL) 
{
    if (length(path) != 1 || is.na(path) || !nzchar(path) || 
        !file.exists(path)) {
        stop("Result file does not exist: ", path, call. = FALSE)
    }
    ext <- tolower(tools::file_ext(path))
    if (!is.null(reader)) {
        return(reader(path, stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (ext %in% c("tsv", "txt")) {
        return(read.delim(path, sep = "\t", stringsAsFactors = FALSE, 
            check.names = FALSE))
    }
    if (ext == "csv") {
        return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
    }
    if (ext == "rds") {
        return(readRDS(path))
    }
    stop("Unsupported result-file type: ", ext, call. = FALSE)
}
sparq_assess_precomputed_files <-
function (full_file, perturbed_file, comparator = sparq_compare_scalar, 
    full_result_id_col = "result_id", perturbed_result_id_col = "result_id", 
    full_value_col = NULL, perturbed_value_col = NULL, reader = NULL, 
    ...) 
{
    full_results <- sparq_read_result_file(full_file, reader = reader)
    perturbed_results <- sparq_read_result_file(perturbed_file, 
        reader = reader)
    if (is.list(full_results) && is.list(perturbed_results) && 
        !is.data.frame(full_results) && !is.data.frame(perturbed_results)) {
        return(sparq_assess_precomputed(full_results = full_results, 
            perturbed_results = perturbed_results, comparator = comparator, 
            ...))
    }
    if (!is.data.frame(full_results) || !is.data.frame(perturbed_results)) {
        stop("Result files must contain data frames or named RDS lists.", 
            call. = FALSE)
    }
    if (!all(c(full_result_id_col) %in% names(full_results))) {
        stop("Full-result file is missing column: ", full_result_id_col, 
            call. = FALSE)
    }
    if (!all(c(perturbed_result_id_col) %in% names(perturbed_results))) {
        stop("Perturbed-result file is missing column: ", perturbed_result_id_col, 
            call. = FALSE)
    }
    if (is.null(full_value_col)) {
        candidates <- c("value", "estimate", "statistic", "score", 
            "theta_full")
        hit <- candidates[candidates %in% names(full_results)]
        if (!length(hit)) {
            stop("Could not identify the full-result value column. ", 
                "Supply full_value_col explicitly.", call. = FALSE)
        }
        full_value_col <- hit[1]
    }
    if (is.null(perturbed_value_col)) {
        candidates <- c("value", "estimate", "statistic", "score", 
            "theta_pert")
        hit <- candidates[candidates %in% names(perturbed_results)]
        if (!length(hit)) {
            stop("Could not identify the perturbed-result value column. ", 
                "Supply perturbed_value_col explicitly.", call. = FALSE)
        }
        perturbed_value_col <- hit[1]
    }
    full_results[[full_result_id_col]] <- as.character(full_results[[full_result_id_col]])
    perturbed_results[[perturbed_result_id_col]] <- as.character(perturbed_results[[perturbed_result_id_col]])
    full_results[[full_value_col]] <- suppressWarnings(as.numeric(full_results[[full_value_col]]))
    perturbed_results[[perturbed_value_col]] <- suppressWarnings(as.numeric(perturbed_results[[perturbed_value_col]]))
    full_list <- lapply(split(full_results[[full_value_col]], 
        full_results[[full_result_id_col]]), function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) 
            NA_real_
        else median(x)
    })
    perturbed_list <- lapply(split(perturbed_results[[perturbed_value_col]], 
        perturbed_results[[perturbed_result_id_col]]), function(x) {
        x[is.finite(x)]
    })
    sparq_assess_precomputed(full_results = full_list, perturbed_results = perturbed_list, 
        comparator = comparator, ...)
}
sparq_assess_cohort_parallel <-
function (sample_data, analysis_function, comparator, stress_model = "uniform_random", 
    retention = 0.75, n_iterations = 50, workers = 1, checkpoint_dir = NULL, 
    resume = TRUE, error_log = NULL, seed = 1, ...) 
{
    if (is.null(names(sample_data))) {
        names(sample_data) <- paste0("sample_", seq_along(sample_data))
    }
    sample_ids <- names(sample_data)
    if (!is.null(checkpoint_dir)) {
        dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (is.null(error_log) && !is.null(checkpoint_dir)) {
        error_log <- file.path(checkpoint_dir, "SPARQ_error_log.tsv")
    }
    run_one <- function(sample_id) {
        checkpoint_file <- NULL
        if (!is.null(checkpoint_dir)) {
            checkpoint_file <- file.path(checkpoint_dir, paste0(sample_id, 
                ".rds"))
        }
        if (resume && !is.null(checkpoint_file) && file.exists(checkpoint_file)) {
            return(list(result = readRDS(checkpoint_file), error = NULL))
        }
        message("Processing sample: ", sample_id)
        ans <- tryCatch({
            result <- sparq_run_with_stress_model(data = sample_data[[sample_id]], 
                analysis_function = analysis_function, comparator = comparator, 
                stress_model = stress_model, retention = retention, 
                n_iterations = n_iterations, result_id = sample_id, 
                seed = seed, ...)
            if (!is.null(checkpoint_file)) {
                saveRDS(result, checkpoint_file)
            }
            list(result = result, error = NULL)
        }, error = function(e) {
            list(result = NULL, error = data.frame(result_id = sample_id, 
                error_message = conditionMessage(e), stringsAsFactors = FALSE))
        })
        ans
    }
    if (workers > 1) {
        workers <- min(as.integer(workers), length(sample_ids))
        cl <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterExport(cl, varlist = c("sample_data", 
            "analysis_function", "comparator", "stress_model", 
            "retention", "n_iterations", "checkpoint_dir", "resume", 
            "seed", "run_one"), envir = environment())
        outputs <- parallel::parLapply(cl, sample_ids, run_one)
    }
    else {
        outputs <- lapply(sample_ids, run_one)
    }
    successful <- outputs[vapply(outputs, function(x) !is.null(x$result), 
        logical(1))]
    failed <- outputs[vapply(outputs, function(x) !is.null(x$error), 
        logical(1))]
    result_table <- NULL
    if (length(successful)) {
        result_table <- do.call(rbind, lapply(successful, function(x) x$result))
        rownames(result_table) <- NULL
    }
    error_table <- NULL
    if (length(failed)) {
        error_table <- do.call(rbind, lapply(failed, function(x) x$error))
        if (!is.null(error_log)) {
            write.table(error_table, error_log, sep = "\t", row.names = FALSE, 
                quote = FALSE)
        }
    }
    list(results = result_table, errors = error_table, expected_samples = sample_ids, 
        completed_samples = if (is.null(result_table)) character() else result_table$result_id, 
        failed_samples = if (is.null(error_table)) character() else error_table$result_id)
}
sparq_calibrate_threshold <-
function (support_quality, truth, held_out_support_quality = NULL, 
    held_out_truth = NULL, folds = 5, B = 2000, seed = 1, min_class_n = 5) 
{
    set.seed(seed)
    ok <- is.finite(support_quality) & !is.na(truth)
    support_quality <- sparq_bound01(as.numeric(support_quality[ok]))
    truth <- as.integer(truth[ok])
    if (length(unique(truth)) != 2 || sum(truth == 1) < min_class_n || 
        sum(truth == 0) < min_class_n) {
        stop("Both truth classes require at least ", min_class_n, 
            " samples.", call. = FALSE)
    }
    metric_at_threshold <- function(threshold, support, truth) {
        predicted <- support >= threshold
        stable <- truth == 1
        unstable <- truth == 0
        sensitivity <- mean(predicted[stable])
        specificity <- mean(!predicted[unstable])
        data.frame(threshold = threshold, sensitivity = sensitivity, 
            specificity = specificity, balanced_accuracy = (sensitivity + 
                specificity)/2)
    }
    threshold_grid <- sort(unique(c(seq(0, 1, length.out = 1001), 
        support_quality)))
    train_table <- do.call(rbind, lapply(threshold_grid, metric_at_threshold, 
        support = support_quality, truth = truth))
    best <- train_table[which.max(train_table$balanced_accuracy), 
        , drop = FALSE]
    calibrated_threshold <- best$threshold
    bootstrap_thresholds <- replicate(B, {
        index <- sample(seq_along(support_quality), replace = TRUE)
        boot_support <- support_quality[index]
        boot_truth <- truth[index]
        boot_table <- do.call(rbind, lapply(threshold_grid, metric_at_threshold, 
            support = boot_support, truth = boot_truth))
        boot_table$threshold[which.max(boot_table$balanced_accuracy)]
    })
    threshold_interval <- as.numeric(quantile(bootstrap_thresholds, 
        probs = c(0.025000000000000001, 0.97499999999999998), 
        na.rm = TRUE))
    threshold_interval <- data.frame(threshold = calibrated_threshold, 
        lower = threshold_interval[1], upper = threshold_interval[2])
    fold_id <- sample(rep(seq_len(folds), length.out = length(support_quality)))
    cv_results <- do.call(rbind, lapply(seq_len(folds), function(fold) {
        train_index <- fold_id != fold
        test_index <- fold_id == fold
        train_support <- support_quality[train_index]
        train_truth <- truth[train_index]
        test_support <- support_quality[test_index]
        test_truth <- truth[test_index]
        fold_table <- do.call(rbind, lapply(threshold_grid, metric_at_threshold, 
            support = train_support, truth = train_truth))
        fold_threshold <- fold_table$threshold[which.max(fold_table$balanced_accuracy)]
        test_metrics <- metric_at_threshold(fold_threshold, test_support, 
            test_truth)
        test_metrics$fold <- fold
        test_metrics
    }))
    held_out_metrics <- NULL
    if (!is.null(held_out_support_quality) && !is.null(held_out_truth)) {
        held_out_ok <- is.finite(held_out_support_quality) & 
            !is.na(held_out_truth)
        held_out_metrics <- metric_at_threshold(calibrated_threshold, 
            sparq_bound01(as.numeric(held_out_support_quality[held_out_ok])), 
            as.integer(held_out_truth[held_out_ok]))
    }
    list(threshold = calibrated_threshold, threshold_interval = threshold_interval, 
        training_metrics = best, cross_validation = cv_results, 
        held_out_metrics = held_out_metrics, bootstrap_thresholds = bootstrap_thresholds)
}
sparq_reliability_curve <-
function (support_quality, true_error, bins = 8, min_observations_per_bin = 10, 
    B = 2000, seed = 1) 
{
    set.seed(seed)
    ok <- is.finite(support_quality) & is.finite(true_error)
    support_quality <- sparq_bound01(as.numeric(support_quality[ok]))
    true_error <- as.numeric(true_error[ok])
    if (length(support_quality) < min_observations_per_bin) {
        stop("Not enough observations for reliability analysis.", 
            call. = FALSE)
    }
    bins <- min(bins, floor(length(support_quality)/min_observations_per_bin))
    rank_order <- order(support_quality)
    ordered_support <- support_quality[rank_order]
    ordered_error <- true_error[rank_order]
    bin_id <- cut(seq_along(ordered_support), breaks = bins, 
        labels = FALSE)
    observed <- lapply(seq_len(bins), function(i) {
        index <- bin_id == i
        support <- ordered_support[index]
        error <- ordered_error[index]
        if (length(support) < min_observations_per_bin) {
            return(NULL)
        }
        bootstrap_error <- replicate(B, median(sample(error, 
            replace = TRUE)))
        error_ci <- as.numeric(quantile(bootstrap_error, probs = c(0.025000000000000001, 
            0.97499999999999998), na.rm = TRUE))
        data.frame(bin = i, n = length(support), support_mean = mean(support), 
            support_median = median(support), true_error_mean = mean(error), 
            true_error_median = median(error), true_error_ci_low = error_ci[1], 
            true_error_ci_high = error_ci[2], calibration_error = mean(abs(support - 
                sparq_bound01(1 - error))), stringsAsFactors = FALSE)
    })
    observed <- observed[!vapply(observed, is.null, logical(1))]
    if (!length(observed)) {
        stop("No reliability bins met the minimum observation requirement.", 
            call. = FALSE)
    }
    curve <- do.call(rbind, observed)
    rownames(curve) <- NULL
    monotonicity <- NA_real_
    monotonicity_p <- NA_real_
    if (nrow(curve) >= 3) {
        test <- suppressWarnings(cor.test(curve$support_mean, 
            curve$true_error_median, method = "spearman", exact = FALSE))
        monotonicity <- as.numeric(test$estimate)
        monotonicity_p <- test$p.value
    }
    overall_calibration_error <- weighted.mean(curve$calibration_error, 
        curve$n)
    list(curve = curve, monotonicity_spearman = monotonicity, 
        monotonicity_p_value = monotonicity_p, overall_calibration_error = overall_calibration_error, 
        n_observations = length(support_quality), n_bins = nrow(curve))
}
