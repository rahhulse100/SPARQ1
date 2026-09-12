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
