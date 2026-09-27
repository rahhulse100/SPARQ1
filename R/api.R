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
