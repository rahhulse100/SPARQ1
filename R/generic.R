sparq_effect_scale <-
function (reference, floor = 9.9999999999999995e-07) 
{
    reference <- reference[is.finite(reference)]
    if (!length(reference)) {
        return(floor)
    }
    max(stats::mad(reference, na.rm = TRUE), floor)
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
