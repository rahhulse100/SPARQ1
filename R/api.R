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
