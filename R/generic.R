sparq_bound01 <-
function (x) 
{
    pmax(0, pmin(1, x))
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
sparq_compare_scalar <-
function (full_result, perturbed_result) 
{
    delta <- as.numeric(perturbed_result) - as.numeric(full_result)
    data.frame(delta = delta, instability = abs(delta), stringsAsFactors = FALSE)
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
