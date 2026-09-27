sparq_ablation_score <-
function (rows, mode = c("full", "no_uncertainty", "no_reproducibility")) 
{
    mode <- match.arg(mode)
    if (mode == "full") {
        return(rows$support_quality)
    }
    if (mode == "no_uncertainty") {
        return(rows$iteration_support * rows$bias_support * rows$reproducibility_support)
    }
    rows$iteration_support * rows$bias_support * rows$uncertainty_support
}
sparq_ablation_continuous_metrics <-
function (truth, score) 
{
    x <- sparq_binary_metrics(truth, score, threshold = 0.5)
    data.frame(sensitivity = NA_real_, specificity = NA_real_, 
        balanced_accuracy = NA_real_, auroc = x$auroc, auprc = x$auprc, 
        unstable_prevalence = x$unstable_prevalence, tp = NA_integer_, 
        fn = NA_integer_, tn = NA_integer_, fp = NA_integer_)
}
sparq_ablation_metric_summary <-
function (metrics, confidence = 0.94999999999999996) 
{
    measures <- c("sensitivity", "specificity", "balanced_accuracy", 
        "auroc", "auprc")
    out <- lapply(split(metrics, metrics$ablation_mode), function(d) {
        sparq_bind_rows(lapply(measures, function(measure) {
            x <- d[[measure]]
            x <- x[is.finite(x)]
            if (!length(x)) {
                return(NULL)
            }
            q <- stats::quantile(x, c((1 - confidence)/2, 1 - 
                (1 - confidence)/2))
            data.frame(ablation_mode = d$ablation_mode[1], score_definition = d$score_definition[1], 
                threshold_mode = d$threshold_mode[1], metric = measure, 
                n_replicates = length(x), mean = mean(x), sd = if (length(x) > 
                  1L) {
                  stats::sd(x)
                }
                else {
                  NA_real_
                }, median = stats::median(x), replicate_interval_low = unname(q[1]), 
                replicate_interval_high = unname(q[2]), interval_type = "central empirical replicate interval; not a confidence interval for the mean")
        }))
    })
    sparq_bind_rows(out)
}
sparq_known_truth_ablation <-
function (benchmark_result, fixed_thresholds = c(0.25, 0.5, 0.75), 
    threshold_bootstrap_B = NULL, confidence = NULL, seed = NULL) 
{
    if (is.list(benchmark_result) && !is.null(benchmark_result$all_results)) {
        rows <- benchmark_result$all_results
        design <- benchmark_result$design
        if (is.null(threshold_bootstrap_B)) {
            threshold_bootstrap_B <- design$threshold_bootstrap_B
        }
        if (is.null(confidence)) {
            confidence <- design$confidence
        }
        if (is.null(seed)) {
            seed <- design$seed
        }
    }
    else if (is.data.frame(benchmark_result)) {
        rows <- benchmark_result
        design <- NULL
        if (is.null(threshold_bootstrap_B)) {
            threshold_bootstrap_B <- 2000L
        }
        if (is.null(confidence)) {
            confidence <- 0.94999999999999996
        }
        if (is.null(seed)) {
            seed <- 1L
        }
    }
    else {
        stop("benchmark_result must be a Figure 3C-E result object or its all_results table.", 
            call. = FALSE)
    }
    required <- c("benchmark_id", "replicate", "seed", "split", 
        "truth", "support_quality", "iteration_support", "bias_support", 
        "uncertainty_support", "reproducibility_support")
    if (!is.data.frame(rows) || !all(required %in% names(rows))) {
        stop("Benchmark results lack required Figure 3C-E columns.", 
            call. = FALSE)
    }
    if (anyNA(rows$benchmark_id) || anyNA(rows$replicate) || 
        anyDuplicated(paste(rows$benchmark_id, rows$replicate, 
            sep = "\r"))) {
        stop("Benchmark results must contain one row per sample and replicate.", 
            call. = FALSE)
    }
    if (!all(rows$split %in% c("training", "held_out")) || !all(rows$truth %in% 
        c(0, 1))) {
        stop("Benchmark results have invalid split or truth values.", 
            call. = FALSE)
    }
    components <- c("support_quality", "iteration_support", "bias_support", 
        "uncertainty_support", "reproducibility_support")
    component_matrix <- as.matrix(rows[, components, drop = FALSE])
    if (any(!is.finite(component_matrix)) || any(component_matrix < 
        0 | component_matrix > 1)) {
        stop("Support components must be finite values in [0, 1].", 
            call. = FALSE)
    }
    fixed_thresholds <- as.numeric(fixed_thresholds)
    if (!length(fixed_thresholds) || any(!is.finite(fixed_thresholds)) || 
        any(fixed_thresholds < 0 | fixed_thresholds > 1) || anyDuplicated(fixed_thresholds)) {
        stop("fixed_thresholds must be distinct finite values in [0, 1].", 
            call. = FALSE)
    }
    sparq_check_number(threshold_bootstrap_B, "threshold_bootstrap_B", 
        1, integer = TRUE)
    sparq_check_number(confidence, "confidence", .Machine$double.eps, 
        1 - .Machine$double.eps)
    sparq_check_number(seed, "seed", 0, .Machine$integer.max, 
        integer = TRUE)
    replicates <- sort(unique(rows$replicate))
    specs <- rbind(data.frame(ablation_mode = c("full_sparq", 
        "no_uncertainty_component", "no_reproducibility_component", 
        "support_score_only_no_calibration"), score_mode = c("full", 
        "no_uncertainty", "no_reproducibility", "full"), score_definition = c("Recorded four-component support quality", 
        "Iteration support x bias support x reproducibility support", 
        "Iteration support x bias support x uncertainty support", 
        "Recorded four-component support quality"), threshold_mode = c("training_calibrated", 
        "training_calibrated", "training_calibrated", "none"), 
        fixed_threshold = NA_real_, stringsAsFactors = FALSE), 
        data.frame(ablation_mode = paste0("fixed_threshold_", 
            formatC(fixed_thresholds, format = "f", digits = 2)), 
            score_mode = "full", score_definition = "Recorded four-component support quality", 
            threshold_mode = "fixed", fixed_threshold = fixed_thresholds, 
            stringsAsFactors = FALSE))
    held_rows <- list()
    metrics <- list()
    thresholds <- list()
    curves <- list()
    item <- 0L
    for (replicate_id in replicates) {
        block <- rows[rows$replicate == replicate_id, , drop = FALSE]
        training <- block[block$split == "training", , drop = FALSE]
        held <- block[block$split == "held_out", , drop = FALSE]
        if (length(unique(training$truth)) != 2L || length(unique(held$truth)) != 
            2L) {
            stop("Every replicate requires both classes in training and held-out rows.", 
                call. = FALSE)
        }
        for (j in seq_len(nrow(specs))) {
            spec <- specs[j, , drop = FALSE]
            train_score <- sparq_ablation_score(training, spec$score_mode)
            held_score <- sparq_ablation_score(held, spec$score_mode)
            if (spec$threshold_mode == "training_calibrated") {
                calibration <- sparq_calibrate_threshold(train_score, 
                  training$truth, held_score, held$truth, B = threshold_bootstrap_B, 
                  confidence = confidence, seed = sparq_seed(seed + 
                    replicate_id, j, 30))
                threshold <- calibration$threshold
                threshold_low <- calibration$threshold_interval$lower
                threshold_high <- calibration$threshold_interval$upper
                metric <- calibration$held_out_metrics
                thresholds[[length(thresholds) + 1L]] <- data.frame(ablation_mode = spec$ablation_mode, 
                  replicate = replicate_id, draw = seq_along(calibration$bootstrap_thresholds), 
                  threshold = calibration$bootstrap_thresholds)
            }
            else if (spec$threshold_mode == "fixed") {
                threshold <- spec$fixed_threshold
                threshold_low <- NA_real_
                threshold_high <- NA_real_
                metric <- sparq_binary_metrics(held$truth, held_score, 
                  threshold)
            }
            else {
                threshold <- NA_real_
                threshold_low <- NA_real_
                threshold_high <- NA_real_
                metric <- sparq_ablation_continuous_metrics(held$truth, 
                  held_score)
            }
            item <- item + 1L
            prediction <- if (is.finite(threshold)) {
                as.integer(held_score >= threshold)
            }
            else {
                rep(NA_integer_, nrow(held))
            }
            keep <- intersect(c("benchmark_id", "difficulty", 
                "true_instability", "truth", "replicate", "seed", 
                "split"), names(held))
            z <- held[, keep, drop = FALSE]
            z$ablation_mode <- spec$ablation_mode
            z$score_definition <- spec$score_definition
            z$threshold_mode <- spec$threshold_mode
            z$ablation_support_quality <- held_score
            z$threshold <- threshold
            z$threshold_ci_low <- threshold_low
            z$threshold_ci_high <- threshold_high
            z$predicted_stable <- prediction
            held_rows[[item]] <- z
            metric$ablation_mode <- spec$ablation_mode
            metric$score_definition <- spec$score_definition
            metric$threshold_mode <- spec$threshold_mode
            metric$threshold <- threshold
            metric$threshold_ci_low <- threshold_low
            metric$threshold_ci_high <- threshold_high
            metric$replicate <- replicate_id
            metric$seed <- held$seed[1]
            metric$n_training <- nrow(training)
            metric$n_held_out <- nrow(held)
            metrics[[item]] <- metric
            curve <- sparq_curves(held$truth, held_score)
            curve$roc$ablation_mode <- spec$ablation_mode
            curve$roc$replicate <- replicate_id
            curve$precision_recall$ablation_mode <- spec$ablation_mode
            curve$precision_recall$replicate <- replicate_id
            curves[[item]] <- curve
        }
    }
    metrics <- sparq_bind_rows(metrics)
    list(held_out_results = sparq_bind_rows(held_rows), replicate_metrics = metrics, 
        metric_summary = sparq_ablation_metric_summary(metrics, 
            confidence), threshold_bootstrap = sparq_bind_rows(thresholds), 
        roc = sparq_bind_rows(lapply(curves, `[[`, "roc")), precision_recall = sparq_bind_rows(lapply(curves, 
            `[[`, "precision_recall")), ablation_specification = specs, 
        fixed_thresholds = fixed_thresholds, design = design, 
        source = "Figure 3C-E all_results; no perturbations were rerun")
}
sparq_export_ablation <-
function (result, output_dir) 
{
    required <- c("held_out_results", "replicate_metrics", "metric_summary", 
        "threshold_bootstrap", "roc", "precision_recall", "ablation_specification")
    if (!is.list(result) || !all(required %in% names(result))) {
        stop("Invalid ablation result.", call. = FALSE)
    }
    if (dir.exists(output_dir) && length(list.files(output_dir, 
        all.files = TRUE, no.. = TRUE))) {
        stop("Use a new, empty output directory; existing results will not be overwritten.", 
            call. = FALSE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (name in required) {
        utils::write.table(result[[name]], file.path(output_dir, 
            paste0(name, ".tsv")), sep = "\t", row.names = FALSE, 
            quote = FALSE)
    }
    if (!is.null(result$design)) {
        sparq_write_design(result$design, file.path(output_dir, 
            "design.yml"))
    }
    saveRDS(result, file.path(output_dir, "figure3_ablation.rds"))
    writeLines(result$source, file.path(output_dir, "source.txt"))
    invisible(normalizePath(output_dir))
}
