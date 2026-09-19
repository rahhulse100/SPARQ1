sparq_known_truth_benchmark <- function (n_samples = 400, n_iterations = 100, 
    train_fraction = 0.7, n_replicates = 20, seed = 20260912, 
    difficulty_parameters = data.frame(difficulty = c("Stable", 
        "Mild instability", "Moderate instability", "Severe instability"), 
        mean_instability = c(0.02, 0.07, 0.13, 0.24), sd_instability = c(0.03, 
            0.06, 0.08, 0.12)), threshold_bootstrap_B = 500, 
    threshold_confidence = 0.95) 
{
    if (n_samples < 20) {
        stop("n_samples must be at least 20.", call. = FALSE)
    }
    if (n_iterations < 10) {
        stop("n_iterations must be at least 10.", call. = FALSE)
    }
    if (train_fraction <= 0 || train_fraction >= 1) {
        stop("train_fraction must be between 0 and 1.", call. = FALSE)
    }
    if (!all(c("difficulty", "mean_instability", "sd_instability") %in% 
        names(difficulty_parameters))) {
        stop("difficulty_parameters must contain difficulty, mean_instability, and sd_instability.", 
            call. = FALSE)
    }
    auc_rank <- function(truth, score) {
        positive <- score[truth == 1]
        negative <- score[truth == 0]
        if (!length(positive) || !length(negative)) {
            return(NA_real_)
        }
        mean(vapply(positive, function(x) mean(x > negative) + 
            0.5 * mean(x == negative), numeric(1)))
    }
    average_precision <- function(truth, score) {
        ord <- order(score, decreasing = TRUE)
        y <- truth[ord]
        if (!any(y == 1)) {
            return(NA_real_)
        }
        precision <- cumsum(y)/seq_along(y)
        recall_increment <- y/sum(y)
        sum(precision * recall_increment)
    }
    metric_at_threshold <- function(truth, score, threshold) {
        predicted <- as.integer(score >= threshold)
        stable <- truth == 1
        unstable <- truth == 0
        tp <- sum(stable & predicted == 1)
        fn <- sum(stable & predicted == 0)
        tn <- sum(unstable & predicted == 0)
        fp <- sum(unstable & predicted == 1)
        sensitivity <- tp/max(1, tp + fn)
        specificity <- tn/max(1, tn + fp)
        data.frame(sensitivity = sensitivity, specificity = specificity, 
            balanced_accuracy = mean(c(sensitivity, specificity)), 
            auroc = auc_rank(truth, score), auprc = average_precision(truth, 
                score), stringsAsFactors = FALSE)
    }
    choose_threshold <- function(truth, score) {
        grid <- sort(unique(c(seq(0, 1, length.out = 501), score)))
        candidates <- lapply(grid, function(threshold) {
            metric <- metric_at_threshold(truth, score, threshold)
            metric$threshold <- threshold
            metric
        })
        candidates <- do.call(rbind, candidates)
        candidates[which.max(candidates$balanced_accuracy), , 
            drop = FALSE]
    }
    bootstrap_threshold_ci <- function(truth, score, B, confidence, 
        seed) {
        set.seed(seed)
        thresholds <- replicate(B, {
            index <- sample(seq_along(truth), size = length(truth), 
                replace = TRUE)
            choose_threshold(truth[index], score[index])$threshold
        })
        alpha <- 1 - confidence
        c(threshold_ci_low = unname(quantile(thresholds, alpha/2)), 
            threshold_ci_high = unname(quantile(thresholds, 1 - 
                alpha/2)))
    }
    all_results <- vector("list", n_replicates)
    all_metrics <- vector("list", n_replicates)
    for (r in seq_len(n_replicates)) {
        set.seed(seed + r)
        difficulty <- sample(difficulty_parameters$difficulty, 
            size = n_samples, replace = TRUE, prob = c(0.5, 0.2, 
                0.2, 0.1))
        parameters <- difficulty_parameters[match(difficulty, 
            difficulty_parameters$difficulty), , drop = FALSE]
        true_instability <- pmax(0, rnorm(n_samples, mean = parameters$mean_instability, 
            sd = parameters$sd_instability))
        truth <- as.integer(difficulty == "Stable")
        sample_rows <- vector("list", n_samples)
        for (i in seq_len(n_samples)) {
            comparison_table <- data.frame(delta = rnorm(n_iterations, 
                mean = true_instability[i], sd = max(0.005, true_instability[i] * 
                  0.35)), instability = true_instability[i], 
                iteration = seq_len(n_iterations))
            support <- sparq_support_from_comparisons(comparison_table = comparison_table, 
                n_iterations = n_iterations, reference_scale = 1, 
                min_iterations = min(10, n_iterations))
            sample_rows[[i]] <- data.frame(benchmark_id = paste0("rep", 
                r, "_sample", i), truth = truth[i], difficulty = difficulty[i], 
                support_quality = as.numeric(support$support_quality), 
                instability_median = as.numeric(support$instability_median), 
                instability_ci_low = as.numeric(support$instability_ci_low), 
                instability_ci_high = as.numeric(support$instability_ci_high), 
                n_iterations = n_iterations, stringsAsFactors = FALSE)
        }
        samples <- do.call(rbind, sample_rows)
        train_index <- unlist(lapply(unique(samples$truth), function(class) {
            ids <- which(samples$truth == class)
            sample(ids, size = max(1, floor(length(ids) * train_fraction)))
        }))
        held_index <- setdiff(seq_len(nrow(samples)), train_index)
        training <- samples[train_index, , drop = FALSE]
        held_out <- samples[held_index, , drop = FALSE]
        threshold_fit <- choose_threshold(training$truth, training$support_quality)
        threshold <- threshold_fit$threshold
        threshold_ci <- bootstrap_threshold_ci(truth = training$truth, 
            score = training$support_quality, B = threshold_bootstrap_B, 
            confidence = threshold_confidence, seed = seed + 
                r)
        held_out$predicted_stable <- as.integer(held_out$support_quality >= 
            threshold)
        held_metrics <- metric_at_threshold(held_out$truth, held_out$support_quality, 
            threshold)
        held_out$seed <- seed + r
        held_out$replicate <- r
        held_out$split <- "held_out"
        held_out$threshold <- threshold
        held_out$threshold_ci_low <- threshold_ci[["threshold_ci_low"]]
        held_out$threshold_ci_high <- threshold_ci[["threshold_ci_high"]]
        held_metrics$seed <- seed + r
        held_metrics$replicate <- r
        held_metrics$threshold <- threshold
        held_metrics$threshold_ci_low <- threshold_ci[["threshold_ci_low"]]
        held_metrics$threshold_ci_high <- threshold_ci[["threshold_ci_high"]]
        held_metrics$n_training <- nrow(training)
        held_metrics$n_held_out <- nrow(held_out)
        all_results[[r]] <- held_out
        all_metrics[[r]] <- held_metrics
    }
    replicate_results <- do.call(rbind, all_results)
    replicate_metrics <- do.call(rbind, all_metrics)
    threshold_summary <- data.frame(threshold_median = median(replicate_metrics$threshold), 
        threshold_ci_low = quantile(replicate_metrics$threshold, 
            (1 - threshold_confidence)/2), threshold_ci_high = quantile(replicate_metrics$threshold, 
            1 - (1 - threshold_confidence)/2), n_replicates = nrow(replicate_metrics), 
        B = threshold_bootstrap_B, confidence = threshold_confidence)
    list(replicate_results = replicate_results, replicate_metrics = replicate_metrics, 
        held_out_results = replicate_results, held_out_metrics = replicate_metrics, 
        threshold_summary = threshold_summary, difficulty_parameters = difficulty_parameters)
}
