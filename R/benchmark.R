sparq_known_truth_benchmark <-
function (n_samples = 200, n_iterations = 100, train_fraction = 0.69999999999999996, 
    seed = 20260912, stable_mean = 1, stable_sd = 0.029999999999999999, 
    unstable_mean = 0.80000000000000004, unstable_sd = 0.10000000000000001) 
{
    if (n_samples < 10) {
        stop("n_samples must be at least 10.", call. = FALSE)
    }
    if (n_iterations < 20) {
        stop("n_iterations must be at least 20.", call. = FALSE)
    }
    set.seed(seed)
    truth <- rep(c(1, 0), length.out = n_samples)
    sample_ids <- paste0("synthetic_", seq_len(n_samples))
    rows <- lapply(seq_len(n_samples), function(i) {
        full_value <- 1
        if (truth[i] == 1) {
            perturbed_values <- rnorm(n_iterations, mean = stable_mean, 
                sd = stable_sd)
        }
        else {
            perturbed_values <- rnorm(n_iterations, mean = unstable_mean, 
                sd = unstable_sd)
        }
        delta <- perturbed_values - full_value
        comparison_table <- data.frame(delta = delta, instability = abs(delta), 
            iteration = seq_len(n_iterations))
        summary <- sparq_support_from_comparisons(comparison_table = comparison_table, 
            n_iterations = n_iterations, reference_scale = full_value, 
            min_iterations = 50)
        data.frame(result_id = sample_ids[i], truth = truth[i], 
            support_quality = summary$support_quality, instability_median = summary$instability_median, 
            instability_ci_low = summary$instability_ci_low, 
            instability_ci_high = summary$instability_ci_high, 
            n_iterations = n_iterations, stringsAsFactors = FALSE)
    })
    results <- do.call(rbind, rows)
    rownames(results) <- NULL
    train_ids <- unlist(lapply(c(0, 1), function(class_value) {
        ids <- results$result_id[results$truth == class_value]
        sample(ids, size = floor(length(ids) * train_fraction))
    }))
    results$split <- ifelse(results$result_id %in% train_ids, 
        "train", "held_out")
    train <- results[results$split == "train", , drop = FALSE]
    held_out <- results[results$split == "held_out", , drop = FALSE]
    threshold_grid <- sort(unique(c(seq(0, 1, length.out = 1001), 
        train$support_quality)))
    threshold_stats <- do.call(rbind, lapply(threshold_grid, 
        function(threshold) {
            predicted <- train$support_quality >= threshold
            stable <- train$truth == 1
            unstable <- train$truth == 0
            sensitivity <- mean(predicted[stable])
            specificity <- mean(!predicted[unstable])
            data.frame(threshold = threshold, sensitivity = sensitivity, 
                specificity = specificity, balanced_accuracy = (sensitivity + 
                  specificity)/2)
        }))
    best_row <- threshold_stats[which.max(threshold_stats$balanced_accuracy), 
        , drop = FALSE]
    calibrated_threshold <- best_row$threshold
    held_out$predicted_stable <- held_out$support_quality >= 
        calibrated_threshold
    held_out_stable <- held_out$truth == 1
    held_out_unstable <- held_out$truth == 0
    held_out_sensitivity <- mean(held_out$predicted_stable[held_out_stable])
    held_out_specificity <- mean(!held_out$predicted_stable[held_out_unstable])
    held_out_balanced_accuracy <- (held_out_sensitivity + held_out_specificity)/2
    held_out$true_error <- ifelse(held_out$truth == 1, held_out$instability_median, 
        held_out$instability_median)
    reliability <- data.frame(support_quality = held_out$support_quality, 
        true_error = held_out$true_error, truth = held_out$truth, 
        predicted_stable = held_out$predicted_stable)
    list(sample_results = results, training_results = train, 
        held_out_results = held_out, threshold_table = threshold_stats, 
        calibrated_threshold = calibrated_threshold, held_out_metrics = data.frame(threshold = calibrated_threshold, 
            sensitivity = held_out_sensitivity, specificity = held_out_specificity, 
            balanced_accuracy = held_out_balanced_accuracy, n_training = nrow(train), 
            n_held_out = nrow(held_out)), reliability = reliability)
}
