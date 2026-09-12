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
