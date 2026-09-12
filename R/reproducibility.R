sparq_reproducibility_feature <- function(data,
                                          result_id_col = "result_id",
                                          theta_full_col = "theta_full",
                                          theta_pert_col = "theta_pert",
                                          scenario_col = "scenario",
                                          iteration_col = "iteration") {
  sparq_validate_results(
    data,
    required_cols = c(result_id_col, theta_full_col, theta_pert_col, scenario_col, iteration_col)
  )

  ids <- unique(as.character(data[[result_id_col]]))

  out <- lapply(ids, function(id) {
    d <- data[as.character(data[[result_id_col]]) == id, , drop = FALSE]
    theta_full <- d[[theta_full_col]]
    theta_pert <- d[[theta_pert_col]]
    delta <- theta_pert - theta_full

    data.frame(
      result_id = id,
      n_perturbations = nrow(d),
      n_scenarios = length(unique(as.character(d[[scenario_col]]))),
      pearson = sparq_safe_cor(theta_full, theta_pert, method = "pearson"),
      spearman = sparq_safe_cor(theta_full, theta_pert, method = "spearman"),
      rmse = sqrt(mean(delta^2, na.rm = TRUE)),
      mae = mean(abs(delta), na.rm = TRUE),
      cv_perturbed = stats::sd(theta_pert, na.rm = TRUE) / max(abs(stats::median(theta_pert, na.rm = TRUE)), .Machine$double.eps),
      sign_agreement = mean(sign(theta_full) == sign(theta_pert), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
