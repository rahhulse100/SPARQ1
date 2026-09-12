sparq_bias <- function(data,
                       result_id_col = "result_id",
                       theta_full_col = "theta_full",
                       theta_pert_col = "theta_pert") {
  sparq_validate_results(data, required_cols = c(result_id_col, theta_full_col, theta_pert_col))

  ids <- unique(as.character(data[[result_id_col]]))

  out <- lapply(ids, function(id) {
    d <- data[as.character(data[[result_id_col]]) == id, , drop = FALSE]
    theta_full <- stats::median(d[[theta_full_col]], na.rm = TRUE)
    delta <- d[[theta_pert_col]] - d[[theta_full_col]]
    delta <- delta[is.finite(delta)]

    bias <- stats::median(delta, na.rm = TRUE)
    uncertainty <- stats::mad(delta, constant = 1.4826, na.rm = TRUE)
    direction_consistency <- sparq_mode_sign_fraction(delta)

    data.frame(
      result_id = id,
      theta_full = theta_full,
      n_perturbations = length(delta),
      bias = bias,
      uncertainty = uncertainty,
      direction_consistency = direction_consistency,
      sign_flip_rate = 1 - direction_consistency,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}
