sparq_validate_results <- function(data,
                                   required_cols = c("result_id", "theta_full", "theta_pert", "scenario", "iteration")) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing)) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(data$theta_full)) {
    stop("theta_full must be numeric.", call. = FALSE)
  }

  if (!is.numeric(data$theta_pert)) {
    stop("theta_pert must be numeric.", call. = FALSE)
  }

  invisible(TRUE)
}

sparq_safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

sparq_bound01 <- function(x) {
  pmax(0, pmin(1, x))
}

sparq_mode_sign_fraction <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)

  med <- stats::median(x, na.rm = TRUE)
  if (!is.finite(med) || med == 0) return(NA_real_)

  mean(sign(x) == sign(med))
}


sparq_validate_columns <- function(data, required_cols) {
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(
      "SPARQ input missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sparq_validate_ranked_feature <- function(
  data,
  result_id_col = "result_id",
  theta_full_col = "theta_full",
  theta_pert_col = "theta_pert",
  rank_full_col = "rank_full",
  rank_pert_col = "rank_pert",
  scenario_col = "scenario",
  iteration_col = "iteration"
) {
  required <- c(
    result_id_col,
    theta_full_col,
    theta_pert_col,
    rank_full_col,
    rank_pert_col,
    scenario_col,
    iteration_col
  )

  sparq_validate_columns(data, required)

  if (!is.numeric(data[[theta_full_col]])) {
    stop(theta_full_col, " must be numeric.", call. = FALSE)
  }

  if (!is.numeric(data[[theta_pert_col]])) {
    stop(theta_pert_col, " must be numeric.", call. = FALSE)
  }

  if (!is.numeric(data[[rank_full_col]]) && !is.integer(data[[rank_full_col]])) {
    stop(rank_full_col, " must be numeric or integer.", call. = FALSE)
  }

  if (!is.numeric(data[[rank_pert_col]]) && !is.integer(data[[rank_pert_col]])) {
    stop(rank_pert_col, " must be numeric or integer.", call. = FALSE)
  }

  if (nrow(data) == 0) {
    stop("SPARQ input has zero rows.", call. = FALSE)
  }

  invisible(TRUE)
}

sparq_bound01 <- function(x) {
  pmax(0, pmin(1, x))
}

sparq_safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (length(unique(x[ok])) < 2) return(NA_real_)
  if (length(unique(y[ok])) < 2) return(NA_real_)

  suppressWarnings(
    stats::cor(x[ok], y[ok], method = method)
  )
}

sparq_mode_sign_fraction <- function(x) {
  x <- x[is.finite(x)]
  x <- sign(x[x != 0])
  if (length(x) == 0) return(0)
  max(mean(x > 0), mean(x < 0))
}


sparq_topk_jaccard <- function(full_ids, pert_ids, k = 100) {
  full_ids <- as.character(full_ids)
  pert_ids <- as.character(pert_ids)

  full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
  pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])

  if (length(full_top) == 0 && length(pert_top) == 0) return(NA_real_)

  length(intersect(full_top, pert_top)) / length(union(full_top, pert_top))
}



sparq_topk_replacement <- function(full_ids, pert_ids, k = 100) {
  full_ids <- as.character(full_ids)
  pert_ids <- as.character(pert_ids)

  full_top <- unique(full_ids[seq_len(min(k, length(full_ids)))])
  pert_top <- unique(pert_ids[seq_len(min(k, length(pert_ids)))])

  if (length(full_top) == 0) return(NA_real_)

  1 - (length(intersect(full_top, pert_top)) / length(full_top))
}


