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

