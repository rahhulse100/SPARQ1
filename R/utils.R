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

# User-facing progress and result-display helpers. These deliberately use base
# R so a standard assessment remains lightweight and usable on a cluster.

sparq_inform <- function(message, verbose = TRUE) {
  if (isTRUE(verbose)) {
    message("SPARQ | ", message)
  }
  invisible(NULL)
}

sparq_progress <- function(completed, total, verbose = TRUE,
                           prefix = "Perturbations") {
  if (!isTRUE(verbose) || total < 1L) {
    return(invisible(NULL))
  }

  width <- 24L
  filled <- min(width, floor(width * completed / total))
  bar <- paste0(
    "[",
    paste(rep("=", filled), collapse = ""),
    paste(rep(".", width - filled), collapse = ""),
    "]"
  )
  message(
    "SPARQ | ", prefix, " ", bar, " ",
    completed, "/", total
  )
  invisible(NULL)
}

sparq_assessment_status <- function(summary_table) {
  if (!is.data.frame(summary_table) || nrow(summary_table) != 1L) {
    return("unavailable")
  }
  if (isTRUE(summary_table$insufficient_support[1])) {
    return("insufficient support")
  }
  "adequate support"
}

sparq_print_assessment <- function(x, label = "assessment") {
  summary_table <- x$summary
  if (!is.data.frame(summary_table) || nrow(summary_table) != 1L) {
    cat("<SPARQ ", label, ">\n", sep = "")
    return(invisible(x))
  }

  result_id <- if (!is.null(x$result_id)) as.character(x$result_id) else "result"
  output_type <- if (!is.null(x$output_type)) as.character(x$output_type) else "custom"
  cat("<SPARQ ", label, ": ", result_id, ">\n", sep = "")
  cat("  Output type: ", output_type, "\n", sep = "")
  cat(
    "  Support quality: ", format(summary_table$support_quality[1], digits = 4),
    " (", sparq_assessment_status(summary_table), ")\n", sep = ""
  )
  cat(
    "  Completed perturbations: ", summary_table$n_iterations[1],
    "/", summary_table$n_requested[1],
    "; failed: ", summary_table$n_failed[1], "\n", sep = ""
  )
  cat(
    "  Median instability: ", format(summary_table$instability_median[1], digits = 4),
    "; reproducibility: ", format(summary_table$reproducibility_score[1], digits = 4),
    "\n", sep = ""
  )
  if (isTRUE(summary_table$insufficient_support[1])) {
    cat("  Reason: ", summary_table$support_reason[1], "\n", sep = "")
  }
  invisible(x)
}

print.sparq_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "assessment")
}

print.sparq_workflow_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "workflow assessment")
}

print.sparq_script_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "script assessment")
}

print.sparq_adaptive_assessment <- function(x, ...) {
  sparq_print_assessment(x, label = "adaptive assessment")
  if (is.data.frame(x$adaptive_decision) && nrow(x$adaptive_decision) == 1L) {
    cat(
      "  Adaptive decision: ", x$adaptive_decision$decision[1],
      " (", x$adaptive_decision$stop_reason[1], ")\n", sep = ""
    )
  }
  invisible(x)
}

print.sparq_report <- function(x, ...) {
  cat("<SPARQ report: ", x$result_id, ">\n", sep = "")
  if (is.data.frame(x$support_summary) && nrow(x$support_summary) == 1L) {
    cat(
      "  Support quality: ", format(x$support_summary$support_quality[1], digits = 4),
      "\n", sep = ""
    )
  }
  if (length(x$output_files)) {
    cat("  Exported files: ", length(x$output_files), "\n", sep = "")
  } else {
    cat("  Exported files: none\n")
  }
  invisible(x)
}
