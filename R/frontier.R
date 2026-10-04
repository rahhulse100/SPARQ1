# Reliability frontiers across explicitly chosen retained fractions.

#' Estimate a spatial reliability frontier
#'
#' Runs the same analysis under a prespecified grid of retained spot fractions.
#' Each grid point receives its own perturbation assessment and reliability
#' certificate. SPARQ never interpolates a critical retention value between
#' evaluated fractions.
#'
#' @param data Spot-level data frame.
#' @param analysis_function Function that reruns the analysis on retained data.
#' @param output_type Type of returned analysis result.
#' @param retentions Unique retained fractions in `(0, 1]`.
#' @param reliability_tolerance Maximum normalized comparator distance defining
#'   preservation.
#' @param required_reliability Minimum preservation probability.
#' @param confidence Family-wise confidence level across the prespecified
#'   retention grid. SPARQ allocates this across grid points with a Bonferroni
#'   correction before defining the critical retention.
#' @param loss_model Optional explicit `sparq_loss_model`.
#' @param seed Master random seed.
#' @param ... Arguments passed to `sparq_assess_workflow()`.
#'
#' @return A `sparq_reliability_frontier` with one reliability assessment per
#'   retention, a tabular frontier, and the lowest evaluated retention that is
#'   supported at every higher evaluated retention.
#' @export
sparq_assess_frontier <- function(
    data, analysis_function,
    output_type = c("scalar", "ranked", "partition", "graph"),
    retentions = c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
    stress_model = c("uniform_random", "contiguous_hole", "none", "custom"),
    reliability_tolerance, required_reliability = 0.9, confidence = 0.95,
    top_k = 100, x_col = NULL, y_col = NULL, custom_function = NULL,
    loss_model = NULL, result_id = "result", reference_scale = 1,
    seed = 1, verbose = TRUE, ...) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a nonempty data.frame.", call. = FALSE)
  }
  if (!is.function(analysis_function)) {
    stop("analysis_function must be a function.", call. = FALSE)
  }
  output_type <- match.arg(output_type)
  stress_model <- match.arg(stress_model)
  if (!is.numeric(retentions) || !length(retentions) || any(!is.finite(retentions)) ||
      any(retentions <= 0 | retentions > 1) || anyDuplicated(retentions)) {
    stop("retentions must be unique finite values in (0, 1].", call. = FALSE)
  }
  sparq_check_number(reliability_tolerance, "reliability_tolerance", 0)
  sparq_check_number(required_reliability, "required_reliability",
    .Machine$double.eps, 1 - .Machine$double.eps
  )
  sparq_check_number(confidence, "confidence", .Machine$double.eps,
    1 - .Machine$double.eps
  )
  sparq_check_number(seed, "seed", 0, integer = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(loss_model) && !inherits(loss_model, "sparq_loss_model")) {
    stop("loss_model must be created by sparq_define_loss_model().", call. = FALSE)
  }

  retentions <- sort(as.numeric(retentions), decreasing = FALSE)
  per_retention_confidence <- 1 - (1 - confidence) / length(retentions)
  assessments <- vector("list", length(retentions))
  names(assessments) <- format(retentions, trim = TRUE, scientific = FALSE)
  rows <- vector("list", length(retentions))

  for (index in seq_along(retentions)) {
    retention <- retentions[index]
    sparq_inform(
      paste0("Frontier ", index, "/", length(retentions), ": ",
        format(retention * 100, trim = TRUE), "% retention."),
      verbose
    )
    assessment <- sparq_assess_workflow(
      data = data,
      analysis_function = analysis_function,
      output_type = output_type,
      stress_model = stress_model,
      retention = retention,
      top_k = top_k,
      x_col = x_col,
      y_col = y_col,
      custom_function = custom_function,
      loss_model = loss_model,
      result_id = paste0(result_id, "_retention_", format(retention, trim = TRUE)),
      reference_scale = reference_scale,
      seed = sparq_seed(seed, index, 71L),
      reliability_tolerance = reliability_tolerance,
      required_reliability = required_reliability,
      reliability_confidence = per_retention_confidence,
      verbose = verbose,
      ...
    )
    certificate <- assessment$reliability_decision
    if (!is.data.frame(certificate) || nrow(certificate) != 1L) {
      stop("Frontier assessment did not return one reliability certificate.",
        call. = FALSE
      )
    }
    assessments[[index]] <- assessment
    rows[[index]] <- data.frame(
      retention = retention,
      support_quality = assessment$summary$support_quality[1],
      n_attempted = certificate$n_attempted[1],
      n_completed = certificate$n_completed[1],
      n_failed = certificate$n_failed[1],
      n_preserved = certificate$n_preserved[1],
      preservation_rate = certificate$preservation_rate[1],
      lower_confidence_bound = certificate$lower_confidence_bound[1],
      upper_confidence_bound = certificate$upper_confidence_bound[1],
      per_retention_confidence = certificate$confidence[1],
      decision = certificate$decision[1],
      stringsAsFactors = FALSE
    )
  }

  frontier <- sparq_bind_rows(rows)
  frontier$supported_at_or_above <- vapply(seq_len(nrow(frontier)), function(index) {
    all(frontier$decision[seq.int(index, nrow(frontier))] == "supported")
  }, logical(1))
  eligible <- frontier$retention[frontier$supported_at_or_above]
  critical_retention <- if (length(eligible)) min(eligible) else NA_real_
  observed_monotonicity <- all(diff(frontier$preservation_rate) >= -1e-12)

  out <- list(
    frontier = frontier,
    assessments = assessments,
    critical_retention = critical_retention,
    reliability_tolerance = reliability_tolerance,
    required_reliability = required_reliability,
    confidence = confidence,
    per_retention_confidence = per_retention_confidence,
    output_type = output_type,
    stress_model = stress_model,
    loss_model = loss_model,
    observed_monotonicity = observed_monotonicity,
    seed = seed
  )
  class(out) <- "sparq_reliability_frontier"
  out
}

print.sparq_reliability_frontier <- function(x, ...) {
  cat("<SPARQ reliability frontier>\n")
  cat("  Evaluated retentions: ", paste(format(x$frontier$retention), collapse = ", "), "\n", sep = "")
  if (is.na(x$critical_retention)) {
    cat("  Critical retention: not certified on this evaluated grid\n")
  } else {
    cat("  Critical retention: ", format(x$critical_retention), "\n", sep = "")
  }
  cat("  Tolerance: ", x$reliability_tolerance,
    "; required reliability: ", x$required_reliability, "\n", sep = "")
  cat("  Family-wise confidence: ", x$confidence,
    "; per-retention confidence: ", x$per_retention_confidence, "\n", sep = "")
  invisible(x)
}

summary.sparq_reliability_frontier <- function(object, ...) object$frontier

plot.sparq_reliability_frontier <- function(x, ...) {
  frontier <- x$frontier
  graphics::plot(
    frontier$retention,
    frontier$preservation_rate,
    ylim = c(0, 1), xlab = "Retained spot fraction",
    ylab = "Preservation probability", pch = 16,
    main = "SPARQ reliability frontier", ...
  )
  graphics::arrows(
    frontier$retention, frontier$lower_confidence_bound,
    frontier$retention, frontier$upper_confidence_bound,
    angle = 90, code = 3, length = 0.04
  )
  graphics::abline(h = x$required_reliability, lty = 2, col = "#B13A3A")
  invisible(x)
}
