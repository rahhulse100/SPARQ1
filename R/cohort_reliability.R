# Hierarchical inference across independent sections.

sparq_extract_cohort_certificates <- function(cohort) {
  assessments <- if (inherits(cohort, "sparq_cohort_assessment")) {
    cohort$assessments
  } else if (is.list(cohort) && !is.null(cohort$assessments)) {
    cohort$assessments
  } else if (is.list(cohort)) {
    cohort
  } else {
    stop("cohort must be a SPARQ cohort assessment or a named list of assessments.",
      call. = FALSE
    )
  }
  if (!length(assessments) || is.null(names(assessments))) {
    stop("cohort assessments must be a nonempty named list.", call. = FALSE)
  }
  rows <- lapply(names(assessments), function(sample_id) {
    assessment <- assessments[[sample_id]]
    certificate <- assessment$reliability_decision
    if (!is.data.frame(certificate) || nrow(certificate) != 1L) return(NULL)
    data.frame(
      sample_id = sample_id,
      n_attempted = certificate$n_attempted[1],
      n_preserved = certificate$n_preserved[1],
      n_completed = certificate$n_completed[1],
      n_failed = certificate$n_failed[1],
      preservation_rate = certificate$preservation_rate[1],
      tolerance = certificate$tolerance[1],
      required_reliability = certificate$required_reliability[1],
      confidence = certificate$confidence[1],
      reference_scale = certificate$reference_scale[1],
      decision = certificate$decision[1],
      stringsAsFactors = FALSE
    )
  })
  certificates <- sparq_bind_rows(rows)
  if (!nrow(certificates)) {
    stop("No assessments contain reliability certificates. Supply reliability_tolerance when assessing the cohort.",
      call. = FALSE
    )
  }
  certificates
}

sparq_beta_binomial_moments <- function(successes, trials) {
  rate <- successes / trials
  observed_mean_rate <- sum(successes) / sum(trials)
  mean_rate <- min(max(observed_mean_rate, 1e-6), 1 - 1e-6)
  if (length(rate) < 2L) {
    return(list(alpha = NA_real_, beta = NA_real_, concentration = NA_real_,
      between_section_variance = NA_real_, mean_reliability = observed_mean_rate,
      estimation_status = "one_section_no_between_section_estimate"))
  }
  observed_variance <- stats::var(rate)
  expected_binomial_variance <- mean(rate * (1 - rate) / trials)
  between_section_variance <- max(0, observed_variance - expected_binomial_variance)
  maximum_variance <- mean_rate * (1 - mean_rate)
  if (!is.finite(between_section_variance) || between_section_variance <= 1e-10) {
    concentration <- 1e6
    status <- "near_zero_between_section_variance"
  } else {
    concentration <- maximum_variance / between_section_variance - 1
    if (!is.finite(concentration) || concentration <= 0) {
      concentration <- 1e-3
      status <- "high_between_section_variance_boundary"
    } else {
      status <- "method_of_moments"
    }
  }
  list(
    alpha = mean_rate * concentration,
    beta = (1 - mean_rate) * concentration,
    concentration = concentration,
    between_section_variance = between_section_variance,
    mean_reliability = observed_mean_rate,
    estimation_status = status
  )
}

#' Infer reliability across independent spatial sections
#'
#' Fits an empirical-Bayes beta-binomial model to section-level preservation
#' counts. The model separates binomial uncertainty from finite perturbation
#' runs within a section from observed variation in reliability between sections.
#' It does not treat spots as biological replicates.
#'
#' @param cohort A `sparq_cohort_assessment` or named list of assessments that
#'   include `reliability_decision` results.
#' @param confidence Posterior interval level for section-level reliability.
#' @param section_metadata Optional data frame with one row per sample.
#' @param sample_id_col Name of the sample-ID column in `section_metadata`.
#'
#' @return A `sparq_cohort_reliability` object containing raw section results,
#'   empirical-Bayes section estimates, and a cohort summary. Posterior intervals
#'   are explicitly labeled as empirical-Bayes intervals.
#' @export
sparq_infer_cohort_reliability <- function(
    cohort, confidence = 0.95,
    section_metadata = NULL, sample_id_col = "sample_id") {
  sparq_check_number(confidence, "confidence", .Machine$double.eps,
    1 - .Machine$double.eps
  )
  certificates <- sparq_extract_cohort_certificates(cohort)
  if (length(unique(certificates$tolerance)) != 1L ||
      length(unique(certificates$required_reliability)) != 1L) {
    stop("All cohort assessments must use the same tolerance and required_reliability.",
      call. = FALSE
    )
  }
  if (any(certificates$n_attempted < 1L) ||
      any(certificates$n_preserved < 0L) ||
      any(certificates$n_preserved > certificates$n_attempted)) {
    stop("Cohort reliability counts are invalid.", call. = FALSE)
  }
  moments <- sparq_beta_binomial_moments(
    successes = certificates$n_preserved,
    trials = certificates$n_attempted
  )
  alpha_tail <- (1 - confidence) / 2
  section_results <- certificates
  if (is.finite(moments$alpha) && is.finite(moments$beta)) {
    posterior_alpha <- section_results$n_preserved + moments$alpha
    posterior_beta <- section_results$n_attempted - section_results$n_preserved + moments$beta
    section_results$posterior_reliability <- posterior_alpha / (posterior_alpha + posterior_beta)
    section_results$posterior_ci_low <- stats::qbeta(alpha_tail, posterior_alpha, posterior_beta)
    section_results$posterior_ci_high <- stats::qbeta(1 - alpha_tail, posterior_alpha, posterior_beta)
    section_results$posterior_decision <- ifelse(
      section_results$posterior_ci_low >= section_results$required_reliability,
      "supported", ifelse(
        section_results$posterior_ci_high < section_results$required_reliability,
        "not_supported", "inconclusive"
      )
    )
  } else {
    section_results$posterior_reliability <- NA_real_
    section_results$posterior_ci_low <- NA_real_
    section_results$posterior_ci_high <- NA_real_
    section_results$posterior_decision <- "not_estimated"
  }

  if (!is.null(section_metadata)) {
    if (!is.data.frame(section_metadata) || !sample_id_col %in% names(section_metadata) ||
        anyDuplicated(section_metadata[[sample_id_col]]) || anyNA(section_metadata[[sample_id_col]])) {
      stop("section_metadata must contain one unique, nonmissing sample-ID column.",
        call. = FALSE
      )
    }
    metadata <- section_metadata
    names(metadata)[names(metadata) == sample_id_col] <- "sample_id"
    section_results <- merge(section_results, metadata,
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    section_results <- section_results[match(certificates$sample_id, section_results$sample_id), , drop = FALSE]
  }

  cohort_summary <- data.frame(
    n_sections = nrow(section_results),
    total_attempted = sum(section_results$n_attempted),
    total_preserved = sum(section_results$n_preserved),
    pooled_preservation_rate = sum(section_results$n_preserved) / sum(section_results$n_attempted),
    mean_section_preservation_rate = mean(section_results$preservation_rate),
    between_section_variance = moments$between_section_variance,
    beta_alpha = moments$alpha,
    beta_beta = moments$beta,
    beta_concentration = moments$concentration,
    estimation_status = moments$estimation_status,
    tolerance = unique(certificates$tolerance),
    required_reliability = unique(certificates$required_reliability),
    confidence = confidence,
    stringsAsFactors = FALSE
  )
  out <- list(
    section_reliability = section_results,
    cohort_summary = cohort_summary,
    method = "empirical_bayes_beta_binomial",
    interval_interpretation = "Empirical-Bayes posterior interval for section-level preservation probability",
    section_metadata_included = !is.null(section_metadata)
  )
  class(out) <- "sparq_cohort_reliability"
  out
}

print.sparq_cohort_reliability <- function(x, ...) {
  summary <- x$cohort_summary
  cat("<SPARQ cohort reliability>\n")
  cat("  Sections: ", summary$n_sections[1], "; pooled preservation: ",
    format(summary$pooled_preservation_rate[1], digits = 4), "\n", sep = "")
  cat("  Between-section variance: ",
    format(summary$between_section_variance[1], digits = 4),
    " (", summary$estimation_status[1], ")\n", sep = "")
  invisible(x)
}

summary.sparq_cohort_reliability <- function(object, ...) object$section_reliability
as.data.frame.sparq_cohort_reliability <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
  x$section_reliability
}
