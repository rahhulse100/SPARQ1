# Reliability inference that remains valid under repeated looks.

sparq_reliability_events <- function(comparison_table, tolerance,
                                     reference_scale = 1) {
  if (!is.data.frame(comparison_table) || !nrow(comparison_table) ||
      !all(c("status", "instability") %in% names(comparison_table))) {
    stop(
      "comparison_table must contain attempted perturbations with status and instability.",
      call. = FALSE
    )
  }
  sparq_check_number(tolerance, "tolerance", 0)
  sparq_check_number(reference_scale, "reference_scale")

  scale <- max(abs(reference_scale), 1e-6)
  distance <- suppressWarnings(as.numeric(comparison_table$instability))
  completed <- !is.na(comparison_table$status) &
    comparison_table$status == "ok" &
    is.finite(distance) &
    distance >= 0
  preserved <- completed & distance / scale <= tolerance

  list(
    n_attempted = nrow(comparison_table),
    completed = completed,
    preserved = preserved,
    reference_scale = scale
  )
}

sparq_beta_mixture_bounds <- function(successes, trials, confidence,
                                      prior_success = 0.5,
                                      prior_failure = 0.5) {
  sparq_check_number(successes, "successes", 0, trials, integer = TRUE)
  sparq_check_number(trials, "trials", 1, integer = TRUE)
  sparq_check_number(confidence, "confidence", .Machine$double.eps,
    1 - .Machine$double.eps
  )
  sparq_check_number(prior_success, "prior_success", .Machine$double.eps)
  sparq_check_number(prior_failure, "prior_failure", .Machine$double.eps)

  log_threshold <- -log(1 - confidence)
  log_mixture_ratio <- function(probability) {
    probability <- min(max(probability, .Machine$double.eps),
      1 - .Machine$double.eps
    )
    lbeta(successes + prior_success,
      trials - successes + prior_failure
    ) - lbeta(prior_success, prior_failure) -
      successes * log(probability) -
      (trials - successes) * log1p(-probability)
  }
  objective <- function(probability) log_mixture_ratio(probability) - log_threshold
  epsilon <- .Machine$double.eps^0.5

  if (successes == 0L) {
    return(c(
      lower = 0,
      upper = stats::uniroot(objective, c(epsilon, 1 - epsilon), tol = 1e-10)$root
    ))
  }
  if (successes == trials) {
    return(c(
      lower = stats::uniroot(objective, c(epsilon, 1 - epsilon), tol = 1e-10)$root,
      upper = 1
    ))
  }

  maximum_likelihood <- successes / trials
  c(
    lower = stats::uniroot(
      objective,
      c(epsilon, maximum_likelihood),
      tol = 1e-10
    )$root,
    upper = stats::uniroot(
      objective,
      c(maximum_likelihood, 1 - epsilon),
      tol = 1e-10
    )$root
  )
}

#' Compute an anytime-valid reliability confidence sequence
#'
#' Converts completed perturbation comparisons into binary preservation events
#' using a prespecified normalized distance tolerance. The returned beta-binomial
#' mixture confidence sequence remains valid across repeated looks at the same
#' accumulating perturbation stream, so it can support adaptive stopping.
#'
#' @param comparison_table A data frame with one row per attempted perturbation
#'   and `status` and `instability` columns.
#' @param tolerance Maximum normalized distance that still preserves the result.
#' @param required_reliability Minimum preservation probability required for a
#'   supported decision.
#' @param confidence Simultaneous confidence level across all sequential looks.
#' @param reference_scale Positive scale used to normalize comparator distance.
#' @param prior_success,prior_failure Beta-mixture prior parameters. The default
#'   Jeffreys mixture is fixed and recorded for reproducibility.
#'
#' @return A data frame with one row per sequential look and an anytime-valid
#'   lower and upper confidence bound. Failures count as non-preserving attempts.
#' @export
sparq_reliability_confidence_sequence <- function(
    comparison_table, tolerance, required_reliability = 0.9,
    confidence = 0.95, reference_scale = 1,
    prior_success = 0.5, prior_failure = 0.5) {
  sparq_check_number(required_reliability, "required_reliability",
    .Machine$double.eps, 1 - .Machine$double.eps
  )
  sparq_check_number(confidence, "confidence", .Machine$double.eps,
    1 - .Machine$double.eps
  )
  events <- sparq_reliability_events(
    comparison_table = comparison_table,
    tolerance = tolerance,
    reference_scale = reference_scale
  )

  rows <- lapply(seq_len(events$n_attempted), function(look) {
    n_preserved <- sum(events$preserved[seq_len(look)])
    bounds <- sparq_beta_mixture_bounds(
      successes = n_preserved,
      trials = look,
      confidence = confidence,
      prior_success = prior_success,
      prior_failure = prior_failure
    )
    decision <- if (bounds[["lower"]] >= required_reliability) {
      "supported"
    } else if (bounds[["upper"]] < required_reliability) {
      "not_supported"
    } else {
      "inconclusive"
    }
    data.frame(
      n_attempted = look,
      n_completed = sum(events$completed[seq_len(look)]),
      n_failed = look - sum(events$completed[seq_len(look)]),
      n_preserved = n_preserved,
      preservation_rate = n_preserved / look,
      lower_confidence_bound = bounds[["lower"]],
      upper_confidence_bound = bounds[["upper"]],
      tolerance = tolerance,
      required_reliability = required_reliability,
      confidence = confidence,
      reference_scale = events$reference_scale,
      prior_success = prior_success,
      prior_failure = prior_failure,
      method = "beta_binomial_mixture_confidence_sequence",
      decision = decision,
      stringsAsFactors = FALSE
    )
  })
  sparq_bind_rows(rows)
}
