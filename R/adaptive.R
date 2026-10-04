# Adaptive perturbation assessment.

sparq_assess_adaptive <-
function(data, analysis_function, comparator,
    stress_model = "uniform_random", retention = NULL,
    support_threshold = NULL, batch_size = NULL, min_iterations = NULL,
    max_iterations = NULL, confidence = 0.95,
    support_bootstrap_B = NULL, result_id = "result",
    x_col = NULL, y_col = NULL, custom_function = NULL,
    reference_scale = 1, seed = 1,
    failure_action = c("record", "stop"), full_result = NULL,
    keep_perturbed_results = FALSE,
    min_reproducibility = 0.5, max_relative_bias = 0.5,
    max_relative_uncertainty = 0.25,
    instability_bootstrap_B = NULL, verbose = NULL, preset = "standard",
    reliability_tolerance = NULL, required_reliability = 0.9,
    reliability_confidence = NULL, loss_model = NULL) {

    failure_action <- match.arg(failure_action)
    preset_values <- sparq_resolve_preset(
        preset = preset, retention = retention, n_iterations = max_iterations,
        min_iterations = min_iterations,
        instability_bootstrap_B = instability_bootstrap_B,
        support_bootstrap_B = support_bootstrap_B, batch_size = batch_size,
        verbose = verbose
    )
    retention <- preset_values$retention
    batch_size <- preset_values$batch_size
    min_iterations <- preset_values$min_iterations
    max_iterations <- preset_values$n_iterations
    support_bootstrap_B <- preset_values$support_bootstrap_B
    instability_bootstrap_B <- preset_values$instability_bootstrap_B
    verbose <- preset_values$verbose
    sparq_inform_preset(preset_values)

    if (!is.function(analysis_function) || !is.function(comparator)) {
        stop("analysis_function and comparator must be functions.", call. = FALSE)
    }
    if (!is.data.frame(data) || !nrow(data)) {
        stop("data must be a nonempty data.frame.", call. = FALSE)
    }
    if (!is.null(loss_model) && !inherits(loss_model, "sparq_loss_model")) {
        stop("loss_model must be created by sparq_define_loss_model().",
            call. = FALSE)
    }
    if (is.null(support_threshold) && is.null(reliability_tolerance)) {
        stop("Supply support_threshold for calibrated-score stopping or reliability_tolerance for tolerance-based stopping.",
            call. = FALSE)
    }
    if (!is.null(support_threshold) && !is.null(reliability_tolerance)) {
        stop("Supply either support_threshold or reliability_tolerance, not both.",
            call. = FALSE)
    }
    decision_mode <- if (is.null(reliability_tolerance)) "calibrated_support" else "reliability_certificate"
    if (decision_mode == "calibrated_support") {
        sparq_check_number(support_threshold, "support_threshold", 0, 1)
    } else {
        sparq_check_number(reliability_tolerance, "reliability_tolerance", 0)
        sparq_check_number(required_reliability, "required_reliability",
            .Machine$double.eps, 1 - .Machine$double.eps)
        if (is.null(reliability_confidence)) reliability_confidence <- confidence
        sparq_check_number(reliability_confidence, "reliability_confidence",
            .Machine$double.eps, 1 - .Machine$double.eps)
    }
    sparq_check_number(batch_size, "batch_size", 1, integer = TRUE)
    sparq_check_number(min_iterations, "min_iterations", 2, integer = TRUE)
    sparq_check_number(max_iterations, "max_iterations", min_iterations, integer = TRUE)
    sparq_check_number(confidence, "confidence", .Machine$double.eps, 1)
    sparq_check_number(support_bootstrap_B, "support_bootstrap_B", 1, integer = TRUE)
    sparq_check_number(reference_scale, "reference_scale")
    sparq_check_number(seed, "seed", 0, integer = TRUE)
    sparq_check_number(min_reproducibility, "min_reproducibility",
        .Machine$double.eps, 1)
    sparq_check_number(max_relative_bias, "max_relative_bias", 0)
    sparq_check_number(max_relative_uncertainty, "max_relative_uncertainty", 0)
    sparq_check_number(instability_bootstrap_B, "instability_bootstrap_B", 1,
        integer = TRUE)

    if (is.null(full_result)) {
        sparq_inform("Running reference analysis for adaptive assessment.", verbose)
        full_result <- sparq_with_seed(
            sparq_seed(seed, 0, 1),
            analysis_function(data)
        )
    }

    summarize <- function(comparison_table, n_requested, summary_seed,
        bootstrap_B = instability_bootstrap_B) {
        sparq_support_from_comparisons(
            comparison_table = comparison_table,
            n_iterations = n_requested,
            reference_scale = reference_scale,
            min_iterations = min_iterations,
            min_reproducibility = min_reproducibility,
            max_relative_bias = max_relative_bias,
            max_relative_uncertainty = max_relative_uncertainty,
            instability_bootstrap_B = bootstrap_B,
            confidence = confidence,
            seed = summary_seed
        )
    }

    bootstrap_support_interval <- function(comparison_table, n_requested,
        interval_seed) {
        valid <- comparison_table[
            comparison_table$status == "ok" &
                is.finite(comparison_table$instability),
            , drop = FALSE
        ]

        if ("delta" %in% names(valid)) {
            valid <- valid[is.finite(valid$delta), , drop = FALSE]
        }

        if (nrow(valid) < 2L) {
            return(c(lower = NA_real_, upper = NA_real_))
        }

        draws <- sparq_with_seed(interval_seed, vapply(
            seq_len(support_bootstrap_B),
            function(draw_id) {
                draw <- valid[sample.int(nrow(valid), nrow(valid), replace = TRUE),
                    , drop = FALSE]
                draw$iteration <- seq_len(nrow(draw))
                summarize(
                    comparison_table = draw,
                    n_requested = n_requested,
                    summary_seed = sparq_seed(seed, n_requested, draw_id + 100),
                    bootstrap_B = 1
                )$support_quality[1]
            },
            numeric(1)
        ))

        alpha <- 1 - confidence
        c(
            lower = unname(stats::quantile(draws, alpha / 2, names = FALSE)),
            upper = unname(stats::quantile(draws, 1 - alpha / 2, names = FALSE))
        )
    }

    outputs <- vector("list", max_iterations)
    history <- list()
    stop_reason <- NULL
    decision <- "inconclusive"

    for (iteration in seq_len(max_iterations)) {
        output <- tryCatch(
            sparq_with_seed(sparq_seed(seed, iteration, 2), {
                perturbed_data <- sparq_stress_model(
                    data = data,
                    retention = retention,
                    model = stress_model,
                    x_col = x_col,
                    y_col = y_col,
                    custom_function = custom_function,
                    iteration = iteration,
                    seed = seed,
                    loss_model = loss_model
                )
                perturbed_result <- analysis_function(perturbed_data)
                comparison <- sparq_comparison_row(
                    comparator = comparator,
                    full_result = full_result,
                    perturbed_result = perturbed_result,
                    iteration = iteration
                )
                comparison$n_retained <- nrow(perturbed_data)
                comparison$loss_model <- if (is.null(loss_model)) stress_model else loss_model$type
                list(
                    comparison = comparison,
                    perturbed_result = if (keep_perturbed_results) perturbed_result else NULL,
                    loss_audit = attr(perturbed_data, "sparq_loss_audit")
                )
            }),
            error = function(error) {
                if (failure_action == "stop") {
                    stop(
                        paste0("Iteration ", iteration, ": ", conditionMessage(error)),
                        call. = FALSE
                    )
                }
                list(
                    comparison = data.frame(
                        instability = NA_real_,
                        iteration = iteration,
                        status = "failed",
                        error_message = conditionMessage(error),
                        n_retained = NA_integer_,
                        stringsAsFactors = FALSE
                    ),
                    perturbed_result = NULL
                )
            }
        )

        outputs[[iteration]] <- output

        evaluate_now <- iteration >= min_iterations &&
            (iteration %% batch_size == 0L || iteration == max_iterations)

        if (!evaluate_now) {
            next
        }

        comparisons <- sparq_bind_rows(lapply(outputs[seq_len(iteration)], `[[`, "comparison"))
        summary <- summarize(
            comparison_table = comparisons,
            n_requested = iteration,
            summary_seed = sparq_seed(seed, iteration, 3)
        )
        if (decision_mode == "reliability_certificate") {
            sequence <- sparq_reliability_confidence_sequence(
                comparison_table = comparisons,
                tolerance = reliability_tolerance,
                required_reliability = required_reliability,
                confidence = reliability_confidence,
                reference_scale = reference_scale
            )
            reliability <- sequence[nrow(sequence), , drop = FALSE]
            interval <- c(lower = reliability$lower_confidence_bound[1],
                upper = reliability$upper_confidence_bound[1])
            batch_decision <- reliability$decision[1]
        } else {
            reliability <- NULL
            interval <- bootstrap_support_interval(
                comparison_table = comparisons,
                n_requested = iteration,
                interval_seed = sparq_seed(seed, iteration, 4)
            )
            batch_decision <- if (is.finite(interval["lower"]) &&
                interval["lower"] > support_threshold) {
                "above_threshold"
            } else if (is.finite(interval["upper"]) &&
                interval["upper"] < support_threshold) {
                "below_threshold"
            } else {
                "inconclusive"
            }
        }

        history[[length(history) + 1L]] <- data.frame(
            n_attempted = iteration,
            n_completed = summary$n_iterations[1],
            n_failed = summary$n_failed[1],
            support_quality = summary$support_quality[1],
            support_ci_low = interval["lower"],
            support_ci_high = interval["upper"],
            support_threshold = if (is.null(support_threshold)) NA_real_ else support_threshold,
            decision_mode = decision_mode,
            reliability_tolerance = if (is.null(reliability_tolerance)) NA_real_ else reliability_tolerance,
            required_reliability = if (is.null(reliability_tolerance)) NA_real_ else required_reliability,
            decision = batch_decision,
            stringsAsFactors = FALSE
        )
        if (decision_mode == "reliability_certificate") {
            sparq_inform(
                paste0("Adaptive check after ", iteration,
                    " perturbations: preservation ",
                    format(reliability$preservation_rate[1], digits = 4),
                    "; confidence sequence [", format(interval["lower"], digits = 4),
                    ", ", format(interval["upper"], digits = 4),
                    "]; ", batch_decision, "."), verbose)
        } else {
            sparq_inform(
                paste0(
                    "Adaptive check after ", iteration, " perturbations: support ",
                    format(summary$support_quality[1], digits = 4), "; interval [",
                    format(interval["lower"], digits = 4), ", ",
                    format(interval["upper"], digits = 4), "]; ", batch_decision, "."
                ), verbose)
        }

        if (batch_decision != "inconclusive") {
            decision <- batch_decision
            stop_reason <- if (decision_mode == "reliability_certificate") {
                if (batch_decision == "supported") "reliability_lower_bound_met" else "reliability_upper_bound_below_requirement"
            } else if (batch_decision == "above_threshold") {
                "support_interval_above_threshold"
            } else {
                "support_interval_below_threshold"
            }
            break
        }
    }

    n_attempted <- which(vapply(outputs, Negate(is.null), logical(1)))[1]
    if (is.na(n_attempted)) {
        stop("No perturbation iterations completed.", call. = FALSE)
    }
    n_attempted <- max(which(vapply(outputs, Negate(is.null), logical(1))))
    outputs <- outputs[seq_len(n_attempted)]
    comparisons <- sparq_bind_rows(lapply(outputs, `[[`, "comparison"))
    summary <- summarize(
        comparison_table = comparisons,
        n_requested = n_attempted,
        summary_seed = sparq_seed(seed, n_attempted, 3)
    )

    final_reliability <- NULL
    reliability_sequence <- NULL
    if (decision_mode == "reliability_certificate") {
        reliability_sequence <- sparq_reliability_confidence_sequence(
            comparison_table = comparisons,
            tolerance = reliability_tolerance,
            required_reliability = required_reliability,
            confidence = reliability_confidence,
            reference_scale = reference_scale
        )
        final_reliability <- reliability_sequence[nrow(reliability_sequence), , drop = FALSE]
        for (name in names(final_reliability)) {
            summary[[paste0("reliability_", name)]] <- final_reliability[[name]][1]
        }
    }

    if (is.null(stop_reason)) {
        stop_reason <- "maximum_iterations_reached"
    }

    perturbed_results <- NULL
    if (keep_perturbed_results) {
        perturbed_results <- lapply(outputs, `[[`, "perturbed_result")
        names(perturbed_results) <- paste0("iteration_", seq_len(n_attempted))
    }
    loss_audits <- lapply(outputs, function(output) output$loss_audit)
    if (all(vapply(loss_audits, is.null, logical(1)))) {
        loss_audits <- NULL
    } else {
        names(loss_audits) <- paste0("iteration_", seq_len(n_attempted))
    }

    out <- list(
        result_id = result_id,
        full_result = full_result,
        perturbed_results = perturbed_results,
        loss_audits = loss_audits,
        comparisons = comparisons,
        summary = summary,
        reliability_decision = final_reliability,
        reliability_confidence_sequence = reliability_sequence,
        adaptive_history = sparq_bind_rows(history),
        adaptive_decision = data.frame(
            decision = decision,
            stop_reason = stop_reason,
            support_threshold = if (is.null(support_threshold)) NA_real_ else support_threshold,
            decision_mode = decision_mode,
            reliability_tolerance = if (is.null(reliability_tolerance)) NA_real_ else reliability_tolerance,
            required_reliability = if (is.null(reliability_tolerance)) NA_real_ else required_reliability,
            n_attempted = n_attempted,
            stringsAsFactors = FALSE
        ),
        settings = list(
            stress_model = stress_model,
            retention = retention,
            min_iterations = min_iterations,
            max_iterations = max_iterations,
            batch_size = batch_size,
            support_threshold = support_threshold,
            decision_mode = decision_mode,
            reliability_tolerance = reliability_tolerance,
            required_reliability = required_reliability,
            reliability_confidence = reliability_confidence,
            loss_model = loss_model,
            confidence = confidence,
            support_bootstrap_B = support_bootstrap_B,
            instability_bootstrap_B = instability_bootstrap_B,
            preset = preset_values$preset,
            reference_scale = reference_scale,
            seed = seed,
            keep_perturbed_results = keep_perturbed_results
        ),
        failures = comparisons[comparisons$status != "ok", , drop = FALSE]
    )
    class(out) <- unique(c("sparq_adaptive_assessment", "sparq_assessment", class(out)))
    sparq_inform(paste0("Adaptive assessment finished: ", stop_reason, "."), verbose)
    out
}
