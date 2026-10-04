# Presets, conventional result methods, and cohort-facing helpers.

sparq_preset <- function(preset = c("quick", "standard", "thorough")) {
  preset <- match.arg(preset)
  values <- switch(
    preset,
    quick = list(
      retention = 0.75,
      n_iterations = 25L,
      min_iterations = 15L,
      instability_bootstrap_B = 200L,
      support_bootstrap_B = 200L,
      batch_size = 5L,
      progress_every = 5L,
      verbose = TRUE
    ),
    standard = list(
      retention = 0.75,
      n_iterations = 100L,
      min_iterations = 50L,
      instability_bootstrap_B = 2000L,
      support_bootstrap_B = 1000L,
      batch_size = 10L,
      progress_every = 10L,
      verbose = TRUE
    ),
    thorough = list(
      retention = 0.75,
      n_iterations = 500L,
      min_iterations = 100L,
      instability_bootstrap_B = 10000L,
      support_bootstrap_B = 5000L,
      batch_size = 25L,
      progress_every = 25L,
      verbose = TRUE
    )
  )
  c(list(name = preset), values)
}

sparq_resolve_preset <- function(preset = "standard", retention = NULL,
                                 n_iterations = NULL, min_iterations = NULL,
                                 instability_bootstrap_B = NULL,
                                 support_bootstrap_B = NULL, batch_size = NULL,
                                 progress_every = NULL, verbose = NULL) {
  values <- sparq_preset(preset)
  resolved <- list(
    preset = values$name,
    retention = if (is.null(retention)) values$retention else retention,
    n_iterations = if (is.null(n_iterations)) values$n_iterations else n_iterations,
    min_iterations = if (is.null(min_iterations)) values$min_iterations else min_iterations,
    instability_bootstrap_B = if (is.null(instability_bootstrap_B)) values$instability_bootstrap_B else instability_bootstrap_B,
    support_bootstrap_B = if (is.null(support_bootstrap_B)) values$support_bootstrap_B else support_bootstrap_B,
    batch_size = if (is.null(batch_size)) values$batch_size else batch_size,
    progress_every = if (is.null(progress_every)) values$progress_every else progress_every,
    verbose = if (is.null(verbose)) values$verbose else verbose
  )
  sparq_check_number(resolved$retention, "retention", .Machine$double.eps, 1)
  sparq_check_number(resolved$n_iterations, "n_iterations", 1, integer = TRUE)
  sparq_check_number(resolved$min_iterations, "min_iterations", 1, integer = TRUE)
  sparq_check_number(resolved$instability_bootstrap_B, "instability_bootstrap_B", 1, integer = TRUE)
  sparq_check_number(resolved$support_bootstrap_B, "support_bootstrap_B", 1, integer = TRUE)
  sparq_check_number(resolved$batch_size, "batch_size", 1, integer = TRUE)
  sparq_check_number(resolved$progress_every, "progress_every", 1, integer = TRUE)
  if (!is.logical(resolved$verbose) || length(resolved$verbose) != 1L || is.na(resolved$verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  resolved
}

sparq_inform_preset <- function(values) {
  sparq_inform(
    paste0(
      "Preset `", values$preset, "`: ", values$n_iterations,
      " perturbations; ", round(values$retention * 100), "% retention; minimum ",
      values$min_iterations, " completed; ", format(values$instability_bootstrap_B, big.mark = ","),
      " instability bootstrap draws; progress every ", values$progress_every, "."
    ),
    values$verbose
  )
}

summary.sparq_assessment <- function(object, ...) {
  out <- object$summary
  if (!is.data.frame(out)) {
    stop("This assessment has no summary table.", call. = FALSE)
  }
  attr(out, "result_id") <- object$result_id
  attr(out, "settings") <- object$settings
  class(out) <- unique(c("summary.sparq_assessment", class(out)))
  out
}

summary.sparq_workflow_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_script_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_spatial_assessment <- function(object, ...) {
  NextMethod("summary")
}

summary.sparq_adaptive_assessment <- function(object, ...) {
  out <- NextMethod("summary")
  if (is.data.frame(object$adaptive_decision) && nrow(object$adaptive_decision) == 1L) {
    for (name in names(object$adaptive_decision)) {
      out[[paste0("adaptive_", name)]] <- object$adaptive_decision[[name]][1]
    }
  }
  out
}

as.data.frame.sparq_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(summary(x), row.names = row.names, optional = optional, ...)
}

as.data.frame.sparq_workflow_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_script_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_spatial_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

as.data.frame.sparq_adaptive_assessment <- function(x, row.names = NULL, optional = FALSE, ...) {
  NextMethod("as.data.frame")
}

sparq_plot_assessment <- function(x, ...) {
  comparisons <- x$comparisons
  summary_table <- x$summary
  if (!is.data.frame(comparisons) || !is.data.frame(summary_table) || !nrow(summary_table)) {
    stop("This assessment does not contain plottable comparison and summary tables.", call. = FALSE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1L))
  valid <- comparisons[comparisons$status == "ok" & is.finite(comparisons$instability), , drop = FALSE]
  if (nrow(valid)) {
    graphics::plot(valid$iteration, valid$instability, pch = 16, col = "#2F6B9A",
      xlab = "Perturbation iteration", ylab = "Instability",
      main = "Instability across perturbations", ...)
    graphics::abline(h = stats::median(valid$instability), lty = 2, col = "#B13A3A")
  } else {
    graphics::plot.new()
    graphics::title(main = "Instability across perturbations")
    graphics::text(0.5, 0.5, "No successful perturbation comparisons")
  }
  components <- c("iteration_support", "bias_support", "uncertainty_support", "reproducibility_support")
  components <- components[components %in% names(summary_table)]
  if (length(components)) {
    graphics::barplot(as.numeric(summary_table[1, components, drop = TRUE]), ylim = c(0, 1),
      col = c("#5B8DB8", "#F2A65A", "#73A857", "#9A6FB0")[seq_along(components)],
      names.arg = gsub("_support", "", components), ylab = "Component value",
      main = "Support-quality components")
  } else {
    graphics::plot.new()
    graphics::title(main = "Support-quality components")
  }
  invisible(x)
}

plot.sparq_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_workflow_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_script_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_spatial_assessment <- function(x, ...) sparq_plot_assessment(x, ...)
plot.sparq_adaptive_assessment <- function(x, ...) sparq_plot_assessment(x, ...)

sparq_format_elapsed <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("unknown")
  if (seconds < 60) return(paste0(round(seconds), "s"))
  if (seconds < 3600) return(paste0(round(seconds / 60, 1), "m"))
  paste0(round(seconds / 3600, 1), "h")
}

sparq_assess_cohort_workflow <- function(sample_data, analysis_function,
                                         output_type = c("scalar", "ranked", "partition", "graph"),
                                         preset = "standard", output_dir = NULL,
                                         progress_every = 1L, seed = 1, verbose = TRUE, ...) {
  if (!is.list(sample_data) || !length(sample_data)) {
    stop("sample_data must be a nonempty named list of sample-level data frames.", call. = FALSE)
  }
  if (is.null(names(sample_data)) || anyNA(names(sample_data)) || any(!nzchar(names(sample_data))) || anyDuplicated(names(sample_data))) {
    stop("sample_data must have unique, nonmissing sample IDs.", call. = FALSE)
  }
  if (!is.function(analysis_function)) stop("analysis_function must be a function.", call. = FALSE)
  output_type <- match.arg(output_type)
  sparq_check_number(progress_every, "progress_every", 1, integer = TRUE)
  sparq_check_number(seed, "seed", 0, integer = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) stop("verbose must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(output_dir) && (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))) {
    stop("output_dir must be NULL or one nonempty path.", call. = FALSE)
  }
  if (!is.null(output_dir) && dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir already exists and is not empty.", call. = FALSE)
  }
  if (!is.null(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  n_samples <- length(sample_data)
  started <- Sys.time()
  manifest <- vector("list", n_samples)
  assessments <- vector("list", n_samples)
  names(assessments) <- names(sample_data)
  sparq_inform(paste0("Cohort assessment started: 0/", n_samples, " sections complete."), verbose)

  for (i in seq_along(sample_data)) {
    sample_id <- names(sample_data)[i]
    sample_started <- Sys.time()
    fit <- tryCatch(
      sparq_assess_workflow(
        data = sample_data[[i]], analysis_function = analysis_function,
        output_type = output_type, preset = preset, result_id = sample_id,
        seed = sparq_seed(seed, i), verbose = FALSE, ...
      ),
      error = function(error) error
    )
    elapsed <- as.numeric(difftime(Sys.time(), sample_started, units = "secs"))
    if (inherits(fit, "error")) {
      manifest[[i]] <- data.frame(
        sample_id = sample_id, status = "failed", error_message = conditionMessage(fit),
        elapsed_seconds = elapsed, seed = NA_integer_, stringsAsFactors = FALSE
      )
    } else {
      assessments[[i]] <- fit
      manifest[[i]] <- data.frame(
        sample_id = sample_id, status = "completed", error_message = "",
        elapsed_seconds = elapsed, seed = fit$settings$seed,
        support_quality = fit$summary$support_quality[1],
        n_completed = fit$summary$n_iterations[1], n_failed = fit$summary$n_failed[1],
        stringsAsFactors = FALSE
      )
    }
    if (verbose && (i %% progress_every == 0L || i == n_samples)) {
      elapsed_total <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      failed <- sum(vapply(manifest[seq_len(i)], function(x) identical(x$status, "failed"), logical(1)))
      remaining <- if (i) elapsed_total / i * (n_samples - i) else NA_real_
      sparq_inform(
        paste0(i, "/", n_samples, " sections complete; ", failed,
          " failed; estimated remaining ", sparq_format_elapsed(remaining), "."),
        TRUE
      )
    }
  }
  manifest <- sparq_bind_rows(manifest)
  summaries <- sparq_bind_rows(lapply(assessments, function(x) if (is.null(x)) NULL else x$summary))
  failures <- manifest[manifest$status == "failed", c("sample_id", "error_message"), drop = FALSE]
  package_version <- tryCatch(as.character(utils::packageVersion("SPARQ")),
    error = function(error) "source checkout")
  manifest$preset <- preset
  manifest$output_type <- output_type
  manifest$cohort_seed <- seed
  manifest$R_version <- R.version.string
  manifest$SPARQ_version <- package_version
  settings <- list(
    preset = preset, output_type = output_type, cohort_seed = seed, started_at = as.character(started),
    completed_at = as.character(Sys.time()), R_version = R.version.string,
    package_version = package_version
  )
  files <- list()
  if (!is.null(output_dir)) {
    files$manifest <- file.path(output_dir, "cohort_manifest.tsv")
    files$summaries <- file.path(output_dir, "cohort_support_summary.tsv")
    files$failures <- file.path(output_dir, "cohort_failures.tsv")
    files$settings <- file.path(output_dir, "cohort_settings.rds")
    utils::write.table(manifest, files$manifest, sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(summaries, files$summaries, sep = "\t", row.names = FALSE, quote = FALSE)
    utils::write.table(failures, files$failures, sep = "\t", row.names = FALSE, quote = FALSE)
    saveRDS(settings, files$settings)
  }
  out <- list(
    summaries = summaries, failures = failures, manifest = manifest,
    settings = settings, assessments = assessments, output_files = files
  )
  class(out) <- "sparq_cohort_assessment"
  out
}

print.sparq_cohort_assessment <- function(x, ...) {
  cat("<SPARQ cohort assessment>\n")
  cat("  Sections completed: ", sum(x$manifest$status == "completed"), "/", nrow(x$manifest), "\n", sep = "")
  cat("  Sections failed: ", sum(x$manifest$status == "failed"), "\n", sep = "")
  if (length(x$output_files)) cat("  Manifest: ", x$output_files$manifest, "\n", sep = "")
  invisible(x)
}

summary.sparq_cohort_assessment <- function(object, ...) object$summaries
as.data.frame.sparq_cohort_assessment <- function(x, row.names = NULL, optional = FALSE, ...) x$summaries
