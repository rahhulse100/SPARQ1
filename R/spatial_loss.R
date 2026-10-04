# Explicit spatial sampling-loss distributions.

#' Define a spatial spot-loss mechanism
#'
#' A loss model specifies how a fixed retained fraction is selected from one
#' spatial section. It makes the perturbation distribution explicit and records
#' the mechanism in every downstream assessment. It does not infer a missingness
#' process from a single section.
#'
#' @param type Loss mechanism: `uniform_random`, `contiguous_hole`,
#'   `quality_weighted`, `edge_weighted`, or `custom`.
#' @param quality_col Numeric quality column required by `quality_weighted`.
#' @param quality_direction Whether larger quality values make a spot more or
#'   less likely to be retained.
#' @param quality_strength Strength of quality weighting; 0 gives uniform
#'   retention and larger positive values increase differential retention.
#' @param edge_strength Strength of radial edge-biased removal; 0 gives uniform
#'   retention.
#' @param custom_function A function returning retained row indices when
#'   `type = "custom"`.
#'
#' @return A `sparq_loss_model` object used by `sparq_apply_loss_model()` or an
#'   assessment's `loss_model` argument.
#' @export
sparq_define_loss_model <- function(
    type = c("uniform_random", "contiguous_hole", "quality_weighted",
      "edge_weighted", "custom"),
    quality_col = NULL,
    quality_direction = c("higher_quality_retained", "lower_quality_retained"),
    quality_strength = 1,
    edge_strength = 1,
    custom_function = NULL) {
  type <- match.arg(type)
  quality_direction <- match.arg(quality_direction)
  sparq_check_number(quality_strength, "quality_strength", 0)
  sparq_check_number(edge_strength, "edge_strength", 0)
  if (type == "quality_weighted" &&
      (!is.character(quality_col) || length(quality_col) != 1L ||
        is.na(quality_col) || !nzchar(quality_col))) {
    stop("quality_col must name one numeric quality column for quality_weighted loss.",
      call. = FALSE
    )
  }
  if (type == "custom" && !is.function(custom_function)) {
    stop("custom_function must be a function for custom loss.", call. = FALSE)
  }
  structure(
    list(
      type = type,
      quality_col = quality_col,
      quality_direction = quality_direction,
      quality_strength = quality_strength,
      edge_strength = edge_strength,
      custom_function = custom_function
    ),
    class = "sparq_loss_model"
  )
}

print.sparq_loss_model <- function(x, ...) {
  cat("<SPARQ spatial loss model>\n")
  cat("  Type: ", x$type, "\n", sep = "")
  if (identical(x$type, "quality_weighted")) {
    cat("  Quality column: ", x$quality_col, "; direction: ",
      x$quality_direction, "; strength: ", x$quality_strength, "\n", sep = "")
  }
  if (identical(x$type, "edge_weighted")) {
    cat("  Edge-removal strength: ", x$edge_strength, "\n", sep = "")
  }
  invisible(x)
}

sparq_loss_probabilities <- function(data, loss_model,
                                     x_col = NULL, y_col = NULL) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a nonempty data.frame.", call. = FALSE)
  }
  if (!inherits(loss_model, "sparq_loss_model")) {
    stop("loss_model must be created by sparq_define_loss_model().", call. = FALSE)
  }
  n <- nrow(data)
  type <- loss_model$type
  weights <- rep(1, n)

  if (identical(type, "quality_weighted")) {
    if (!loss_model$quality_col %in% names(data)) {
      stop("quality_col is not present in data.", call. = FALSE)
    }
    quality <- data[[loss_model$quality_col]]
    if (!is.numeric(quality) || any(!is.finite(quality))) {
      stop("quality_col must contain finite numeric values for every spot.",
        call. = FALSE
      )
    }
    quality_rank <- rank(quality, ties.method = "average") / (n + 1)
    if (identical(loss_model$quality_direction, "lower_quality_retained")) {
      quality_rank <- 1 - quality_rank
    }
    weights <- pmax(quality_rank, .Machine$double.eps)^loss_model$quality_strength
  } else if (identical(type, "edge_weighted")) {
    if (length(x_col) != 1L || length(y_col) != 1L ||
        !all(c(x_col, y_col) %in% names(data))) {
      stop("edge_weighted loss requires valid x_col and y_col.", call. = FALSE)
    }
    xy <- data[, c(x_col, y_col), drop = FALSE]
    if (!all(vapply(xy, is.numeric, logical(1))) ||
        any(!is.finite(as.matrix(xy)))) {
      stop("x_col and y_col must contain finite numeric coordinates.",
        call. = FALSE
      )
    }
    center_x <- stats::median(xy[[1]])
    center_y <- stats::median(xy[[2]])
    radius <- sqrt((xy[[1]] - center_x)^2 + (xy[[2]] - center_y)^2)
    scaled_radius <- if (max(radius) > 0) radius / max(radius) else rep(0, n)
    weights <- 1 / (1 + loss_model$edge_strength * scaled_radius)
  }

  data.frame(
    spot_index = seq_len(n),
    retention_weight = as.numeric(weights),
    relative_loss_weight = as.numeric(1 / weights),
    stringsAsFactors = FALSE
  )
}

#' Apply a defined spatial loss model
#'
#' @param data Spot-level data frame.
#' @param loss_model A `sparq_loss_model` object.
#' @param retention Fraction of spots retained.
#' @param x_col,y_col Coordinate columns when required by the loss model.
#' @param iteration,seed Reproducibility controls.
#'
#' @return A list with perturbed `data`, retained and removed row indices, and a
#'   per-spot loss audit. The retained count is exactly determined by retention.
#' @export
sparq_apply_loss_model <- function(data, loss_model, retention = 0.75,
                                   x_col = NULL, y_col = NULL,
                                   iteration = 1, seed = NULL) {
  if (!is.data.frame(data) || !nrow(data)) {
    stop("data must be a nonempty data.frame.", call. = FALSE)
  }
  if (!inherits(loss_model, "sparq_loss_model")) {
    stop("loss_model must be created by sparq_define_loss_model().", call. = FALSE)
  }
  sparq_check_number(retention, "retention", .Machine$double.eps, 1)
  sparq_check_number(iteration, "iteration", 1, integer = TRUE)
  n <- nrow(data)
  n_keep <- min(n, max(2L, floor(n * retention)))

  select_indices <- function() {
    if (n_keep == n || identical(loss_model$type, "none")) return(seq_len(n))
    if (identical(loss_model$type, "contiguous_hole")) {
      if (length(x_col) != 1L || length(y_col) != 1L ||
          !all(c(x_col, y_col) %in% names(data))) {
        stop("contiguous_hole loss requires valid x_col and y_col.", call. = FALSE)
      }
      xy <- data[, c(x_col, y_col), drop = FALSE]
      if (!all(vapply(xy, is.numeric, logical(1))) ||
          any(!is.finite(as.matrix(xy)))) {
        stop("x_col and y_col must contain finite numeric coordinates.",
          call. = FALSE
        )
      }
      center <- sample.int(n, 1L)
      distance_squared <- (xy[[1]] - xy[[1]][center])^2 +
        (xy[[2]] - xy[[2]][center])^2
      return(order(-distance_squared, seq_len(n))[seq_len(n_keep)])
    }
    if (identical(loss_model$type, "custom")) {
      selected <- loss_model$custom_function(
        data = data,
        retention = retention,
        iteration = iteration
      )
      if (!is.numeric(selected) || length(selected) != n_keep ||
          anyNA(selected) || anyDuplicated(selected) ||
          any(selected != floor(selected)) || any(selected < 1L | selected > n)) {
        stop("custom loss must return exactly n_keep unique integer row indices.",
          call. = FALSE
        )
      }
      return(as.integer(selected))
    }
    probabilities <- sparq_loss_probabilities(data, loss_model, x_col, y_col)
    sample.int(n, n_keep, replace = FALSE,
      prob = probabilities$retention_weight
    )
  }

  retained_index <- if (is.null(seed)) select_indices() else {
    sparq_with_seed(sparq_seed(seed, iteration, 61L), select_indices())
  }
  audit <- sparq_loss_probabilities(data, loss_model, x_col, y_col)
  audit$retained <- audit$spot_index %in% retained_index
  audit$loss_model <- loss_model$type
  audit$retention <- retention
  out <- list(
    data = data[retained_index, , drop = FALSE],
    retained_index = retained_index,
    removed_index = setdiff(seq_len(n), retained_index),
    audit = audit,
    loss_model = loss_model
  )
  class(out) <- "sparq_loss_application"
  out
}
