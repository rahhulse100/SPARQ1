# Spatial-object adapters. Optional upstream packages are intentionally not
# imported: users only need the package that owns the object they provide.

sparq_detect_spatial_object <- function(object) {
  if (is.data.frame(object)) return("data.frame")
  if (inherits(object, "Seurat")) return("seurat")
  if (inherits(object, "SpatialExperiment")) return("spatialexperiment")
  if (inherits(object, "giotto")) return("giotto")
  "unknown"
}

sparq_optional_export <- function(package, function_name) {
  if (!package %in% loadedNamespaces()) {
    loaded <- tryCatch({
      loadNamespace(package)
      TRUE
    }, error = function(error) FALSE)
    if (!loaded) {
      stop(paste0(package, " must be installed to use this object type."), call. = FALSE)
    }
  }
  getExportedValue(package, function_name)
}

sparq_spatial_coordinate_columns <- function(data, coordinate_columns = NULL) {
  if (!is.null(coordinate_columns)) {
    if (!is.character(coordinate_columns) || length(coordinate_columns) != 2L ||
        anyNA(coordinate_columns) || !all(coordinate_columns %in% names(data))) {
      stop("coordinate_columns must name exactly two columns in the spatial coordinates.",
        call. = FALSE)
    }
    return(coordinate_columns)
  }

  x_candidates <- c("x", "X", "imagecol", "pxl_col_in_fullres", "col", "array_col")
  y_candidates <- c("y", "Y", "imagerow", "pxl_row_in_fullres", "row", "array_row")
  x <- x_candidates[x_candidates %in% names(data)][1]
  y <- y_candidates[y_candidates %in% names(data)][1]
  if (is.na(x) || is.na(y)) {
    stop(
      "SPARQ could not identify x/y coordinates. Supply coordinate_columns = c('x_column', 'y_column').",
      call. = FALSE
    )
  }
  c(x, y)
}

sparq_make_spot_table <- function(spot_ids, coordinates, coordinate_columns = NULL) {
  spot_ids <- as.character(spot_ids)
  if (!length(spot_ids) || anyNA(spot_ids) || any(!nzchar(spot_ids)) || anyDuplicated(spot_ids)) {
    stop("Spatial spot IDs must be nonmissing and unique.", call. = FALSE)
  }
  coordinates <- as.data.frame(coordinates, stringsAsFactors = FALSE)
  coordinate_columns <- sparq_spatial_coordinate_columns(coordinates, coordinate_columns)
  if (is.null(rownames(coordinates)) || anyNA(rownames(coordinates)) ||
      anyDuplicated(rownames(coordinates))) {
    stop("Spatial coordinate rows must be named by unique spot IDs.", call. = FALSE)
  }
  coordinates <- coordinates[match(spot_ids, rownames(coordinates)), , drop = FALSE]
  if (anyNA(coordinates[[coordinate_columns[1]]]) || anyNA(coordinates[[coordinate_columns[2]]]) ||
      !is.numeric(coordinates[[coordinate_columns[1]]]) || !is.numeric(coordinates[[coordinate_columns[2]]]) ||
      any(!is.finite(coordinates[[coordinate_columns[1]]])) ||
      any(!is.finite(coordinates[[coordinate_columns[2]]]))) {
    stop("Spatial coordinates must be finite numeric values for every retained spot.",
      call. = FALSE)
  }
  data.frame(
    .sparq_spot_id = spot_ids,
    .sparq_x = as.numeric(coordinates[[coordinate_columns[1]]]),
    .sparq_y = as.numeric(coordinates[[coordinate_columns[2]]]),
    stringsAsFactors = FALSE
  )
}

sparq_prepare_spatial <- function(object,
    object_type = c("auto", "data.frame", "seurat", "spatialexperiment", "giotto"),
    coordinate_columns = NULL, spot_id_column = NULL, image = NULL,
    adapter = NULL) {

  object_type <- match.arg(object_type)
  if (!is.null(adapter) && (!is.list(adapter) || !is.function(adapter$extract) ||
      !is.function(adapter$subset))) {
    stop("adapter must be a list with extract(object) and subset(object, spot_ids) functions.",
      call. = FALSE)
  }
  detected_type <- sparq_detect_spatial_object(object)
  if (object_type == "auto") object_type <- detected_type
  if (!is.null(adapter)) object_type <- "custom"
  if (object_type == "unknown") {
    stop("Unsupported spatial object. Supply a custom adapter or use a data.frame, Seurat, SpatialExperiment, or Giotto object.",
      call. = FALSE)
  }

  if (object_type == "custom") {
    extracted <- adapter$extract(object)
    if (!is.list(extracted) || is.null(extracted$spot_ids) || is.null(extracted$coordinates)) {
      stop("adapter$extract(object) must return a list with spot_ids and coordinates.",
        call. = FALSE)
    }
    spot_table <- sparq_make_spot_table(extracted$spot_ids, extracted$coordinates,
      coordinate_columns)
    subsetter <- adapter$subset
  } else if (object_type == "data.frame") {
    ids <- if (!is.null(spot_id_column)) {
      if (!spot_id_column %in% names(object)) stop("spot_id_column is not in object.", call. = FALSE)
      object[[spot_id_column]]
    } else if (!is.null(rownames(object)) && !anyDuplicated(rownames(object)) &&
      all(nzchar(rownames(object)))) {
      rownames(object)
    } else {
      paste0("spot_", seq_len(nrow(object)))
    }
    coordinates <- object
    rownames(coordinates) <- as.character(ids)
    spot_table <- sparq_make_spot_table(ids, coordinates, coordinate_columns)
    subsetter <- function(x, spot_ids) x[match(spot_ids, as.character(ids)), , drop = FALSE]
  } else if (object_type == "seurat") {
    seurat_images <- sparq_optional_export("SeuratObject", "Images")
    seurat_coordinates <- sparq_optional_export("SeuratObject", "GetTissueCoordinates")
    seurat_cells <- sparq_optional_export("SeuratObject", "Cells")
    image_names <- seurat_images(object)
    if (is.null(image)) image <- image_names[1]
    if (length(image) != 1L || is.na(image) || !image %in% image_names) {
      stop("image must name one spatial image stored in the Seurat object.", call. = FALSE)
    }
    coordinates <- seurat_coordinates(object[[image]])
    ids <- intersect(seurat_cells(object), rownames(coordinates))
    if (!length(ids)) stop("No Seurat cell IDs match the selected image coordinates.", call. = FALSE)
    spot_table <- sparq_make_spot_table(ids, coordinates, coordinate_columns)
    subsetter <- function(x, spot_ids) x[, spot_ids, drop = FALSE]
  } else if (object_type == "spatialexperiment") {
    spatial_coordinates <- sparq_optional_export("SpatialExperiment", "spatialCoords")
    coordinates <- as.data.frame(spatial_coordinates(object))
    ids <- colnames(object)
    rownames(coordinates) <- ids
    spot_table <- sparq_make_spot_table(ids, coordinates, coordinate_columns)
    subsetter <- function(x, spot_ids) x[, spot_ids, drop = FALSE]
  } else if (object_type == "giotto") {
    candidates <- c("GiottoClass", "Giotto")
    available <- vapply(candidates, function(package) {
      package %in% loadedNamespaces() || isTRUE(tryCatch({
        loadNamespace(package)
        TRUE
      }, error = function(error) FALSE))
    }, logical(1))
    if (!any(available)) {
      stop("GiottoClass or Giotto must be installed to use a giotto object.", call. = FALSE)
    }
    giotto_package <- candidates[which(available)[1]]
    spatial_locations <- getExportedValue(giotto_package, "spatLocs")(object)
    coordinates <- as.data.frame(spatial_locations, stringsAsFactors = FALSE)
    id_column <- c("cell_ID", "cell_id", "spat_ID", "spat_id")[c("cell_ID", "cell_id", "spat_ID", "spat_id") %in% names(coordinates)][1]
    if (is.na(id_column)) stop("Giotto spatial locations must include a cell_ID column.", call. = FALSE)
    ids <- coordinates[[id_column]]
    rownames(coordinates) <- as.character(ids)
    spot_table <- sparq_make_spot_table(ids, coordinates, coordinate_columns)
    subsetter <- function(x, spot_ids) {
      getExportedValue(giotto_package, "subsetGiotto")(x, cell_ids = spot_ids, verbose = FALSE)
    }
  } else {
    stop("object_type does not match a supported object.", call. = FALSE)
  }

  out <- list(
    object = object,
    object_type = object_type,
    spot_data = spot_table,
    n_spots = nrow(spot_table),
    coordinate_columns = c(x = ".sparq_x", y = ".sparq_y"),
    subset_function = subsetter
  )
  class(out) <- "sparq_spatial_input"
  out
}

sparq_subset_spatial <- function(spatial_input, spot_ids) {
  if (!inherits(spatial_input, "sparq_spatial_input")) {
    stop("spatial_input must be created by sparq_prepare_spatial().", call. = FALSE)
  }
  spot_ids <- as.character(spot_ids)
  available <- spatial_input$spot_data$.sparq_spot_id
  if (!length(spot_ids) || anyNA(spot_ids) || anyDuplicated(spot_ids) ||
      !all(spot_ids %in% available)) {
    stop("spot_ids must be unique IDs retained in spatial_input.", call. = FALSE)
  }
  spatial_input$subset_function(spatial_input$object, spot_ids)
}

sparq_assess_spatial <- function(object, analysis_function,
    output_type = c("scalar", "ranked", "partition", "graph"),
    stress_model = c("uniform_random", "contiguous_hole", "none", "custom"),
    retention = NULL, n_iterations = NULL, top_k = 100,
    result_id = "result", reference_scale = 1, seed = 1,
    object_type = c("auto", "data.frame", "seurat", "spatialexperiment", "giotto"),
    coordinate_columns = NULL, spot_id_column = NULL, image = NULL,
    adapter = NULL, verbose = NULL, preset = "standard", progress_every = NULL,
    min_iterations = NULL, instability_bootstrap_B = NULL, ...) {

  if (!is.function(analysis_function)) {
    stop("analysis_function must accept one perturbed spatial object.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }
  prepared <- sparq_prepare_spatial(
    object = object,
    object_type = object_type,
    coordinate_columns = coordinate_columns,
    spot_id_column = spot_id_column,
    image = image,
    adapter = adapter
  )
  sparq_inform(
    paste0("Prepared ", prepared$object_type, " input with ", prepared$n_spots, " spatial observations."),
    verbose
  )
  rerun <- function(perturbed_spots) {
    perturbed_object <- sparq_subset_spatial(prepared, perturbed_spots$.sparq_spot_id)
    analysis_function(perturbed_object)
  }
  fit <- sparq_assess_workflow(
    data = prepared$spot_data,
    analysis_function = rerun,
    output_type = output_type,
    stress_model = stress_model,
    retention = retention,
    n_iterations = n_iterations,
    top_k = top_k,
    x_col = ".sparq_x",
    y_col = ".sparq_y",
    result_id = result_id,
    reference_scale = reference_scale,
    seed = seed,
    verbose = verbose,
    preset = preset,
    progress_every = progress_every,
    min_iterations = min_iterations,
    instability_bootstrap_B = instability_bootstrap_B,
    ...
  )
  fit$spatial_input <- prepared
  fit$spatial_object_type <- prepared$object_type
  class(fit) <- unique(c("sparq_spatial_assessment", class(fit)))
  fit
}
