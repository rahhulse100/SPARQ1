test_that("data-frame spatial inputs preserve the source object through assessment", {
  source_object <- data.frame(
    value = c(1, 2, 3, 4, 5, 6),
    x = 1:6,
    y = 6:1,
    row.names = paste0("spot_", 1:6)
  )

  prepared <- sparq_prepare_spatial(
    source_object,
    coordinate_columns = c("x", "y")
  )

  expect_s3_class(prepared, "sparq_spatial_input")
  expect_equal(prepared$n_spots, 6)
  expect_equal(
    rownames(sparq_subset_spatial(prepared, c("spot_1", "spot_3"))),
    c("spot_1", "spot_3")
  )

  fit <- sparq_assess_spatial(
    object = source_object,
    analysis_function = function(x) mean(x$value),
    output_type = "scalar",
    retention = 0.75,
    n_iterations = 4,
    min_iterations = 3,
    instability_bootstrap_B = 20,
    seed = 5,
    verbose = FALSE
  )

  expect_s3_class(fit, "sparq_spatial_assessment")
  expect_equal(fit$spatial_object_type, "data.frame")
  expect_equal(nrow(fit$comparisons), 4)
})

test_that("custom adapters declare extraction and subsetting contracts", {
  object <- list(
    values = setNames(c(1, 2, 3, 4), paste0("s", 1:4)),
    coordinates = data.frame(x = 1:4, y = 4:1, row.names = paste0("s", 1:4))
  )
  adapter <- list(
    extract = function(x) list(spot_ids = names(x$values), coordinates = x$coordinates),
    subset = function(x, spot_ids) x$values[spot_ids]
  )

  prepared <- sparq_prepare_spatial(object, adapter = adapter)
  expect_equal(sparq_subset_spatial(prepared, c("s2", "s4")), c(s2 = 2, s4 = 4))
})
