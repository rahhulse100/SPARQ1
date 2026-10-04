test_that("legacy methods validation remains executable as a package test", {
  source(testthat::test_path("..", "methods-validation.R"), local = environment())
  expect_silent(sparq_test_methods(full_size = FALSE))
})
