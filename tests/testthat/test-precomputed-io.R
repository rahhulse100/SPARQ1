test_that("precomputed results can be assessed from lists and files", {
  full <- list(sample_a = 1, sample_b = 2)
  perturbed <- list(
    sample_a = c(0.9, 1.1, 1.0),
    sample_b = c(1.9, 2.1, 2.0)
  )
  assessment <- sparq_assess_precomputed(
    full_results = full,
    perturbed_results = perturbed,
    comparator = sparq_compare_scalar,
    min_iterations = 3,
    instability_bootstrap_B = 20
  )
  expect_equal(nrow(assessment), 2)

  full_file <- tempfile(fileext = ".tsv")
  perturbed_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(full_file, perturbed_file)), add = TRUE)
  write.table(
    data.frame(result_id = c("sample_a", "sample_b"), value = c(1, 2)),
    full_file, sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    data.frame(
      result_id = rep(c("sample_a", "sample_b"), each = 3),
      value = c(0.9, 1.1, 1.0, 1.9, 2.1, 2.0)
    ),
    perturbed_file, sep = "\t", row.names = FALSE, quote = FALSE
  )
  from_files <- sparq_assess_precomputed_files(
    full_file, perturbed_file,
    min_iterations = 3,
    instability_bootstrap_B = 20
  )
  expect_equal(nrow(from_files), 2)
})
