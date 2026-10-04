test_that("every NAMESPACE export resolves after package load", {
  exported <- getNamespaceExports("SPARQ")
  available <- vapply(exported, exists, logical(1),
    envir = asNamespace("SPARQ"), inherits = FALSE)

  expect_true(all(available))
})
