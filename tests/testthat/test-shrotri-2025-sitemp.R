test_that("Shrotri SItemp returns pairwise and community synchrony", {
  data(example_wide)
  res <- shrotri_SI_temp(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("pairwise_SItemp", "individual_SItemp", "community_mean"))
  expect_true(is.matrix(res$pairwise_SItemp))
  expect_true(is.numeric(res$community_mean) || is.na(res$community_mean))
})
