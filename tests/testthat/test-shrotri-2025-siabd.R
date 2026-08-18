test_that("Shrotri SIabd returns pairwise and community synchrony", {
  data(example_wide)
  res <- shrotri_SI_abd(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("pairwise_SIabd", "individual_SIabd", "community_mean"))
  expect_true(is.matrix(res$pairwise_SIabd))
  expect_true(is.numeric(res$community_mean) || is.na(res$community_mean))
})
