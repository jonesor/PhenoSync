test_that("Koenig 2003 returns pairwise correlations and population rp", {
  data(example_wide)
  res <- koenig_rp(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("pairwise_correlation", "pairwise_values", "individual_mean_r", "population_rp"))
  expect_true(is.matrix(res$pairwise_correlation))
  expect_true(is.numeric(res$population_rp) || is.na(res$population_rp))
})
