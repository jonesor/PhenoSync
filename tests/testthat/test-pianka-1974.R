test_that("Pianka 1974 returns pairwise overlap matrix", {
  data(example_wide)
  res <- pianka_overlap(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("proportional_matrix", "pairwise_overlap", "individual_mean_overlap", "population_mean"))
  expect_true(is.matrix(res$pairwise_overlap))
  expect_true(is.numeric(res$population_mean) || is.na(res$population_mean))
})
