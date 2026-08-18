test_that("Primack 1980 returns pairwise table and Z", {
  data(example_wide)
  res <- primack_Z(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("pairwise_table", "Z", "strong_matrix"))
  expect_true(is.data.frame(res$pairwise_table))
  expect_true(is.numeric(res$Z) || is.na(res$Z))
})
