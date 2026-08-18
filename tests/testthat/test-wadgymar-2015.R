test_that("Wadgymar 2015 returns individual results and mean Si", {
  data(example_wide)
  res <- wadgymar_Si(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("individual_results", "Si_individual", "mean_Si"))
  expect_length(res$Si_individual, ncol(example_wide))
  expect_true(is.numeric(res$mean_Si) || is.na(res$mean_Si))
})
