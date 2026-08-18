test_that("Marquis 1988 returns individual and population SM", {
  data(example_wide)
  res <- marquis_SM(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("SM_individual", "SM_population_mean", "p_t", "flowering_presence", "total_output"))
  expect_length(res$SM_individual, ncol(example_wide))
  expect_true(is.numeric(res$SM_population_mean) || is.na(res$SM_population_mean))
})
