test_that("Albert 2001 returns individual and population synchrony", {
  data(example_wide)
  res <- albert_Si(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("individual_Si", "population_mean"))
  expect_length(res$individual_Si, ncol(example_wide))
  expect_true(is.numeric(res$population_mean))
  expect_true(res$population_mean >= 0 && res$population_mean <= 1)
})
