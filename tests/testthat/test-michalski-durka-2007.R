test_that("Michalski and Durka 2007 returns CV summaries", {
  data(example_wide)
  res <- michalski_durka_CV(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("population_profile", "population_CV", "individual_CV", "mean_individual_CV"))
  expect_length(res$individual_CV, ncol(example_wide))
  expect_true(is.numeric(res$population_CV) || is.na(res$population_CV))
})
