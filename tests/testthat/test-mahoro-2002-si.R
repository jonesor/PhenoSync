test_that("Mahoro 2002 Si returns one value per individual", {
  data(example_wide)
  res <- mahoro_Si(example_wide, verbose = FALSE)

  expect_type(res, "double")
  expect_length(res, ncol(example_wide))
})
