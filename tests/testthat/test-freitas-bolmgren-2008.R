test_that("Freitas and Bolmgren 2008 returns population S", {
  data(example_wide)
  res <- freitas_bolmgren_S(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("relative_intensity", "Ti", "di", "population_S"))
  expect_length(res$di, ncol(example_wide))
  expect_true(is.numeric(res$population_S) || is.na(res$population_S))
})
