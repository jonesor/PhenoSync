test_that("Dainou 2012 returns male synchrony and population S", {
  data(example_dainou)
  res <- dainou_SI(example_dainou, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("male_dj", "female_prop_by_day", "population_S"))
  expect_true(is.numeric(res$population_S))
  expect_true(res$population_S >= 0 && res$population_S <= 1)
})
