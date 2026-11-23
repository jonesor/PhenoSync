test_that("Dainou 2012 returns male synchrony and population S", {
  data(example_dainou)
  res <- dainou_SI(example_dainou, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("male_dj", "female_prop_by_day", "population_S"))
  expect_true(is.numeric(res$population_S))
  expect_true(res$population_S >= 0 && res$population_S <= 1)
})

test_that("Dainou 2012 accepts id_col and sex_col as column names", {
  data(example_dainou)
  by_number <- dainou_SI(example_dainou, id_col = 1, sex_col = 2, verbose = FALSE)
  by_name <- dainou_SI(example_dainou,
                       id_col = names(example_dainou)[1],
                       sex_col = names(example_dainou)[2],
                       verbose = FALSE)
  expect_equal(by_name, by_number)
})

test_that("Dainou 2012 gives a clear error for a column that does not exist", {
  data(example_dainou)
  expect_error(dainou_SI(example_dainou, id_col = "no_such_column", verbose = FALSE),
               "not a column")
})
