test_that("wide index functions preserve explicit census column", {
  x <- make_pheno_dummy(format = "wide", n_species = 4, n_days = 10, seed = 1)
  df <- cbind(census_day = rownames(x), as.data.frame(x, check.names = FALSE))

  res <- albert_Si(df, census_col = "census_day", verbose = FALSE)

  expect_type(res, "list")
  expect_true(is.numeric(res$population_mean))
})

test_that("calculate_indices preserves explicit census column as row names", {
  x <- make_pheno_dummy(format = "wide", n_species = 4, n_days = 10, seed = 1)
  df <- cbind(census_day = rownames(x), as.data.frame(x, check.names = FALSE))

  res <- calculate_indices(
    df,
    columns_are = "species",
    binary = FALSE,
    census_col = "census_day",
    indices = c("albert_Si", "pianka_overlap"),
    verbose = FALSE,
    make_plot = FALSE
  )

  expect_equal(rownames(res$prepared_data), rownames(x))
  expect_equal(colnames(res$prepared_data), colnames(x))
})
