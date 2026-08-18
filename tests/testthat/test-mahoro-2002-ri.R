test_that("Mahoro 2002 Ri returns rank-order output", {
  data(example_wide)
  res <- mahoro_Ri(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("Ri", "peak_day", "daily_ranks", "cumulative_proportion"))
  expect_length(res$Ri, ncol(example_wide))

  # Mahoro Ri ranks cumulative flowering only up to the population peak day,
  # so daily_ranks has peak_day rows, not the full number of census days.
  expect_equal(nrow(res$daily_ranks), unname(res$peak_day))
  expect_equal(ncol(res$daily_ranks), ncol(example_wide))

  # The cumulative proportion matrix still tracks the full phenological season.
  expect_equal(dim(res$cumulative_proportion), dim(example_wide))
})
