test_that("Bolmgren 1998 returns individual results and mean S", {
  data(example_bolmgren_long)
  res <- bolmgren_S(example_bolmgren_long, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("individual_results", "mean_S"))
  expect_true(is.data.frame(res$individual_results))
  expect_true(is.numeric(res$mean_S) || is.na(res$mean_S))
})
