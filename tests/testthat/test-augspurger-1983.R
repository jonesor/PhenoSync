test_that("Augspurger 1983 returns Xi table and Z", {
  data(example_wide)
  res <- augspurger_Xi(example_wide, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c("Xi_table", "Z", "presence_matrix", "overlap_matrix"))
  expect_equal(nrow(res$Xi_table), ncol(example_wide))
  expect_true(is.numeric(res$Z))
  expect_true(res$Z >= 0 && res$Z <= 1)
})
