test_that("Weis 2014 returns S, rho, and matrices", {
  data(example_weiss)
  res <- weis_sync(example_weiss, n_boot = 0, verbose = FALSE)

  expect_type(res, "list")
  expect_named(res, c(
    "input_info", "trait_vector", "flowering_matrix", "M_matrix", "P_matrix",
    "PHI_matrix", "eigenvalues", "lambda1", "eigen_sum", "S", "rho", "bootstrap"
  ))
  expect_true(is.matrix(res$PHI_matrix))
  expect_true(is.numeric(res$S))
  expect_true(is.numeric(res$rho))
})
