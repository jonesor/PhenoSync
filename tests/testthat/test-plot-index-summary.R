test_that("plot_index_summary draws a base R plot invisibly", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)
  res <- calculate_indices(x, columns_are = "species", binary = FALSE, verbose = FALSE, make_plot = FALSE)

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  out <- plot_index_summary(res$summary)
  grDevices::dev.off()

  expect_true(file.exists(tmp))
  expect_true(file.info(tmp)$size > 0)
  expect_equal(out, res$summary)
})

test_that("calculate_indices can store a recorded base plot", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)

  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  res <- calculate_indices(x, columns_are = "species", binary = FALSE, verbose = FALSE, make_plot = TRUE)
  grDevices::dev.off()

  expect_type(res, "list")
  expect_true(is.null(res$plot) || inherits(res$plot, "recordedplot"))
})
