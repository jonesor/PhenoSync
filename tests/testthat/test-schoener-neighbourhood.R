test_that("Schoener neighbourhood synchrony works with x/y coordinates", {
  data(example_neighbourhood)
  res <- schoener_SI(
    metadata = example_neighbourhood$metadata,
    flower_mat = example_neighbourhood$flower_mat,
    id_col = "id",
    x_col = "x",
    y_col = "y",
    neighbor_distance = 50,
    verbose = FALSE
  )

  expect_type(res, "list")
  expect_named(res, c("summary", "SI_individual", "mean_SI", "distance_matrix", "coordinates_xy", "details"))
  expect_equal(nrow(res$summary), ncol(example_neighbourhood$flower_mat))
  expect_true(is.numeric(res$mean_SI) || is.na(res$mean_SI))
})
