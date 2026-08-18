test_that("calculate_indices calculates selected wide indices", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)
  res <- calculate_indices(
    x,
    columns_are = "species",
    binary = FALSE,
    indices = c("albert_Si", "pianka_overlap"),
    verbose = FALSE,
    make_plot = FALSE
  )

  expect_type(res, "list")
  expect_equal(res$detected_type, "wide")
  expect_equal(res$columns_are, "species")
  expect_false(res$binary)
  expect_true(all(c("albert_Si", "pianka_overlap") %in% res$calculated_indices))
  expect_equal(nrow(res$summary), 2)
  expect_true(all(c("range_class", "range_panel") %in% names(res$summary)))
})

test_that("calculate_indices filters auto indices by species vs individuals", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, seed = 1)

  species_res <- calculate_indices(x, columns_are = "species", binary = FALSE, verbose = FALSE, make_plot = FALSE)
  individual_res <- calculate_indices(x, columns_are = "individuals", binary = FALSE, verbose = FALSE, make_plot = FALSE)

  expect_true("poole_rathcke" %in% species_res$possible_indices)
  expect_false("poole_rathcke" %in% individual_res$possible_indices)
  expect_true("mahoro_Si" %in% individual_res$possible_indices)
  expect_false("mahoro_Si" %in% species_res$possible_indices)
})

test_that("calculate_indices filters auto indices for binary data", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, binary = TRUE, seed = 1)
  res <- calculate_indices(x, columns_are = "species", binary = TRUE, verbose = FALSE, make_plot = FALSE)

  expect_true(res$binary)
  expect_true(all(res$possible_indices %in% c("albert_Si", "pianka_overlap", "shrotri_SI_temp")))
  expect_false("shrotri_SI_abd" %in% res$possible_indices)
})

test_that("calculate_indices cleans raw wide tables with first row and first column labels", {
  x <- make_pheno_dummy(format = "wide", n_species = 4, n_days = 10, seed = 1)
  raw <- cbind(day = rownames(x), as.data.frame(x, check.names = FALSE))
  raw <- rbind(c("day", colnames(x)), raw)
  rownames(raw) <- NULL
  names(raw) <- paste0("V", seq_len(ncol(raw)))

  res <- calculate_indices(raw, columns_are = "species", binary = FALSE, census_col = "day", verbose = FALSE, make_plot = FALSE)

  expect_equal(dim(res$prepared_data), dim(x))
  expect_equal(colnames(res$prepared_data), colnames(x))
  expect_equal(rownames(res$prepared_data), rownames(x))
})

test_that("calculate_indices no longer exposes data_type or specialised input arguments", {
  fmls <- names(formals(calculate_indices))
  expect_false("data_type" %in% fmls)
  expect_false("trait_col" %in% fmls)
  expect_false("metadata" %in% fmls)
})
