test_that("make_pheno_dummy creates expected wide data", {
  x <- make_pheno_dummy(format = "wide", n_species = 5, n_days = 30, binary = FALSE, unit = "species", seed = 1)

  expect_true(is.matrix(x))
  expect_equal(dim(x), c(30, 5))
  expect_true(is.numeric(x))
  expect_true(all(x >= 0))
})

test_that("make_pheno_dummy creates expected specialised data", {
  bolmgren <- make_pheno_dummy(format = "bolmgren_long", n_individuals = 5, n_days = 30, seed = 1)
  dainou <- make_pheno_dummy(format = "dainou", n_individuals = 6, n_days = 30, seed = 1)
  weis <- make_pheno_dummy(format = "weis", n_individuals = 6, n_days = 30, seed = 1)
  nb <- make_pheno_dummy(format = "neighbourhood", n_individuals = 6, n_days = 30, seed = 1)

  expect_true(all(c("plant_id", "onset_day", "end_day") %in% names(bolmgren)))
  expect_true(all(c("id", "sex") %in% names(dainou)))
  expect_true("f.date" %in% names(weis))
  expect_true(is.list(nb))
  expect_true(all(c("metadata", "flower_mat") %in% names(nb)))
})
