data(example_wide)

wide_indices <- list(
  albert_Si = albert_Si,
  augspurger_Xi = augspurger_Xi,
  freitas_bolmgren_S = freitas_bolmgren_S,
  koenig_rp = koenig_rp,
  mahoro_Ri = mahoro_Ri,
  mahoro_Si = mahoro_Si,
  marquis_SM = marquis_SM,
  michalski_durka_CV = michalski_durka_CV,
  pianka_overlap = pianka_overlap,
  poole_rathcke = poole_rathcke,
  primack_Z = primack_Z,
  shrotri_SI_abd = shrotri_SI_abd,
  shrotri_SI_temp = shrotri_SI_temp,
  wadgymar_Si = wadgymar_Si
)

run_wide <- function(fn, d) fn(d, verbose = FALSE)

test_that("wide-format indices stop on negative values", {
  d <- example_wide
  d[20, 1] <- -3
  for (nm in names(wide_indices)) {
    expect_error(run_wide(wide_indices[[nm]], d), "negative or infinite", info = nm)
  }
})

test_that("wide-format indices stop on infinite values", {
  d <- example_wide
  d[20, 1] <- Inf
  for (nm in names(wide_indices)) {
    expect_error(run_wide(wide_indices[[nm]], d), "negative or infinite", info = nm)
  }
})

test_that("wide-format indices stop on repeated column names", {
  d <- example_wide
  colnames(d)[2] <- colnames(d)[1]
  for (nm in names(wide_indices)) {
    expect_error(run_wide(wide_indices[[nm]], d), "same name", info = nm)
  }
})

test_that("wide-format indices stop on empty data or a single census", {
  for (nm in names(wide_indices)) {
    expect_error(run_wide(wide_indices[[nm]], example_wide[0, , drop = FALSE]),
                 "empty", info = nm)
    expect_error(run_wide(wide_indices[[nm]], example_wide[20, , drop = FALSE]),
                 "At least two censuses", info = nm)
  }
})

test_that("calculate_indices() stops on bad inputs", {
  d <- example_wide
  d[20, 1] <- -3
  expect_error(calculate_indices(d, make_plot = FALSE, verbose = FALSE),
               "negative or infinite")
})

test_that("augspurger_Xi() needs at least two individuals", {
  expect_error(augspurger_Xi(example_wide[, 1, drop = FALSE], verbose = FALSE),
               "At least two individuals")
})

test_that("dainou_SI() stops on bad flowering values and repeated IDs", {
  data(example_dainou)
  d <- example_dainou
  d[3, 4] <- -2
  expect_error(dainou_SI(d, verbose = FALSE), "negative or infinite")
  d <- example_dainou
  d[3, 4] <- Inf
  expect_error(dainou_SI(d, verbose = FALSE), "negative or infinite")
  d <- example_dainou
  d$id[2] <- d$id[1]
  expect_error(dainou_SI(d, verbose = FALSE), "repeated IDs")
})

test_that("weis_sync() stops on infinite values", {
  data(example_weiss)
  d <- example_weiss
  d[2, 3] <- Inf
  expect_error(weis_sync(d, verbose = FALSE), "negative or infinite")
  d <- example_weiss
  d[2, 1] <- Inf
  expect_error(weis_sync(d, verbose = FALSE), "Trait column contains infinite")
})

test_that("bolmgren_S() stops on infinite days but accepts negative days", {
  data(example_bolmgren_long)
  d <- example_bolmgren_long
  d$end_day[2] <- Inf
  expect_error(bolmgren_S(d, verbose = FALSE), "infinite values")
  d <- example_bolmgren_long
  d$onset_day <- d$onset_day - 30
  d$end_day <- d$end_day - 30
  expect_no_error(bolmgren_S(d, verbose = FALSE))
})
