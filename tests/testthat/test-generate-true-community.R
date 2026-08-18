test_that("phenology archetype defaults contain all seven presets", {
  d <- phenology_archetype_defaults()

  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 7L)
  expect_setequal(
    d$community_type,
    c(
      "even_weak_seasonal",
      "strongly_seasonal_even",
      "seasonal_herbaceous",
      "phenologically_staggered",
      "bimodal_seasonal",
      "mass_flowering_dominated",
      "ephemeral_pulse_flowering"
    )
  )
  expect_equal(
    d$rank_abundance_alpha[d$community_type == "mass_flowering_dominated"],
    1.50
  )
  expect_equal(
    d$target_top_species_output_share[
      d$community_type == "mass_flowering_dominated"
    ],
    0.80
  )
})


test_that("generate_true_community returns documented matrices and parameters", {
  x <- generate_true_community(
    community_type = "seasonal_herbaceous",
    n_species = 7,
    season_length = 60,
    seed = 149,
    control = list(min_individuals = 3, max_individuals = 8)
  )

  expect_s3_class(x, "true_phenology_community")
  expect_equal(nrow(x$species_parameters), 7L)
  expect_equal(nrow(x$new_flowers_species), 60L)
  expect_equal(ncol(x$new_flowers_species), 7L)
  expect_equal(nrow(x$open_flowers_species), 60L)
  expect_equal(dim(x$presence_species), dim(x$open_flowers_species))
  expect_true(is.integer(x$new_flowers_individual))
  expect_true(is.integer(x$open_flowers_individual))
  expect_true(is.integer(x$new_flowers_species))
  expect_true(is.integer(x$open_flowers_species))
  expect_true(is.logical(x$presence_individual))
  expect_true(is.logical(x$presence_species))
  expect_true(all(c(
    "realized_total_flowers",
    "realized_peak_day",
    "realized_onset_day",
    "realized_offset_day",
    "realized_duration_days"
  ) %in% names(x$species_parameters)))
})


test_that("all seven archetypes can be generated", {
  types <- phenology_archetype_defaults()$community_type

  for (type in types) {
    x <- generate_true_community(
      community_type = type,
      n_species = 7,
      season_length = 60,
      seed = 10,
      control = list(min_individuals = 3, max_individuals = 8)
    )

    expect_identical(x$community_type, type)
    expect_equal(ncol(x$new_flowers_species), 7L)
  }
})


test_that("true community generation is reproducible for a fixed seed", {
  args <- list(
    community_type = "bimodal_seasonal",
    n_species = 7,
    season_length = 60,
    seed = 151,
    control = list(min_individuals = 3, max_individuals = 8)
  )

  x1 <- do.call(generate_true_community, args)
  x2 <- do.call(generate_true_community, args)

  expect_identical(x1$species_parameters, x2$species_parameters)
  expect_identical(x1$individual_parameters, x2$individual_parameters)
  expect_identical(x1$new_flowers_individual, x2$new_flowers_individual)
})


test_that("generated community converts to PhenoSync phenology tables", {
  x <- generate_true_community(
    community_type = "even_weak_seasonal",
    n_species = 7,
    season_length = 60,
    seed = 107,
    control = list(min_individuals = 3, max_individuals = 8)
  )

  wide <- as_phenology_table(
    x,
    level = "species",
    measure = "open",
    census_column = "census_day",
    census_interval = 1
  )

  expect_s3_class(wide, "data.frame")
  expect_identical(names(wide)[1], "census_day")
  expect_equal(nrow(wide), 60L)
  expect_equal(ncol(wide), 8L)
  expect_equal(wide$census_day, seq_len(60L))
})



test_that("generate_true_community defaults are 50 species and 90 days", {
  expect_identical(formals(generate_true_community)$n_species, 50)
  expect_identical(formals(generate_true_community)$season_length, 90)

  d <- phenology_archetype_defaults()
  expect_true(all(d$species_default == 50L))
  expect_true(all(d$season_length_days_default == 90L))
})


test_that("as_phenology_table samples true phenology at the requested interval", {
  x <- generate_true_community(
    community_type = "even_weak_seasonal",
    n_species = 7,
    season_length = 60,
    seed = 107,
    control = list(min_individuals = 3, max_individuals = 8)
  )

  weekly <- as_phenology_table(
    x,
    level = "species",
    measure = "open",
    census_column = "day",
    census_interval = 7
  )

  expected_days <- seq.int(1L, 60L, by = 7L)
  expect_equal(weekly$day, expected_days)
  expect_equal(nrow(weekly), length(expected_days))
  expect_equal(ncol(weekly), 8L)
  expect_equal(
    unname(as.matrix(weekly[, -1, drop = FALSE])),
    unname(x$open_flowers_species[expected_days, , drop = FALSE])
  )
})


test_that("as_phenology_table validates census_interval", {
  x <- generate_true_community(
    n_species = 5,
    season_length = 30,
    seed = 1,
    control = list(min_individuals = 2, max_individuals = 4)
  )

  expect_error(as_phenology_table(x, census_interval = 0), "positive whole number")
  expect_error(as_phenology_table(x, census_interval = 2.5), "positive whole number")
})
