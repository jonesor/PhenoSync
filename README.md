# PhenoSync

`PhenoSync` provides R functions for calculating plant phenological synchrony indices from wide phenology matrices and selected index-specific data structures, including flower-level, inter-sex, trait-schedule, and spatial-neighbourhood data.


## Main wide-data assumption

For `calculate_indices()`, the standard input is a wide phenology table:

```text
one row    = one census day or census date
one column = one individual or one species, specified by `columns_are`
one cell   = flowering state, count, abundance, or intensity
```

If census names or dates are stored in a data column, pass that column explicitly with `census_col`.

## Main functions

Core index functions include:

- `albert_Si()`
- `augspurger_Xi()`
- `bolmgren_S()`
- `dainou_SI()`
- `freitas_bolmgren_S()`
- `koenig_rp()`
- `mahoro_Si()`
- `mahoro_Ri()`
- `marquis_SM()`
- `michalski_durka_CV()`
- `pianka_overlap()`
- `poole_rathcke()`
- `primack_Z()`
- `schoener_SI()`
- `shrotri_SI_temp()`
- `shrotri_SI_abd()`
- `wadgymar_Si()`
- `weis_sync()`

Helper functions:

- `make_pheno_dummy()` generates trial phenology datasets.
- `calculate_indices()` calculates several compatible indices from one dataset.
- `define_index()` explains an index, including definition, formula, citation, input, and output structure.

## Local development workflow

From inside the package project:

```r
source("data-raw/create_all_example_data.R")
devtools::document()
devtools::test()
devtools::install()
```

Then restart R and use:

```r
library(PhenoSync)

data(example_wide)
albert_Si(example_wide)
```

## Trial examples

```r
x <- make_pheno_dummy(format = "wide", n_species = 6, n_days = 40, seed = 1)
calculate_indices(x, indices = c("albert_Si", "pianka_overlap"))

y <- make_pheno_dummy(format = "weis", n_individuals = 10, n_days = 30, seed = 1)
weis_sync(y, trait_col = "f.date")

define_index("pianka")
```


## Update v0.0.0.9003

`calculate_indices()` supports the standard wide workflow only: one row is one census day, one column is one individual or species, and each cell is flowering state/count/abundance/intensity. Use `columns_are = "individuals"` or `"species"`, `binary = TRUE/FALSE/NULL`, and `census_col` when the census labels are stored in a column. The base R summary plot uses a fixed 0–1 x-axis for bounded indices in the left panel and a data-driven x-axis for other-scale indices in the right panel.


## Census-column handling

For ordinary CSV/Excel-style wide phenology data, keep the census names or dates
in a named column and pass that column explicitly with `census_col`. For example:

```r
wide_df <- data.frame(
  census_day = paste0("day", 1:5),
  sp1 = c(0, 1, 3, 0, 0),
  sp2 = c(0, 0, 2, 4, 1)
)

res <- calculate_indices(
  wide_df,
  census_col = "census_day",
  columns_are = "species",
  binary = FALSE
)

albert_Si(wide_df, census_col = "census_day")
```

The specified census column is removed from the numeric phenology matrix and
preserved as row names. If your data are already a numeric matrix with meaningful
row names, `census_col = NULL` is still acceptable.


## Update v0.0.0.9008

`calculate_indices()` no longer has a `data_type` argument and no longer exposes `detect_pheno_data()`. The function is now deliberately simple: rows are census days, columns are individuals or species as specified by `columns_are`, and cells are flowering state/count/abundance/intensity. For Bolmgren, Dainou, Weis, and Schoener special data structures, call the corresponding single-index function directly.

## True phenological community generator

`generate_true_community()` generates completely known daily flowering
phenologies for seven contrasting plant-community archetypes. By default, a
community contains 50 species, a 90-day season, and 20--500 individuals per
species. The returned object contains individual- and species-level parameter
tables plus daily matrices for newly opened flowers, currently open flowers,
and flowering presence/absence.

Available archetypes are:

- `even_weak_seasonal`
- `strongly_seasonal_even`
- `seasonal_herbaceous`
- `phenologically_staggered`
- `bimodal_seasonal`
- `mass_flowering_dominated`
- `ephemeral_pulse_flowering`

Inspect the complete preset table with:

```r
phenology_archetype_defaults()
```

Generate and inspect a community with:

```r
community <- generate_true_community(
  community_type = "seasonal_herbaceous",
  n_species = 50,
  season_length = 90,
  seed = 149
)

community
summarise_true_community(community)
validate_true_community(community)
plot_true_community(community)
```

Convert a generated matrix to the standard PhenoSync wide-data structure:

```r
species_open <- as_phenology_table(
  community,
  level = "species",
  measure = "open",
  census_column = "day",
  census_interval = 7
)

res <- calculate_indices(
  species_open,
  census_col = "day",
  columns_are = "species",
  binary = FALSE
)
```

Preset values can be modified without rewriting the simulator. For example:

```r
modified <- generate_true_community(
  community_type = "seasonal_herbaceous",
  seed = 202,
  control = list(
    rank_abundance_alpha = 1.20,
    duration_probabilities = c(
      short = 0.60,
      medium = 0.35,
      long = 0.05
    )
  )
)
```

## Authors

Saket S. Shrotri (`shrotri@biology.sdu.dk`) and Owen R. Jones
(`jones@biology.sdu.dk`), Faculty of Science, Department of Biology, PopBio,
University of Southern Denmark, Campusvej 55, DK-5230, Odense M.

## Update v0.0.0.9009

Added the seven-archetype true-community simulator and its conversion,
validation, summary, and plotting helpers. Package author metadata and package
description were also updated.


## Update v0.0.0.9010

`generate_true_community()` now defaults to 50 species and a 90-day season.
`as_phenology_wide()` was renamed to `as_phenology_table()`. The conversion
helper can return species- or individual-level wide phenology tables and can
sample the daily true phenology using `census_interval`; `census_interval = 1`
retains every day.
