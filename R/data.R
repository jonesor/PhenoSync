#' Shared Wide-Format Example Phenology Dataset
#'
#' A numeric matrix with census days in rows and individuals/species in columns.
#' Values are simulated flowering counts.
#'
#' @format A numeric matrix with 60 rows and 12 columns.
#' @name example_wide
NULL

#' Example Long Flower-Level Dataset for Bolmgren (1998)
#'
#' A simulated long-format data frame with one row per flower.
#'
#' @format A data frame with columns `plant_id`, `onset_day`, and `end_day`.
#' @name example_bolmgren_long
NULL

#' Example Inter-Sex Phenology Dataset for Dainou et al. (2012)
#'
#' A simulated data frame with individual ID and sex columns followed by
#' census-day flowering counts.
#'
#' @format A data frame with columns `id`, `sex`, and `day1` to `day24`.
#' @name example_dainou
NULL

#' Example Weis et al. (2014) Dataset
#'
#' A simulated data frame with individuals in rows, a trait column, and flowering
#' schedule columns.
#'
#' @format A data frame with `f.date` and flowering census columns.
#' @name example_weiss
NULL

#' Example Neighbourhood Phenology Dataset
#'
#' A list containing metadata and a matching flowering matrix for neighbourhood
#' synchrony calculations.
#'
#' @format A list with `metadata` and `flower_mat`.
#' @name example_neighbourhood
NULL
