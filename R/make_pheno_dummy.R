#' Generate dummy phenology datasets
#'
#' Creates small simulated phenology datasets for testing and examples. The
#' generated data are intended for trial runs, demonstrations, and package tests,
#' not for biological inference.
#'
#' @param format Character. One of `"wide"`, `"bolmgren_long"`, `"dainou"`,
#'   `"weis"`, or `"neighbourhood"`. For backward compatibility,
#'   `"weiss"` is also accepted and treated as `"weis"`.
#' @param n_individuals Number of individuals to simulate. Used unless
#'   `n_species` is supplied for a wide species-level dataset.
#' @param n_species Optional number of species/columns to simulate.
#' @param n_days Number of census days.
#' @param binary Logical. If `TRUE`, values are returned as presence/absence.
#'   If `FALSE`, values are simulated as counts.
#' @param unit Character. For `format = "wide"`, should columns represent
#'   `"individual"` or `"species"`?
#' @param mean_duration Approximate mean flowering duration. If `NULL`, a
#'   default based on `n_days` is used.
#' @param mean_peak_count Approximate mean peak count for numeric data.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A matrix, data frame, or list depending on `format`:
#' \describe{
#'   \item{`wide`}{A matrix with census days in rows and individuals/species in columns.}
#'   \item{`bolmgren_long`}{A flower-level data frame with plant IDs, onset days, and end days.}
#'   \item{`dainou`}{An individual-by-day data frame with `id` and `sex` columns.}
#'   \item{`weis`}{An individual-by-day data frame with the trait column first.}
#'   \item{`neighbourhood`}{A list with `metadata` and `flower_mat`.}
#' }
#'
#' @examples
#' wide_counts <- make_pheno_dummy(format = "wide", n_species = 6, n_days = 40, seed = 1)
#' dim(wide_counts)
#'
#' bolmgren_dat <- make_pheno_dummy(format = "bolmgren_long", seed = 1)
#' head(bolmgren_dat)
#'
#' @export
make_pheno_dummy <- function(format = c("wide", "bolmgren_long", "dainou", "weis", "weiss", "neighbourhood"),
                             n_individuals = 12,
                             n_species = NULL,
                             n_days = 60,
                             binary = FALSE,
                             unit = c("individual", "species"),
                             mean_duration = NULL,
                             mean_peak_count = 25,
                             seed = NULL) {
  format <- match.arg(format)
  if (identical(format, "weiss")) format <- "weis"
  unit <- match.arg(unit)

  if (!is.null(seed)) set.seed(seed)

  if (!is.numeric(n_days) || length(n_days) != 1 || is.na(n_days) || n_days < 2) {
    stop("n_days must be a single number >= 2.", call. = FALSE)
  }
  n_days <- as.integer(n_days)

  n_units <- if (!is.null(n_species)) n_species else n_individuals
  if (!is.numeric(n_units) || length(n_units) != 1 || is.na(n_units) || n_units < 2) {
    stop("n_individuals or n_species must be a single number >= 2.", call. = FALSE)
  }
  n_units <- as.integer(n_units)

  if (is.null(mean_duration)) {
    mean_duration <- max(5, round(n_days / 4))
  }

  make_wide_matrix <- function(n_units, n_days, binary, unit, mean_duration, mean_peak_count) {
    ids <- if (unit == "species") {
      paste0("sp", sprintf("%02d", seq_len(n_units)))
    } else {
      paste0("ind", sprintf("%02d", seq_len(n_units)))
    }

    day_names <- paste0("day", seq_len(n_days))
    mat <- matrix(0, nrow = n_days, ncol = n_units,
                  dimnames = list(day_names, ids))

    peaks <- sample(seq_len(n_days), n_units, replace = TRUE)
    durations <- pmax(2, round(stats::rnorm(n_units, mean_duration, mean_duration / 4)))
    peak_counts <- pmax(1, round(stats::rpois(n_units, lambda = mean_peak_count)))

    for (i in seq_len(n_units)) {
      sigma <- max(1, durations[i] / 4)
      curve <- peak_counts[i] * exp(-0.5 * ((seq_len(n_days) - peaks[i]) / sigma)^2)
      curve[curve < 0.05 * max(curve)] <- 0

      if (binary) {
        mat[, i] <- as.integer(curve > 0)
      } else {
        mat[, i] <- stats::rpois(n_days, lambda = pmax(curve, 0))
        mat[curve == 0, i] <- 0
      }
    }

    storage.mode(mat) <- "numeric"
    mat
  }

  if (format == "wide") {
    return(make_wide_matrix(n_units, n_days, binary, unit, mean_duration, mean_peak_count))
  }

  if (format == "bolmgren_long") {
    plant_ids <- paste0("plant", sprintf("%02d", seq_len(n_units)))
    n_flowers <- pmax(2, stats::rpois(n_units, lambda = 8))

    out <- do.call(rbind, lapply(seq_along(plant_ids), function(i) {
      peak <- sample(seq_len(n_days), 1)
      onset <- round(stats::rnorm(n_flowers[i], mean = peak, sd = max(2, n_days / 12)))
      onset <- pmin(pmax(onset, 1), n_days)
      persistence <- sample(2:8, n_flowers[i], replace = TRUE)
      end <- pmin(onset + persistence - 1, n_days)
      data.frame(
        plant_id = plant_ids[i],
        onset_day = onset,
        end_day = end,
        stringsAsFactors = FALSE
      )
    }))

    rownames(out) <- NULL
    return(out)
  }

  if (format == "dainou") {
    ids <- paste0("ind", sprintf("%02d", seq_len(n_units)))
    sex <- rep(c("F", "M"), length.out = n_units)
    sex <- sample(sex)
    flower_mat <- make_wide_matrix(n_units, n_days, binary, "individual", mean_duration, mean_peak_count)
    out <- data.frame(id = ids, sex = sex, as.data.frame(t(flower_mat)),
                      check.names = FALSE, stringsAsFactors = FALSE)
    names(out) <- c("id", "sex", paste0("day", seq_len(n_days)))
    rownames(out) <- NULL
    return(out)
  }

  if (format == "weis") {
    ids <- paste0("ind", sprintf("%02d", seq_len(n_units)))
    flower_mat <- make_wide_matrix(n_units, n_days, binary, "individual", mean_duration, mean_peak_count)
    flower_mat <- t(flower_mat)
    trait <- apply(flower_mat, 1, function(x) {
      if (all(x == 0)) NA_real_ else stats::weighted.mean(seq_len(n_days), w = x + 1e-6)
    })
    out <- data.frame(f.date = trait, flower_mat, check.names = FALSE)
    rownames(out) <- ids
    names(out)[-1] <- paste0("day", seq_len(n_days))
    return(out)
  }

  if (format == "neighbourhood") {
    ids <- paste0("ind", sprintf("%02d", seq_len(n_units)))
    flower_mat <- make_wide_matrix(n_units, n_days, binary, "individual", mean_duration, mean_peak_count)
    colnames(flower_mat) <- ids

    metadata <- data.frame(
      id = ids,
      x = stats::runif(n_units, 0, 100),
      y = stats::runif(n_units, 0, 100),
      species = sample(paste0("sp", seq_len(max(2, min(4, n_units)))), n_units, replace = TRUE),
      stringsAsFactors = FALSE
    )

    return(list(metadata = metadata, flower_mat = flower_mat))
  }
}
