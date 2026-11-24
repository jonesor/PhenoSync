#' Koenig et al. (2003) Mean Pairwise Correlation Synchrony
#'
#' Calculates synchrony as the mean of all pairwise Pearson correlations among
#' flowering or reproductive time series.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   dates or years, columns are individuals or species, and values are flowering
#'   or reproductive intensities.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param use Missing-value handling argument passed to [stats::cor()].
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `pairwise_correlation`, `pairwise_values`,
#'   `individual_mean_r`, and `population_rp`.
#'
#' @examples
#' data(example_wide)
#' koenig_rp(example_wide, verbose = FALSE)$population_rp
#'
#' @export
koenig_rp <- function(data, census_col = NULL, use = "pairwise.complete.obs", verbose = TRUE) {
  
  # STEP 0: Input checks
 
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
  if (!is.numeric(flower_mat)) {
    stop("Input data must be a numeric matrix or data frame in wide format.")
  }
  
  n <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  if (n < 2) {
    stop("At least two species/individuals are required.")
  }
  
  # warn for zero-variance series, because Pearson r is undefined
  sd_vec <- apply(flower_mat, 2, stats::sd, na.rm = TRUE)
  if (any(sd_vec == 0 | is.na(sd_vec))) {
    warning(
      "Some columns have zero variance; pairwise correlations involving them will be NA: ",
      paste(ids[sd_vec == 0 | is.na(sd_vec)], collapse = ", ")
    )
  }
  
 
  # STEP 1: Pairwise Pearson correlation matrix
  cor_mat <- stats::cor(flower_mat, method = "pearson", use = use)
  
  
  # STEP 2: Extract unique pairwise correlations
  pair_vals <- cor_mat[upper.tri(cor_mat)]
  
  # remove NA correlations before averaging
  valid_pair_vals <- pair_vals[!is.na(pair_vals)]
  
  
  # STEP 3: Population synchrony
  if (length(valid_pair_vals) == 0) {
    population_rp <- NA_real_
  } else {
    population_rp <- mean(valid_pair_vals)
  }
  
  # STEP 4: Species-level mean correlation with all others
  individual_mean_r <- numeric(n)
  names(individual_mean_r) <- ids
  
  for (i in seq_len(n)) {
    individual_mean_r[i] <- mean(cor_mat[i, -i], na.rm = TRUE)
  }
  
  individual_mean_r[is.nan(individual_mean_r)] <- NA_real_
  
  if (verbose) {
    message("Koenig population synchrony calculated for ", n, " species/individuals.")
    message("Number of unique pairs = ", n * (n - 1) / 2)
    message("Population mean r_p = ", round(population_rp, 4))
  }
  
  return(list(
    pairwise_correlation = cor_mat,
    pairwise_values = valid_pair_vals,
    individual_mean_r = individual_mean_r,
    population_rp = population_rp
  ))
}

