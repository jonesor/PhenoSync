#' Shrotri et al. (2025) Abundance-Weighted Synchrony Index
#'
#' Calculates abundance-weighted pairwise temporal synchrony among species or
#' individuals from wide flowering abundance data.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are species or individuals, and values are abundance or
#'   flowering intensity.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `pairwise_SIabd`, `individual_SIabd`, and
#'   `community_mean`.
#'
#' @examples
#' data(example_wide)
#' shrotri_SI_abd(example_wide, verbose = FALSE)$community_mean
#'
#' @export
shrotri_SI_abd <- function(data, census_col = NULL, verbose = TRUE) {
  
 
  # STEP 0: Input checks
  
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
  if (!is.numeric(flower_mat)) {
    stop("Input data must be a numeric matrix or data frame in wide format.")
  }
  
  if (any(is.na(flower_mat))) {
    stop("Input contains NA values. Please remove or replace missing values before running the function.")
  }
  
  if (any(flower_mat < 0)) {
    stop("Negative values are not allowed.")
  }
  
  n_sp <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(n_sp))
    colnames(flower_mat) <- ids
  }
  
  if (n_sp < 2) {
    stop("At least two species are required.")
  }
  
 
  # STEP 1: Pairwise SI_abd
  
  pairwise_mat <- matrix(NA_real_, nrow = n_sp, ncol = n_sp)
  rownames(pairwise_mat) <- ids
  colnames(pairwise_mat) <- ids
  
  for (j in 1:n_sp) {
    for (k in 1:n_sp) {
      
      if (j == k) {
        pairwise_mat[j, k] <- 1
      } else {
        
        x <- flower_mat[, j]
        y <- flower_mat[, k]
        
        # co-flowering period: both species present in flowering phenophase
        valid_days <- which(x > 0 & y > 0)
        
        if (length(valid_days) == 0) {
          pairwise_mat[j, k] <- NA
        } else {
          x_valid <- x[valid_days]
          y_valid <- y[valid_days]
          
          a_i <- pmin(x_valid, y_valid)   # less common species each day
          b_i <- pmax(x_valid, y_valid)   # more common species each day
          
          pairwise_mat[j, k] <- mean(a_i / b_i)
        }
      }
    }
  }
  
 
  # STEP 2: Species-level mean SI_abd
 
  individual_SIabd <- numeric(n_sp)
  names(individual_SIabd) <- ids
  
  for (j in 1:n_sp) {
    individual_SIabd[j] <- mean(pairwise_mat[j, -j], na.rm = TRUE)
  }
  
  # if a species has no valid pairwise comparisons, mean(...) gives NaN
  individual_SIabd[is.nan(individual_SIabd)] <- NA_real_
  
 
  # STEP 3: Community mean SI_abd
 
  community_mean <- mean(pairwise_mat[upper.tri(pairwise_mat)], na.rm = TRUE)
  
  if (is.nan(community_mean)) {
    community_mean <- NA_real_
  }
  
  if (verbose) {
    message("SI_abd calculated for ", n_sp, " species.")
  }
  
  return(list(
    pairwise_SIabd = pairwise_mat,
    individual_SIabd = individual_SIabd,
    community_mean = community_mean
  ))
}

