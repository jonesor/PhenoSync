#' Shrotri et al. (2025) Temporal Synchrony Index
#'
#' Calculates pairwise temporal synchrony from presence/absence overlap in wide
#' phenological data.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are species, and values are phenophase counts or abundances.
#'   Values greater than zero are treated as presence.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `pairwise_SItemp`, `individual_SItemp`, and
#'   `community_mean`.
#'
#' @examples
#' data(example_wide)
#' shrotri_SI_temp(example_wide, verbose = FALSE)$community_mean
#'
#' @export
shrotri_SI_temp <- function(data, census_col = NULL, verbose = TRUE) {
  
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
  if (!is.numeric(flower_mat)) {
    stop("Input data must be numeric.")
  }
  
  n <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  if (n < 2) {
    stop("At least two species are required.")
  }
  
  # Convert to presence/absence
  presence_mat <- (flower_mat > 0) * 1
  
  # Empty matrix to store pairwise SI_temp values
  pairwise_mat <- matrix(NA, nrow = n, ncol = n)
  rownames(pairwise_mat) <- ids
  colnames(pairwise_mat) <- ids
  
  # Calculate pairwise temporal synchrony
  for (i in 1:n) {
    for (j in 1:n) {
      
      if (i == j) {
        pairwise_mat[i, j] <- 1
      } else {
        F1 <- sum(presence_mat[, i])                         # flowering duration of species i
        F2 <- sum(presence_mat[, j])                         # flowering duration of species j
        F0 <- sum(presence_mat[, i] & presence_mat[, j])     # shared flowering duration
        
        denom <- (F1 + F2) - F0
        
        pairwise_mat[i, j] <- ifelse(denom == 0, 0, F0 / denom)
      }
    }
  }
  
  # Mean synchrony of each species with all others
  individual_SItemp <- numeric(n)
  names(individual_SItemp) <- ids
  
  for (i in 1:n) {
    individual_SItemp[i] <- mean(pairwise_mat[i, -i], na.rm = TRUE)
  }
  
  # Overall mean synchrony across all unique pairs
  community_mean <- mean(pairwise_mat[upper.tri(pairwise_mat)], na.rm = TRUE)
  
  if (verbose) {
    message("SI temp calculated for ", n, " species.")
  }
  
  return(list(
    pairwise_SItemp = pairwise_mat,
    individual_SItemp = individual_SItemp,
    community_mean = community_mean
  ))
}

