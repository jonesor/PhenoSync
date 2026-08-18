#' Pianka (1974) Niche-Overlap Synchrony Index
#'
#' Calculates pairwise temporal overlap among flowering schedules using Pianka's
#' proportional overlap formulation.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days or time intervals, columns are individuals or species, and values are
#'   flowering intensity or presence/absence.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `proportional_matrix`, `pairwise_overlap`,
#'   `individual_mean_overlap`, and `population_mean`.
#'
#' @examples
#' data(example_wide)
#' pianka_overlap(example_wide, verbose = FALSE)$population_mean
#'
#' @export
pianka_overlap <- function(data, census_col = NULL, verbose = TRUE) {
  
  phen_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
  # Input check
  if (!is.numeric(phen_mat)) {
    stop("Input data must be numeric.")
  }
  
  n <- ncol(phen_mat)
  ids <- colnames(phen_mat)
  
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(n))
    colnames(phen_mat) <- ids
  }
  
  if (n < 2) {
    stop("At least two species/individuals are required.")
  }
  
  
  # STEP 1: Convert each column to proportional temporal use
  # Pij = proportion of total flowering/fruiting of species j occurring in census day i
  
  col_totals <- colSums(phen_mat, na.rm = TRUE)
  
  if (any(col_totals < 0)) {
    stop("Negative values are not allowed.")
  }
  
  # Species with total = 0 cannot have a proportional distribution
  zero_cols <- which(col_totals == 0)
  
  if (length(zero_cols) > 0) {
    warning(
      "The following columns have total abundance = 0 and will return NA in pairwise overlaps: ",
      paste(ids[zero_cols], collapse = ", ")
    )
  }
  
  prop_mat <- matrix(NA_real_, nrow = nrow(phen_mat), ncol = ncol(phen_mat))
  rownames(prop_mat) <- rownames(phen_mat)
  colnames(prop_mat) <- ids
  
  for (j in 1:n) {
    if (col_totals[j] > 0) {
      prop_mat[, j] <- phen_mat[, j] / col_totals[j]
    }
  }
  
  # STEP 2: Calculate pairwise Pianka overlap
  # Ojk = sum(Pij * Pik) / sqrt(sum(Pij^2) * sum(Pik^2))
  
  pairwise_mat <- matrix(NA_real_, nrow = n, ncol = n)
  rownames(pairwise_mat) <- ids
  colnames(pairwise_mat) <- ids
  
  for (j in 1:n) {
    for (k in 1:n) {
      
      pj <- prop_mat[, j]
      pk <- prop_mat[, k]
      
      if (all(is.na(pj)) || all(is.na(pk))) {
        pairwise_mat[j, k] <- NA
      } else {
        numerator <- sum(pj * pk, na.rm = TRUE)
        denominator <- sqrt(sum(pj^2, na.rm = TRUE) * sum(pk^2, na.rm = TRUE))
        
        pairwise_mat[j, k] <- ifelse(denominator == 0, NA, numerator / denominator)
      }
    }
  }
  
 
  # STEP 3: Species-level mean overlap
  # mean overlap of each species with all others
 
  individual_Oj <- numeric(n)
  names(individual_Oj) <- ids
  
  for (j in 1:n) {
    individual_Oj[j] <- mean(pairwise_mat[j, -j], na.rm = TRUE)
  }
  
 
  # STEP 4: Overall mean pairwise overlap
 
  
  population_mean <- mean(pairwise_mat[upper.tri(pairwise_mat)], na.rm = TRUE)
  
  if (verbose) {
    message("Pianka overlap calculated for ", n, " species/individuals.")
  }
  
  return(list(
    proportional_matrix = prop_mat,
    pairwise_overlap = pairwise_mat,
    individual_mean_overlap = individual_Oj,
    population_mean = population_mean
  ))
}

