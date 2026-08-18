#' Michalski and Durka (2007) Coefficient of Variation Index
#'
#' Calculates population and individual coefficients of variation in flowering
#' intensity over census dates.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts or flowering intensity.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `population_profile`, `population_CV`, `individual_CV`,
#'   and `mean_individual_CV`.
#'
#' @examples
#' data(example_wide)
#' michalski_durka_CV(example_wide, verbose = FALSE)$population_CV
#'
#' @export
michalski_durka_CV <- function(data, census_col = NULL, verbose = TRUE) {
  
 
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
  
  n <- ncol(flower_mat)
  ids <- colnames(flower_mat)
  
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  if (n < 1) {
    stop("At least one species/individual is required.")
  }
  
  if (nrow(flower_mat) < 2) {
    stop("At least two census days are required to calculate standard deviation.")
  }
  
 
  # STEP 1: Population-level flowering profile through time
  pop_profile <- rowSums(flower_mat)
  pop_mean <- mean(pop_profile)
  pop_sd <- stats::sd(pop_profile)
  
  population_CV <- if (pop_mean == 0) NA_real_ else pop_sd / pop_mean
  
  
  # STEP 2: Individual/species-level CV through time
  
  individual_CV <- numeric(n)
  names(individual_CV) <- ids
  
  for (i in seq_len(n)) {
    x <- flower_mat[, i]
    x_mean <- mean(x)
    x_sd <- stats::sd(x)
    
    individual_CV[i] <- if (x_mean == 0) NA_real_ else x_sd / x_mean
  }
  
  
  # STEP 3: Mean individual CV
  mean_individual_CV <- mean(individual_CV, na.rm = TRUE)
  if (is.nan(mean_individual_CV)) mean_individual_CV <- NA_real_
  
  if (verbose) {
    message("Michalski & Durka CV calculated for ", n, " species/individuals.")
    message("Population CV = ", round(population_CV, 4))
  }
  
  return(list(
    population_profile = pop_profile,
    population_CV = population_CV,
    individual_CV = individual_CV,
    mean_individual_CV = mean_individual_CV
  ))
}

