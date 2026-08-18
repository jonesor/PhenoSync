#' Marquis (1988) Flowering Overlap Index
#'
#' Calculates individual and population-level synchrony based on overlap between
#' an individual's flowering days and the population flowering profile.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param na_as_zero Logical. If `TRUE`, missing values are converted to zero.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `SM_individual`, `SM_population_mean`, `p_t`,
#'   `flowering_presence`, and `total_output`.
#'
#' @examples
#' data(example_wide)
#' marquis_SM(example_wide, verbose = FALSE)$SM_population_mean
#'
#' @export
marquis_SM <- function(data, census_col = NULL, na_as_zero = TRUE, verbose = TRUE) {
  
  # STEP 0: Input checks
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  
  if (!is.numeric(flower_mat)) {
    stop("Input data must be a numeric matrix or data frame in wide format.")
  }
  
  if (na_as_zero) {
    flower_mat[is.na(flower_mat)] <- 0
  }
  
  if (any(is.na(flower_mat))) {
    stop("Input contains NA values. Set na_as_zero = TRUE or remove missing values.")
  }
  
  if (any(flower_mat < 0)) {
    stop("Negative values are not allowed.")
  }
  
  n <- ncol(flower_mat)
  T <- nrow(flower_mat)
  
  if (n < 2) {
    stop("At least two individuals are required.")
  }
  
  ids <- colnames(flower_mat)
  if (is.null(ids)) {
    ids <- paste0("ind", seq_len(n))
    colnames(flower_mat) <- ids
  }
  

  # STEP 1: Compute p_t
  # p_t = proportion of individuals flowering at census day t
  # Individual is considered flowering if x_it > 0

  
  flowering_presence <- (flower_mat > 0) * 1
  p_t <- rowMeans(flowering_presence)
  

  # STEP 2: Compute individual Marquis SM
  total_output <- colSums(flower_mat)
  
  if (any(total_output == 0)) {
    warning(
      "Some individuals have total annual output = 0. Their SM will be NA: ",
      paste(ids[total_output == 0], collapse = ", ")
    )
  }
  
  SM_values <- numeric(n)
  names(SM_values) <- ids
  SM_values[] <- NA_real_
  
  for (i in 1:n) {
    
    x_t <- flower_mat[, i]
    total_i <- total_output[i]
    
    if (total_i == 0) {
      next
    }
    
    weight_t <- x_t / total_i
    SM_values[i] <- sum(weight_t * p_t)
  }
  

  # STEP 3: Population mean
  SM_population_mean <- mean(SM_values, na.rm = TRUE)
  
  if (verbose) {
    message("Marquis SM calculated for ", n, " individuals across ", T, " census days.")
    message("Mean population SM = ", round(SM_population_mean, 4))
  }
  
  return(list(
    SM_individual = SM_values,
    SM_population_mean = SM_population_mean,
    p_t = p_t,
    flowering_presence = flowering_presence,
    total_output = total_output
  ))
}

