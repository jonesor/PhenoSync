#' Primack (1980) Strong-Flowering Overlap Index
#'
#' Calculates pairwise and population overlap among strong flowering periods,
#' where strong flowering is defined relative to each individual's maximum
#' flowering intensity.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param strong_threshold Numeric threshold between 0 and 1. A census day is
#'   treated as strong flowering when the count is at least this fraction of the
#'   individual's maximum count.
#' @param na_as_zero Logical. If `TRUE`, missing values are converted to zero.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `pairwise_table`, `Z`, and `strong_matrix`.
#'
#' @examples
#' data(example_wide)
#' primack_Z(example_wide, verbose = FALSE)$Z
#'
#' @export
primack_Z <- function(data,
                      census_col = NULL,
                      strong_threshold = 0.5,
                      na_as_zero = TRUE,
                      verbose = TRUE) {
  
  # STEP 0: Input checks + coercion
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data")
  mode(flower_mat) <- "numeric"
  
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
  if (n < 2) {
    stop("At least two individuals are required.")
  }
  
  ids <- colnames(flower_mat)
  if (is.null(ids)) {
    ids <- paste0("ind", seq_len(n))
    colnames(flower_mat) <- ids
  }
  
  # STEP 1: Identify strong flowering days
  max_per_ind <- apply(flower_mat, 2, max, na.rm = TRUE)
  
  strong_mat <- matrix(
    0,
    nrow = nrow(flower_mat),
    ncol = ncol(flower_mat),
    dimnames = dimnames(flower_mat)
  )
  
  valid_cols <- max_per_ind > 0
  
  if (any(valid_cols)) {
    strong_mat[, valid_cols] <- sweep(
      flower_mat[, valid_cols, drop = FALSE],
      2,
      max_per_ind[valid_cols] * strong_threshold,
      FUN = ">="
    ) * 1
  }
  
  strong_days <- colSums(strong_mat)
  
  # STEP 2: Pairwise overlap calculation
  pairs <- utils::combn(n, 2)
  i <- pairs[1, ]
  j <- pairs[2, ]
  
  a <- as.integer(colSums(strong_mat[, i, drop = FALSE] == 1 & strong_mat[, j, drop = FALSE] == 1))
  b <- unname(pmin(strong_days[i], strong_days[j]))
  
  c_values <- ifelse(b == 0, NA_real_, a / b)
  
  pair_results <- data.frame(
    id_i = ids[i],
    id_j = ids[j],
    a_overlap = a,
    b_shortest = b,
    c_ij = c_values
  )
  
  # STEP 3: Population mean
  Z <- mean(c_values, na.rm = TRUE)
  
  if (verbose) {
    message(sprintf("Number of individuals (n) = %d", n))
    message(sprintf("Number of pairwise comparisons (N) = %d", length(c_values)))
    message(sprintf("Primack population overlap Z = %.4f", Z))
    
    if (any(!valid_cols)) {
      message("Note: some individuals had all-zero flowering schedules and were treated as having no strong flowering days.")
    }
  }
  
  return(list(
    pairwise_table = pair_results,
    Z = Z,
    strong_matrix = strong_mat
  ))
}

