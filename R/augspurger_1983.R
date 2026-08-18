#' Augspurger (1983) Individual and Population Synchrony
#'
#' Calculates Augspurger's individual synchrony index (`Xi`) and population
#' synchrony (`Z`) from a wide flowering matrix. Flower-count values greater
#' than zero are treated as flowering presence.
#'
#' @param flower_mat A numeric matrix or data frame in wide format. Rows are
#'   census days, columns are individuals, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param na_as_zero Logical. If `TRUE`, missing values are converted to zero
#'   before calculating flowering presence/absence.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `Xi_table`, `Z`, `presence_matrix`, and `overlap_matrix`.
#'
#' @examples
#' data(example_wide)
#' augspurger_Xi(example_wide, verbose = FALSE)$Z
#'
#' @export
augspurger_Xi <- function(flower_mat,
                            census_col = NULL,
                            na_as_zero = TRUE, #Note: Individuals with fi == 0 (never flowered) will have Xi = NA (and are included in n)
                            verbose = TRUE) {
  
  flower_mat <- prepare_census_matrix(flower_mat, census_col = census_col, data_arg = "flower_mat")
  
  if (na_as_zero) {
    flower_mat[is.na(flower_mat)] <- 0
  }
  
  if (is.null(colnames(flower_mat))) {
    colnames(flower_mat) <- paste0("ind", seq_len(ncol(flower_mat)))
  }
  
  ids <- colnames(flower_mat)
  n <- ncol(flower_mat)
  
  # STEP 1: Convert flower counts to presence/absence
  
  presence_mat <- (flower_mat > 0) * 1
  
 
  # STEP 2: Flowering duration for each individual
  
  fi <- colSums(presence_mat)
  

  # STEP 3: Pairwise overlap matrix
  overlap_mat <- t(presence_mat) %*% presence_mat #t() is a transpose function, thus the matrix multiplication gives ind x indi matrix. with 1 is overlap and 0 is non overlap. 
  
  # Now let's add the 1s to get total number of overlap days
  # sum of overlaps of each individual with all others
  
  sum_ei <- rowSums(overlap_mat) - diag(overlap_mat) #diag function: Extract or replace the diagonal of a matrix, or construct a diagonal matrix. this removes the ovserlap of an individual with itself.
 
  # STEP 4: Individual synchrony Xi
  
  Xi <- rep(NA_real_, n)
  
  for (i in 1:n) {
    if (fi[i] > 0) {
      Xi[i] <- sum_ei[i] / ((n - 1) * fi[i])
    }
  }
 
  # STEP 5: Population synchrony Z
   Z <- mean(Xi, na.rm = TRUE)
  
  Xi_table <- data.frame(
    id = ids,
    fi = as.numeric(fi),
    sum_ei = as.numeric(sum_ei),
    Xi = as.numeric(Xi)
  )
  
  if (verbose) {
    message(sprintf("Number of individuals (n) = %d", n))
    message(sprintf("Augspurger population synchrony Z = %.4f", Z))
    
    if (any(fi == 0)) {
      message(sprintf(
        "Note: %d individual(s) never flowered (fi = 0); Xi set to NA.",
        sum(fi == 0)
      ))
    }
  }
  
  return(list(
    Xi_table = Xi_table,
    Z = Z,
    presence_matrix = presence_mat,
    overlap_matrix = overlap_mat
  ))
}

