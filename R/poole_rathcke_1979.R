#' Poole and Rathcke (1979) Peak-Spacing Index
#'
#' Extracts species peak flowering dates from a wide matrix and evaluates
#' whether peak dates are regularly spaced, aggregated, or approximately random
#' across the flowering season.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are species, and values are counts of open flowers.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param season_length Numeric season length used to normalize peak dates. If
#'   `NULL`, the number of rows in `data` is used.
#' @param perform_test Logical. If `TRUE`, calculates an approximate chi-square
#'   style test.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with peak dates, ordered peaks, interval distances, observed
#'   and expected spacing variance, ratio, optional test statistics, and an
#'   interpretation.
#'
#' @examples
#' data(example_wide)
#' poole_rathcke(example_wide, verbose = FALSE)$ratio
#'
#' @export
poole_rathcke <- function(data,
                          census_col = NULL,
                          season_length = NULL,
                          perform_test = TRUE,
                          verbose = TRUE) {
  
 
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
  
  k <- ncol(flower_mat)
  n_days <- nrow(flower_mat)
  
  if (k < 2) {
    stop("At least two species are required.")
  }
  
  if (n_days < 1) {
    stop("The dataset must contain at least one census day.")
  }
  
  if (is.null(season_length)) {
    season_length <- n_days
  }
  
  if (!is.numeric(season_length) || length(season_length) != 1 || is.na(season_length) || season_length <= 0) {
    stop("season_length must be a single positive numeric value.")
  }
  
  ids <- colnames(flower_mat)
  if (is.null(ids)) {
    ids <- paste0("sp", seq_len(k))
    colnames(flower_mat) <- ids
  }
  
  
  # STEP 1: Extract peak dates from the wide flowering matrix
  # Peak date = first census day with maximum flowering count
  
  peak_dates <- numeric(k)
  names(peak_dates) <- ids
  
  for (i in seq_len(k)) {
    x <- flower_mat[, i]
    
    if (all(x == 0)) {
      peak_dates[i] <- NA_real_
    } else {
      peak_dates[i] <- which.max(x)
    }
  }
  
  # Species with no flowering cannot have a peak date
  valid <- !is.na(peak_dates)
  
  if (sum(valid) < 2) {
    stop("Fewer than two species have valid peak dates.")
  }
  
  if (any(!valid)) {
    warning(
      "Some species had no flowering and were excluded: ",
      paste(ids[!valid], collapse = ", ")
    )
  }
  
  peak_dates <- peak_dates[valid]
  ids_valid <- ids[valid]
  k_valid <- length(peak_dates)
  
  
  # STEP 2: Normalize peak dates to [0, 1] and order them
 
  
  peak_df <- data.frame(
    species = ids_valid,
    peak_date = peak_dates,
    normalized_peak = peak_dates / season_length,
    stringsAsFactors = FALSE
  )
  
  peak_df <- peak_df[order(peak_df$normalized_peak), , drop = FALSE]
  
  y <- peak_df$normalized_peak
  
 
  # STEP 3: Distances between adjacent peaks, including boundaries
  
  D <- numeric(k_valid + 1)
  D[1] <- y[1]
  
  if (k_valid > 1) {
    D[2:k_valid] <- diff(y)
  }
  
  D[k_valid + 1] <- 1 - y[k_valid]
  
  interval_labels <- c(
    paste0("0_to_", peak_df$species[1]),
    if (k_valid > 1) paste0(peak_df$species[1:(k_valid - 1)], "_to_", peak_df$species[2:k_valid]),
    paste0(peak_df$species[k_valid], "_to_1")
  )
  
  distance_table <- data.frame(
    interval = interval_labels,
    distance = D,
    stringsAsFactors = FALSE
  )
  
  
  # STEP 4: Observed variance of interval distances
  
  
  D_bar <- mean(D)
  P <- sum((D - D_bar)^2) / (k_valid + 1)
  
  
  # STEP 5: Expected variance under random spacing
  
  expected_P <- k_valid / ((k_valid + 1)^2 * (k_valid + 2))
  ratio <- P / expected_P
  
  
  # STEP 6: Optional chi-square style test
  
  chi_square <- NA_real_
  df <- NA_integer_
  p_value <- NA_real_
  
  if (perform_test) {
    chi_square <- k_valid * ratio
    df <- k_valid
    p_value <- stats::pchisq(chi_square, df = df, lower.tail = FALSE)
  }
  
  
  # STEP 7: Interpretation
  
  tol <- 0.05
  
  interpretation <- if (ratio < (1 - tol)) {
    "Regular spacing"
  } else if (ratio > (1 + tol)) {
    "Aggregated peaks"
  } else {
    "Approximately random spacing"
  }
  
  if (verbose) {
    message("Poole & Rathcke index calculated for ", k_valid, " species with valid peak dates.")
    message("Observed P = ", round(P, 5))
    message("Expected E(P) = ", round(expected_P, 5))
    message("P / E(P) ratio = ", round(ratio, 3))
    message("Interpretation: ", interpretation)
  }
  
  return(list(
    peak_dates = peak_dates,
    ordered_peaks = peak_df,
    distances = D,
    distance_table = distance_table,
    mean_distance = D_bar,
    P = P,
    expected_P = expected_P,
    ratio = ratio,
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    interpretation = interpretation
  ))
}

