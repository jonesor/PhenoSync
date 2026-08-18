#' Bolmgren (1998) Within-Individual Synchrony Index
#'
#' Calculates within-individual synchrony from flower-level records using mean
#' floral persistence divided by the standard deviation of flower onset dates.
#'
#' @param data A data frame in long format, with one row per flower.
#' @param census_col Optional. Not used by this flower-level long-format
#'   function; included for interface consistency with other PhenoSync indices.
#' @param plant_col Name of the column containing plant or individual IDs.
#' @param onset_col Name of the column containing flower onset day.
#' @param end_col Name of the column containing flower end day.
#' @param inclusive Logical. If `TRUE`, persistence is calculated as
#'   `end_day - onset_day + 1`; otherwise as `end_day - onset_day`.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `individual_results` and `mean_S`.
#'
#' @examples
#' data(example_bolmgren_long)
#' bolmgren_S(example_bolmgren_long, verbose = FALSE)$mean_S
#'
#' @export
bolmgren_S <- function(data,
                       census_col = NULL,
                       plant_col = "plant_id",
                       onset_col = "onset_day",
                       end_col = "end_day",
                       inclusive = TRUE,
                       verbose = TRUE) {
  
  df <- as.data.frame(data)
  
  # STEP 0: Input checks
  required_cols <- c(plant_col, onset_col, end_col)
  
  if (!all(required_cols %in% names(df))) {
    stop("Data must contain the columns: ",
         paste(required_cols, collapse = ", "))
  }
  
  df[[plant_col]] <- as.character(df[[plant_col]])
  df[[onset_col]] <- as.numeric(df[[onset_col]])
  df[[end_col]]   <- as.numeric(df[[end_col]])
  
  # remove rows with missing essential information
  keep <- !(is.na(df[[plant_col]]) | is.na(df[[onset_col]]) | is.na(df[[end_col]]))
  df <- df[keep, , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("No complete flower records found after removing missing values.")
  }
  
  # invalid records
  bad_rows <- which(df[[end_col]] < df[[onset_col]])
  if (length(bad_rows) > 0) {
    stop("Some rows have end day earlier than onset day.")
  }
  
  plants <- unique(df[[plant_col]])
  
  results <- data.frame(
    plant_id = plants,
    n_flowers = NA_integer_,
    mean_persistence = NA_real_,
    SD_onset = NA_real_,
    S = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # STEP 1: Calculate S for each plant

  for (i in seq_along(plants)) {
    
    plant_data <- df[df[[plant_col]] == plants[i], , drop = FALSE]
    n_flowers <- nrow(plant_data)
    
    results$n_flowers[i] <- n_flowers
    
    # Need at least 2 flowers to compute SD(onset)
    if (n_flowers < 2) {
      next
    }
    
    # Persistence per flower
    if (inclusive) {
      persistence <- plant_data[[end_col]] - plant_data[[onset_col]] + 1
    } else {
      persistence <- plant_data[[end_col]] - plant_data[[onset_col]]
    }
    
    P <- mean(persistence, na.rm = TRUE)
    SD_onset <- stats::sd(plant_data[[onset_col]], na.rm = TRUE)
    
    if (is.na(SD_onset)) {
      S_value <- NA_real_
    } else if (SD_onset == 0) {
      S_value <- Inf
    } else {
      S_value <- P / SD_onset
    }
    
    results$mean_persistence[i] <- P
    results$SD_onset[i] <- SD_onset
    results$S[i] <- S_value
  }
  
 # STEP 2: Population summary

  finite_S <- results$S[is.finite(results$S)]
  mean_S <- if (length(finite_S) == 0) NA_real_ else mean(finite_S, na.rm = TRUE)
  
  if (verbose) {
    message("Bolmgren synchrony calculated for ", length(plants), " plants.")
    message("Plants with valid S values: ", sum(!is.na(results$S)))
    message("Mean finite S across plants = ", round(mean_S, 4))
  }
  
  return(list(
    individual_results = results,
    mean_S = mean_S
  ))
}

