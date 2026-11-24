#' Dainou et al. (2012) Inter-Sex Synchrony Index
#'
#' Calculates male synchrony with female flowering by averaging the proportion
#' of females flowering on the days when each male is flowering.
#'
#' @param data A data frame with individual ID and sex columns followed by
#'   numeric census-day columns. Values greater than zero are treated as
#'   flowering presence.
#' @param census_col Optional. Not used by this Dainou-style data format
#'   because census days are represented by the numeric day columns. Included
#'   for interface consistency with other PhenoSync indices.
#' @param id_col Column index or name identifying individual IDs.
#' @param sex_col Column index or name identifying sex.
#' @param female_label Value in `sex_col` used for females.
#' @param male_label Value in `sex_col` used for males.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `male_dj`, `female_prop_by_day`, and `population_S`.
#'
#' @examples
#' data(example_dainou)
#' dainou_SI(example_dainou, verbose = FALSE)$population_S
#'
#' @export
dainou_SI <- function(data,
                        census_col = NULL,
                        id_col = 1,
                        sex_col = 2,
                        female_label = "F",
                        male_label = "M",
                        verbose = TRUE) {
  
  df <- as.data.frame(data)
  
  # Work out which columns hold the IDs and sexes. id_col and sex_col can be
  # given as column numbers or column names.
  find_col <- function(col, arg_name) {
    if (length(col) != 1 || is.na(col)) {
      stop(arg_name, " must be a single column number or column name.")
    }
    if (is.character(col)) {
      pos <- match(col, names(df))
      if (is.na(pos)) {
        stop(arg_name, " '", col, "' is not a column in data.")
      }
    } else if (is.numeric(col) && col == round(col) && col >= 1 && col <= ncol(df)) {
      pos <- as.integer(col)
    } else {
      stop(arg_name, " must be a column number between 1 and ", ncol(df),
           ", or a column name.")
    }
    pos
  }
  
  id_pos  <- find_col(id_col, "id_col")
  sex_pos <- find_col(sex_col, "sex_col")
  
  if (id_pos == sex_pos) {
    stop("id_col and sex_col must be different columns.")
  }
  
  # Extract metadata
  ids  <- df[[id_pos]]
  sex  <- df[[sex_pos]]
  
  # Extract census matrix
  flower_mat <- as.matrix(df[, -c(id_pos, sex_pos), drop = FALSE])
  
  if (!is.numeric(flower_mat)) {
    stop("Census-day columns must be numeric.")
  }

  check_flowering_values(flower_mat, "data")

  if (anyDuplicated(ids) > 0) {
    stop("The ID column has repeated IDs: ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "),
         ". Each individual needs its own row.")
  }
  
  # Presence / absence
  presence_mat <- (flower_mat > 0) * 1
  
  rownames(presence_mat) <- ids
  
  # Split males and females
  male_idx <- which(sex == male_label)
  female_idx <- which(sex == female_label)
  
  if (length(male_idx) == 0) {
    stop("No males found in the dataset.")
  }
  
  if (length(female_idx) == 0) {
    stop("No females found in the dataset.")
  }
  
  male_ids <- ids[male_idx]
  female_ids <- ids[female_idx]
  
  male_mat <- presence_mat[male_idx, , drop = FALSE]
  female_mat <- presence_mat[female_idx, , drop = FALSE]
  
  Nm <- nrow(male_mat) #calculating number of males in a population
  
  # Proportion of females flowering at each census day
  female_prop_by_day <- colMeans(female_mat)
  
  # Store d_j for each male
  dj <- numeric(Nm)
  names(dj) <- male_ids
  
  for (j in 1:Nm) {
    
    male_flowering_days <- which(male_mat[j, ] == 1)
    Tj <- length(male_flowering_days)
    
    if (Tj == 0) {
      dj[j] <- NA
    } else {
      dj[j] <- mean(female_prop_by_day[male_flowering_days])
    }
  }
  
  # Population synchrony S
  S <- mean(dj, na.rm = TRUE)
  
  if (verbose) {
    message("Dainou inter-sex synchrony calculated for ", Nm, " males and ", nrow(female_mat), " females.")
  }
  
  return(list(
    male_dj = dj,
    female_prop_by_day = female_prop_by_day,
    population_S = S
  ))
}

