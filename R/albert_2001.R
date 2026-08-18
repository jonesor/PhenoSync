#' Albert (2001) Pairwise Flowering Overlap Index
#'
#' Calculates individual-level and population-level flowering synchrony as
#' pairwise overlap in flowering presence/absence. Flower-count values greater
#' than zero are treated as flowering presence.
#'
#' @param data A numeric matrix or data frame in wide format. Rows are census
#'   days, columns are individuals or species, and values are counts of open
#'   flowers or other flowering-intensity measures.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param verbose Logical. If `TRUE`, prints a short progress message.
#'
#' @return A list with `individual_Si`, a named numeric vector of individual
#'   synchrony values, and `population_mean`, the mean synchrony across
#'   individuals.
#'
#' @examples
#' data(example_wide)
#' albert_Si(example_wide, verbose = FALSE)
#'
#' @export
albert_Si <- function(data, census_col = NULL, verbose = TRUE) {
  
  flower_mat <- prepare_census_matrix(data, census_col = census_col, data_arg = "data") #there is no need to convert the table into presence absence matrix because codes takes care ot it. 
  
  n <- ncol(flower_mat) #how many individuals in a population
  ids <- colnames(flower_mat) #identity of each individual
  
  if (n < 2) {
    stop("At least two individuals required.")
  }
  
  Si <- numeric(n) #creating a numeric empty vector to store average synchrony of individual i with all others
  names(Si) <- ids 
  
  for (i in 1:n) {
    
    pair_values <- numeric(n - 1)
    counter <- 1 #start from species pair 1 or keeps track of the next free position in pair_values
    
    for (j in 1:n) {
      
      if (i != j) {    #to skip self comparison
        
        # Intersection (simultaneous flowering)
        aij <- sum(flower_mat[, i] & flower_mat[, j]) #converts data to logical
        
        # Union (at least one flowering)
        bij <- sum(flower_mat[, i] | flower_mat[, j]) #converts data to logical
        
        pair_values[counter] <- ifelse(bij == 0, 0, aij / bij)
        counter <- counter + 1 #move to next species pair
      }
    }
    
    Si[i] <- mean(pair_values)
  }
  
  population_mean <- mean(Si)
  
  if (verbose) {
    message("Albert synchrony calculated for ", n, " individuals.")
  }
  
  return(list(
    individual_Si = Si,
    population_mean = population_mean
  ))
}

