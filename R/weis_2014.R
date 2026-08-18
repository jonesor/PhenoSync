#' Weis et al. (2014) Phenotypic Assortment and Synchrony
#'
#' Calculates the mating matrix, synchrony `S`, and assortative mating
#' correlation `rho` using the flowering-schedule approach of Weis et al. The
#' first column, or the column selected by `trait_col`, contains the phenotypic
#' trait used for assortment; remaining columns contain flowering schedules.
#'
#' @param data A data frame or matrix with individuals in rows. One column is a
#'   numeric trait, usually flowering date, and all other columns are numeric
#'   flowering counts across census dates.
#' @param census_col Optional. Not used by this Weis-style data format
#'   because census days are represented by flowering-schedule column names.
#'   Included for interface consistency with other PhenoSync indices.
#' @param trait_col Numeric column index or character column name identifying the trait column.
#' @param n_boot Non-negative integer number of bootstrap replicates. Use zero
#'   to skip bootstrapping.
#' @param ci_probs Numeric vector of probabilities used for bootstrap intervals.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list containing input metadata, trait vector, flowering matrix,
#'   `M_matrix`, `P_matrix`, `PHI_matrix`, eigenvalues, `S`, `rho`, and optional
#'   bootstrap results.
#'
#' @examples
#' data(example_weiss)
#' weis_sync(example_weiss, n_boot = 0, verbose = FALSE)$rho
#'
#' @export
weis_sync <- function(data,
                      census_col = NULL,
                      trait_col = 1,
                      n_boot = 0,
                      ci_probs = c(0.05, 0.5, 0.95),
                      verbose = TRUE) {
  
  
  # STEP 0: Input checks

  df <- as.data.frame(data)
  
  if (ncol(df) < 2) {
    stop("Input must have at least two columns: one trait column and at least one flowering-schedule column.")
  }
  
  if (length(trait_col) != 1 || is.na(trait_col)) {
    stop("trait_col must identify one valid column in the dataset.")
  }

  if (is.character(trait_col)) {
    if (!(trait_col %in% names(df))) {
      stop("trait_col was supplied as a name but was not found in the dataset.")
    }
    trait_col <- match(trait_col, names(df))
  }

  if (!is.numeric(trait_col) || trait_col != as.integer(trait_col) ||
      trait_col < 1 || trait_col > ncol(df)) {
    stop("trait_col must identify one valid column in the dataset.")
  }
  trait_col <- as.integer(trait_col)
  
  # Reorder so trait is first, schedule columns after
  if (trait_col != 1) {
    df <- df[, c(trait_col, setdiff(seq_len(ncol(df)), trait_col)), drop = FALSE]
  }
  
  trait_name <- names(df)[1]
  flr_names <- names(df)[-1]
  
  # Numeric coercion checks
  df[[1]] <- as.numeric(df[[1]])
  for (j in 2:ncol(df)) {
    df[[j]] <- as.numeric(df[[j]])
  }
  
  if (any(is.na(df[[1]]))) {
    stop("Trait column contains NA or non-numeric values after coercion.")
  }
  
  if (any(is.na(as.matrix(df[, -1, drop = FALSE])))) {
    stop("Flowering-schedule columns contain NA or non-numeric values after coercion.")
  }
  
  flr.matrix <- as.matrix(df[, -1, drop = FALSE])
  mode(flr.matrix) <- "numeric"
  
  if (any(flr.matrix < 0)) {
    stop("Flowering-schedule values cannot be negative.")
  }
  
  L <- nrow(df)                 # number of individuals
  W <- ncol(df)                 # total number of columns
  n_census <- ncol(flr.matrix)  # number of census dates
  
  if (L < 2) {
    stop("At least two individuals are required.")
  }
  
  if (n_census < 1) {
    stop("At least one flowering-schedule column is required.")
  }
  
  total_flowers <- sum(flr.matrix)
  if (total_flowers == 0) {
    stop("All flowering-schedule values are zero. Cannot construct mating matrix.")
  }
  
  row_flower_totals <- rowSums(flr.matrix)
  if (any(row_flower_totals == 0)) {
    warning(
      "Some individuals have zero flowers across all census dates. ",
      "They remain in the trait vector but contribute no mating opportunities."
    )
  }
  
  col_flower_totals <- colSums(flr.matrix)
  if (any(col_flower_totals == 0)) {
    warning(
      "Some census dates have zero flowers across all individuals. ",
      "Division-by-zero protection will be applied in the paternal matrix."
    )
  }
  
  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) || n_boot < 0) {
    stop("n_boot must be a single non-negative integer.")
  }
  n_boot <- as.integer(n_boot)
  
  if (any(ci_probs < 0 | ci_probs > 1)) {
    stop("ci_probs must contain probabilities between 0 and 1.")
  }
  
 
  # Helper function: compute PHI, S, rho
 
  compute_weis <- function(input_df) {
    
    # trait vector, standardized exactly as in supplied code
    f.date <- scale(input_df[, 1])
    
    # flowering schedule matrix
    flr.matrix <- as.matrix(input_df[, -1, drop = FALSE])
    mode(flr.matrix) <- "numeric"
    
    # M matrix: proportion of all open flowers over the season
    M <- flr.matrix / sum(flr.matrix)
    
    # P matrix: proportion of flowers on each census date that were on each donor
    dad.colsumx <- apply(flr.matrix, 2, sum)
    dad.colsum  <- replace(dad.colsumx, dad.colsumx == 0, 1)
    
    # Preserve the original logic from the supplied code
    P <- apply(flr.matrix, 1, function(d) d / dad.colsum)
    
    # Mating matrix
    PHI <- M %*% P
    
    # Eigenvalues
    EV <- eigen(PHI, sym = FALSE)
    ev.vector <- EV$values
    ev.1 <- ev.vector[1]
    ev.sum <- sum(ev.vector)
    S <- ev.1 / ev.sum
    
    # Correlation between mates
    rho <- t(f.date) %*% PHI %*% f.date
    
    list(
      PHI = PHI,
      M = M,
      P = P,
      eigenvalues = ev.vector,
      lambda1 = ev.1,
      eigen_sum = ev.sum,
      S = as.numeric(Re(S)),
      rho = as.numeric(Re(rho))
    )
  }
  
 
  # STEP 1: Observed values
   observed <- compute_weis(df)
  
  
  # STEP 2: Optional bootstrap confidence intervals
   boot_results <- NULL
  
  if (n_boot > 0) {
    index <- seq_len(L)
    SYNC <- numeric(n_boot)
    RHO <- numeric(n_boot)
    
    for (i in seq_len(n_boot)) {
      bootstrapped.index <- sample(index, L, replace = TRUE)
      boot.sample <- df[bootstrapped.index, , drop = FALSE]
      
      tmp <- compute_weis(boot.sample)
      SYNC[i] <- tmp$S
      RHO[i] <- tmp$rho
    }
    
    S.CL <- quantile(SYNC, probs = ci_probs, na.rm = TRUE, type = 7)
    rho.CL <- quantile(RHO, probs = ci_probs, na.rm = TRUE, type = 7)
    
    boot_results <- list(
      n_boot = n_boot,
      ci_probs = ci_probs,
      boot_S = SYNC,
      boot_rho = RHO,
      S_CI = S.CL,
      rho_CI = rho.CL
    )
  }
  
  
  # STEP 3: Output
  if (verbose) {
    message("Weis et al. synchrony calculated for ", L, " individuals.")
    message("Trait column used: ", trait_name)
    message("Number of census dates: ", n_census)
    message("Observed S = ", round(observed$S, 4))
    message("Observed rho = ", round(observed$rho, 4))
    if (n_boot > 0) {
      message("Bootstrap replicates: ", n_boot)
    }
  }
  
  return(list(
    input_info = list(
      n_individuals = L,
      n_columns_total = W,
      n_census_dates = n_census,
      trait_column_name = trait_name,
      flowering_columns = flr_names
    ),
    trait_vector = df[, 1],
    flowering_matrix = flr.matrix,
    M_matrix = observed$M,
    P_matrix = observed$P,
    PHI_matrix = observed$PHI,
    eigenvalues = observed$eigenvalues,
    lambda1 = observed$lambda1,
    eigen_sum = observed$eigen_sum,
    S = observed$S,
    rho = observed$rho,
    bootstrap = boot_results
  ))
}

