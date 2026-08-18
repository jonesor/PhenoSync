#' Schoener-Style Neighbourhood Synchrony Index
#'
#' Calculates the synchrony between each focal individual and the summed
#' flowering profile of its local neighbourhood using a Schoener-style temporal
#' profile overlap index.
#'
#' @param metadata A data frame with one row per individual and columns for IDs
#'   and coordinates.
#' @param flower_mat A numeric matrix or data frame in wide format. Rows are
#'   census days, columns are individuals, and column names must match the IDs in
#'   `metadata`.
#' @param census_col Optional column name or column index containing census
#'   names/dates in wide-format data. When supplied, this column is removed
#'   from the phenology matrix and used as row names. If `NULL`, existing row
#'   names are preserved.
#' @param id_col Name of the metadata column containing individual IDs.
#' @param x_col,y_col Names of projected coordinate columns in meters.
#' @param lon_col,lat_col Names of longitude and latitude columns in decimal
#'   degrees. Used only when `x_col` and `y_col` are not supplied.
#' @param target_crs Optional projected CRS used for lon/lat conversion. If
#'   `NULL`, a UTM CRS is selected automatically.
#' @param neighbor_distance Positive numeric neighbourhood radius in meters.
#' @param grouping_factor Optional metadata column used to restrict neighbours
#'   to the same group as the focal individual.
#' @param exclude_focal Logical. If `TRUE`, the focal individual is excluded
#'   from its own neighbourhood.
#' @param na_as_zero Logical. If `TRUE`, missing flowering values are converted
#'   to zero.
#' @param verbose Logical. If `TRUE`, prints summary messages.
#'
#' @return A list with `summary`, `SI_individual`, `mean_SI`, `distance_matrix`,
#'   `coordinates_xy`, and detailed per-individual neighbourhood results.
#'
#' @details The function can use projected `x`/`y` coordinates directly. If only
#'   longitude/latitude are supplied, the package `sf` is required.
#'
#' @examples
#' data(example_neighbourhood)
#' schoener_SI(
#'   metadata = example_neighbourhood$metadata,
#'   flower_mat = example_neighbourhood$flower_mat,
#'   id_col = "id", x_col = "x", y_col = "y",
#'   neighbor_distance = 50,
#'   verbose = FALSE
#' )$mean_SI
#'
#' @export
schoener_SI <- function(metadata,
                                   flower_mat,
                                   id_col,
                                   census_col = NULL,
                                   x_col = NULL,
                                   y_col = NULL,
                                   lon_col = NULL,
                                   lat_col = NULL,
                                   target_crs = NULL, #Coordinate reference system, default: 4236, with auto-pick UTM zone from mean lon/lat
                                   neighbor_distance, #no default value, must add input for each run
                                   grouping_factor = NULL,
                                   exclude_focal = TRUE, #exclude the focal individual when calculating average synchrony
                                   na_as_zero = TRUE,
                                   verbose = TRUE) {
  
  
  # STEP 0: Input checks
 
  meta <- as.data.frame(metadata)
  flower_mat <- prepare_census_matrix(flower_mat, census_col = census_col, data_arg = "flower_mat")
  
  if (!is.numeric(flower_mat)) {
    stop("flower_mat must be a numeric matrix or data frame.")
  }
  
  if (!is.numeric(neighbor_distance) || length(neighbor_distance) != 1 || is.na(neighbor_distance) || neighbor_distance <= 0) {
    stop("neighbor_distance must be a single positive numeric value in meters.")
  }
  
  if (!(id_col %in% names(meta))) {
    stop("id_col not found in metadata.")
  }
  
  ids_meta <- as.character(meta[[id_col]])
  
  if (anyDuplicated(ids_meta) > 0) {
    stop("Metadata contains duplicate individual IDs.")
  }
  
  if (is.null(colnames(flower_mat))) {
    stop("flower_mat must have column names corresponding to individual IDs.")
  }
  
  ids_flower <- colnames(flower_mat)
  
  if (anyDuplicated(ids_flower) > 0) {
    stop("flower_mat contains duplicate column names (individual IDs).")
  }
  
  if (!all(ids_meta %in% ids_flower)) {
    missing_ids <- ids_meta[!ids_meta %in% ids_flower]
    stop("These metadata IDs are missing from flower_mat columns: ",
         paste(missing_ids, collapse = ", "))
  }
  
  # Reorder flower_mat columns to match metadata order
  flower_mat <- flower_mat[, ids_meta, drop = FALSE]
  
  if (na_as_zero) {
    flower_mat[is.na(flower_mat)] <- 0
  }
  
  if (any(is.na(flower_mat))) {
    stop("flower_mat contains NA values. Use na_as_zero = TRUE or clean the data.")
  }
  
  if (any(flower_mat < 0, na.rm = TRUE)) {
    stop("Negative flowering values are not allowed.")
  }
  
  n <- ncol(flower_mat)
  T <- nrow(flower_mat)
  
  if (n < 2) {
    stop("At least two individuals are required.")
  }
  
 
  # STEP 1: Obtain x,y coordinates in meters
  
  
  use_xy <- !is.null(x_col) && !is.null(y_col)
  use_lonlat <- !is.null(lon_col) && !is.null(lat_col)
  
  if (!use_xy && !use_lonlat) {
    stop("Provide either x_col and y_col, or lon_col and lat_col.")
  }
  
  if (use_xy) {
    if (!(x_col %in% names(meta)) || !(y_col %in% names(meta))) {
      stop("x_col and/or y_col not found in metadata.")
    }
    
    x <- as.numeric(meta[[x_col]])
    y <- as.numeric(meta[[y_col]])
    
    if (any(is.na(x)) || any(is.na(y))) {
      stop("x/y coordinates contain NA values.")
    }
    
  } else {
    if (!(lon_col %in% names(meta)) || !(lat_col %in% names(meta))) {
      stop("lon_col and/or lat_col not found in metadata.")
    }
    
    lon <- as.numeric(meta[[lon_col]])
    lat <- as.numeric(meta[[lat_col]])
    
    if (any(is.na(lon)) || any(is.na(lat))) {
      stop("lon/lat coordinates contain NA values.")
    }
  
    # if sf package is not installed, this error will appear  
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required when using longitude/latitude coordinates.")
    }
    
    pts <- sf::st_as_sf(meta,
                        coords = c(lon_col, lat_col),
                        crs = 4326,
                        remove = FALSE)
    
    # If target_crs not given, auto-pick UTM zone from mean lon/lat
    if (is.null(target_crs)) {
      mean_lon <- mean(lon, na.rm = TRUE)
      mean_lat <- mean(lat, na.rm = TRUE)
      
      utm_zone <- floor((mean_lon + 180) / 6) + 1
      
      if (mean_lat >= 0) {
        target_crs <- paste0("EPSG:", 32600 + utm_zone)  # Northern hemisphere
      } else {
        target_crs <- paste0("EPSG:", 32700 + utm_zone)  # Southern hemisphere
      }
    }
    
    pts_proj <- sf::st_transform(pts, crs = target_crs)
    coords <- sf::st_coordinates(pts_proj)
    
    x <- coords[, 1]
    y <- coords[, 2]
  }
  
 
  # STEP 2: Build Euclidean distance matrix
  
  
  coord_mat <- cbind(x, y)
  rownames(coord_mat) <- ids_meta
  
  dist_mat <- as.matrix(stats::dist(coord_mat, method = "euclidean"))
  rownames(dist_mat) <- ids_meta
  colnames(dist_mat) <- ids_meta
  
 
  # STEP 3: For each focal individual, assign neighborhoods

  
  results_list <- vector("list", n)
  names(results_list) <- ids_meta
  
  SI_values <- numeric(n)
  names(SI_values) <- ids_meta
  SI_values[] <- NA_real_
  
  n_neighbors <- integer(n)
  names(n_neighbors) <- ids_meta
  
  total_focal_flowers <- colSums(flower_mat)
  
  if (any(total_focal_flowers == 0)) {
    warning("Some individuals have zero total flowers and will return NA for synchrony: ",
            paste(ids_meta[total_focal_flowers == 0], collapse = ", "))
  }
  
  for (i in seq_len(n)) {
    
    focal_id <- ids_meta[i]
    
    # Start with all individuals
    eligible <- rep(TRUE, n)
    names(eligible) <- ids_meta
    
    # Exclude focal if requested
    if (exclude_focal) {
      eligible[focal_id] <- FALSE
    }
    
    # Apply grouping filter if requested
    if (!is.null(grouping_factor)) {
      if (!(grouping_factor %in% names(meta))) {
        stop("grouping_factor not found in metadata.")
      }
      
      focal_group <- meta[[grouping_factor]][i]
      eligible <- eligible & (meta[[grouping_factor]] == focal_group)
    }
    
    # Apply distance rule
    within_radius <- dist_mat[focal_id, ] <= neighbor_distance
    
    neighbor_ids <- ids_meta[eligible & within_radius]
    n_neighbors[i] <- length(neighbor_ids)
    
    # Focal flowering curve
    focal_curve <- flower_mat[, focal_id]
    focal_total <- sum(focal_curve)
    
    # If focal never flowers, SI is undefined
    if (focal_total == 0) {
      SI_values[i] <- NA_real_
      
      results_list[[i]] <- list(
        focal_id = focal_id,
        neighbor_ids = neighbor_ids,
        n_neighbors = length(neighbor_ids),
        focal_total = focal_total,
        neighborhood_total = NA_real_,
        SI = NA_real_
      )
      
      next
    }
    
    # If no neighbors, SI cannot be computed
    if (length(neighbor_ids) == 0) {
      SI_values[i] <- NA_real_
      
      results_list[[i]] <- list(
        focal_id = focal_id,
        neighbor_ids = neighbor_ids,
        n_neighbors = 0,
        focal_total = focal_total,
        neighborhood_total = 0,
        SI = NA_real_
      )
      
      next
    }
    
    # Neighborhood curve = sum of all neighbors by census day
    neighborhood_curve <- rowSums(flower_mat[, neighbor_ids, drop = FALSE])
    neighborhood_total <- sum(neighborhood_curve)
    
    # If neighbors never flower, SI undefined
    if (neighborhood_total == 0) {
      SI_values[i] <- NA_real_
      
      results_list[[i]] <- list(
        focal_id = focal_id,
        neighbor_ids = neighbor_ids,
        n_neighbors = length(neighbor_ids),
        focal_total = focal_total,
        neighborhood_total = neighborhood_total,
        SI = NA_real_
      )
      
      next
    }
    
  
    # STEP 4: Compute Schoener/Hall synchrony
  
    
    p_profile <- focal_curve / focal_total
    q_profile <- neighborhood_curve / neighborhood_total
    
    SI_i <- 1 - 0.5 * sum(abs(p_profile - q_profile))
    
    SI_values[i] <- SI_i
    
    results_list[[i]] <- list(
      focal_id = focal_id,
      neighbor_ids = neighbor_ids,
      n_neighbors = length(neighbor_ids),
      focal_total = focal_total,
      neighborhood_total = neighborhood_total,
      SI = SI_i
    )
  }
  

  # STEP 5: Build output table
 
  
  summary_df <- data.frame(
    focal_id = ids_meta,
    n_neighbors = n_neighbors,
    focal_total = total_focal_flowers,
    SI = SI_values,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(grouping_factor)) {
    summary_df[[grouping_factor]] <- meta[[grouping_factor]]
  }
  
  mean_SI <- mean(SI_values, na.rm = TRUE)
  if (is.nan(mean_SI)) mean_SI <- NA_real_
  
  if (verbose) {
    message("Neighborhood synchrony calculated for ", n, " individuals across ", T, " census days.")
    message("Neighbor radius = ", neighbor_distance, " meters.")
    if (!is.null(grouping_factor)) {
      message("Neighbors restricted within grouping factor: ", grouping_factor)
    }
    message("Mean individual synchrony = ", round(mean_SI, 4))
  }
  
  return(list(
    summary = summary_df,
    SI_individual = SI_values,
    mean_SI = mean_SI,
    distance_matrix = dist_mat,
    coordinates_xy = data.frame(id = ids_meta, x = x, y = y),
    details = results_list
  ))
}

