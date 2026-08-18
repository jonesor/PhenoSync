# Internal helpers ---------------------------------------------------------

.pheno_clamp <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

.pheno_lognormal_parameters_from_cv <- function(median_value, cv) {
  if (length(median_value) != 1L || is.na(median_value) || median_value <= 0) {
    stop("median_value must be a single number > 0.", call. = FALSE)
  }
  if (length(cv) != 1L || is.na(cv) || cv < 0) {
    stop("cv must be a single number >= 0.", call. = FALSE)
  }

  sigma <- sqrt(log(cv^2 + 1))
  mu <- log(median_value)
  list(mu = mu, sigma = sigma)
}

.pheno_deterministic_lognormal <- function(n, median_value, cv) {
  if (n <= 0) {
    return(numeric(0))
  }

  pars <- .pheno_lognormal_parameters_from_cv(median_value, cv)
  if (pars$sigma == 0) {
    return(rep(median_value, n))
  }

  probabilities <- (seq_len(n) - 0.5) / n
  stats::qlnorm(probabilities, meanlog = pars$mu, sdlog = pars$sigma)
}

.pheno_allocate_categories <- function(n, probabilities, labels, seed = NULL) {
  if (length(probabilities) != length(labels)) {
    stop("probabilities and labels must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(probabilities)) || any(probabilities < 0)) {
    stop("probabilities must contain finite, non-negative values.", call. = FALSE)
  }
  if (abs(sum(probabilities) - 1) > 1e-8) {
    stop("probabilities must sum to 1.", call. = FALSE)
  }

  expected <- n * probabilities
  counts <- floor(expected)
  remainder <- n - sum(counts)

  if (remainder > 0) {
    add_to <- order(expected - counts, decreasing = TRUE)[seq_len(remainder)]
    counts[add_to] <- counts[add_to] + 1L
  }

  result <- rep(labels, counts)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample(result, size = length(result), replace = FALSE)
}

.pheno_rank_abundance_sizes <- function(n_species,
                                        min_individuals = 20,
                                        max_individuals = 500,
                                        alpha = 0.9) {
  if (n_species < 2) {
    stop("n_species must be at least 2.", call. = FALSE)
  }
  if (min_individuals < 1 || max_individuals < min_individuals) {
    stop("Check min_individuals and max_individuals.", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0) {
    stop("rank_abundance_alpha must be > 0.", call. = FALSE)
  }

  species_rank <- seq_len(n_species)
  rank_weight <- species_rank^(-alpha)
  scaled_weight <- (rank_weight - min(rank_weight)) /
    (max(rank_weight) - min(rank_weight))

  population_size <- round(
    min_individuals +
      (max_individuals - min_individuals) * scaled_weight
  )

  population_size[1] <- max_individuals
  population_size[n_species] <- min_individuals
  as.integer(population_size)
}

.pheno_assign_duration_values <- function(n_species,
                                          probabilities,
                                          duration_ranges,
                                          seed = NULL) {
  labels <- names(probabilities)

  if (is.null(labels) || any(labels == "")) {
    stop("duration_probabilities must be a named vector.", call. = FALSE)
  }
  if (!all(labels %in% names(duration_ranges))) {
    stop("Each duration class must have a corresponding duration range.", call. = FALSE)
  }

  duration_class <- .pheno_allocate_categories(
    n = n_species,
    probabilities = as.numeric(probabilities),
    labels = labels,
    seed = seed
  )

  duration_days <- numeric(n_species)

  for (class_name in labels) {
    indices <- which(duration_class == class_name)
    n_in_class <- length(indices)
    if (n_in_class == 0) {
      next
    }

    class_range <- duration_ranges[[class_name]]
    if (length(class_range) != 2L || any(!is.finite(class_range)) ||
        class_range[1] <= 0 || class_range[2] < class_range[1]) {
      stop(
        "Each duration range must contain a positive minimum and maximum.",
        call. = FALSE
      )
    }

    values <- seq(
      from = class_range[1],
      to = class_range[2],
      length.out = n_in_class
    )

    if (!is.null(seed)) {
      set.seed(seed + match(class_name, labels) * 101)
    }

    duration_days[indices] <- sample(
      round(values),
      size = n_in_class,
      replace = FALSE
    )
  }

  data.frame(
    duration_class = duration_class,
    duration_days = as.integer(duration_days),
    stringsAsFactors = FALSE
  )
}

.pheno_community_presets <- function() {
  default_duration_ranges <- list(
    short = c(10, 20),
    medium = c(25, 45),
    long = c(50, 80)
  )

  staggered_duration_ranges <- list(
    short = c(5, 10),
    medium = c(12, 20),
    long = c(25, 35)
  )

  pulse_duration_ranges <- list(
    short = c(4, 8),
    medium = c(9, 14),
    long = c(15, 22)
  )

  list(
    even_weak_seasonal = list(
      description = paste(
        "Weak dominance, broadly distributed flowering peaks,",
        "mostly symmetric flowering curves."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.15,
      peak_pattern = "broad",
      peak_beta_shape = c(1.4, 1.4),
      peak_margin = 12,
      duration_probabilities = c(short = 0.20, medium = 0.55, long = 0.25),
      duration_ranges = default_duration_ranges,
      shape_probabilities = c(symmetric = 0.80, skewed = 0.15, bimodal = 0.05),
      peak_sd_ratio = 0.35,
      minimum_peak_sd = 2,
      maximum_peak_sd = 25,
      species_flower_median = 15,
      species_flower_cv = 0.50,
      individual_flower_cv = 0.45,
      individual_duration_cv = 0.15,
      longevity_probabilities = c(`1` = 0.30, `2` = 0.40, `3` = 0.20, `5` = 0.10),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    ),

    strongly_seasonal_even = list(
      description = paste(
        "Weak abundance dominance but strong community-level seasonality,",
        "with species peaks concentrated near the middle of the season."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.20,
      peak_pattern = "central",
      peak_beta_shape = c(10, 10),
      peak_margin = 12,
      duration_probabilities = c(short = 0.50, medium = 0.40, long = 0.10),
      duration_ranges = default_duration_ranges,
      shape_probabilities = c(symmetric = 0.75, skewed = 0.20, bimodal = 0.05),
      peak_sd_ratio = 0.25,
      minimum_peak_sd = 1,
      maximum_peak_sd = 18,
      species_flower_median = 18,
      species_flower_cv = 0.55,
      individual_flower_cv = 0.45,
      individual_duration_cv = 0.15,
      longevity_probabilities = c(`1` = 0.35, `2` = 0.40, `3` = 0.20, `5` = 0.05),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    ),

    seasonal_herbaceous = list(
      description = paste(
        "Moderate dominance, one concentrated flowering season,",
        "mostly short- to medium-duration species."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.90,
      peak_pattern = "central",
      peak_beta_shape = c(4, 4),
      peak_margin = 12,
      duration_probabilities = c(short = 0.35, medium = 0.50, long = 0.15),
      duration_ranges = default_duration_ranges,
      shape_probabilities = c(symmetric = 0.55, skewed = 0.35, bimodal = 0.10),
      peak_sd_ratio = 0.30,
      minimum_peak_sd = 2,
      maximum_peak_sd = 22,
      species_flower_median = 20,
      species_flower_cv = 0.80,
      individual_flower_cv = 0.55,
      individual_duration_cv = 0.18,
      longevity_probabilities = c(`1` = 0.25, `2` = 0.40, `3` = 0.25, `5` = 0.10),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    ),

    phenologically_staggered = list(
      description = paste(
        "Low-to-moderate abundance dominance, with species peaks spread",
        "approximately evenly through the season and relatively short",
        "flowering schedules that promote temporal turnover."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.40,
      peak_pattern = "staggered",
      peak_margin = 12,
      duration_probabilities = c(short = 0.60, medium = 0.35, long = 0.05),
      duration_ranges = staggered_duration_ranges,
      shape_probabilities = c(symmetric = 0.85, skewed = 0.10, bimodal = 0.05),
      peak_sd_ratio = 0.18,
      minimum_peak_sd = 1,
      maximum_peak_sd = 10,
      species_flower_median = 15,
      species_flower_cv = 0.60,
      individual_flower_cv = 0.45,
      individual_duration_cv = 0.12,
      longevity_probabilities = c(`1` = 0.40, `2` = 0.40, `3` = 0.15, `5` = 0.05),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    ),

    bimodal_seasonal = list(
      description = paste(
        "Moderate abundance inequality and two distinct community-level",
        "flowering seasons centered in the early and late parts of",
        "the season."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.60,
      peak_pattern = "bimodal_seasonal",
      peak_margin = 12,
      first_peak_fraction = 0.30,
      second_peak_fraction = 0.70,
      bimodal_cluster_sd_fraction = 0.045,
      first_cluster_species_fraction = 0.50,
      duration_probabilities = c(short = 0.35, medium = 0.50, long = 0.15),
      duration_ranges = default_duration_ranges,
      shape_probabilities = c(symmetric = 0.70, skewed = 0.20, bimodal = 0.10),
      peak_sd_ratio = 0.28,
      minimum_peak_sd = 2,
      maximum_peak_sd = 20,
      species_flower_median = 17,
      species_flower_cv = 0.70,
      individual_flower_cv = 0.50,
      individual_duration_cv = 0.17,
      longevity_probabilities = c(`1` = 0.30, `2` = 0.40, `3` = 0.20, `5` = 0.10),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    ),

    mass_flowering_dominated = list(
      description = paste(
        "Strong dominance, clustered peaks, short/pulsed schedules,",
        "with the top 20% of species contributing about 80% of",
        "seasonal floral output."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 1.50,
      peak_pattern = "clustered",
      peak_margin = 12,
      dominant_peak_fraction = 0.50,
      dominant_peak_sd_fraction = 0.035,
      early_peak_fraction = 0.30,
      late_peak_fraction = 0.70,
      secondary_peak_sd_fraction = 0.055,
      duration_probabilities = c(short = 0.60, medium = 0.30, long = 0.10),
      duration_ranges = default_duration_ranges,
      shape_probabilities = c(symmetric = 0.50, skewed = 0.20, bimodal = 0.30),
      peak_sd_ratio = 0.25,
      dominant_peak_sd_ratio = 0.12,
      minimum_peak_sd = 1,
      maximum_peak_sd = 18,
      species_flower_median = 12,
      species_flower_cv = 0.85,
      individual_flower_cv = 0.60,
      individual_duration_cv = 0.15,
      longevity_probabilities = c(`1` = 0.45, `2` = 0.35, `3` = 0.15, `5` = 0.05),
      top_fraction = 0.20,
      target_top_output_share = 0.80
    ),

    ephemeral_pulse_flowering = list(
      description = paste(
        "A brief, strongly synchronized community-wide flowering pulse,",
        "with predominantly short species and individual flowering",
        "schedules but without enforced floral-output dominance."
      ),
      min_individuals = 20,
      max_individuals = 500,
      rank_abundance_alpha = 0.70,
      peak_pattern = "pulse",
      peak_margin = 12,
      pulse_peak_fraction = 0.50,
      pulse_peak_sd_fraction = 0.020,
      duration_probabilities = c(short = 0.75, medium = 0.20, long = 0.05),
      duration_ranges = pulse_duration_ranges,
      shape_probabilities = c(symmetric = 0.80, skewed = 0.15, bimodal = 0.05),
      peak_sd_ratio = 0.12,
      minimum_peak_sd = 1,
      maximum_peak_sd = 8,
      species_flower_median = 20,
      species_flower_cv = 0.65,
      individual_flower_cv = 0.50,
      individual_duration_cv = 0.10,
      longevity_probabilities = c(`1` = 0.60, `2` = 0.30, `3` = 0.08, `5` = 0.02),
      top_fraction = 0.20,
      target_top_output_share = NA_real_
    )
  )
}

.pheno_assign_species_peaks <- function(n_species,
                                        season_length,
                                        configuration,
                                        dominant_species_indices,
                                        seed) {
  peak_pattern <- configuration$peak_pattern
  peak_margin <- configuration$peak_margin
  minimum_peak <- 1 + peak_margin
  maximum_peak <- season_length - peak_margin

  if (minimum_peak >= maximum_peak) {
    stop(
      "season_length is too short for the selected peak_margin. ",
      "Increase season_length or reduce peak_margin.",
      call. = FALSE
    )
  }

  if (peak_pattern %in% c("broad", "central")) {
    probabilities <- (seq_len(n_species) - 0.5) / n_species
    beta_shape <- configuration$peak_beta_shape
    peak_dates <- minimum_peak +
      (maximum_peak - minimum_peak) *
      stats::qbeta(probabilities, shape1 = beta_shape[1], shape2 = beta_shape[2])

    set.seed(seed)
    peak_dates <- sample(round(peak_dates), n_species, replace = FALSE)

  } else if (peak_pattern == "staggered") {
    peak_dates <- seq(
      from = minimum_peak,
      to = maximum_peak,
      length.out = n_species
    )

    set.seed(seed)
    peak_dates <- sample(round(peak_dates), n_species, replace = FALSE)

  } else if (peak_pattern == "bimodal_seasonal") {
    first_fraction <- configuration$first_cluster_species_fraction
    if (first_fraction <= 0 || first_fraction >= 1) {
      stop("first_cluster_species_fraction must lie between 0 and 1.", call. = FALSE)
    }

    n_first <- as.integer(round(n_species * first_fraction))
    n_first <- max(1L, min(n_species - 1L, n_first))

    set.seed(seed + 1)
    shuffled_ids <- sample(seq_len(n_species), n_species, replace = FALSE)
    first_ids <- shuffled_ids[seq_len(n_first)]
    second_ids <- shuffled_ids[(n_first + 1L):n_species]
    peak_dates <- numeric(n_species)

    p_first <- (seq_along(first_ids) - 0.5) / length(first_ids)
    peak_dates[first_ids] <- round(
      season_length * configuration$first_peak_fraction +
        stats::qnorm(p_first) *
        season_length * configuration$bimodal_cluster_sd_fraction
    )

    p_second <- (seq_along(second_ids) - 0.5) / length(second_ids)
    peak_dates[second_ids] <- round(
      season_length * configuration$second_peak_fraction +
        stats::qnorm(p_second) *
        season_length * configuration$bimodal_cluster_sd_fraction
    )

  } else if (peak_pattern == "pulse") {
    probabilities <- (seq_len(n_species) - 0.5) / n_species
    peak_dates <- round(
      season_length * configuration$pulse_peak_fraction +
        stats::qnorm(probabilities) *
        season_length * configuration$pulse_peak_sd_fraction
    )

    set.seed(seed)
    peak_dates <- sample(peak_dates, n_species, replace = FALSE)

  } else if (peak_pattern == "clustered") {
    peak_dates <- numeric(n_species)
    dominant_ids <- dominant_species_indices
    other_ids <- setdiff(seq_len(n_species), dominant_ids)

    p_dominant <- (seq_along(dominant_ids) - 0.5) / length(dominant_ids)
    peak_dates[dominant_ids] <- round(
      season_length * configuration$dominant_peak_fraction +
        stats::qnorm(p_dominant) *
        season_length * configuration$dominant_peak_sd_fraction
    )

    if (length(other_ids) > 0L) {
      set.seed(seed + 1)
      other_ids <- sample(other_ids, length(other_ids), replace = FALSE)
      n_early <- ceiling(length(other_ids) / 2)
      early_ids <- other_ids[seq_len(n_early)]
      late_ids <- setdiff(other_ids, early_ids)

      if (length(early_ids) > 0L) {
        p_early <- (seq_along(early_ids) - 0.5) / length(early_ids)
        peak_dates[early_ids] <- round(
          season_length * configuration$early_peak_fraction +
            stats::qnorm(p_early) *
            season_length * configuration$secondary_peak_sd_fraction
        )
      }

      if (length(late_ids) > 0L) {
        p_late <- (seq_along(late_ids) - 0.5) / length(late_ids)
        peak_dates[late_ids] <- round(
          season_length * configuration$late_peak_fraction +
            stats::qnorm(p_late) *
            season_length * configuration$secondary_peak_sd_fraction
        )
      }
    }

  } else {
    stop("Unknown peak_pattern: ", peak_pattern, call. = FALSE)
  }

  as.integer(.pheno_clamp(peak_dates, minimum_peak, maximum_peak))
}

.pheno_enforce_top_output_share <- function(population_size,
                                            mean_flowers_per_individual,
                                            top_species_indices,
                                            target_share) {
  if (is.na(target_share)) {
    return(mean_flowers_per_individual)
  }

  if (target_share <= 0 || target_share >= 1) {
    stop("target_top_output_share must lie between 0 and 1.", call. = FALSE)
  }

  expected_output <- population_size * mean_flowers_per_individual
  top_output <- sum(expected_output[top_species_indices])
  lower_output <- sum(expected_output[-top_species_indices])

  if (top_output <= 0 || lower_output <= 0) {
    stop("Cannot enforce target output share because expected output is zero.", call. = FALSE)
  }

  multiplier <- (target_share * lower_output) /
    ((1 - target_share) * top_output)

  adjusted <- mean_flowers_per_individual
  adjusted[top_species_indices] <-
    adjusted[top_species_indices] * multiplier
  adjusted
}

.pheno_build_species_parameters <- function(n_species,
                                            season_length,
                                            configuration,
                                            seed) {
  species_id <- sprintf("sp%03d", seq_len(n_species))
  species_rank <- seq_len(n_species)

  population_size <- .pheno_rank_abundance_sizes(
    n_species = n_species,
    min_individuals = configuration$min_individuals,
    max_individuals = configuration$max_individuals,
    alpha = configuration$rank_abundance_alpha
  )

  n_dominant <- max(1L, ceiling(n_species * configuration$top_fraction))
  if (n_dominant >= n_species && !is.na(configuration$target_top_output_share)) {
    stop(
      "top_fraction must leave at least one non-dominant species when ",
      "target_top_output_share is enforced.",
      call. = FALSE
    )
  }
  dominant_species_indices <- seq_len(n_dominant)

  peak_date <- .pheno_assign_species_peaks(
    n_species = n_species,
    season_length = season_length,
    configuration = configuration,
    dominant_species_indices = dominant_species_indices,
    seed = seed + 100
  )

  duration_info <- .pheno_assign_duration_values(
    n_species = n_species,
    probabilities = configuration$duration_probabilities,
    duration_ranges = configuration$duration_ranges,
    seed = seed + 200
  )

  if (!is.na(configuration$target_top_output_share)) {
    short_range <- configuration$duration_ranges$short
    dominant_durations <- round(seq(
      short_range[1],
      short_range[2],
      length.out = n_dominant
    ))
    duration_info$duration_class[dominant_species_indices] <- "short"
    duration_info$duration_days[dominant_species_indices] <- dominant_durations
  }

  curve_shape <- .pheno_allocate_categories(
    n = n_species,
    probabilities = as.numeric(configuration$shape_probabilities),
    labels = names(configuration$shape_probabilities),
    seed = seed + 300
  )

  skew_direction <- rep(NA_character_, n_species)
  skewed_indices <- which(curve_shape == "skewed")
  if (length(skewed_indices) > 0) {
    skew_direction[skewed_indices] <- rep(
      c("right_tail", "left_tail"),
      length.out = length(skewed_indices)
    )
  }

  species_flower_values <- .pheno_deterministic_lognormal(
    n = n_species,
    median_value = configuration$species_flower_median,
    cv = configuration$species_flower_cv
  )

  set.seed(seed + 400)
  mean_flowers_per_individual <- sample(
    species_flower_values,
    n_species,
    replace = FALSE
  )

  mean_flowers_per_individual <- .pheno_enforce_top_output_share(
    population_size = population_size,
    mean_flowers_per_individual = mean_flowers_per_individual,
    top_species_indices = dominant_species_indices,
    target_share = configuration$target_top_output_share
  )

  peak_sd_ratio <- rep(configuration$peak_sd_ratio, n_species)
  if (!is.null(configuration$dominant_peak_sd_ratio)) {
    peak_sd_ratio[dominant_species_indices] <-
      configuration$dominant_peak_sd_ratio
  }

  individual_peak_sd <- .pheno_clamp(
    duration_info$duration_days * peak_sd_ratio,
    configuration$minimum_peak_sd,
    configuration$maximum_peak_sd
  )

  longevity_values <- as.integer(names(configuration$longevity_probabilities))
  floral_longevity <- as.integer(
    .pheno_allocate_categories(
      n = n_species,
      probabilities = as.numeric(configuration$longevity_probabilities),
      labels = as.character(longevity_values),
      seed = seed + 500
    )
  )

  expected_species_output <- population_size * mean_flowers_per_individual

  data.frame(
    species_id = species_id,
    species_rank = species_rank,
    population_size = population_size,
    dominant_species = species_rank <= n_dominant,
    species_peak_day = peak_date,
    duration_class = duration_info$duration_class,
    species_duration_days = duration_info$duration_days,
    individual_peak_sd_days = individual_peak_sd,
    curve_shape = curve_shape,
    skew_direction = skew_direction,
    mean_flowers_per_individual = mean_flowers_per_individual,
    expected_species_output = expected_species_output,
    floral_longevity_days = floral_longevity,
    stringsAsFactors = FALSE
  )
}

.pheno_build_individual_parameters <- function(species_parameters,
                                               season_length,
                                               configuration,
                                               seed) {
  total_individuals <- sum(species_parameters$population_size)

  individual_id <- character(total_individuals)
  species_id <- character(total_individuals)
  species_rank <- integer(total_individuals)
  individual_peak_day <- integer(total_individuals)
  individual_duration_days <- integer(total_individuals)
  total_flowers <- integer(total_individuals)
  curve_shape <- character(total_individuals)
  skew_direction <- character(total_individuals)
  floral_longevity_days <- integer(total_individuals)

  current_index <- 1L

  for (species_row in seq_len(nrow(species_parameters))) {
    species <- species_parameters[species_row, ]
    n_individuals <- species$population_size
    indices <- current_index:(current_index + n_individuals - 1L)

    individual_id[indices] <- sprintf(
      "%s_ind%04d",
      species$species_id,
      seq_len(n_individuals)
    )
    species_id[indices] <- species$species_id
    species_rank[indices] <- species$species_rank

    peak_quantiles <- stats::qnorm(
      (seq_len(n_individuals) - 0.5) / n_individuals
    )
    set.seed(seed + species_row * 11)
    peak_quantiles <- sample(peak_quantiles, n_individuals, replace = FALSE)

    individual_peak_day[indices] <- round(
      species$species_peak_day +
        peak_quantiles * species$individual_peak_sd_days
    )
    individual_peak_day[indices] <- as.integer(
      .pheno_clamp(individual_peak_day[indices], 1, season_length)
    )

    duration_multipliers <- .pheno_deterministic_lognormal(
      n = n_individuals,
      median_value = 1,
      cv = configuration$individual_duration_cv
    )
    set.seed(seed + species_row * 13)
    duration_multipliers <- sample(
      duration_multipliers,
      n_individuals,
      replace = FALSE
    )
    individual_duration_days[indices] <- as.integer(
      .pheno_clamp(
        round(species$species_duration_days * duration_multipliers),
        3,
        season_length
      )
    )

    flower_multipliers <- .pheno_deterministic_lognormal(
      n = n_individuals,
      median_value = 1,
      cv = configuration$individual_flower_cv
    )
    set.seed(seed + species_row * 17)
    flower_multipliers <- sample(
      flower_multipliers,
      n_individuals,
      replace = FALSE
    )
    total_flowers[indices] <- pmax(
      1L,
      as.integer(
        round(species$mean_flowers_per_individual * flower_multipliers)
      )
    )

    curve_shape[indices] <- species$curve_shape
    skew_direction[indices] <- species$skew_direction
    floral_longevity_days[indices] <- species$floral_longevity_days
    current_index <- current_index + n_individuals
  }

  data.frame(
    individual_id = individual_id,
    species_id = species_id,
    species_rank = species_rank,
    individual_peak_day = individual_peak_day,
    individual_duration_days = individual_duration_days,
    total_flowers = total_flowers,
    curve_shape = curve_shape,
    skew_direction = skew_direction,
    floral_longevity_days = floral_longevity_days,
    stringsAsFactors = FALSE
  )
}

.pheno_make_flowering_weights <- function(days,
                                           peak_day,
                                           duration_days,
                                           curve_shape = "symmetric",
                                           skew_direction = NA_character_) {
  duration_days <- max(3L, as.integer(round(duration_days)))
  peak_day <- as.integer(round(peak_day))
  weights <- numeric(length(days))

  if (curve_shape == "symmetric") {
    left_span <- floor((duration_days - 1) / 2)
    right_span <- duration_days - 1 - left_span
    active <- days >= (peak_day - left_span) &
      days <= (peak_day + right_span)
    sigma <- max(duration_days / 5, 0.8)
    weights[active] <- exp(
      -0.5 * ((days[active] - peak_day) / sigma)^2
    )

  } else if (curve_shape == "skewed") {
    if (is.na(skew_direction)) {
      skew_direction <- "right_tail"
    }

    if (skew_direction == "right_tail") {
      left_span <- max(1L, round(duration_days / 3))
      right_span <- duration_days - 1 - left_span
    } else {
      right_span <- max(1L, round(duration_days / 3))
      left_span <- duration_days - 1 - right_span
    }

    active <- days >= (peak_day - left_span) &
      days <= (peak_day + right_span)
    left_sigma <- max(left_span / 2.2, 0.8)
    right_sigma <- max(right_span / 2.2, 0.8)
    left_days <- active & days <= peak_day
    right_days <- active & days > peak_day

    weights[left_days] <- exp(
      -0.5 * ((days[left_days] - peak_day) / left_sigma)^2
    )
    weights[right_days] <- exp(
      -0.5 * ((days[right_days] - peak_day) / right_sigma)^2
    )

  } else if (curve_shape == "bimodal") {
    left_span <- floor((duration_days - 1) / 2)
    right_span <- duration_days - 1 - left_span
    active <- days >= (peak_day - left_span) &
      days <= (peak_day + right_span)
    peak_gap <- max(4, round(duration_days * 0.35))
    first_peak <- peak_day - peak_gap / 2
    second_peak <- peak_day + peak_gap / 2
    sigma <- max(1, duration_days * 0.10)

    weights[active] <-
      0.55 * stats::dnorm(days[active], mean = first_peak, sd = sigma) +
      0.45 * stats::dnorm(days[active], mean = second_peak, sd = sigma)

  } else {
    stop("Unknown curve_shape: ", curve_shape, call. = FALSE)
  }

  if (sum(weights) <= 0) {
    weights[which.min(abs(days - peak_day))] <- 1
  }
  weights / sum(weights)
}

.pheno_generate_new_flower_matrix <- function(individual_parameters,
                                               season_length,
                                               seed = 1) {
  n_individuals <- nrow(individual_parameters)
  days <- seq_len(season_length)
  new_flowers <- matrix(
    0L,
    nrow = season_length,
    ncol = n_individuals
  )
  rownames(new_flowers) <- as.character(days)
  colnames(new_flowers) <- individual_parameters$individual_id

  set.seed(seed)
  for (individual_index in seq_len(n_individuals)) {
    individual <- individual_parameters[individual_index, ]
    weights <- .pheno_make_flowering_weights(
      days = days,
      peak_day = individual$individual_peak_day,
      duration_days = individual$individual_duration_days,
      curve_shape = individual$curve_shape,
      skew_direction = individual$skew_direction
    )

    new_flowers[, individual_index] <- as.integer(
      stats::rmultinom(
        n = 1,
        size = individual$total_flowers,
        prob = weights
      )[, 1]
    )
  }

  new_flowers
}

.pheno_apply_floral_longevity <- function(new_flower_matrix,
                                          longevity_days) {
  if (ncol(new_flower_matrix) != length(longevity_days)) {
    stop(
      "The length of longevity_days must equal the number of columns ",
      "in new_flower_matrix.",
      call. = FALSE
    )
  }
  if (!is.integer(new_flower_matrix) || any(new_flower_matrix < 0L)) {
    stop(
      "new_flower_matrix must contain non-negative integer flower counts.",
      call. = FALSE
    )
  }

  open_flowers <- matrix(
    0L,
    nrow = nrow(new_flower_matrix),
    ncol = ncol(new_flower_matrix),
    dimnames = dimnames(new_flower_matrix)
  )

  unique_longevities <- sort(unique(longevity_days))
  for (longevity in unique_longevities) {
    column_indices <- which(longevity_days == longevity)
    block <- new_flower_matrix[, column_indices, drop = FALSE]
    accumulated <- block

    if (longevity > 1L && nrow(block) > 1L) {
      maximum_lag <- min(longevity - 1L, nrow(block) - 1L)
      for (lag in seq_len(maximum_lag)) {
        shifted <- rbind(
          matrix(0L, nrow = lag, ncol = ncol(block)),
          block[seq_len(nrow(block) - lag), , drop = FALSE]
        )
        accumulated <- accumulated + shifted
      }
    }

    storage.mode(accumulated) <- "integer"
    open_flowers[, column_indices] <- accumulated
  }

  storage.mode(open_flowers) <- "integer"
  open_flowers
}

.pheno_aggregate_individual_matrix_by_species <- function(individual_matrix,
                                                          individual_species) {
  if (!is.integer(individual_matrix) || any(individual_matrix < 0L)) {
    stop(
      "individual_matrix must contain non-negative integer flower counts.",
      call. = FALSE
    )
  }

  aggregated <- rowsum(
    x = t(individual_matrix),
    group = individual_species,
    reorder = FALSE
  )
  aggregated <- t(aggregated)
  storage.mode(aggregated) <- "integer"
  aggregated
}

.pheno_check_configuration <- function(configuration) {
  if (configuration$top_fraction <= 0 || configuration$top_fraction > 1) {
    stop("top_fraction must lie in (0, 1].", call. = FALSE)
  }
  if (configuration$peak_margin < 0) {
    stop("peak_margin must be >= 0.", call. = FALSE)
  }

  invisible(.pheno_allocate_categories(
    n = 10,
    probabilities = as.numeric(configuration$duration_probabilities),
    labels = names(configuration$duration_probabilities)
  ))
  invisible(.pheno_allocate_categories(
    n = 10,
    probabilities = as.numeric(configuration$shape_probabilities),
    labels = names(configuration$shape_probabilities)
  ))
  invisible(.pheno_allocate_categories(
    n = 10,
    probabilities = as.numeric(configuration$longevity_probabilities),
    labels = names(configuration$longevity_probabilities)
  ))

  invisible(TRUE)
}

# Public functions --------------------------------------------------------

#' Generate a true phenological community
#'
#' @description
#' Generates daily flowering phenologies for a synthetic plant community. Seven
#' plant community archetypes are provided to model the phenologies. These
#' archetypes vary in their species abundance, dominance, seasonal
#' concentration, temporal turnover, multimodality, flowering duration, daily
#' flower production, within-species variation, and floral longevity. The
#' function produces individual-level and species-level daily flower-count
#' matrices, together with presence/absence matrices. Biological parameters used
#' to generate them are also provided either as default values or specified by
#' the user.
#'
#' Reproducibility: the same seed number and the same combination of parameter
#' values produce the same simulated community.
#'
#' @param community_type Character string selecting one of the seven preset
#'   archetypes: `"even_weak_seasonal"`, `"strongly_seasonal_even"`,
#'   `"seasonal_herbaceous"`, `"phenologically_staggered"`,
#'   `"bimodal_seasonal"`, `"mass_flowering_dominated"`, or
#'   `"ephemeral_pulse_flowering"`. The default is
#'   `"even_weak_seasonal"`.
#' @param n_species Integer number of species to generate. Default is `50`.
#'   Must be at least 2.
#' @param season_length Integer number of consecutive daily censuses in the
#'   simulated season. Default is `90` days.
#' @param seed Integer random seed controlling species assignments,
#'   individual-level permutations, and daily multinomial flower allocation.
#'   Default is `123`.
#' @param control Named list of preset values to override. Default is
#'   `list()`, which uses the selected archetype unchanged. Common controls
#'   include `min_individuals`, `max_individuals`, `rank_abundance_alpha`,
#'   `peak_margin`, `peak_beta_shape`, `duration_probabilities`,
#'   `duration_ranges`, `shape_probabilities`, `peak_sd_ratio`,
#'   `dominant_peak_sd_ratio`, `minimum_peak_sd`, `maximum_peak_sd`,
#'   `species_flower_median`, `species_flower_cv`, `individual_flower_cv`,
#'   `individual_duration_cv`, `longevity_probabilities`, `top_fraction`, and
#'   `target_top_output_share`. Archetype-specific timing controls include
#'   `first_peak_fraction`, `second_peak_fraction`,
#'   `bimodal_cluster_sd_fraction`, `first_cluster_species_fraction`,
#'   `dominant_peak_fraction`, `dominant_peak_sd_fraction`,
#'   `early_peak_fraction`, `late_peak_fraction`,
#'   `secondary_peak_sd_fraction`, `pulse_peak_fraction`, and
#'   `pulse_peak_sd_fraction`.
#'
#' @details
#' The default simulation uses 50 species, 90 daily censuses, a 12-day peak
#' margin, and 20--500 individuals per species. Species abundance follows a
#' truncated rank-abundance rule. Species and individual flowering schedules
#' are then generated from preset peak structures, flowering-duration classes,
#' curve-shape classes, flower-output distributions, and floral-longevity
#' mixtures. The six main daily matrices are retained at both individual and
#' species levels so that users can distinguish newly opened flowers, standing
#' floral display, and presence/absence.
#'
#' @section Principal archetype defaults:
#' The seven presets are reference simulation environments rather than a formal
#' classification of all natural plant communities. Their principal defaults
#' are:
#'
#' * `even_weak_seasonal`: rank-abundance exponent 0.15; broad Beta(1.4, 1.4)
#'   peak structure; duration proportions 0.20/0.55/0.25; curve-shape
#'   proportions 0.80/0.15/0.05; peak-SD ratio 0.35 (2--25 days).
#' * `strongly_seasonal_even`: exponent 0.20; central Beta(10, 10) peaks;
#'   duration proportions 0.50/0.40/0.10; curve-shape proportions
#'   0.75/0.20/0.05; peak-SD ratio 0.25 (1--18 days).
#' * `seasonal_herbaceous`: exponent 0.90; central Beta(4, 4) peaks; duration
#'   proportions 0.35/0.50/0.15; curve-shape proportions 0.55/0.35/0.10;
#'   peak-SD ratio 0.30 (2--22 days).
#' * `phenologically_staggered`: exponent 0.40; peaks spaced across the season;
#'   shorter duration ranges 5--10/12--20/25--35 days; duration proportions
#'   0.60/0.35/0.05; peak-SD ratio 0.18 (1--10 days).
#' * `bimodal_seasonal`: exponent 0.60; peak centers at 0.30 and 0.70 of the
#'   season with SD 0.045 of season length and a 50:50 species split; peak-SD
#'   ratio 0.28 (2--20 days).
#' * `mass_flowering_dominated`: exponent 1.50; early/central/late peak centers
#'   at 0.30/0.50/0.70 of season length; the top 20% of species target 80% of
#'   expected floral output; peak-SD ratio 0.25 (1--18 days), with 0.12 for
#'   dominant species.
#' * `ephemeral_pulse_flowering`: exponent 0.70; one pulse centered at 0.50 of
#'   season length with SD 0.020 of season length; duration ranges
#'   4--8/9--14/15--22 days; peak-SD ratio 0.12 (1--8 days).
#'
#' Use `phenology_archetype_defaults()` to inspect the complete preset table,
#' including flower-output CVs, duration and curve-shape proportions, and
#' floral-longevity probabilities.
#'
#' @section Output matrices:
#' Every generated matrix has one row per day. Individual matrices have one
#' column per simulated individual and species matrices have one column per
#' species. The matrices are:
#'
#' * `new_flowers_individual` and `new_flowers_species`: numbers of flowers
#'   opening on each day.
#' * `open_flowers_individual` and `open_flowers_species`: numbers of flowers
#'   currently open after floral longevity is applied.
#' * `presence_individual` and `presence_species`: logical presence/absence
#'   versions of the open-flower matrices.
#'
#' @return An object of class `true_phenology_community`, implemented as a list
#'   containing:
#'
#' * `community_type`, `description`, `season_length`, and `seed`.
#' * `configuration`: the full preset after any values in `control` are applied.
#' * `species_parameters`: species ID, abundance rank, population size,
#'   dominance status, assigned peak and duration, within-species peak SD,
#'   curve shape, skew direction, flower production, floral longevity, and
#'   realized total flowers, peak, onset, offset, and duration.
#' * `individual_parameters`: individual ID, species ID/rank, individual peak,
#'   individual duration, total flower production, curve shape, skew direction,
#'   and floral longevity.
#' * the six daily matrices described in **Output matrices** above.
#'
#' @examples
#' community <- generate_true_community(
#'   community_type = "seasonal_herbaceous",
#'   n_species = 12,
#'   season_length = 60,
#'   seed = 149,
#'   control = list(min_individuals = 3, max_individuals = 8)
#' )
#'
#' community
#' head(community$species_parameters)
#'
#' # Modify selected archetype defaults with custom values.
#' modified_community <- generate_true_community(
#'   community_type = "seasonal_herbaceous",
#'   n_species = 12,
#'   season_length = 60,
#'   seed = 202,
#'   control = list(
#'     min_individuals = 3,
#'     max_individuals = 8,
#'     rank_abundance_alpha = 1.20,
#'     duration_probabilities = c(
#'       short = 0.60,
#'       medium = 0.35,
#'       long = 0.05
#'     )
#'   )
#' )
#'
#' # Convert the true daily phenology to a species-level weekly census table.
#' weekly_species <- as_phenology_table(
#'   community,
#'   level = "species",
#'   measure = "open",
#'   census_column = "day",
#'   census_interval = 7
#' )
#' head(weekly_species)
#'
#' \dontrun{
#' # Export daily species- and individual-level phenology tables as CSV files.
#' species_daily <- as_phenology_table(
#'   community,
#'   level = "species",
#'   measure = "open",
#'   census_column = "day",
#'   census_interval = 1
#' )
#' individual_daily <- as_phenology_table(
#'   community,
#'   level = "individual",
#'   measure = "open",
#'   census_column = "day",
#'   census_interval = 1
#' )
#' write.csv(species_daily, "species_phenology.csv", row.names = FALSE)
#' write.csv(individual_daily, "individual_phenology.csv", row.names = FALSE)
#' }
#'
#' @export
generate_true_community <- function(
    community_type = c(
      "even_weak_seasonal",
      "strongly_seasonal_even",
      "seasonal_herbaceous",
      "phenologically_staggered",
      "bimodal_seasonal",
      "mass_flowering_dominated",
      "ephemeral_pulse_flowering"
    ),
    n_species = 50,
    season_length = 90,
    seed = 123,
    control = list()) {

  community_type <- match.arg(community_type)

  if (length(n_species) != 1L || is.na(n_species) || n_species < 2 ||
      n_species != as.integer(n_species)) {
    stop("n_species must be a single integer >= 2.", call. = FALSE)
  }
  if (length(season_length) != 1L || is.na(season_length) ||
      season_length < 3 || season_length != as.integer(season_length)) {
    stop("season_length must be a single integer >= 3.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("seed must be one finite numeric value.", call. = FALSE)
  }
  if (!is.list(control)) {
    stop("control must be a named list.", call. = FALSE)
  }
  if (length(control) > 0L &&
      (is.null(names(control)) || any(names(control) == ""))) {
    stop("Every element of control must be named.", call. = FALSE)
  }

  n_species <- as.integer(n_species)
  season_length <- as.integer(season_length)
  seed <- as.integer(seed)

  presets <- .pheno_community_presets()
  configuration <- presets[[community_type]]

  if (length(control) > 0L) {
    valid_control_names <- unique(unlist(lapply(presets, names), use.names = FALSE))
    unknown <- setdiff(names(control), valid_control_names)
    if (length(unknown) > 0L) {
      stop(
        "Unknown control parameter(s): ",
        paste(unknown, collapse = ", "),
        ". See ?generate_true_community and phenology_archetype_defaults().",
        call. = FALSE
      )
    }
    configuration <- utils::modifyList(configuration, control)
  }

  .pheno_check_configuration(configuration)

  species_parameters <- .pheno_build_species_parameters(
    n_species = n_species,
    season_length = season_length,
    configuration = configuration,
    seed = seed
  )

  individual_parameters <- .pheno_build_individual_parameters(
    species_parameters = species_parameters,
    season_length = season_length,
    configuration = configuration,
    seed = seed + 1000L
  )

  new_flowers_individual <- .pheno_generate_new_flower_matrix(
    individual_parameters = individual_parameters,
    season_length = season_length,
    seed = seed + 2000L
  )

  open_flowers_individual <- .pheno_apply_floral_longevity(
    new_flower_matrix = new_flowers_individual,
    longevity_days = individual_parameters$floral_longevity_days
  )
  presence_individual <- open_flowers_individual > 0L

  new_flowers_species <- .pheno_aggregate_individual_matrix_by_species(
    individual_matrix = new_flowers_individual,
    individual_species = individual_parameters$species_id
  )
  open_flowers_species <- .pheno_aggregate_individual_matrix_by_species(
    individual_matrix = open_flowers_individual,
    individual_species = individual_parameters$species_id
  )
  presence_species <- open_flowers_species > 0L

  realized_total_flowers <- as.integer(colSums(new_flowers_species))
  realized_peak_day <- apply(
    new_flowers_species,
    2,
    function(x) which.max(x)
  )
  realized_onset_day <- apply(
    presence_species,
    2,
    function(x) {
      days_present <- which(x)
      if (length(days_present) == 0L) NA_integer_ else min(days_present)
    }
  )
  realized_offset_day <- apply(
    presence_species,
    2,
    function(x) {
      days_present <- which(x)
      if (length(days_present) == 0L) NA_integer_ else max(days_present)
    }
  )
  realized_duration_days <- realized_offset_day - realized_onset_day + 1L

  species_order <- match(
    species_parameters$species_id,
    colnames(new_flowers_species)
  )
  species_parameters$realized_total_flowers <-
    realized_total_flowers[species_order]
  species_parameters$realized_peak_day <-
    realized_peak_day[species_order]
  species_parameters$realized_onset_day <-
    realized_onset_day[species_order]
  species_parameters$realized_offset_day <-
    realized_offset_day[species_order]
  species_parameters$realized_duration_days <-
    realized_duration_days[species_order]

  output <- list(
    community_type = community_type,
    description = configuration$description,
    season_length = season_length,
    seed = seed,
    configuration = configuration,
    species_parameters = species_parameters,
    individual_parameters = individual_parameters,
    new_flowers_individual = new_flowers_individual,
    open_flowers_individual = open_flowers_individual,
    presence_individual = presence_individual,
    new_flowers_species = new_flowers_species,
    open_flowers_species = open_flowers_species,
    presence_species = presence_species
  )

  class(output) <- "true_phenology_community"
  output
}

#' Inspect default parameters for the seven phenological archetypes
#'
#' @description
#' Returns a compact table of the principal default values used by
#' `generate_true_community()`. The table is intended as an easily inspectable
#' version of the archetype structure described in that function's help page.
#'
#' @return A data frame with one row per community archetype and columns for
#'   default richness, season length, abundance structure, peak structure,
#'   duration ranges and proportions, curve-shape proportions, peak-date
#'   variation, flower-output parameters, floral longevity, and the optional
#'   dominant-species floral-output target.
#'
#' @examples
#' phenology_archetype_defaults()
#'
#' @export
phenology_archetype_defaults <- function() {
  presets <- .pheno_community_presets()
  ids <- names(presets)

  display_name <- c(
    even_weak_seasonal = "Even, weakly seasonal",
    strongly_seasonal_even = "Strongly seasonal, even",
    seasonal_herbaceous = "Seasonal herbaceous",
    phenologically_staggered = "Phenologically staggered",
    bimodal_seasonal = "Bimodal seasonal",
    mass_flowering_dominated = "Mass-flowering dominated",
    ephemeral_pulse_flowering = "Ephemeral pulse-flowering"
  )

  fmt_named <- function(x) {
    paste(
      paste0(names(x), "=", format(as.numeric(x), trim = TRUE, scientific = FALSE)),
      collapse = " / "
    )
  }

  fmt_ranges <- function(x) {
    paste(
      vapply(
        names(x),
        function(nm) paste0(nm, "=", x[[nm]][1], "-", x[[nm]][2]),
        character(1)
      ),
      collapse = " / "
    )
  }

  peak_structure_label <- c(
    even_weak_seasonal = "Broad beta",
    strongly_seasonal_even = "Narrow central beta",
    seasonal_herbaceous = "Central beta",
    phenologically_staggered = "Evenly staggered",
    bimodal_seasonal = "Two seasonal clusters",
    mass_flowering_dominated = "Three clusters",
    ephemeral_pulse_flowering = "Single narrow pulse"
  )

  peak_parameter_label <- c(
    even_weak_seasonal = "Beta(1.4, 1.4)",
    strongly_seasonal_even = "Beta(10, 10)",
    seasonal_herbaceous = "Beta(4, 4)",
    phenologically_staggered = "Even spacing within peak margins",
    bimodal_seasonal = "Centers 0.30L / 0.70L; SD 0.045L; 50:50 species",
    mass_flowering_dominated = paste0(
      "Centers 0.30L / 0.50L / 0.70L; SD 0.055L secondary, ",
      "0.035L dominant"
    ),
    ephemeral_pulse_flowering = "Center 0.50L; SD 0.020L"
  )

  data.frame(
    community_type = ids,
    archetype = unname(display_name[ids]),
    species_default = rep(50L, length(ids)),
    season_length_days_default = rep(90L, length(ids)),
    individuals_per_species_default = vapply(
      presets,
      function(p) paste0(p$min_individuals, "-", p$max_individuals),
      character(1)
    ),
    rank_abundance_alpha = vapply(presets, `[[`, numeric(1), "rank_abundance_alpha"),
    peak_structure = unname(peak_structure_label[ids]),
    peak_parameters = unname(peak_parameter_label[ids]),
    peak_margin_days = vapply(presets, `[[`, numeric(1), "peak_margin"),
    duration_ranges_days = vapply(
      presets,
      function(p) fmt_ranges(p$duration_ranges),
      character(1)
    ),
    duration_proportions_short_medium_long = vapply(
      presets,
      function(p) fmt_named(p$duration_probabilities),
      character(1)
    ),
    curve_shape_proportions = vapply(
      presets,
      function(p) fmt_named(p$shape_probabilities),
      character(1)
    ),
    peak_sd_ratio = vapply(presets, `[[`, numeric(1), "peak_sd_ratio"),
    peak_sd_bounds_days = vapply(
      presets,
      function(p) paste0(p$minimum_peak_sd, "-", p$maximum_peak_sd),
      character(1)
    ),
    dominant_peak_sd_ratio = vapply(
      presets,
      function(p) if (is.null(p$dominant_peak_sd_ratio)) NA_real_ else p$dominant_peak_sd_ratio,
      numeric(1)
    ),
    median_flowers_per_individual = vapply(
      presets,
      `[[`,
      numeric(1),
      "species_flower_median"
    ),
    species_flower_output_cv = vapply(presets, `[[`, numeric(1), "species_flower_cv"),
    individual_flower_output_cv = vapply(presets, `[[`, numeric(1), "individual_flower_cv"),
    individual_duration_cv = vapply(presets, `[[`, numeric(1), "individual_duration_cv"),
    floral_longevity_proportions_1_2_3_5_days = vapply(
      presets,
      function(p) fmt_named(p$longevity_probabilities),
      character(1)
    ),
    top_species_fraction = vapply(presets, `[[`, numeric(1), "top_fraction"),
    target_top_species_output_share = vapply(
      presets,
      `[[`,
      numeric(1),
      "target_top_output_share"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Convert a simulated community to a phenology table
#'
#' @description
#' Converts a `true_phenology_community` object returned by
#' `generate_true_community()` to a wide phenology table suitable for direct use
#' with `calculate_indices()`. The user can return species-level or
#' individual-level phenology and can retain the complete daily true phenology
#' or sample it at a regular census interval.
#'
#' @param simulation A `true_phenology_community` object returned by
#'   `generate_true_community()`.
#' @param level Character. `"species"` (default) or `"individual"`.
#' @param measure Character. `"open"` (default) for currently open flowers,
#'   `"new"` for newly opened flowers, or `"presence"` for presence/absence.
#' @param census_column Character name for the first column containing the
#'   retained census-day numbers. Default is `"day"`.
#' @param census_interval Positive whole number giving the number of days between
#'   retained censuses. Default is `1`, which retains the complete daily true
#'   phenology. For example, `census_interval = 7` retains days 1, 8, 15, and so
#'   on until the end of the simulated season.
#'
#' @return A wide data frame with the census column first, followed by one column
#'   per species or individual. Rows correspond to the retained census days and
#'   cells contain the selected flowering measure. The returned table can be
#'   passed directly to `calculate_indices()` by supplying the same
#'   `census_column` name to its `census_col` argument.
#'
#' @examples
#' community <- generate_true_community(
#'   n_species = 8,
#'   season_length = 60,
#'   seed = 1,
#'   control = list(min_individuals = 3, max_individuals = 8)
#' )
#'
#' # Complete daily species-level standing floral display.
#' species_daily <- as_phenology_table(
#'   community,
#'   level = "species",
#'   measure = "open",
#'   census_column = "day"
#' )
#' head(species_daily)
#'
#' # Sample the same true phenology every 7 days.
#' species_weekly <- as_phenology_table(
#'   community,
#'   level = "species",
#'   measure = "open",
#'   census_column = "day",
#'   census_interval = 7
#' )
#' species_weekly
#'
#' # Individual-level presence/absence sampled every 5 days.
#' individual_presence <- as_phenology_table(
#'   community,
#'   level = "individual",
#'   measure = "presence",
#'   census_column = "day",
#'   census_interval = 5
#' )
#'
#' @export
as_phenology_table <- function(
    simulation,
    level = c("species", "individual"),
    measure = c("open", "new", "presence"),
    census_column = "day",
    census_interval = 1) {
  if (!inherits(simulation, "true_phenology_community")) {
    stop("simulation must be a true_phenology_community object.", call. = FALSE)
  }

  level <- match.arg(level)
  measure <- match.arg(measure)

  if (length(census_column) != 1L || is.na(census_column) ||
      !is.character(census_column) || !nzchar(census_column)) {
    stop("census_column must be one non-empty character string.", call. = FALSE)
  }
  if (length(census_interval) != 1L || is.na(census_interval) ||
      !is.finite(census_interval) || census_interval < 1 ||
      census_interval != as.integer(census_interval)) {
    stop("census_interval must be a single positive whole number.", call. = FALSE)
  }
  census_interval <- as.integer(census_interval)

  object_name <- paste0(measure, "_flowers_", level)
  if (measure == "presence") {
    object_name <- paste0("presence_", level)
  }

  matrix_object <- simulation[[object_name]]
  if (is.null(matrix_object)) {
    stop("Requested matrix is not available: ", object_name, call. = FALSE)
  }

  census_days <- seq.int(
    from = 1L,
    to = simulation$season_length,
    by = census_interval
  )
  sampled_matrix <- matrix_object[census_days, , drop = FALSE]

  output <- as.data.frame(sampled_matrix, check.names = FALSE)
  output <- cbind(
    stats::setNames(
      data.frame(census_days),
      census_column
    ),
    output
  )
  rownames(output) <- NULL
  output
}

#' Summarise a generated true phenological community
#'
#' @param simulation A `true_phenology_community` object returned by
#'   `generate_true_community()`.
#'
#' @return A one-row data frame reporting community type, numbers of species and
#'   individuals, season length, population-size range, the population and
#'   flower-output share of the configured top species fraction, and summary
#'   statistics for species peak dates and flowering durations.
#'
#' @export
summarise_true_community <- function(simulation) {
  if (!inherits(simulation, "true_phenology_community")) {
    stop("simulation must be a true_phenology_community object.", call. = FALSE)
  }

  species_parameters <- simulation$species_parameters
  output_by_species <- colSums(simulation$new_flowers_species)

  top_n <- max(
    1L,
    ceiling(
      nrow(species_parameters) * simulation$configuration$top_fraction
    )
  )
  top_species <- species_parameters$species_id[seq_len(top_n)]

  top_output_share <- sum(output_by_species[top_species]) / sum(output_by_species)
  top_population_share <- sum(
    species_parameters$population_size[seq_len(top_n)]
  ) / sum(species_parameters$population_size)

  data.frame(
    community_type = simulation$community_type,
    species = nrow(species_parameters),
    individuals = nrow(simulation$individual_parameters),
    season_length = simulation$season_length,
    minimum_population_size = min(species_parameters$population_size),
    median_population_size = stats::median(species_parameters$population_size),
    maximum_population_size = max(species_parameters$population_size),
    top_20_percent_population_share = top_population_share,
    top_20_percent_flower_share = top_output_share,
    median_species_peak_day = stats::median(species_parameters$species_peak_day),
    peak_day_range = paste(range(species_parameters$species_peak_day), collapse = "-"),
    median_assigned_duration = stats::median(species_parameters$species_duration_days),
    median_realized_duration = stats::median(species_parameters$realized_duration_days),
    stringsAsFactors = FALSE
  )
}

#' Validate a generated true phenological community
#'
#' @param simulation A `true_phenology_community` object returned by
#'   `generate_true_community()`.
#'
#' @return A data frame with one row per validation check and logical column
#'   `passed`. Checks cover matrix count types, conservation of individual and
#'   species flower totals, presence/absence types, peak bounds, population-size
#'   bounds, and the optional dominant-species flower-share target.
#'
#' @examples
#' community <- generate_true_community(n_species = 8, season_length = 60, seed = 1)
#' validate_true_community(community)
#'
#' @export
validate_true_community <- function(simulation) {
  if (!inherits(simulation, "true_phenology_community")) {
    stop("simulation must be a true_phenology_community object.", call. = FALSE)
  }

  species_parameters <- simulation$species_parameters
  individual_parameters <- simulation$individual_parameters

  individual_total_generated <- colSums(simulation$new_flowers_individual)
  species_total_from_individuals <- sum(simulation$new_flowers_individual)
  species_total_from_aggregate <- sum(simulation$new_flowers_species)
  open_total_from_individuals <- sum(simulation$open_flowers_individual)
  open_total_from_aggregate <- sum(simulation$open_flowers_species)
  top_share <- summarise_true_community(simulation)$top_20_percent_flower_share
  target_share <- simulation$configuration$target_top_output_share

  individual_new_are_integer <-
    is.integer(simulation$new_flowers_individual) &&
    all(simulation$new_flowers_individual >= 0L)
  individual_open_are_integer <-
    is.integer(simulation$open_flowers_individual) &&
    all(simulation$open_flowers_individual >= 0L)
  species_new_are_integer <-
    is.integer(simulation$new_flowers_species) &&
    all(simulation$new_flowers_species >= 0L)
  species_open_are_integer <-
    is.integer(simulation$open_flowers_species) &&
    all(simulation$open_flowers_species >= 0L)

  data.frame(
    test = c(
      "Number of species is correct",
      "Number of individuals equals summed population sizes",
      "Individual newly opened flower values are non-negative integers",
      "Individual open-flower values are non-negative integers",
      "Species newly opened flower values are non-negative integers",
      "Species open-flower values are non-negative integers",
      "Every individual retains its assigned total flower output exactly",
      "Species aggregation conserves total newly opened flowers exactly",
      "Species aggregation conserves total open flowers exactly",
      "Individual presence matrix is logical",
      "Species presence matrix is logical",
      "All individual peaks fall inside the season",
      "Population sizes remain inside configured bounds",
      "Dominated-community flower-share target is approximately met"
    ),
    passed = c(
      nrow(species_parameters) == ncol(simulation$new_flowers_species),
      nrow(individual_parameters) == sum(species_parameters$population_size),
      individual_new_are_integer,
      individual_open_are_integer,
      species_new_are_integer,
      species_open_are_integer,
      all(individual_total_generated == individual_parameters$total_flowers),
      species_total_from_individuals == species_total_from_aggregate,
      open_total_from_individuals == open_total_from_aggregate,
      is.logical(simulation$presence_individual),
      is.logical(simulation$presence_species),
      all(
        individual_parameters$individual_peak_day >= 1L &
          individual_parameters$individual_peak_day <= simulation$season_length
      ),
      min(species_parameters$population_size) >=
        simulation$configuration$min_individuals &&
        max(species_parameters$population_size) <=
        simulation$configuration$max_individuals,
      if (is.na(target_share)) TRUE else abs(top_share - target_share) <= 0.05
    ),
    stringsAsFactors = FALSE
  )
}

#' Plot structural diagnostics for a true phenological community
#'
#' @param simulation A `true_phenology_community` object returned by
#'   `generate_true_community()`.
#' @param top_fraction Proportion of the most seasonally productive species to
#'   highlight in the fourth panel. Default is `0.20`.
#'
#' @return Invisibly returns `NULL`. Produces a four-panel base-R figure showing
#'   rank abundance, a species flowering heat map, total community flowering,
#'   and the daily flower-output share of the top species.
#'
#' @export
plot_true_community <- function(simulation, top_fraction = 0.20) {
  if (!inherits(simulation, "true_phenology_community")) {
    stop("simulation must be a true_phenology_community object.", call. = FALSE)
  }
  if (length(top_fraction) != 1L || is.na(top_fraction) ||
      top_fraction <= 0 || top_fraction > 1) {
    stop("top_fraction must lie in (0, 1].", call. = FALSE)
  }

  species_parameters <- simulation$species_parameters
  species_open <- simulation$open_flowers_species
  species_new <- simulation$new_flowers_species
  days <- seq_len(simulation$season_length)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  graphics::plot(
    species_parameters$species_rank,
    species_parameters$population_size,
    type = "b",
    pch = 16,
    xlab = "Species abundance rank",
    ylab = "Number of individuals",
    main = "Population-size structure"
  )

  species_max <- apply(species_open, 2, max)
  species_max[species_max == 0] <- 1
  normalized_open <- sweep(species_open, 2, species_max, FUN = "/")
  peak_order <- order(species_parameters$species_peak_day)

  graphics::image(
    x = days,
    y = seq_len(ncol(normalized_open)),
    z = normalized_open[, peak_order, drop = FALSE],
    xlab = "Day",
    ylab = "Species ordered by peak",
    main = "Relative flowering intensity",
    useRaster = TRUE
  )

  community_output <- rowSums(species_new)
  graphics::plot(
    days,
    community_output,
    type = "l",
    lwd = 2,
    xlab = "Day",
    ylab = "Newly opened flowers",
    main = "Community flowering curve"
  )

  n_top <- max(1L, ceiling(nrow(species_parameters) * top_fraction))
  seasonal_output <- colSums(species_new)
  top_species <- names(sort(seasonal_output, decreasing = TRUE))[seq_len(n_top)]
  daily_total <- rowSums(species_new)
  daily_top <- rowSums(species_new[, top_species, drop = FALSE])
  daily_share <- ifelse(daily_total > 0, daily_top / daily_total, NA_real_)

  graphics::plot(
    days,
    daily_share,
    type = "l",
    lwd = 2,
    ylim = c(0, 1),
    xlab = "Day",
    ylab = "Top-species share",
    main = paste0(
      "Daily output share: top ",
      round(top_fraction * 100),
      "% species"
    )
  )

  if (!is.na(simulation$configuration$target_top_output_share)) {
    graphics::abline(
      h = simulation$configuration$target_top_output_share,
      lty = 2
    )
  }

  invisible(NULL)
}

#' Plot flowering curves for selected species
#'
#' @param simulation A `true_phenology_community` object returned by
#'   `generate_true_community()`.
#' @param species_ids Optional character vector of species IDs to plot. If
#'   `NULL`, species are selected across the range of assigned peak dates.
#' @param n_species Number of species to select when `species_ids = NULL`.
#'   Default is `9`.
#'
#' @return Invisibly returns the species IDs that were plotted. Produces base-R
#'   line plots of daily open flowers; dashed vertical lines show assigned
#'   species peak days.
#'
#' @export
plot_species_examples <- function(simulation,
                                  species_ids = NULL,
                                  n_species = 9) {
  if (!inherits(simulation, "true_phenology_community")) {
    stop("simulation must be a true_phenology_community object.", call. = FALSE)
  }

  species_parameters <- simulation$species_parameters
  species_open <- simulation$open_flowers_species
  days <- seq_len(simulation$season_length)

  if (is.null(species_ids)) {
    n_species <- max(1L, as.integer(n_species))
    ordered <- order(species_parameters$species_peak_day)
    selected_positions <- unique(
      round(
        seq(
          1,
          length(ordered),
          length.out = min(n_species, length(ordered))
        )
      )
    )
    species_ids <- species_parameters$species_id[ordered[selected_positions]]
  }

  missing_species <- setdiff(species_ids, colnames(species_open))
  if (length(missing_species) > 0L) {
    stop(
      "Unknown species IDs: ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  n_panels <- length(species_ids)
  n_columns <- ceiling(sqrt(n_panels))
  n_rows <- ceiling(n_panels / n_columns)
  graphics::par(mfrow = c(n_rows, n_columns), mar = c(3, 3, 2, 1))

  for (species_id in species_ids) {
    species_row <- species_parameters[
      species_parameters$species_id == species_id,
      ,
      drop = FALSE
    ]

    graphics::plot(
      days,
      species_open[, species_id],
      type = "l",
      lwd = 2,
      xlab = "Day",
      ylab = "Open flowers",
      main = paste0(
        species_id,
        "\n",
        species_row$curve_shape,
        ", duration = ",
        species_row$species_duration_days
      )
    )
    graphics::abline(v = species_row$species_peak_day, lty = 2)
  }

  invisible(species_ids)
}

#' @export
print.true_phenology_community <- function(x, ...) {
  cat("\nTrue phenological community\n")
  cat("----------------------------\n")
  cat("Type:        ", x$community_type, "\n", sep = "")
  cat("Description: ", x$description, "\n", sep = "")
  cat("Species:     ", nrow(x$species_parameters), "\n", sep = "")
  cat("Individuals: ", nrow(x$individual_parameters), "\n", sep = "")
  cat("Season:      ", x$season_length, " days\n", sep = "")
  cat("Flower data: non-negative integer daily counts\n\n")

  print(summarise_true_community(x), row.names = FALSE)
  invisible(x)
}
