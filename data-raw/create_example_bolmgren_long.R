# Run this script from the package root.
# It creates one example dataset and saves it to data/ using usethis::use_data().

set.seed(1028)

n_plants <- 8
plant_ids <- paste0("plant_", sprintf("%02d", seq_len(n_plants)))

example_bolmgren_long <- data.frame()

for (p in plant_ids) {
  n_flowers <- sample(3:8, 1)
  onset_vals <- sort(sample(1:20, n_flowers, replace = FALSE))
  durations <- sample(2:6, n_flowers, replace = TRUE)
  end_vals <- onset_vals + durations

  example_bolmgren_long <- rbind(
    example_bolmgren_long,
    data.frame(
      plant_id = p,
      onset_day = onset_vals,
      end_day = end_vals,
      stringsAsFactors = FALSE
    )
  )
}

usethis::use_data(example_bolmgren_long, overwrite = TRUE)
