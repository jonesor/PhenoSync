# Run this script from the package root.
# It creates one example dataset and saves it to data/ using usethis::use_data().

set.seed(1028)

n_ind <- 16
n_days <- 40
ids <- paste0("ind", sprintf("%02d", seq_len(n_ind)))

example_metadata <- data.frame(
  id = ids,
  species = rep(c("spA", "spB"), length.out = n_ind),
  site = rep(c("site1", "site2"), each = n_ind / 2),
  x = rep(seq(0, 45, length.out = 4), times = 4),
  y = rep(seq(0, 45, length.out = 4), each = 4),
  stringsAsFactors = FALSE
)

flower_mat <- matrix(0, nrow = n_days, ncol = n_ind)
colnames(flower_mat) <- ids
rownames(flower_mat) <- paste0("day", seq_len(n_days))

starts <- sample(1:24, n_ind, replace = TRUE)
durations <- sample(8:16, n_ind, replace = TRUE)

for (i in seq_along(ids)) {
  s <- starts[i]
  e <- min(n_days, s + durations[i] - 1)
  flower_mat[s:e, i] <- sample(1:20, e - s + 1, replace = TRUE)
}

example_neighbourhood <- list(
  metadata = example_metadata,
  flower_mat = flower_mat
)

usethis::use_data(example_neighbourhood, overwrite = TRUE)
