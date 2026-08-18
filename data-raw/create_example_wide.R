# Run this script from the package root.
# It creates one example dataset and saves it to data/ using usethis::use_data().

set.seed(1028)

n_ind <- 12
n_days <- 60
ids <- paste0("ind", sprintf("%02d", seq_len(n_ind)))

example_wide <- matrix(0, nrow = n_days, ncol = n_ind)
colnames(example_wide) <- ids
rownames(example_wide) <- paste0("day", seq_len(n_days))

starts <- sample(1:42, n_ind, replace = TRUE)
durations <- sample(8:20, n_ind, replace = TRUE)

for (i in seq_len(n_ind)) {
  s <- starts[i]
  e <- min(n_days, s + durations[i] - 1)
  example_wide[s:e, i] <- sample(1:25, e - s + 1, replace = TRUE)
}

usethis::use_data(example_wide, overwrite = TRUE)
