# Run this script from the package root.
# It creates one example dataset and saves it to data/ using usethis::use_data().

set.seed(407)

n_female <- 8
n_male <- 7
n_days <- 24

female_ids <- paste0("F", sprintf("%02d", seq_len(n_female)))
male_ids <- paste0("M", sprintf("%02d", seq_len(n_male)))

ids <- c(female_ids, male_ids)
sex <- c(rep("F", n_female), rep("M", n_male))

mat <- matrix(0, nrow = length(ids), ncol = n_days)
colnames(mat) <- paste0("day", seq_len(n_days))
rownames(mat) <- ids

starts <- sample(1:15, length(ids), replace = TRUE)
durations <- sample(5:12, length(ids), replace = TRUE)

for (i in seq_along(ids)) {
  s <- starts[i]
  e <- min(n_days, s + durations[i] - 1)
  mat[i, s:e] <- sample(1:15, e - s + 1, replace = TRUE)
}

example_dainou <- data.frame(
  id = ids,
  sex = sex,
  mat,
  check.names = FALSE
)

usethis::use_data(example_dainou, overwrite = TRUE)
