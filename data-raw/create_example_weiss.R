# Run this script from the package root.
# It creates one example dataset and saves it to data/ using usethis::use_data().

set.seed(1028)

n_ind <- 15
n_days <- 12

f.date <- sample(1:6, n_ind, replace = TRUE)
flr.mat <- matrix(0, nrow = n_ind, ncol = n_days)
colnames(flr.mat) <- paste0("flrs.", seq_len(n_days))

for (i in seq_len(n_ind)) {
  peak <- f.date[i] + sample(2:5, 1)
  duration <- sample(4:7, 1)

  start <- max(1, peak - floor(duration / 2))
  end <- min(n_days, peak + floor(duration / 2))

  flr.mat[i, start:end] <- sample(1:8, end - start + 1, replace = TRUE)
}

# Ensure every census date has at least one open flower.
for (d in seq_len(n_days)) {
  if (sum(flr.mat[, d]) == 0) {
    flr.mat[sample(seq_len(n_ind), 1), d] <- 1
  }
}

example_weiss <- data.frame(
  f.date = f.date,
  flr.mat,
  check.names = FALSE
)

usethis::use_data(example_weiss, overwrite = TRUE)
