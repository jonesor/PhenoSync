# Run this script from the package root to regenerate all example datasets.
#
# This script removes legacy data/*.R source files if present. Package datasets
# should be saved as .rda files in data/, while the generation code belongs in
# data-raw/.

if (dir.exists("data")) {
  unlink(list.files("data", pattern = "\\.R$", full.names = TRUE))
}

source("data-raw/create_example_wide.R", local = TRUE)
source("data-raw/create_example_bolmgren_long.R", local = TRUE)
source("data-raw/create_example_dainou.R", local = TRUE)
source("data-raw/create_example_weiss.R", local = TRUE)
source("data-raw/create_example_neighbourhood.R", local = TRUE)
