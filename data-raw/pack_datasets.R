## code to prepare `matching` dataset
matching <- read.table("../../matching/data/matching_data.txt", header = TRUE,
                  sep = ";")

usethis::use_data(matching, overwrite = TRUE)

## code to prepare 'default_gpm_priors' dataset
default_gpm_priors <- read.csv("../data-raw/default_gpm_priors.csv", sep = ";") |>
  (\(x) setNames(as.list(x$value), x$parname))()

usethis::use_data(default_gpm_priors, overwrite = TRUE)
