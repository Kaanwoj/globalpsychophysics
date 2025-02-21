## code to prepare `matching` dataset
matching <- read.table("../../matching/data/matching_data.txt", header = TRUE,
                  sep = ";")

usethis::use_data(matching, overwrite = TRUE)
