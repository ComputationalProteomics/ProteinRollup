library(tidyverse)

print(getwd())

source("../../ProteinRollup.R")

df <- read_tsv("../../testdata/dummy_df.tsv")


pr$protein_rollup(df$Protein, df[, -c(1,2)] %>% as.matrix())


# is_output_identical <- function() {
#     
# }

# expect_true(1 == 1)