#' Read-in and prepare transcriptomics data
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-05 - Initial commit

# import
library(tidyverse)

# read in counts
data <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_RNA_COUNTS_V1V2.csv')

# make matrix
data <- data %>%
  separate(sample_id, c(NA, NA, 'vid')) %>%
  column_to_rownames('vid') %>%
  as.matrix()

# save
saveRDS(data, 'data/processed/png/import_transcriptomics.rds', compress = F)
