#' Read-in and prepare metabolomics data with additional batch correction
#' efforts:
#'
#' @author: CShannon
#'
#' Changelog:
#'
#' 2021-11-22 - Initial commit
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# read-in
df <- readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_MET_COUNTS_WITHXENO.csv')

# fix names and make matrix
df <- df %>%
  rename(sample_id = '...1') %>%
  column_to_rownames('sample_id') %>%
  as.matrix()

saveRDS(df, 'data/processed/gam/import_metabolomics.rds', compress = F)
