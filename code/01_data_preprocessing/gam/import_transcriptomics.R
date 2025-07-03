#' Read-in and prepare transcriptomics data
#' @author: CShannon
#'
#' Changelog:
#'
#' 2021-06-17 - Initial commit
#' 2021-06-29 - Use cDNA batch to correct (per A.Lee)
#' 2021-07-06 - Use edgeR filterByExpr to exclude low count genes
#' 2021-11-17 - Can't use 'submission-ready' data from DMC - older version, not batch-corrected using ComBat-seq...
#' 2021-11-23 - Use a Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# read in counts
data <- readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_RNA_COUNTS.csv')

# make matrix
data <- data %>%
  column_to_rownames('...1') %>%
  as.matrix()

# save
saveRDS(data, 'data/processed/gam/import_transcriptomics.rds', compress = F)
