#' Read-in and prepare luminex (cytokine) data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-05 - Initial commit
#' 2023-12-15 - Fixed colnames to match gambia
#' 2024-01-18 - Removed profiles from V3 or V4 (only ran for cytokines)

# import
library(tidyverse)

# read in - cytokines
data <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_CYT_invivo_updated_15Dec2023.csv')

# make matrix and name rows as Unique Identifier + Visit Identifier
data <- data %>%
  column_to_rownames('vid') %>%
  select(-flagged, -comment) %>%
  as.matrix()

# remove columns with too many nas
data <- data[ , !matrixStats::colAnyNAs(data)]

# finally, remove profiles from V3 or V4 (only ran for cytokines...)
ids_v3_v4 <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_CLINICAL.csv') %>%
  select(Visit_ID_V3, Visit_ID_V4) %>%
  as.list() %>%
  reduce(union) %>%
  na.omit()

data <- data[!rownames(data) %in% ids_v3_v4, ]

saveRDS(data, 'data/processed/png/import_luminex_cytokines.rds')
