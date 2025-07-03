#' Read-in and prepare luminex (cytokine) data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2020-02-19 - Initial commit
#' 2021-11-17 - Using 'submission-ready' data from DMC
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# source
meta <- readRDS('data/processed/gam/import_metadata.rds')

# read in - cytokines
data <- readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_CYT.csv')

# subset to Visit Identifiers in metadata
table(data$vid %in% meta$`Visit Identifier`) # all present

# make matrix and name rows as Unique Identifier + Visit Identifier
data <- data %>%
  column_to_rownames('vid') %>%
  select(-flagged, -comment) %>%
  as.matrix()

data <- data[intersect(meta$`Visit Identifier`, rownames(data)), ]

saveRDS(data, 'data/processed/gam/import_luminex_cytokines.rds')
