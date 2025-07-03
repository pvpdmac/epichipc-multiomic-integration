#' Read-in and prepare epigenetic data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2021-11-01 - Initial commit
#' 2021-11-17 - Using 'submission-ready' data from DMC
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# load prepared object to temporary environment
env <- new.env()
load(
  envir = env,
  file = 'data/raw/gam/GAMMAIN_SINGLEOMICS_EPI.RData'
)

data <- env$EPI_M_clean_2 %>% as.matrix() %>% t()

saveRDS(data, 'data/processed/gam/import_epigenetics.rds', compress = F)
