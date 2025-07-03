#' Read-in and prepare HepB antibody response data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2020-02-19 - Initial commit
#' 2020-03-25 - Updated to reflect renaming of files on S3
#' 2021-11-17 - Using 'submission-ready' data from DMC
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# read-in and make it long
titers <- readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_TITERS.csv')

titers <- titers %>%
  select(`Subject ID`, `Visit Num`, ConcentrationImputed) %>% # we're only using imputed titer values
  spread(`Visit Num`, ConcentrationImputed) %>%
  mutate(across(c(MAT, V1, V3, V4), as.numeric))

saveRDS(titers, 'data/processed/gam/import_responses.rds')
