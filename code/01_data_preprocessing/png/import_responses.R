#' Read-in and prepare HepB antibody response data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-28 - Initial commit

# import
library(tidyverse)

# read-in and make it long
titers <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_TITERS.csv')

# fix names
titers <- titers %>%
  mutate(Visit_Num = factor(Visit_Num, levels = c('MAT_V1', 'INF_V1', 'INF_V3', 'INF_V4'), labels = c('MAT', 'V1', 'V3', 'V4'))) %>%
  rename('Visit Num' = Visit_Num)

titers <- titers %>%
  select(`Subject ID`, `Visit Num`, ConcentrationImputed) %>% # we're only using imputed titer values
  spread(`Visit Num`, ConcentrationImputed) %>%
  mutate(across(c(MAT, V1, V3, V4), as.numeric))

saveRDS(titers, 'data/processed/png/import_responses.rds')
