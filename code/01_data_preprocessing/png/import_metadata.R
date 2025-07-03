#' Read-in and prepare metadata:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-05 - Initial commit

# import
library(tidyverse)

# read-in and make it long
meta <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_CLINICAL.csv')

meta <- meta %>%
  dplyr::select(Unique_Identifier, Visit_ID_V1, Visit_ID_V2, Sex, VaccineGroup) %>%
  dplyr::rename(V1 = Visit_ID_V1, V2 = Visit_ID_V2) %>%
  tidyr::gather(Visit, VID, V1, V2) %>%
  dplyr::mutate(DOL = Visit %>% factor(levels = c('V1', 'V2'), labels = c(0, 7)) %>% as.character() %>% as.numeric()) %>%
  dplyr::select(VID, Unique_Identifier, Sex, VaccineGroup, DOL) %>%
  dplyr::arrange(Unique_Identifier, DOL)

saveRDS(meta, 'data/processed/png/import_metadata.rds')
