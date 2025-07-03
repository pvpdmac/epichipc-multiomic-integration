#' Read-in and prepare metadata:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2020-02-19 - Initial commit
#' 2020-03-25 - Updated to reflect renaming of files on S3
#' 2020-12-14 - Fungal contamination list added
#' 2021-11-17 - Using 'submission-ready' data from DMC
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# read-in and make it long
meta <-readr::read_csv('data/raw/gam/GAM_Main_2020-05-01_EssentialClinicalData_Clean.csv') %>%
  gather(Visit, 'Visit Identifier', contains('Visit')) %>%
  mutate(Visit = gsub('^(.+)_Visit ID$', '\\1', Visit)) %>%
  select(`Visit Identifier`, Visit, `Unique Identifier`,  `Randomization Group`, Sex) %>%
  arrange(`Unique Identifier`)

# remove cord blood and maternal sample info
meta <- meta %>% filter(Visit %in% c('V1', 'V2'))

# add corresponding vaccine and DOL
randomization <- readr::read_csv('data/raw/gam/GAM_RandomizationGroupDefinition.csv')

# split to sample- and subject-specific metadata
meta_samp <- meta %>% distinct(`Unique Identifier`, `Visit Identifier`, Visit)
meta_subj <- meta %>% distinct(`Unique Identifier`, `Randomization Group`, Sex) %>% left_join(randomization)

# fix day column to reflect DOL0 - DOLX format
meta <- full_join(meta_subj, meta_samp)
meta <- meta %>%
  arrange(`Unique Identifier`, Visit) %>%
  group_by(`Unique Identifier`) %>%
  mutate(DOL = ifelse(Visit == 'V1', 0, Day)) %>%
  ungroup() %>%
  select(`Unique Identifier`, `Randomization Group`, VaccineGrp = Vaccine, DayGrp = Day, Sex, DOL, `Visit Identifier`) %>%
  mutate(DayGrp = paste0('DOL', DayGrp))

rm(meta_subj, meta_samp, randomization)

saveRDS(meta, 'data/processed/gam/import_metadata.rds')
