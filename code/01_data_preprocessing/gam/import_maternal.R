#' Read-in and prepare HepB antibody response data:
#' @author: CShannon
#'
#' Changelog:
#'

# import
library(tidyverse)

# read-in maternal metadata
maternal <- aws.s3::s3read_using(readr::read_csv, object = 's3://epichipc-main/Clean_Data/Clinical/GAM_Main_Tier1_Clean.csv')

# clean-up
maternal <- maternal %>%
  select(
    subject_id         = 'Unique Identifier',
    birth_weight       = 'Weight (grams) V1',
    maternal_age       = 'Age (years) Maternal',
    maternal_ethnicity = 'Ethnicity Neonatal',
    maternal_race      = 'Race Neonatal'
  )

saveRDS(maternal, 'data/processed/gam/import_maternal.rds')
