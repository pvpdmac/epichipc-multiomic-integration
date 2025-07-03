#' Read-in and prepare proteomics data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2020-02-19 - Initial commit
#' 2020-03-25 - Updated to reflect renaming of files on S3
#' 2020-03-25 - Updated to use normalized/imputed data provided by B.Fatou
#' 2021-11-17 - Using 'submission-ready' data from DMC
#' 2021-11-23 - Use Makefile
#' 2022-03-06 - Additional metabolomics exclusions

# import
library(tidyverse)

# source
meta <- readRDS('data/processed/gam/import_metadata.rds')

data <- readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_PRT.csv', col_names = F, skip = 2)

colnames <- readr::read_lines('data/raw/gam/GAMMAIN_SINGLEOMICS_PRT.csv', n_max = 1, skip = 1) %>%
  strsplit(',') %>%
  dplyr::first()

# save plate info for later
plates <- data %>% select(X1) %>% separate(X1, into = c('Plate', 'Unique Identifier', 'Visit Identifier'))

# drop all but the sample ids and protein data
data <- data %>%
  select(-X1) %>%
  column_to_rownames('X2') %>%
  as.matrix()

colnames(data) <- colnames %>% tail(-2) %>% gsub('.+(sp.+HUMAN).+', '\\1', .)

data <- data[intersect(meta$`Visit Identifier`, rownames(data)), ]

# impute ------------------------------------------------------------------

# missing value filter
df_long <- data %>%
  as_tibble(rownames = 'Visit Identifier') %>%
  gather(protein, conc, starts_with('sp')) %>%
  inner_join(plates, .)

df_missing <- df_long %>%
  group_by(Plate, protein) %>%
  summarize(n_missing = sum(is.nan(conc))/n()) %>%
  ungroup() %>%
  mutate(
    Plate = gsub('^Plate(.+)$', '\\1', Plate),
    Plate = sprintf('Plate %02d', as.numeric(Plate))
  )

keep <- df_missing %>%
  group_by(protein) %>%
  summarise(missing = median(n_missing)) %>%
  arrange(missing) %>%
  filter(missing <= 0.1) %>%
  pull(protein) %>%
  unique()

data_sub <- data[ , keep]

# impute remaining missing values
data_sub_impute <- apply(data_sub, 2, function(x) {
  imp <- min(x, na.rm = T)
  x[is.nan(x)] <- imp/2
  x
})

# write out ---------------------------------------------------------------

saveRDS(data_sub_impute, 'data/processed/gam/import_proteomics.rds')
