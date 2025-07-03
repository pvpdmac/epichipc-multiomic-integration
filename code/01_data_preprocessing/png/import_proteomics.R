#' Read-in and prepare proteomics data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-05 - Initial commit
#' 2023-12-15 - Fixed colnames to match gambia

# import
library(tidyverse)

# read-in
data <- readr::read_csv('data/raw/png/PNGVAL_SINGLEOMICS_PRT_COUNTS_V1V2.csv')

data <- data %>%
  separate(sample_id, c(NA, 'VID', NA)) %>%
  head(-4) %>%
  column_to_rownames('VID') %>%
  as.matrix()

# impute ------------------------------------------------------------------

# missing value filter
df_long <- data %>%
  as_tibble(rownames = 'VID') %>%
  gather(protein, conc, everything(), -VID)

df_missing <- df_long %>%
  group_by(protein) %>%
  summarize(n_missing = sum(is.na(conc))/n()) %>%
  ungroup()

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
  x[is.na(x)] <- imp/2
  x
})

data_sub_impute <- log(data_sub_impute)

# write out ---------------------------------------------------------------

saveRDS(data_sub_impute, 'data/processed/png/import_proteomics.rds')
