#' Prepare omics profile to train on:
#'
#' * Only include sample with complete/paired omics
#' * Create both paired and unpaired data matrices
#'
#' @author C.Shannon

# imports ---
library(tidyverse)

suffix <- 'top20p'

# omics ----

# define core set of samples to use in integration

# read-in metadata
meta <- readRDS('data/processed/gam/import_metadata.rds') %>% select(-`Randomization Group`)

# omic-complete samples
df <- tribble(
  ~block, ~data,
  "flow_type",         readRDS("data/processed/gam/import_flow_flowtype.rds"),
  "luminex_cytokines", readRDS("data/processed/gam/import_luminex_cytokines.rds"),
  "proteomics",        readRDS("data/processed/gam/import_proteomics.rds"),
  "metabolomics",      readRDS("data/processed/gam/import_metabolomics.rds"),
  "epigenetics",       readRDS("data/processed/gam/import_epigenetics.rds"),
  "transcriptomics",   readRDS("data/processed/gam/import_transcriptomics.rds")
)

names(df$data) <- df$block

# any NAs
df$data %>% map(anyNA) %>% unlist() %>% all() %>% '!'() %>% stopifnot()

# define inner product of the rownames across all data
inner <- df$data %>%
  map(rownames) %>%
  reduce(intersect)

# intersection of samples
meta_inner <- meta %>%
  arrange(`Unique Identifier`, DOL) %>%
  filter(`Visit Identifier` %in% inner)

# get rid of delayed group
meta_inner <- filter(meta_inner, VaccineGrp != 'Delayed')

# define ab naive subset
df_responses <- readRDS('data/processed/gam/import_responses.rds')

meta_inner <- df_responses %>%
  mutate(abstatus = factor(V1 <= 2.5, levels = c(F, T), labels = c('non-naive', 'naive'))) %>%
  select(`Unique Identifier` = `Subject ID`, abstatus) %>%
  inner_join(meta_inner, .)

# cleanup
rm(meta, inner, df_responses)

# # print size
# meta_inner %>%
#   count(VaccineGrp, DOL, abstatus) %>%
#   spread(DOL, n) %>%
#   split(.$abstatus)

# finally, subset omics
df$data <- map(df$data, ~ .[meta_inner$`Visit Identifier`, ])

# sanity - check order of samples
df$data %>%
  map(~ all(rownames(.) == meta_inner$`Visit Identifier`)) %>%
  unlist() %>%
  all()

# add withinVar adjusted
df$data_within <- df$data

# only pairwise complete samples
pairwise_complete <- meta_inner %>%
  group_by(`Unique Identifier`) %>%
  filter(n() > 1) %>%
  ungroup()

# withinVar
df$data_within <- map(
  df$data_within,
  ~ mixOmics::withinVariation(X = .[pairwise_complete$`Visit Identifier`, ], design = data.frame(pairwise_complete$`Unique Identifier`))
)

# collapse to merge later
df <- df %>%
  gather(paired, data, starts_with('data')) %>%
  mutate(paired = grepl('within$', paired))

# simple variance filter - top 20%
df$data <- map(df$data, ~ {

  if(ncol(.) > 2000) {

    var <- matrixStats::colMads(.)
    .[ , var >= quantile(var, probs = 0.8)] # top20p

  } else {
    .
  }

})

# save data grid
file <- sprintf('data/processed/cross_validation/data_master_set_%s.rds', suffix)
saveRDS(df, file = file, compress = F)
