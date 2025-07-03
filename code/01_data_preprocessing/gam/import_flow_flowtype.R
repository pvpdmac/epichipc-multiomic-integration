#' Read-in and prepare auto *high-dimensional* flow cytometry data:
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

# read in - flow
f <- list(
  bcells  = readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_BCELLFREQS.csv'),
  myeloid = readr::read_csv('data/raw/gam/GAMMAIN_SINGLEOMICS_MYELOIDFREQS.csv')
)

# hard coded
f$myeloid <- f$myeloid %>% select(-root_live)

# transform to matrix
g <- f %>%
  map(~ {
    .$samples_name <- gsub('^FCT_(G.+)_(.+)\\.fcs$', '\\2', .$samples_name)
    column_to_rownames(., 'samples_name') %>%
      as.matrix()
  })

# check visit ids
g %>%
  map(rownames) %>%
  map(~ all(. %in% meta$`Visit Identifier`))

# re-order
g <- g %>%
  map(~ .[intersect(meta$`Visit Identifier`, rownames(.)), ])

# distinct feature names?
g %>% map(colnames) %>% reduce(intersect)

# common observations?
common <- g %>% map(rownames) %>% reduce(intersect)

# combine
data <- g %>% map(~ .x[common, ]) %>% reduce(cbind)

rm(f, g)

# CLR transform
data <- mixOmics::logratio.transfo(data, logratio = 'CLR', offset = min(data[data > 0])/2)

# save
saveRDS(data, 'data/processed/gam/import_flow_flowtype.rds')
