#' Prepare data subsets and hyperparameter grid to train on:
#'
#' * Subset metadata to only include sample with complete/paired omics
#'
#' @author C.Shannon

# imports ---
library(tidyverse)

df <- readr::read_rds('data/processed/cross_validation/data_master_set_top20p.rds')

# metadata: subsets ----

# create subsets by query... in two-steps because day-groups at DOL0
queries <- tibble(
  daygrp = c('ALL', 'DOL1', 'DOL3', 'DOL7', 'DOL1', 'DOL3', 'DOL7'),
  dol = c(0, 0, 0, 0, 1, 3, 7)
)

# cross-reference rownames from df
paired_f <- df %>% filter(!paired) %>% slice(2) %>% pull(data) %>% dplyr::first() %>% rownames()
paired_t <- df %>% filter(paired) %>% slice(2) %>% pull(data) %>% dplyr::first() %>% rownames()
queries <- tibble(
  paired = c(F, T),
  vids = list(paired_f, paired_t)
) %>% crossing(queries, .)

queries <- queries %>%
  crossing(
    grp = c('ALL', 'HBV', 'BCG', 'HBV+BCG'),
    status = c('ALL', 'NAIVE', 'NON-NAIVE')
  )

# ...and subset metadata - dol, grp, status
queries$meta <- pmap(queries, function(daygrp, dol, vids, grp, status, ...) {

  r <- meta_inner %>% filter(`Visit Identifier` %in% vids, DOL == dol)
  if (daygrp != 'ALL')
    r <- filter(r, DayGrp == daygrp)

  if (grp != 'ALL')
    r <- filter(r, VaccineGrp == grp)

  if (status != 'ALL')
    r <- filter(r, abstatus == tolower(status))

  r
})

queries <- queries %>% select(paired, dol, daygrp, grp, status, meta)
queries

rm(meta_inner, pairwise_complete)

# metadata: responses ----

response <- readRDS('data/processed/gam/import_responses.rds') %>% dplyr::rename('Unique Identifier' = `Subject ID`)

responses <- list(
  'log(V4/V1)'    = transmute(response, `Unique Identifier`, y = log(V4) - log(V1)),
  'log(V3/V1)' = transmute(response, `Unique Identifier`, y = log(V3) - log(V1))
)

responses <- tibble(response = names(responses), y = responses)

# remove missing responses
responses <- responses %>% mutate(y = map(y, na.omit))

rm(response)

# hyperparameters ----

# designs
blocks <- unique(df$block)
strength <- c(0, 0.01, 0.1, 1)

designs <- map(strength, ~ {
  d <- matrix(., nrow = length(blocks), ncol = length(blocks))
  diag(d) <- 0

  # ada unconnected
  i <- which(blocks %in% c('ada', 'flow_freq', 'luminex_cytokines', 'luminex_complement'))
  d[i, ] <- 0
  d[, i] <- 0

  d
})

names(designs) <- designs %>% map(as.vector) %>% map(unique) %>% map(last) %>% unlist()

# complexity
# ~~max ncomp should be < min nrow for smallest subset~~
# instead adjust this dynamically when running cross-validation
ncomp <- 20
nkeep <- ncomp * c(5, 10, 50, 100)

hyperparameters <- crossing(
  # ncomp
  ncomp = ncomp,
  # nkeep
  nkeep = nkeep,
  # design
  design = designs
)

# cleanup
rm(blocks, strength, designs)

# finally, combine omics and queries ----

# fix names
names(df$data) <- df$block

# metadata: combine ----

# get rid of subset of queries for dol == 0 and paired == TRUE
queries <- queries %>% filter(!(paired & dol == 0))

experiments <- crossing(queries, responses) %>%
  mutate(meta = map2(meta, y, ~ inner_join(.x, .y))) %>%
  crossing(hyperparameters)

# create friendly 'y' object for diablo
experiments$outcome <- experiments$meta %>%
  map(column_to_rownames, 'Visit Identifier') %>%
  map(select, y) %>%
  map(as.matrix)

# # finally re-order
experiments <- experiments %>% select(paired:status, response, meta, outcome, ncomp:design)

# final cleanup
rm(queries, responses, hyperparameters)

# save subsets grid
file <- sprintf('data/processed/cross_validation/diablo_queries_master_set_%s.rds', suffix)
saveRDS(experiments, file = file, compress = F)
