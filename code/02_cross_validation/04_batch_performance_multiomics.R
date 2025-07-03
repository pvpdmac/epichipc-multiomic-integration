#' Obtain performance estimates and save
#'
#' @author C.Shannon


# performance summary
library(tidyverse)

# load cv
experiments <- readr::read_rds('results/other_outputs/diablo_variablecv_top20p_results.rds')

# some CVs did not produce results
i <- experiments$results %>% map(ncol) %>% flatten_dbl()

# only keep splits that produced performance result
experiments <- experiments[i == 3, ]

# only consider cv data
experiments <- experiments %>% unnest(results) %>% filter(type == 'cv') %>% select(-type, results = assess)

# multicore
future::plan(future::multisession, workers = 20)

experiments$tab <- experiments %>%
  furrr::future_pmap(
    ., .progress = T,
    .options = furrr::furrr_options(seed = TRUE),
    function(outcome, results, ...) {
      source('code/utils/helpers_performance_multiomics.R')
      summarize_performance(outcome = outcome, cv = results, alpha = 1-0.95, times = 100) %>%
        mutate(type = 'rho') %>%
        select(type, cvm, cvlo, cvup)
    })

future::plan(future::sequential)
experiments <- experiments %>% select(-results, -meta, -outcome) %>% unnest(tab)

# save
saveRDS(experiments, file = 'results/other_outputs/diablo_variablecv_top20p_results_summarized_rho.rds', compress = F)
