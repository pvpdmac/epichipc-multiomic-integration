# Performance Summary
library(tidyverse)

# load cv
experiments <- readr::read_rds('results/other_outputs/glmnet_variablecv_top20p_results.rds')

experiments <- experiments[1:20, ]

# summarize to spearman's rho
pb <- progress::progress_bar$new(
  format = "  computing [:bar] :percent eta: :eta",
  total = nrow(experiments),
  width = 150,
  clear = FALSE
)

experiments$tab <- experiments %>%
  purrr::pmap(
    function(outcome, results, ...) {
      pb$tick()
      source('code/utils/helpers_performance_singleomics.R')
      summarize_performance(outcome = outcome, cv = results, alpha = 1-0.95, times = 100) %>%
        mutate(type = 'rho') %>%
        select(type, cvm, cvlo, cvup)
    })

experiments %>%
  select(-results, -meta, -outcome) %>%
  unnest(tab)
  readr::write_rds('results/other_outputs/glmnet_variablecv_top20p_results_summarized_rho.rds')
