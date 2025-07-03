#' Run cross-validation:
#'
#' @author C.Shannon

params_cores   <- 24
params_nfold   <- 5
params_nrepeat <- 20

# imports
library(tidyverse)
source('code/utils/helpers_diablo.R')

# plots
ggplot2::theme_set(cowplot::theme_cowplot() + theme(strip.background = element_blank()))

# helper
scor_vec <- function(truth, estimate, na_rm = TRUE, ...) {
  scor_impl <- function(truth, estimate) {
    cor(truth, estimate, method='spearman')
  }
  yardstick::metric_vec_template(
    metric_impl = scor_impl,
    truth = truth,
    estimate = estimate,
    na_rm = na_rm,
    cls = "numeric",
    ...
  )
}
scor <- function(data, ...) {
  UseMethod("scor")
}
scor.data.frame <- function(data, truth, estimate, na_rm = TRUE, ...) {
  yardstick::metric_summarizer(
    metric_nm = "scor",
    metric_fn = scor_vec,
    data = data,
    truth = !! dplyr::enquo(truth),
    estimate = !! dplyr::enquo(estimate),
    na_rm = na_rm,
    ...
  )
}
class(scor) <- c("numeric_metric", class(scor))
attr(scor, "direction") <- "maximize"

# START ----
experiments <- readRDS('data/processed/cross_validation/diablo_queries_master_set_top20p.rds')

# SOME FIXES

# define max comp based on subset n
experiments <- experiments %>%
  mutate(
    n = map(meta, nrow) %>% unlist(),
    # +2 to account for ill-defined ncomp - 1 and loocv - 1
    ncomp = ifelse(n <= max(ncomp) + 2, n - 3, ncomp),
    nkeep = ifelse(n <= max(ncomp) + 2, nkeep / max(ncomp) * (n -3), nkeep)
  )

# define type of CV based on subset n
experiments <- experiments %>%
  mutate(
    n = map(meta, nrow) %>% unlist(),
    cv_type = ifelse(n < 40, 'loocv', 'kfold')
  )

# DOL0 -> unpaired
experiments <- experiments %>% filter(!(paired & dol == 0))

# sanity
experiments %>% distinct(dol, daygrp, grp, paired, response, status)

# shuffle to run multicore
experiments <- experiments %>% sample_frac()

# Run integration ----

# finally, load globals
omics <- readRDS('data/processed/cross_validation/data_master_set_top20p.rds')

# multicore
future::plan(future::multisession, workers = params_cores)
furrr::furrr_options(scheduling = Inf)
# allow for omics to be shared per futures
# omics %>% object.size() %>% format(unit = 'B')
options(future.globals.maxSize = 2.1e9)

# final model + cross-validation performance
experiments$results <- experiments %>%
  furrr::future_pmap(
    .,
    .progress = T,
    .options = furrr::furrr_options(seed = TRUE),

    function(data, outcome, nkeep, ncomp, design, cv_type, ...) {

      RhpcBLASctl::blas_set_num_threads(1)

      # pull in correct data and subset
      data <- omics %>%
        dplyr::select(block, paired, mat = data) %>%
        tidyr::nest(data = c(block, mat)) %>%
        dplyr::filter(paired == list(...)$paired)%>%
        dplyr::pull(data) %>%
        dplyr::first() %>%
        dplyr::pull(mat) %>%
        purrr::map(~ .[rownames(outcome), ])

      # make outcome dataframe
      df <- outcome %>%
        tidyr::as_tibble(rownames = 'rownames') %>%
        dplyr::rename(actual = y)

      # fit model
      diablo <- purrr::possibly(make_fit, otherwise = NA)(
        X = data,
        Y = outcome,
        ncomp = ncomp,
        nkeep = nkeep,
        design = design
      )

      # training performance
      library(mixOmics)
      predictions <- make_predict(fit = diablo, newx = data)

      # per component
      perf_training <- predictions %>%
        dplyr::right_join(df, .) %>%
        tidyr::nest(assess = c(rownames, actual, prediction))

      # cross-validation perfomance
      perf_cv <- purrr::possibly(make_kfoldcv, otherwise = NA)(
        X = data,
        Y = outcome,
        ncomp = ncomp,
        nkeep = nkeep,
        design = design,
        nfold = ifelse(cv_type == 'loocv', nrow(outcome), params_nfold),
        nrepeat = ifelse(cv_type == 'loocv', 1, params_nrepeat)
      )

      list(training = perf_training, cv = perf_cv) %>% dplyr::bind_rows(.id = 'type')

    })

future::plan(future::sequential)

# SAVE ----

saveRDS(experiments, file = 'results/other_outputs/diablo_variablecv_top20p_results.rds', compress = FALSE)

print("Done!")
