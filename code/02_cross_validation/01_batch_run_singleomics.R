#' Run cross-validation:
#'
#' @author C.Shannon

params_cores   <- 24
params_nfold   <- 10
params_nrepeat <- 10

# imports
library(tidyverse)

# START ----
experiments <- readRDS('data/processed/cross_validation/glmnet_queries_master_set_top20p.rds')

# SOME FIXES

# define type of CV based on subset n
experiments <- experiments %>%
  mutate(
    n = map(meta, nrow) %>% unlist(),
    cv_type = ifelse(n < 40, 'loocv', 'kfold')
  )

# DOL0 -> unpaired
experiments <- experiments %>% filter(!(paired & dol == 0))

# sanity
# experiments %>% distinct(dol, daygrp, grp, paired, response, status)

# shuffle to run multicore
experiments <- experiments %>% sample_frac()

# run cross-validation ----

# single-omics
omics <- readRDS('data/processed/cross_validation/data_master_set_top20p.rds')
experiments <- experiments %>% crossing(block = omics$block)

# multicore
future::plan(future::multisession, workers = params_cores)
furrr::furrr_options(scheduling = Inf)

# allow for omics to be shared per futures
options(future.globals.maxSize = 4e9)

# helper
make_stratfolds <- function(y, nfold = 10, nrepeat = 10) {

  folds <- y %>%
    tibble::as_tibble(rownames = 'rownames') %>%
    rsample::vfold_cv(v = nfold, repeats = nrepeat, strata = y) %>%
    tidyr::nest(repeats = -id)

  folds$vector <- folds$repeats %>%
    map(~ {
      vec <- .x$splits %>%
        map(rsample::assessment) %>%
        map(pull, 1) %>%
        set_names(1:nrow(.x)) %>%
        map(tibble) %>%
        bind_rows(.id = 'fold')

      ret <- as.numeric(vec[[1]])
      names(ret) <- vec[[2]]
      ret
      ret <- ret[rownames(y)]
    })

  folds %>% select(id, vector)

}

# final model + cross-validation performance
experiments$results <- experiments %>%
  furrr::future_pmap(
    .,
    .progress = T,
    .options = furrr::furrr_options(seed = TRUE),

    function(outcome, alpha, cv_type, paired, block, ...) {

      RhpcBLASctl::blas_set_num_threads(1)
      nfold <- params_nfold
      nrepeat <- params_nrepeat

      # subset data
      x <- tibble(paired = paired, block = block) %>%
        left_join(omics, by = c('paired', 'block')) %>%
        pull(data) %>%
        first() %>%
        '['(rownames(outcome), )

      # create folds and fit; loocv is special case of (repeated) kfold
      if(cv_type == 'loocv') {

        print('loocv')
        cv <- glmnet::cv.glmnet(
          x = x, y = outcome[[1]],
          alpha = alpha,
          nfolds = nrow(outcome),
          alignment = 'fraction',
          grouped = F,
          keep = T
        )

        # container
        cv <- tibble(id = 'Repeat00', cv = list(cv))

      } else {

        print('kfoldcv')
        cv <- outcome %>% make_stratfolds(nfold = params_nfold, nrepeat = params_nrepeat)
        cv$cv <- cv$vector %>%
          map(~ {
            glmnet::cv.glmnet(
              x = x, y = outcome[[1]],
              alpha = alpha,
              nfolds = params_nfold,
              foldid = .,
              alignment = 'fraction',
              grouped = F,
              keep = T
            )
          })

        cv <- cv %>% select(id, cv)
      }

      return(cv)
})

future::plan(future::sequential)

# SAVE ----
saveRDS(experiments, file = 'results/other_outputs/glmnet_variablecv_top20p_results.rds', compress = FALSE)

print("Done!")
