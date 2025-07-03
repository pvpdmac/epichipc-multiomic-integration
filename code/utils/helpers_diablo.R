#' Helpers for running DIABLO integration analysis
#' @author C.Shannon
#'
#' Changelog:
#' 2021/07/07 - re-write k-fold CV to use rsample (support for stratified cv)
#' 2021/01/04 - update to LOOCV (10x10-fold = 100 models... we often have fewer samples than that...)
#' 2020/12/24 - update to evaluate performance at all ncomps, not just max ncomp...
#' 2020/12/12 - generalize to allow specification of design up front.

# wrappers ----------------------------------------------------------------

#' create all designs
#' hint: produce all possible graphs given n vertices...
#' to create random weights could mat multiply by matrix of weights size blocks x blocks, over range of values
#' e.g. sample(seq(0.1, 0.9, by = 0.1), size = 25, replace = T)
make_designs <- function(X) {
  nvertices <- length(X)
  graphs <- 2:nvertices %>%
    purrr::map(~ {
      utils::combn(nvertices, ., simplify = F) %>%
        purrr::map(~ utils::combn(., 2, simplify = T)) %>%
        purrr::map(igraph::graph, n = nvertices, directed = F)
    }) %>%
    purrr::flatten()

  adj <- graphs %>%
    purrr::map(igraph::as_adjacency_matrix) %>%
    purrr::map(as.matrix)

  c(list(matrix(0, nrow = nvertices, ncol = nvertices)), adj)
}

#' create keepX
make_nkeep <- function(X, ncomp, nkeep) {
  map(X, ~ {
    upper <- floor(min(nkeep, ncol(.))/ncomp)
    rep(upper, ncomp)
  })
}

#' fit - dispatching on type of Y
make_fit <- function(X, Y, ncomp, nkeep, design) {

  # dispatch on type of response
  if(is.numeric(Y))
    fit <- purrr::quietly(mixOmics::block.spls)(
      X = X,
      Y = Y,
      ncomp = ncomp,
      keepX = make_nkeep(X, ncomp, nkeep),
      design = design
    )
  else
    fit <- purrr::quietly(mixOmics::block.splsda)(
      X = X,
      Y = Y[ , 1],
      ncomp = ncomp,
      keepX = make_nkeep(X, ncomp, nkeep),
      design = design
    )

  # return (quiet) result
  fit$result
}

#' predict
make_predict <- function(fit, newx) {

  # call mixOmics predict function quietly!
  suppressPackageStartupMessages(library(mixOmics))
  pred <- predict(object = fit, newdata = newx)
  detach("package:mixOmics", unload = TRUE)
  detach("package:MASS", unload = TRUE)
  ncomp <- fit$ncomp[[1]]

  # dispatch on response - use WeightedPredict for regression; use WeightedVote
  # and Mahalanobis for classification
  if(is.null(fit$Y)) {
    x <- pred$WeightedPredict[ , , , drop = F] %>%
      as_tibble(rownames = 'rownames')
  } else {
    x <- pred$WeightedVote$mahalanobis.dist[ , , drop = F] %>%
      as_tibble(rownames = 'rownames')
  }

  # cleanup and gather
  colnames(x) <- c('rownames', paste0('comp', 1:ncomp))
  x <- x %>% gather(comp, prediction, -rownames)

}

# make cross-validation

#' repeated k-fold cv
make_kfoldcv <- function(X, Y, ncomp, nkeep, design, nfold, nrepeat = NULL) {

  # sanity
  if(any(is.na(Y))) {
    warning(sprintf('response (Y) vector includes NAs; observations (n = %d) will be removed.', length(is.na(Y))))
    X <- map(X, ~ .[!is.na(Y), ])
    Y <- Y[!is.na(Y), , drop = F]
    nfold <- min(nfold, nrow(Y))
  }

  if(nfold == length(Y) & !is.null(nrepeat)) {
    warning('loocv scheme selected; nrepeat will be ignored.')
  }

  if(nfold > length(Y)) {
    stop('nfold must be less than number of available samples.')
  }

  # construct folds

  # create folds; loocv is special case of (repeated) kfold
  if(nfold == length(Y)) {
    folds <- Y %>%
      tibble::as_tibble(rownames = 'rownames') %>%
      rsample::vfold_cv(v = nfold, repeats = 1)
    class(folds)[1] <- 'loo_cv'
  } else {
    folds <- Y %>%
      tibble::as_tibble(rownames = 'rownames') %>%
      rsample::vfold_cv(v = nfold, repeats = nrepeat, strata = y)
  }

  # train in 'analysis', test in 'assessment'
  folds$assess <- purrr::map(folds$splits, ~ {

    # train
    train <- rsample::analysis(.)
    X_p <- purrr::map(X, ~ .[train$rownames, , drop = F])
    Y_p <- Y[train$rownames, , drop = F]
    fit <- make_fit(X = X_p, Y = Y_p, ncomp, nkeep, design)

    # test - at each comp
    test <- rsample::assessment(.) %>% dplyr::select(rownames, actual = y)
    pred <- make_predict(fit = fit, newx = X) %>% inner_join(test, by = 'rownames')

  })

  # reformat and return
  folds %>%
    dplyr::select(starts_with('id'), assess) %>%
    tidyr::unnest(assess) %>%
    tidyr::nest(assess = c(starts_with('id'), 'rownames', 'actual', 'prediction')) %>%
    dplyr::mutate(assess = map(assess, group_by, id)) %>%
    dplyr::mutate_if(.predicate = is.character, .funs = tolower)
}

#' get feature loadings per component, return tibble of feature scores, by component and block
extract_loadings <- function(fit, return_symbols = T) {

  # create iterator over ncomps
  comps <- 1:dplyr::first(fit$ncomp)
  # names(comps) <- paste0('comp', comps)

  # map over and pull out variables
  vars <- comps %>%
    purrr::map(~ {
      mixOmics::selectVar(fit, comp = .) %>%
        purrr::map(~ {
          purrr::safely(tibble, otherwise = NA)(feature_name = .$name, coefficient = .$value$value.var)$result
        }) %>%
        head(-2) %>%
        dplyr::bind_rows(.id = 'block')
    })

  # combine and return
  vars <- vars %>%
    dplyr::bind_rows(.id = 'comp') %>%
    dplyr::mutate(comp = as.integer(comp)) %>%
    dplyr::arrange(comp, block) %>%
    dplyr::select(block, feature_name, comp, coefficient)

  if(return_symbols) {
    vars <- readRDS(here::here('/home/cshannon/EPIC-HIPC P001 v2.0/code/misc/features_to_hgnc_symbols.rds')) %>%
      unnest(map) %>%
      dplyr::left_join(vars, .) %>%
      dplyr::mutate(hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == '', sprintf('(%s)', feature_name), hgnc_symbol)) %>%
      dplyr::select(block, feature_name, hgnc_symbol, comp, coefficient)
  }

  vars
}

#' get features per component, return as a tibble of components, blocks, and character vectors of feature names
extract_features <- function(fit) {

  # create iterator over ncomps
  comps <- 1:dplyr::first(fit$ncomp)
  # names(comps) <- paste0('comp', comps)

  # map over and pull out variables
  vars <- comps %>%
    purrr::map(~ {
      mixOmics::selectVar(fit, comp = .) %>%
        purrr::map('name') %>%
        head(-1) %>%
        purrr::map(~ dplyr::tibble(feature = list(.))) %>%
        dplyr::bind_rows(.id = 'block')
    })

  # combine and return
  vars %>%
    dplyr::bind_rows(.id = 'component') %>%
    dplyr::mutate(component = as.integer(component)) %>%
    dplyr::arrange(component, block)
}

#' map features to hgnc_symbol
map_features <- function(vars) {

  tabs <- readRDS(here::here('code/misc/features_to_hgnc_symbols.rds'))
  return <- tabs %>% dplyr::left_join(vars, .)

  return$hgnc_symbol <- pmap(return, function(feature, map, ...) {
    if(is.null(map))
      feature
    else {
      map %>%
        dplyr::filter(feature_name %in% feature) %>%
        dplyr::filter(hgnc_symbol != '') %>%
        dplyr::pull(hgnc_symbol)
    }
  })

  return %>% dplyr::select(component, block, feature, hgnc_symbol)
}
