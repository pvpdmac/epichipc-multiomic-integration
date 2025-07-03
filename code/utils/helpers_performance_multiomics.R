#' Helpers for evaluating CV performance for the single-omic models
#' @author C.Shannon
#'
#' Changelog:
#' 2022/11/15 - initial commit

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

summarize_performance <- function(outcome, cv, alpha = 0.05, times = 1000) {

  # bootstrap
  lim <- alpha/2
  cv %>%
    dplyr::rename(truth = actual, estimate = prediction) %>%
    rsample::bootstraps(strata = truth, times = times) %>%
    dplyr::mutate(
      .estimate = splits %>%
        purrr::map(dplyr::as_tibble) %>%
        purrr::map(scor, truth = truth, estimate = estimate) %>%
        purrr::map_dbl('.estimate')
    ) %>%
    dplyr::summarise(
      cvm  = mean(.estimate, na.rm = T),
      cvlo = quantile(.estimate, na.rm = T, probs = lim),
      cvup = quantile(.estimate, na.rm = T, probs = 1-lim)
    )

}

summarize_performance_auc <- function(cv, alpha = 0.05, times = 1000) {

  # bootstrap
  lim <- alpha/2
  cv %>%
    rename(truth = actual, estimate = prediction) %>%
    mutate(truth = factor(truth > median(truth))) %>%
    rsample::bootstraps(strata = truth, times = times) %>%
    dplyr::mutate(
      .estimate = splits %>%
        map(as_tibble) %>%
        map_dbl(~ precrec::evalmod(scores = .x$estimate, labels = .x$truth, mode = 'aucroc')$uaucs$aucs)
    ) %>%
    dplyr::summarise(
      cvm  = mean(.estimate, na.rm = T),
      cvlo = quantile(.estimate, na.rm = T, probs = lim),
      cvup = quantile(.estimate, na.rm = T, probs = 1-lim)
    )
}

safe_rho <- purrr::possibly(summarize_performance, otherwise = NA_real_)

safe_auc <- purrr::possibly(summarize_performance_auc, otherwise = NA_real_)
