#' GAM Validation using samples with incomplete omic profiles
#'
#' @author C.Shannon

# imports
library(tidyverse)
ggplot2::theme_set(cowplot::theme_cowplot() + ggplot2::theme(strip.background = ggplot2::element_blank()))

# identify samples with incomplete omic profiles
df <- tribble(
  ~block, ~data,
  'flow_type',         readRDS('data/processed/gam/import_flow_flowtype.rds'),
  'luminex_cytokines', readRDS('data/processed/gam/import_luminex_cytokines.rds'),
  'proteomics',        readRDS('data/processed/gam/import_proteomics.rds'),
  'metabolomics',      readRDS('data/processed/gam/import_metabolomics.rds'),
  'epigenetics',       readRDS('data/processed/gam/import_epigenetics.rds'),
  'transcriptomics',   readRDS('data/processed/gam/import_transcriptomics.rds')
)

names(df$data) <- df$block
all <- df$data %>% map(rownames) %>% reduce(intersect)

# only considering transcriptomics - too few available samples otherwise
df <- df %>% filter(block %in% c('transcriptomics'))

# compute set difference
keepers <- df$data %>% map(rownames) %>% reduce(intersect) %>% setdiff(all)

# add metadata
df$meta <- list(
  readr::read_rds('data/processed/gam/import_metadata.rds')  %>%
    select(uid = `Unique Identifier`, vid = `Visit Identifier`, sex = Sex, grp = VaccineGrp, dol = DOL) %>%
    filter(vid %in% keepers) %>%
    na.omit()
)

# subset
df$data <- df %>% pmap(function(data, meta, ...) data[meta$vid, ])

# split by dol
df_unpaired <- df %>%
  crossing(tibble(dol = c(0, 1, 3, 7))) %>%
  select(block, dol, meta, data) %>%
  pmap(function(block, meta, data, dol) {

    meta <- meta %>% filter(dol == !!dol)
    data <- data[meta$vid, ]
    tibble(paired = FALSE, block = block, dol = dol, meta = list(meta), data = list(data))
  }) %>% bind_rows()

# paired
df_paired <- df %>%
  crossing(tibble(dol = c(1, 3, 7))) %>%
  pmap(function(block, dol, meta, data, ...) {
    # withinvar
    y <- meta %>% filter(vid %in% rownames(data), dol == !!dol | dol == 0) %>% group_by(uid) %>% filter(n() > 1) %>% ungroup()
    x <- data[y$vid, ]
    x <- mixOmics::withinVariation(x, design = data.frame(y$uid))
    y <- y %>% filter(dol != 0)
    tibble(paired = TRUE, dol = !!dol, block = block, meta = list(y), data = list(x[y$vid, ]))
  }) %>% bind_rows()

# write-out
saveRDS(bind_rows(df_unpaired, df_paired),  file = 'data/processed/process_validation/data_gam_incomplete_omic_profiles.rds', compress = F)
