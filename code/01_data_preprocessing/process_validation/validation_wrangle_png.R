#' PNG Validation including batch-correction
#'
#' @author C.Shannon

# imports
library(tidyverse)

# load omics and compute feature overlaps
df <- tribble(
  ~block,              ~source,
  'epigenetics',       'gam',
  'luminex_cytokines', 'gam',
  'proteomics',        'gam',
  'transcriptomics',   'gam',
  'epigenetics',       'png',
  'luminex_cytokines', 'png',
  'proteomics',        'png',
  'transcriptomics',   'png'
)

# construct path from vars
df$path <- map2_chr(df$block, df$source, ~ file.path('data/processed', .y, sprintf('import_%s.rds', .x)))
df$path %>% map_lgl(file.exists) %>% all()

# read in data
df$data <- df$path %>% map(readr::read_rds)

# NA cols?
df$data %>% map(matrixStats::colAnyNAs) %>% map(table)

# resize
df <- df %>%
  select(block, source, data) %>%
  spread(source, data)

# overlap with model features
source('code/utils/helpers_diablo.R')

# read in model tibble
mods <- readr::read_rds('data/processed/our_favourite_models.rds')

# simplify to features per block
mods <- mods %>%
  pmap(function(comp, diablo, ...) {
    lim <- gsub('comp', '', comp) %>% as.numeric()
    extract_loadings(diablo, return_symbols = F) %>%
      filter(comp <= lim)
  }) %>%
  set_names(paste(mods$status, mods$response)) %>%
  bind_rows(.id = 'model') %>%
  select(-comp, -coefficient) %>%
  distinct() %>%
  nest(features = -block)

# identify overlapping samples between gam and png
df$overlap <- df %>%
  inner_join(mods) %>%
  pmap(function(gam, png, features, ...) {
    l <- list()
    l$model <- features %>% split(.$model) %>% map(pull, feature_name) %>% reduce(union)
    l$gam <- colnames(gam)
    l$png <- colnames(png)
    l$intersection <- intersect(colnames(gam), colnames(png))
    l
  })

# gam - subset to common vars
df$a <- df %>% pmap(function(gam, overlap, ...) gam[ , overlap$intersection])

# png - subset to common vars
df$b <- df %>% pmap(function(png, overlap, ...) png[ , overlap$intersection])

# combine
df$data <- df %>% pmap(function(a, b, ...) rbind(a, b))
df <- df %>% select(block, gam, png, overlap, data)

# add metadata
df$meta <- list(
  bind_rows(
    readr::read_rds('data/processed/gam/import_metadata.rds')  %>%
      transmute(source = 'gam',  uid = `Unique Identifier`, vid = `Visit Identifier`, sex = Sex, grp = VaccineGrp, dol = DOL) %>%
      na.omit(),

    readr::read_rds('data/processed/png/import_metadata.rds') %>%
      transmute(source = 'png', uid = Unique_Identifier, vid = VID, sex = Sex, grp = VaccineGroup, dol = DOL) %>%
      na.omit()
  )
)

# batch correction
#
# 1. Carried out using the common set of variables only.
# 2. Removing the effect of `source` (GAM/PNG).
# 3. Additionally, retaining variance related to: randomization groups (vaccine grp, dol grp), time (dol), sex
df$data_combat <- df %>%
  pmap(function(data, meta, ...) {
    mod <- meta %>% filter(vid %in% rownames(data))
    dat <- t(data[mod$vid, ])
    combat <- sva::ComBat(dat, batch = mod$source, ref.batch = 'gam', mod = model.matrix(~ dol + grp + sex, data = mod)) %>% t()
  })

# only validation samples
df <- df %>%
  pmap(function(block, meta, data_combat, ...) {
    m <- meta %>% filter(source == 'png') %>% filter(vid %in% rownames(data_combat))
    tibble(block = block, meta = list(m), data = list(data_combat[m$vid, ]))
  }) %>% bind_rows()

# unpaired
df_unpaired <- df %>%
  crossing(tibble(dol = c(0, 7))) %>%
  select(block, dol, meta, data) %>%
  pmap(function(block, meta, data, dol) {
    # subset to dol
    meta <- meta %>% filter(dol == !!dol)
    data <- data[meta$vid, ]
    tibble(paired = FALSE, block = block, dol = dol, meta = list(meta), data = list(data))
  }) %>% bind_rows()

# paired
df_paired <- df %>%
  crossing(tibble(dol = c(7))) %>%
  pmap(function(block, dol, meta, data, ...) {
    # withinvar
    y <- meta %>% filter(vid %in% rownames(data)) %>% group_by(uid) %>% filter(n() > 1) %>% ungroup()
    x <- data[y$vid, ]
    x <- mixOmics::withinVariation(x, design = data.frame(y$uid))
    y <- y %>% filter(dol != 0)
    tibble(paired = TRUE, dol = !!dol, block = block, meta = list(y), data = list(x[y$vid, ]))
  }) %>% bind_rows()

# write-out
saveRDS(bind_rows(df_unpaired, df_paired), file = 'data/processed/process_validation/data_png_combat_corrected.rds', compress = F)
