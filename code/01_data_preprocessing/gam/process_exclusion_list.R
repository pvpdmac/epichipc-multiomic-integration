#' Read-in and prepare metadata:
#' @author: CShannon
#'
#'exclude samples from participants who were hospitalized w/in 30 days or had
#'evidence of sex-mismatch in epi data QC

# metadata
meta <- readr::read_rds('data/processed/gam/import_metadata.rds')

# exclusion list
excl <- readr::read_csv('data/raw/gam/exclusion_list.csv') %>%
  inner_join(meta) %>%
  drop_na() %>%
  select(`Unique Identifier`, `Visit Identifier`)

# titers ----

# special case - titers are identified at the participant level
readr::read_rds('data/processed/gam/import_responses.rds') %>%
  filter(!`Subject ID` %in% excl$`Unique Identifier`) %>%
  saveRDS('data/processed/gam/import_responses.rds')

# omics ----

# target processed files
tbl <- tidyr::tribble(
  ~path, ~id,
  'data/processed/gam/import_epigenetics.rds',       'Visit Identifier',
  'data/processed/gam/import_transcriptomics.rds',   'Visit Identifier',
  'data/processed/gam/import_proteomics.rds',        'Visit Identifier',
  'data/processed/gam/import_luminex_cytokines.rds', 'Visit Identifier',
  'data/processed/gam/import_metabolomics.rds',      'Visit Identifier',
  'data/processed/gam/import_flow_flowtype.rds',     'Visit Identifier',
)

# read in processed data sets
tbl <- tbl %>% mutate(data = map(path, readr::read_rds))

# remove excl
tbl$data_sub <- tbl %>%
  pmap(\(data, id, ...) {
    keep <- setdiff(rownames(data), pull(excl, !!id))
    data[keep, ]
  })

# write out
pwalk(tbl, \(path, data_sub, ...) saveRDS(data_sub, file = path))
