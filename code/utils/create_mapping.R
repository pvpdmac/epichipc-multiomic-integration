#' Create mapping tables to hgnc_symbol for reporting/interpretation:
#' @author C.Shannon
#'

# imports ----
library(biomaRt)
library(stringdist)
library(tidyverse)

# helpers ----

# get feature names ----
features <- list(
  'data/import_luminex_cytokines.rds',
  'data/import_proteomics.rds',
  'data/import_transcriptomics.rds',
  'data/import_epigenetics.rds',
  'data/import_metabolomics.rds'
)

features <- features %>% map(readRDS) %>% map(colnames)
names(features) <- c('luminex_cytokines', 'proteomics', 'transcriptomics', 'epigenetics', 'metabolomics')

# transcriptomics

# map features to hgnc_symbol for various data types using biomaRt
ensembl <- biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
table_transcriptomics <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)

features$transcriptomics <- tibble(ensembl_gene_id = features$transcriptomics) %>% left_join(table_transcriptomics)

# luminex
table_luminex <- tribble(
  ~luminex_id, ~hgnc_symbol,
  'EGF',          'EGF',
  'Eotaxin',      'CCL11',
  'FGF-2',        'FGF2',
  'Flt-3L',       'FLT3LG',
  'Fractalkine',  'CX3CL1',
  'G-CSF',        'CSF3',
  'GM-CSF',       'CSF2',
  'GRO',          'CXCL1',
  'IFNa2',        'IFNA2',
  'IFNg',         'IFNG',
  'IL-10',        'IL10',
  'IL-12P40',     'IL12B',
  'IL-12P70',     'IL12A',
  'IL-13',        'IL13',
  'IL-15',        'IL15',
  'IL-17A',       'IL17A',
  'IL-1a',        'IL1A',
  'IL-1b',        'IL1B',
  'IL-1RA',       'IL1RN',
  'IL-2',         'IL2',
  'IL-3',         'IL3',
  'IL-4',         'IL4',
  'IL-5',         'IL5',
  'IL-6',         'IL6',
  'IL-7',         'IL7',
  'IL-8',         'CXCL8',
  'IL-9',         'IL9',
  'IP-10',        'CXCL10',
  'MCP-1',        'CCL2',
  'MCP-3',        'CCL7',
  'MDC',          'ADAM11',
  'MDC',          'CCL22',
  'MIP-1a',       'CCL3',
  'MIP-1b',       'CCL4',
  'PDGF-AA',      'PDGFA',
  'PDGF-AB/BB',   'PDGFB',
  'RANTES',       'CCL5',
  'sCD40L',       'CD40LG',
  'TGF-a',        'TGFA',
  'TNFa',         'TNF',
  'TNFb',         'LTA',
  'VEGF',         'VEGFA'
)

features$luminex_cytokines <- tibble(luminex_id = features$luminex) %>% left_join(table_luminex)

# proteomics
table_uniprot <- biomaRt::getBM(attributes = c('uniprot_gn_id', 'hgnc_symbol'), mart = ensembl)

features$proteomics <- features$proteomics %>%
  tibble(feature_name = .) %>%
  mutate(uniprot_gn_id = gsub('^sp\\|(.+)\\|(.+)_HUMAN$', '\\1', feature_name)) %>%
  left_join(table_uniprot) %>%
  dplyr::select(feature_name, hgnc_symbol)

# fix a few nas manually
features$proteomics <- features$proteomics %>%
  mutate(
    hgnc_symbol = ifelse(
      is.na(hgnc_symbol),
      gsub('^sp\\|(.+)\\|(.+)_HUMAN$', '\\2', feature_name),
      hgnc_symbol
    )
  )

# epigenetics
env <- new.env()
load(
  envir = env,
  file = 'data/raw/gam/GAMMAIN_MULTIOMICS_EPI.RData'
)

features$epigenetics <- env$map2genome %>%
  as_tibble() %>%
  dplyr::select(Name, UCSC_RefGene_Name)

# metabolomics
metamap <- readr::read_delim('data/raw/gam/GAMMAIN_MULTIOMICS_MET_ROWFEATURE_WITHXENO.txt')

features$metabolomics <- metamap %>%
  as_tibble() %>%
  dplyr::select(MET_ID, PLOT_NAME)

# colnames
features <- features %>% purrr::map(dplyr::select, feature_name = 1, hgnc_symbol = 2)

# container
features <- bind_rows(features, .id = 'block') %>% nest(map = c(feature_name, hgnc_symbol))

saveRDS(features, 'code/misc/features_to_hgnc_symbols.rds')
