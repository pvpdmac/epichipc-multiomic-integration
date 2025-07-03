#' Read-in and prepare epigenetic data:
#' @author: CShannon
#'
#' Changelog:
#'
#' 2023-11-05 - Initial commit

# import
library(tidyverse)

# load prepared object to temporary environment
env <- new.env()
load(
  envir = env,
  file = 'data/raw/png/PNGVAL_SINGLEOMICS_EPI_DATA_V1V2.RData'
)

# relabel using visit identifiers
tmp <- env$metadata %>%
  as_tibble() %>%
  select(Slide, Array, Sample_Name) %>%
  unite(rowid, Slide, Array, sep = '_') %>%
  separate(Sample_Name, c('Unique_Identifier', 'VID', NA), sep = '_') %>%
  select(VID, rowid) %>%
  filter(rowid %in% colnames(env$M))

# make matrix
data <- env$M %>% as.matrix() %>% t()
data <- data[tmp$rowid, ]
rownames(data) <- tmp$VID

# relabel cpg identifiers
# colnames(data) <- data %>% colnames() %>% stringr::str_split(pattern = '_') %>% map_chr(~ .[1])

# relabel cpg identifiers
mapping <- env$map2genome %>%
  as_tibble() %>%
  select(Name, EPICv1_Loci) %>%
  # select(Name, EPICv1_Loci, Manifest_probe_match) %>%
  filter(EPICv1_Loci != '')
  # filter(EPICv1_Loci != '', Manifest_probe_match == TRUE)
map <- mapping$EPICv1_Loci
names(map) <- mapping$Name
rm(mapping)

data <- data[ , names(map)]
colnames(data) <- map[colnames(data)]

saveRDS(data, 'data/processed/png/import_epigenetics.rds', compress = F)
