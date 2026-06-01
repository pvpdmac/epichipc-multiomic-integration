# imports
library(tidyverse)

# obtain participant ids for samples used for multiomic work ----

# for testing association with HBsAb titers, use the same subset of participants used to train models
my_participants <- readr::read_rds('results/other_outputs/trained_models.rds')$x %>%
  map('transcriptomics') %>%
  map(rownames) %>%
  reduce(union)

# map to gene symbols
my_map <- readr::read_rds('code/utils/features_to_hgnc_symbols.rds') %>%
  unnest(map) %>%
  filter(block == 'transcriptomics') %>%
  select(-block) %>%
  group_by(feature_name) %>%
  slice(1) %>%
  ungroup() %>%
  rename(id = feature_name)

# obtain some pathways ----
library(org.Hs.eg.db)

reactome_data <- ReactomePA::gson_Reactome(organism = 'human')
reactome_data <- left_join(
  tibble(reactome_data@gsid2name),
  tibble(reactome_data@gsid2gene)
)

# map to gene symbol
reactome_data$hgnc_symbol <- mapIds(
  org.Hs.eg.db,
  keys = reactome_data$gene, # Keys must be characters
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

detach('package:org.Hs.eg.db', force = T)
detach('package:AnnotationDbi', force = T)
detach('package:IRanges', force = T)
detach('package:S4Vectors', force = T)

reactome_data <- reactome_data %>%
  select(gsid, name, hgnc_symbol) %>%
  nest(data = hgnc_symbol) %>%
  mutate(data = data %>% map('hgnc_symbol') %>% map(unname)) %>%
  rename(geneset = gsid, description = name, members_mrna = data)

# travis modified this one very long name in his table...
reactome_data$description[reactome_data$description == 'Diseases of signal transduction by growth factor receptors and second messengers'] <- 'Diseases of signal transduction'

# identify important subset of pathways implicated by our multiomic models ----
f <- readr::read_tsv('results/tables/sigora_pathways.tsv')

# flag target
reactome_data <- reactome_data %>%
  mutate(target = description %in% c('Neutrophil degranulation', 'Interleukin-1 signaling', 'Interleukin-4 and Interleukin-13 signaling', 'Platelet activation, signaling and aggregation'))

# flag candidates
reactome_data <- reactome_data %>%
  mutate(candidate = geneset %in% f$pathway_id)

# add model modules - these are the model features and interactor
# (NetworkAnalyst) genes that were part of the leadingEdge (GSEA) when the
# association between the 4 focus pathways and HBsAb titers were assessed at
# baseline (DOL0)
shannon <- readr::read_rds('reviews/foremostEdge.rds') %>%
  group_by(geneset = analysis) %>%
  summarize(members_mrna = foremostEdge %>% reduce(union) %>% list()) %>%
  mutate(members_mrna = set_names(members_mrna, geneset), target = T, candidate = T)

# add R.Sekaly 2016 modules
sekaly <- readxl::read_excel('data/external/ncomms10369-s1.xlsx')
sekaly <- sekaly %>% filter(Module %in% c(1, 16), `Gene Symbol` != '---') %>% split(.$Module) %>% map('Gene Symbol')
names(sekaly) <- paste0('mod', names(sekaly), ' (PMID:26742691)')
sekaly <- tibble(geneset = c('mod01', 'mod16'), members_mrna = sekaly, target = T, candidate = T)

# combine with reactome pathways
msigdb <- bind_rows(
  shannon,
  sekaly,
  reactome_data
)

# cleanup
msigdb <- msigdb %>%
  mutate(
    description = ifelse(is.na(description), geneset, description),
    members_mrna = set_names(members_mrna, description)
  )

msigdb <- msigdb %>% group_by(description) %>% filter(n() < 2) %>% ungroup()

msigdb <- msigdb %>% filter(map_dbl(members_mrna, length) >= 5)

rm(f, shannon, sekaly)

# create GAM data subsets ----
x <- readr::read_rds('data/processed/gam/import_transcriptomics.rds')

a <- readr::read_rds('data/processed/gam/import_metadata.rds') %>%
  dplyr::select(uid = `Unique Identifier`, vid = `Visit Identifier`, grp = VaccineGrp, daygrp = DayGrp, sex = Sex, dol = DOL) %>%
  dplyr::filter(vid %in% rownames(x))

b <- readr::read_rds('data/processed/gam/import_responses.rds') %>%
  dplyr::select(uid = `Subject ID`, V1:V4) %>%
  dplyr::mutate(
    status = factor(V1 <= 2.5, levels = c(T, F), labels = c('MA-', 'MA+')),
    DOL30 = log(V3)-log(V1),
    DOL128 = log(V4)-log(V1)
  )

y <- left_join(a, b) %>% select(uid, sex, grp, daygrp, status, DOL30, DOL128, dol, vid)
rm(a, b)

con_gam <- tribble(
  ~site, ~query,        ~formula,   ~y,

  # ab status @dol0
  'GAM', 'baseline_status', '~ status', y %>% filter(dol == 0) %>% droplevels() %>% drop_na(),

  # vax status @dol1
  'GAM', 'dol1',            '~ grp',    y %>% filter(status == 'MA-', dol == 1) %>% mutate(grp = factor(grp, levels = c('Delayed', 'BCG', 'HBV', 'HBV+BCG'))) %>% droplevels(),

  # vax status @dol7
  # 'GAM', 'dol7',            '~ grp',    y %>% filter(status == 'MA-', dol == 7, grp == 'Delayed' | grp == 'BCG') %>% mutate(grp = factor(grp, levels = c('Delayed', 'BCG'))) %>% droplevels(),

  # bcg@fol7 vs. hbv@dol0
  'GAM', 'anyhbv0_bcg7',    '~ dol',    y %>% filter(status == 'MA-', (dol == 0 & grp %in% c('HBV')) | (dol == 7 & grp == 'BCG')) %>% mutate(dol = factor(dol, levels = c(0, 7))) %>% droplevels(),

  # ~HBsAb titers @dol0
  'GAM', 'MA-DOL30',        '~ DOL30',  y %>% filter(vid %in% my_participants, status == 'MA-', dol == 0) %>% droplevels() %>% drop_na(),
  'GAM', 'MA-DOL128',       '~ DOL128', y %>% filter(vid %in% my_participants, status == 'MA-', dol == 0) %>% droplevels() %>% drop_na(),

  'GAM', 'MA+DOL30',        '~ DOL30',  y %>% filter(vid %in% my_participants, status == 'MA+', dol == 0) %>% droplevels() %>% drop_na(),
  'GAM', 'MA+DOL128',       '~ DOL128', y %>% filter(vid %in% my_participants, status == 'MA+', dol == 0) %>% droplevels() %>% drop_na()
)

con_gam <- con_gam %>% mutate(x = purrr::map(y, ~ x[.$vid, ]))

# create in vivo, in vitro data subsets ----
invivo <- tribble(
  ~ site,     ~query,   ~formula,           ~coef,         ~y,
  'invivo',   'status', '~ grp + abstatus', 'abstatusMA+',  y %>% filter(dol == 1) %>% mutate(grp = factor(grp, levels = c('Delayed', 'HBV', 'BCG', 'HBV+BCG'), labels = c('Delayed', 'HBV', 'BCG', 'COM')), abstatus = factor(status, levels = c('MA-', 'MA+'))) %>% select(vid, grp, abstatus) %>% drop_na()
)

invivo <- invivo %>% mutate(x = purrr::map(y, ~ x[.$vid, ]))

rm(x, y)

# invitro
x <- readr::read_rds('data/processed/invitro/RNA_WBA_VEH_adjusted_GAM_PNG.rds')
y <- readr::read_rds('data/processed/invitro/Metadata_WBA_VEH_adjusted_GAM_PNG.rds')

# clean-up
x <- x %>% as.matrix()
y <- y %>% rename(grp = stimulation, vid = sample_id) %>% droplevels()

invitro <- tribble(
  ~ site,     ~query,   ~formula,           ~coef,         ~x, ~y,
  'invitro',  'status', '~ grp + abstatus', 'abstatusMA+',  x,  y
)

con_status <- bind_rows(invitro, invivo)

rm(x, y, invitro, invivo)

# create PNG data subsets ----
x <- readr::read_rds('data/processed/png/import_transcriptomics.rds')

a <- readr::read_rds('data/processed/png/import_metadata.rds') %>%
  dplyr::select(uid = Unique_Identifier, vid = VID, grp = VaccineGroup, sex = Sex, dol = DOL) %>%
  dplyr::filter(vid %in% rownames(x))

b <- readr::read_rds('data/processed/png/import_responses.rds') %>%
  dplyr::select(uid = `Subject ID`, V1:V4) %>%
  dplyr::mutate(
    status = factor(V1 <= 2.5, levels = c(F, T), labels = c('MA+', 'MA-')),
    DOL30 = log(V3)-log(V1),
    DOL128 = log(V4)-log(V1)
  )

y <- left_join(a, b) %>% select(uid, sex, grp, status, DOL30, DOL128, dol, vid) %>% drop_na()
rm(a, b)

con_png <- tribble(
  ~site, ~query,        ~formula,   ~y,

  # ab status @dol0
  'PNG', 'baseline_status', '~ status', y %>% filter(dol == 0) %>% droplevels() %>% drop_na(),

  # vax status @dol1
  # 'PNG', 'dol1'     - no dol1 in png...

  # vax status @dol7
  # 'PNG', 'dol7',            '~ grp',    y %>% filter(status == 'MA-', dol == 7, grp == 'Delayed' | grp == 'BCG') %>% mutate(grp = factor(grp, levels = c('Delayed', 'BCG'))) %>% droplevels(),

  # bcg@fol7 vs. hbv@dol0
  'PNG', 'anyhbv0_bcg7',    '~ dol',    y %>% filter(status == 'MA-', (dol == 0 & grp %in% c('HBV')) | (dol == 7 & grp == 'BCG')) %>% mutate(dol = factor(dol, levels = c(0, 7))) %>% droplevels(),

  # ~HBsAb titers @dol0
  # 'PNG', 'MA-DOL30' - no dol30 response in png
  'PNG', 'MA-DOL128',       '~ DOL128', y %>% filter(status == 'MA-', dol == 0) %>% droplevels() %>% drop_na(),

  # 'PNG', 'MA+DOL30' - no dol30 response in png
  'PNG', 'MA+DOL128',       '~ DOL128', y %>% filter(status == 'MA+', dol == 0) %>% droplevels() %>% drop_na()
)

con_png <- con_png %>% mutate(x = purrr::map(y, ~ x[.$vid, ]))

rm(x, y)

# combine
con <- bind_rows(con_gam, con_png) %>% mutate(coef = NULL)
con <- bind_rows(con_status, con)

rm(con_gam, con_png, con_status)

# create subsets of genesets based on the question ----
con$msigdb <- con %>%
  pmap(\(site, query, formula, ...) {

    if(query =='baseline_status' | site %in% c('invivo', 'invitro'))
      msigdb
    else if(query == 'baseline' | grepl('^MA', query) | query == 'anyhbv0_bcg7')
      msigdb %>% filter(candidate)
    else
      msigdb %>% filter(target)
  })

# fit lm ----
con$fit <- con %>%
  pmap(.progress = T, \(x, y, formula, coef, ...) {
    if(is.na(coef))
      coef <- NULL
    d <- model.matrix(as.formula(formula), data = y)
    x %>%
      t() %>%
      limma::lmFit(design = d) %>%
      limma::eBayes() %>%
      limma::topTable(coef = coef, number = Inf)  %>%
      as_tibble(rownames = 'id')  %>%
      inner_join(my_map, by = 'id')
  })

# fgsea ----
set.seed(3341L)
con$gsea <- con %>%
  pmap(.progress = T, \(fit, msigdb, ...) {

    # only keep the most significant by gene symbol
    fit <- fit %>% group_by(hgnc_symbol) %>% arrange(P.Value) %>% slice(1) %>% ungroup() %>% filter(hgnc_symbol != '')

    # deal with diff. coefs
    fit <- fit %>% select(id, coef = matches('logFC$|grpBCG$'), P.Value, adj.P.Val, hgnc_symbol)

    # create ranked list for gsea
    a <- fit %>% arrange(desc(coef))
    b <- a$coef
    names(b) <- a$hgnc_symbol

    # fgsea
    library(fgsea)
    gsea <- fgsea(
      pathways = msigdb$members_mrna,
      stats    = b,
      nPermSimple = 1e4,
      eps      = 0
    ) %>%
      as_tibble() %>%
      arrange(pval)

    tibble(stats = list(b), gsea = list(gsea))

  })

con$ssgsea <- con %>%
  pmap(\(x, y, msigdb, ...) {

    # map to hgnc_symbol
    xx <- x %>%
      t() %>%
      as_tibble(rownames = NA) %>%
      rownames_to_column('id') %>%
      inner_join(my_map, ., by = c('id')) %>%
      group_by(hgnc_symbol) %>% filter(n() < 2) %>% ungroup() %>%
      filter(hgnc_symbol != '') %>%
      drop_na()

    # make matrix
    xx <- xx %>%
      column_to_rownames('hgnc_symbol') %>%
      select(-id) %>%
      as.matrix()

    # scale input counts
    xx <- t(scale(t(xx)))

    # ssGSEA
    gsva_ssgsea <- GSVA::gsva(GSVA::ssgseaParam(xx, msigdb$members_mrna, minSize = 5), verbose = F)

    # container for results
    gsva_ssgsea %>%
      as_tibble() %>%
      mutate(pathway = gsva_ssgsea %>% attr('geneSet') %>% names()) %>%
      gather(vid, nes, -pathway) %>%
      left_join(y)

  })

# save result ----
con <- con %>% select(site:ssgsea) %>% unnest(gsea)
con <- con %>% mutate(gsea = gsea %>% set_names(query), ssgsea = ssgsea %>% set_names(query))

readr::write_rds(con, 'results/other_outputs/gsea_tables.rds')
