library(fgsea)

#Create a ranked list (query) to be used with fgsea

rankedList <- function(rankedData, effectSize = "logFC", ID = "geneID"){
  rList <- rankedData %>% dplyr::arrange(P.Value)

  rList <- rList[, effectSize] %>% unlist()
  names(rList) <- rankedData[,ID] %>% unlist()
  return(rList)
}

##Create a pathway in fgsea accepted form

pathwayList <- function(msigCollection = "C2", keyType = "human_gene_symbol", species = "Homo sapiens"){

  ##Creating the pathway list to for fgsea

  pathways <- msigdbr::msigdbr(species = species, category = msigCollection)
  newPathway <- pathways %>% dplyr::select(c("gs_name", keyType))
  newPathway <- newPathway %>% nest(key = c(keyType))
  n2 <- lapply(newPathway$key, function(x) unlist(x) %>% unname())
  names(n2) <- newPathway$gs_name
  return(n2)
}


##fGSEA function to be run after a query and pathway is ready

fGSEA <- function(query, pathway){

  fgseaRes <- fgseaMultilevel(pathways = pathway,
                              stats    = query,
                              minSize  = 5,
                              nPermSimple = 1000)


  return(fgseaRes)

}

# Function to create GESECA table (This code is a modified version of fgsea::GESECA() function code)

myTable <- function (gesecaRes, pathways, E, center = TRUE, scale = FALSE,
                     colwidths = c(5, 3, 0.8, 1.2, 1.2), titles = colnames(E),
                     colors = c("blue", "white", "red"), pathwayLabelStyle = NULL,
                     headerLabelStyle = NULL, valueStyle = NULL, axisLabelStyle = NULL,
                     axisLabelHeightScale = NULL)
{
  pathwayLabelStyleDefault <- list(size = 12, hjust = 1, x = 0.95,
                                   vjust = 0)
  pathwayLabelStyle <- modifyList(pathwayLabelStyleDefault,
                                  as.list(pathwayLabelStyle))
  headerLabelStyleDefault <- list(size = 12)
  headerLabelStyle <- modifyList(headerLabelStyleDefault, as.list(headerLabelStyle))
  valueStyleDefault <- list(size = 12, vjust = 0)
  valueStyle <- modifyList(valueStyleDefault, as.list(valueStyle))
  axisLabelStyleDefault <- list(angle = 90, hjust = 1, size = 10)
  axisLabelStyle <- modifyList(axisLabelStyleDefault, as.list(axisLabelStyle))
  if (is.null(axisLabelHeightScale)) {
    axisLabelHeightScale <- max(sapply(titles, nchar))/4 *
      axisLabelStyle$size/pathwayLabelStyle$size
  }
  gesecaRes <- gesecaRes[pathway %in% names(pathways)]
  pathways <- pathways[gesecaRes$pathway]
  E <- t(base::scale(t(E), center = center, scale = scale))
  colnames(E) <- titles
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(fmatch(p, rownames(E)))))
  })
  prjs <- t(do.call(cbind, lapply(pathways, function(p) {
    scale(colSums(E[p, , drop = FALSE]))
  })))
  rownames(prjs) <- names(pathways)
  prjspd <- as.data.table(prjs, keep.rownames = "pathway")

  prjspd <- copy(melt(prjspd, id.vars = "pathway", measure.vars = colnames(prjspd)[2:ncol(prjspd)],
                      variable.name = "sample", variable.factor = FALSE))
  prjspd[, `:=`(pathway, factor(pathway, levels = rev(rownames(prjs))))]
  prjspd[, `:=`(sample, factor(sample, levels = titles))]

  return(prjspd)
}
