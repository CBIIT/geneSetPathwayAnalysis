#--------------------------------------------------------------------------------------------------
# UPDATE GENE (SHORT) ANNOTATIONS TO USE DTB GENE LIST MEMBERSHIP WHERE POSSIBLE
#--------------------------------------------------------------------------------------------------
geneAnnotTab <- geneSetPathwayAnalysis::geneAnnotTab
geneSets <- geneSetPathwayAnalysis::geneSets

dtbGenes <- geneSets[["All Gene Sets"]]
dtbGenesToAdd <- setdiff(dtbGenes, rownames(geneAnnotTab))
if (length(dtbGenesToAdd) > 0){
  tmp <- data.frame(GENE_NAME = dtbGenesToAdd, ANNOT = NA, SHORT_ANNOT = NA,
                    stringsAsFactors = FALSE)
  rownames(tmp) <- tmp$GENE_NAME
  stopifnot(identical(colnames(tmp), colnames(geneAnnotTab)))
  geneAnnotTab <- rbind(geneAnnotTab, tmp)
  geneAnnotTab <- geneAnnotTab[order(geneAnnotTab$GENE_NAME), ]

  stopifnot(all(!duplicated(rownames(geneAnnotTab))))
  stopifnot(all(dtbGenes %in% rownames(geneAnnotTab)))
}

geneSets[["All Gene Sets"]] <- NULL
tmp <- names(geneSets)
tmp[which(tmp == "DDR")] <- "DNA Damage Response (DDR)"
names(geneSets) <- tmp

getDtbAnnot <- function(gene){
  annotStr <- ""
  tmp <- vapply(geneSets, function(x){ gene %in% x}, logical(1))
  if (any(tmp)){
    tmp <- tmp[tmp]
    annotStr <- paste0(names(tmp), collapse = "; ")
  }
  return(annotStr)
}

# If gene is in one or more DTB gene lists, annotate based on these.
# Otherwise retain CellMiner annotation strings.
for (gene in rownames(geneAnnotTab)){
  annotStr <- getDtbAnnot(gene)
  if (annotStr != ""){
    geneAnnotTab[gene, "SHORT_ANNOT"] <- annotStr
  }
}

save(geneAnnotTab, file = "data/geneAnnotTab.RData")
#--------------------------------------------------------------------------------------------------
