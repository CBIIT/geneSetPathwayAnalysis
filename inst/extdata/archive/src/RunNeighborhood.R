library(jsonlite)

x <- fromJSON("~/Downloads/del.json")

library(paxtoolsr)
library(data.table)
library(igraph)

tmp <- downloadPc2("Pathway Commons.7.All.BINARY_SIF.hgnc.sif.gz")
tmp <- downloadPc2("PathwayCommons.8.All.BINARY_SIF.hgnc.txt.sif.gz")
#tmp <- downloadPc2("PathwayCommons.8.All.BINARY_SIF.hgnc.txt.sif.gz")
sif <- tmp
x <- getSifInteractionCategories()
sif <- filterSif(sif, c(x$BetweenProteins, x$BetweenProteinSmallMolecule))

g <- loadSifInIgraph(sif)

#sif <- tmp$edges

if(class(sif) == "data.table") {
    setDF(sif)
}

y <- unique(c(sif[,1], sif[,3]))
genes <- c("CCND1", "CETN2", "SLFN11", "MDM4", "UBR5", "RPAIN", "SMARCC1", "FANCE", "BRIP1", "GEN1", "RAD17", "TERF2")
z <- genes %in% y
names(z) <- genes

z

# Get shortest path and return SIF
inputFile <- file.path(.lmp, "gene_set_pathway_analysis", "data", "genePairs.txt")
dat <- read.table(inputFile, header = FALSE, sep="\t", stringsAsFactors=FALSE)

results <- list()

for(i in 1:nrow(dat)) {
    a <- dat[i, "V1"]
    b <- dat[i, "V2"]

    results[[i]] <- getShortestPathSif(sif, a, b)
}

results

## DIM
library(rcellminerPubchem)

chebiId <- "CHEBI:63637"
getNscsFromChebi(chebiId)

aIdx <- match("RPAIN", V(g)$name)
bIdx <- match("CETN2", V(g)$name)
sort(neighborhood(g, 1, nodes=bIdx)[[1]]$name)
are.connected(g, aIdx, bIdx)

getShortestPathSif(sif, "TCF3", "MAZ")

plot(n1[[1]], layout = layout.fruchterman.reingold)

chebiId <- "CHEBI:64153"
aIdx <- match(chebiId, V(g)$name)
neighborhood(g, 1, nodes=aIdx)[[1]]$name

# #genes <- c("AKT1", "IRS1", "MTOR", "IGF1R")
# t1 <- graphPc(source = genes, kind = "PATHSBETWEEN", format = "BINARY_SIF", 
#               verbose = TRUE)
# 
# #t2 <- t1
# t2 <- t1[which(t1[, 2] %in% x$BetweenProteins), ]
# 
# g <- graph.edgelist(as.matrix(t2[, c(1, 3)]), directed = FALSE)
# 
# plot(g, layout = layout.fruchterman.reingold)


