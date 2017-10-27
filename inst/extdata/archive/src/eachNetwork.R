## Get Network
library(paxtoolsr)
library(rcellminer)
library(igraph)

sif <- downloadPc2("PathwayCommons.8.All.EXTENDED_BINARY_SIF.hgnc.txt.gz")
#sif <- downloadPc2("PathwayCommons.8.drugbank.BINARY_SIF.hgnc.txt.sif.gz")

x <- getSifInteractionCategories()
excludeInteractions <- c("in-complex-with")
sifX <- filterSif(sif$edges, setdiff(c(x$BetweenProteins, x$BetweenProteinSmallMolecule), excludeInteractions))

ignoreDb <- "CTD"
t1 <- lapply(sifX$INTERACTION_DATA_SOURCE, function(x) {
    if(ignoreDb %in% x && length(x) == 1) {
        return(FALSE)
    } else {
        return(TRUE)
    }
})

t2 <- unlist(t1)

sifY <- sifX[t2, ]

g <- loadSifInIgraph(sifY)

## Load CHEBI/NSC Mapping
chebiNscMap <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "chebiNscMap.txt"), header=TRUE, sep="\t", stringsAsFactors = FALSE)

## Create Weights

listGenes <- unique(unlist(lov))

eL <- as_edgelist(g)

weights <- numeric(nrow(eL)) 

idx <- which(!(sifY[, 1] %in% listGenes))
weights[idx] <- weights[idx] + 0.5

idx <- which(!(sifY[, 3] %in% listGenes))
weights[idx] <- weights[idx] + 0.5

idx <- which(!(sifY[, 2] %in% "in-complex-with"))
weights[idx] <- weights[idx] + 0.5

## EXAMPLE

# nsc <- "105014"
# idB <- "YIPF3"
# 
# idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]
# 
# y <- getShortestPathSif(g, idA, idB, "out")
# y[,1:6]


##
nsc <- "701852"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("ARHGAP4", "PHF5FP", "FANCI", "PTPN20A", "GUSB", "GJA10", "KIF11", "RIMBP2", "DEFA3")

netR <- NULL

for(idB in idBs) {
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR <- rbind(netR, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

g1 <- loadSifInIgraph(netR)

plot(g1, layout = layout.fruchterman.reingold, edge.arrow.size=0.35)

##
nsc <- "727989"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("RPL15", "RPL15P14", "TMED9", "PRSS57", "UGT2B26P")

netR2 <- NULL

for(idB in idBs) {
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR2 <- rbind(netR2, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

g1 <- loadSifInIgraph(netR2)

plot(g1, layout = layout.fruchterman.reingold, edge.arrow.size=0.75)

##
nsc <- "743414"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("GAGE1", "RHOXF2B", "PKD2L1", "MCHR1", "LOC388813", "GYPA", "HCAR3", "HEMGN")

netR3 <- NULL

for(idB in idBs) {
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR3 <- rbind(netR3, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

g1 <- loadSifInIgraph(netR3)

plot(g1, layout = layout.fruchterman.reingold, edge.arrow.size=0.75)

## 
nsc <- "750690"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("KCNJ3", "TAS2R1", "RPL27AP5", "SCN4B", "SLC22A14", "INSIG1", "NGRN", "ACER1")

netR4 <- NULL

for(i in 1:length(idBs)) {
    idB <- idBs[i]
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR4 <- rbind(netR4, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

g1 <- loadSifInIgraph(netR4)

plot(g1, layout = layout.fruchterman.reingold, vertex.size=1, edge.arrow.size=0.75)

yT1 <- getShortestPathSif(g, "ACER1", "SCN4B", "all")
yT2 <- getShortestPathSif(g, "ACER1", "TAS2R1", "all")
yT3 <- getShortestPathSif(g, "ACER1", "KCNJ3", "all")
yT4 <- getShortestPathSif(g, "ACER1", "INSIG1", "all")
yT5 <- getShortestPathSif(g, "ACER1", "SLC22A14", "all")
yT6 <- rbind(yT1, yT2, yT3, yT4, yT5, netR4)

write.table(yT6, file="~/Downloads/yT6.txt", sep="\t", col.names = TRUE, row.names = FALSE, quote=FALSE)

g1 <- loadSifInIgraph(yT6)

plot(g1, layout = layout.fruchterman.reingold, vertex.size=1, edge.arrow.size=0.75)

aIdx <- match(idA, V(g1)$name)
length(neighborhood(g, 1, aIdx, mode="out")[[1]])
length(neighborhood(g, 2, aIdx, mode="out")[[1]])
length(neighborhood(g, 3, aIdx, mode="out")[[1]])
length(neighborhood(g, 4, aIdx, mode="out")[[1]])

## 
nsc <- "681239"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("PPP1R14A", "CUTA", "RETNLB", "CD70", "ACIN1", "PPP2R2A", "KRT17", "BSN",
          "DMD")

netR4 <- NULL

for(i in 1:length(idBs)) {
    idB <- idBs[i]
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR4 <- rbind(netR4, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

g1 <- loadSifInIgraph(netR4)

plot(g1, layout = layout.fruchterman.reingold, vertex.size=1, edge.arrow.size=0.5)

aIdx <- match(idA, V(g1)$name)
length(neighborhood(g, 1, aIdx, mode="out")[[1]])
length(neighborhood(g, 2, aIdx, mode="out")[[1]])
length(neighborhood(g, 3, aIdx, mode="out")[[1]])
length(neighborhood(g, 4, aIdx, mode="out")[[1]])

##
nsc <- "759856"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- c("MRPS15",
          "RBBP7",
          "UBE3D",
          "NDUFB7",
          "PTPN11",
          "KTN1",
          "ASB9",
          "EIF4H",
          "MMP14")

netR4 <- NULL

for(i in 1:length(idBs)) {
    idB <- idBs[i]
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out", weights)
        netR4 <- rbind(netR4, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}

netR4 <- unique(netR4)
g1 <- loadSifInIgraph(netR4)

plot(g1, layout = layout.fruchterman.reingold, vertex.size=1, edge.arrow.size=0.5)

aIdx <- match(idA, V(g1)$name)
length(neighborhood(g, 1, aIdx, mode="out")[[1]])
length(neighborhood(g, 2, aIdx, mode="out")[[1]])
length(neighborhood(g, 3, aIdx, mode="out")[[1]])
length(neighborhood(g, 4, aIdx, mode="out")[[1]])



