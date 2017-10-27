library(rcellminerUtils)
library(geneSetPathwayAnalysis)
library(GO.db)
library(ctrpData)
library(rcellminer)

#goTermsPath <- file.path(.lmp, "geneSetPathwayAnalysis", "inst", "extdata", "archive", "data", "all_go_slim_generic.txt")
goTermsPath <- file.path("inst/scripts/apoptosis_go_terms.txt")
selectedGoTerms <- read.table(goTermsPath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#selectedGoTerms <- selectedGoTerms[1:2,]
#selectedGoTerms <- selectedGoTerms[which(selectedGoTerms$ontology == "BP"), ]
#goIdx <- which(!(selectedGoTerms$go_child_term_name %in% c("molecular_function", "protein binding")))
#selectedGoTerms <- selectedGoTerms[goIdx,]

symbols <- geneSets[["All Gene Sets"]]
#symbols <- "BAX"

#goResultsNames <- getGoCnts(symbols, selectedGoTerms, returnType="list", useNames=TRUE)
#goResultsTerms <- getGoCnts(symbols, selectedGoTerms, returnType="list", useNames=FALSE)

dat <- select(org.Hs.eg.db, keys=symbols, columns=c("SYMBOL", "GO"), keytype="SYMBOL")
t2 <- dat[dat$GO %in% selectedGoTerms$go_child_term,]

goids <- dat$GO
dat2 <- select(GO.db, keys=goids, columns=c("GOID", "TERM"), keytype="GOID")
t3 <- unique(dat2)

t4 <- merge(t2, t3, by.x="GO", by.y="GOID")
t5 <- t4[, c("GO", "SYMBOL", "ONTOLOGY", "TERM")]
t6 <- unique(t5)

a <- t6[which(t6$TERM == "negative regulation of apoptotic process"), "SYMBOL"]
b <- t6[which(t6$TERM == "positive regulation of apoptotic process"), "SYMBOL"]
c <- t6[which(t6$TERM == "positive regulation of cell division"), "SYMBOL"]
d <- t6[which(t6$TERM == "negative regulation of cell division"), "SYMBOL"]

length(a)
length(b)
length(c)
length(d)

#saveRDS(goResultsNames, file=goResultsNamesFilename)
#saveRDS(goResultsTerms, file=goResultsTermsFilename)

drugAct <- exprs(getAct(ctrpData::drugData))
geneExp <- getAllFeatureData(ctrpData::molData)[["exp"]]

genes <- rownames(geneExp)
genes <- genes[genes %in% b]

drug <- "ABT-199"

plotData <- as.data.frame(t(rbind(drugAct[drug, , drop=FALSE], geneExp[genes, ])))
c1 <- drugAct[drug, , drop=FALSE]
c2 <- geneExp[genes, ]
c3 <- crossCors(c1, c2, method = "spearman")
sort(as.vector(c3$cor))

# Cross-cors 
ddrGenes <- geneSets$DDR
d1 <- geneExp[genes, ]
d2 <- pairwiseCor(d1)
d2 <- crossCors(d1)

pairwiseCor <- function(dataframe){
  pairs <- combn(rownames(dataframe), 2, simplify=FALSE)
  df <- data.frame(var1=rep(0,length(pairs)), var2=rep(0,length(pairs)), 
                   absCor=rep(0,length(pairs)), cor=rep(0,length(pairs)))
  
  for(i in 1:length(pairs)){
    df[i,1] <- pairs[[i]][1]
    df[i,2] <- pairs[[i]][2]
    df[i,3] <- round(abs(cor(dataframe[pairs[[i]][1], ], dataframe[pairs[[i]][2], ])),3)
    df[i,4] <- round(cor(dataframe[pairs[[i]][1], ], dataframe[pairs[[i]][2], ]), 3)
  }
  
  pairwiseCorDF <- df
  pairwiseCorDF <- pairwiseCorDF[order(pairwiseCorDF$absCor, decreasing=TRUE),]
  row.names(pairwiseCorDF) <- 1:length(pairs)

  pairwiseCorDF
}

