library(rcellminer)
library(rcellminerData)
library(paxtoolsr)

#' # Read predictor table 
#+ readPredictorTable
enResults <- read.table(system.file("extdata", "en_cvresults_table.txt", package="rcellminerElasticNet"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
enNscs <- enResults[,"nsc"]
nscs <- unique(enResults[,"nsc"])

#' # Get Pathway Commons Data 
#+ getPathDb
if(!exists("sif")) {
  if(file.exists("sif.rds")) {
    sif <- readRDS("sif.rds")
  } else {
    sif <- downloadPc2("PathwayCommons.8.All.EXTENDED_BINARY_SIF.hgnc.txt.gz", version="8")
    saveRDS(sif, "sif.rds")
  }    
}

x <- getSifInteractionCategories()
# in-complex-with may include high-throughput interaction assays (e.g. Yeast2Hybrid)
#excludeInteractions <- c("in-complex-with", "interacts-with")
excludeInteractions <- c("neighbor-of", "interacts-with")
sifX <- filterSif(sif$edges, setdiff(c(x$BetweenProteins, x$BetweenProteinSmallMolecule, x$BetweenProteinsOther), excludeInteractions))

# CTD contains many "fuzzy" drug/gene interactions 
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

#' # Load CHEBI/NSC Mapping ----
#+ getChebiNscMapping
chebiNscMap <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "chebiNscMap.txt"), header=TRUE, sep="\t", stringsAsFactors = FALSE)

results <- NULL

for(i in 1:nrow(chebiNscMap)) {
  #i <- 1
  
  drug <- chebiNscMap$NSC[i]
  actIdx <- which(rownames(drugActData) == drug)
  
  tmpSif <- filterSif(sifY, ids=chebiNscMap$CHEBI[i])
  genes <- tmpSif$PARTICIPANT_B
  
  molDB <- getMolDataMatrices()
  drugActData <- exprs(getAct(rcellminerData::drugData))
  
  for(j in 1:length(genes)) {
    #j <- 1
    
    gene <- genes[j]
    
    if(paste0("swa", gene) %in% rownames(molDB[["swa"]])) {
      cat("G: ", gene, "\n")
      
      x <- molDB[["swa"]][paste0("swa", gene), ]
      y <- drugActData[actIdx, ]
      
      val <- cor(x, y, use="pairwise.complete.obs")
      
      if(val >= 0.3) {
        cat("H: ", gene, " HIT\n")
        
        #results[drug] <- c(results[drug], gene)
        
        results <- rbind(results, c(drug, val, gene))
        #results[[drug]] <- c(results[[drug]], gene)  
      }      
    } else {
      next
    }
  }
}


