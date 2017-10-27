#' Filter Binary SIF network by interaction type 
#' 
#' @param networkFile a gzipped binary SIF network (Example: http://www.pathwaycommons.org/pc2/downloads)
#' @param interactionTypes a vector of interaction types 
#'   (List of interaction types: https://docs.google.com/document/d/1coFo66uuPQQ4ZMSHr8IzCV7I2DwXCoDBfZw7Vg4MgUE/edit)
#'   
#' @return a data.frame of filtered interactions with three columns: "PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B"
#' 
#' @examples 
#' networkFile <- file.path(.lmp, "gene_set_pathway_analysis", "data", "Pathway Commons.7.All.BINARY_SIF.hgnc.sif.gz")
#' intTypes <- c("controls-state-change-of", "controls-phosphorylation-of", "controls-transport-of", "controls-expression-of", "catalysis-precedes", "in-complex-with")
#' filteredNetwork <- filterBinarySif(networkFile, intTypes)
#' 
filterBinarySif <- function(networkFile, interactionTypes) {
    network <- read.table(gzfile(networkFile), sep="\t", quote="")
    colnames(network) <- c("PARTICIPANT_A", "INTERACTION_TYPE", "PARTICIPANT_B")

    idx <- which(network$INTERACTION_TYPE %in% interactionTypes) 

    filteredNetwork <- network[idx, ]
    
    return(filteredNetwork)
}
