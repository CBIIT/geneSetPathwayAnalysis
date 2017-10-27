library(paxtoolsr)

#interactions <- read.table("C:/Users/cannin/Downloads/pc_sif_101813/test_result/test_result.txt", sep=" ", as.is=TRUE)
#interactions <- read.table("../data/pc_sif_102813.gz", sep=" ", as.is=TRUE, quote="")
interactions <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "pc_sif_050814.gz"), sep="\t", as.is=TRUE, quote="")

load(file.path(.lmp, "gene_set_pathway_analysis", "data", "ddr_genelist.Rdata"))
genes <- sample(rownames(ddr.list), 15)

#pathways <- readGmtPathwayLists("../data/c2.cp.v4.0.symbols.gmt")
#genes <- pathways[["REACTOME_SIGNALLING_TO_ERKS"]]

#genes <- c("TP53", "MDM2", "MDMX") 

# TCA Cycle
genes <- c("IDH3B", "DLST", "PCK2", "CS", "PDHB", "PCK1", "PDHA1", "LOC642502", "PDHA2", "LOC283398", "FH", "SDHD", "OGDH", "SDHB", "IDH3A", "SDHC", "IDH2", "IDH1", "ACO1", "ACLY", "MDH2", "DLD", "MDH1", "DLAT", "OGDHL", "PC", "SDHA", "SUCLG1", "SUCLA2", "SUCLG2", "IDH3G", "ACO2")

#genes <- c("FRS2", "THEM4", "KLB", "DOK1", "EIF4B", "EIF4E", "EIF4EBP1", "EIF4G1", "AKT2", "FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", "FGF10", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "MTOR", "GAB1", "FGF20", "FGF22", "GRB2", "EEF2K", "PIK3R4", "INS", "INSR", "IRS1", "PDE3B", "PRKAG2", "PDPK1", "CAB39", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R2", "PRKAG3", "TLR9", "PPM1A", "STRADB", "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "RPTOR", "TRIB3", "RHEB", "RPS6", "RPS6KB1", "MLST8", "LOC644462", "STK11", "TSC1", "TSC2", "LOC729120", "LOC730244", "FGF23", "CAB39L", "IRS2", "FGF18", "FGF17", "STRADA", "KL", "FGF19")

tmp <- interactions[which(interactions$V1 %in% genes), ]
tmp2 <- tmp[which(tmp$V3 %in% genes), ]

# The interaction types to remove 
#notTheseIntTypes <- c("GENERIC_OF", "IN_SAME_COMPONENT")
#tmp2 <- tmp2[!(tmp2$V2 %in% notTheseIntTypes),]

# Remove self interactions 
tmp2 <- tmp2[tmp2$V1 != tmp2$V3,]

tmp3 <- unique(tmp2[,c(1,3)])

write.table(tmp3, file="filteredNetwork.txt", row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
#tmp2 <- tmp

genesNotInPathway <- setdiff(unique(c(as.vector(tmp2$V1), as.vector(tmp2$V3))), genes)


require(RCytoscape)

g <- new("graphNEL", edgemode="directed")

genesToBeAdded <- unique(c(as.vector(tmp2$V1), as.vector(tmp2$V3)))

for(gene in genesToBeAdded) {
	g <- graph::addNode(gene, g)
}

#g <- graph::addNode("A", g)
#g <- graph::addNode("B", g)
#g <- graph::addNode("C", g)

for(i in 1:nrow(tmp2)) {
	#if((tmp2[i,1] != tmp2[i,3]) && (tmp2[i,2] != "GENERIC_OF") && (tmp2[i,2] != "IN_SAME_COMPONENT") && (tmp2[i,2] != "SEQUENTIAL_CATALYSIS")) {
	if((tmp2[i,1] != tmp2[i,3]) && (tmp2[i,2] != "GENERIC_OF") && (tmp2[i,2] != "IN_SAME_COMPONENT")) {
		g <- graph::addEdge(tmp2[i,1], tmp2[i,3], g)
		#edgeData(g, tmp2[i,1], tmp2[i,3], 'edgeType') <- tmp2[i,2]
		cat(tmp2[i,1], " ", tmp2[i,2], " ", tmp2[i,3], "\n") 
	}
} 

#g <- graph::addEdge('A', 'B', g)
#g <- graph::addEdge('B', 'C', g)
#g <- graph::addEdge('C', 'A', g)

#edgeData(g, 'A', 'B', 'edgeType') <- 'phosphorylates'
#edgeData(g, 'B', 'C', 'edgeType') <- 'promotes'
#edgeData(g, 'C', 'A', 'edgeType') <- 'indirectly activates'

cw <- new.CytoscapeWindow("vignette", graph=g, overwriteWindow=TRUE)
setDefaultBackgroundColor (cw, "#0000FF")
setDefaultNodeColor(cw, "#FFFFFF")
setDefaultNodeBorderColor(cw, "#000000")
setDefaultNodeBorderWidth(cw, 5)
setDefaultEdgeColor (cw, "#00FF00")

g <- cw@graph 
g <- initNodeAttribute(graph=g, attribute.name="moleculeType", attribute.type="char", default.value="undefined")
g <- initEdgeAttribute(graph=g, attribute.name="edgeType", attribute.type="char", default.value="unspecified")

for(i in 1:length(genesToBeAdded)) {	
	#if(i %% 2 == 0) {
	if(genesToBeAdded[i] %in% genes) {
		nodeData(g, genesToBeAdded[i], "moleculeType") <- "inPathway"
	} else {
		nodeData(g, genesToBeAdded[i], "moleculeType") <- "outPathway"	
	}
}

for(i in 1:nrow(tmp2)) {
	#if((tmp2[i,1] != tmp2[i,3]) && (tmp2[i,2] != "GENERIC_OF") && (tmp2[i,2] != "IN_SAME_COMPONENT") && (tmp2[i,2] != "SEQUENTIAL_CATALYSIS")) {
	if((tmp2[i,1] != tmp2[i,3]) && (tmp2[i,2] != "GENERIC_OF") && (tmp2[i,2] != "IN_SAME_COMPONENT")) {
		edgeData(g, tmp2[i,1], tmp2[i,3], 'edgeType') <- tmp2[i,2]
	}
} 

cw <- setGraph(cw, g)

#nodeAttributeValues <- c("kinase", "transcription factor", "glycoprotein")
#lineWidths <- c(0, 8, 16)
#setNodeBorderWidthRule(cw, "moleculeType", nodeAttributeValues, lineWidths)

edgeTypeValues <- c("REACTS_WITH", "INTERACTS_WITH", "STATE_CHANGE", "CO_CONTROL", "SEQUENTIAL_CATALYSIS")
edgeColors <- c("#FF0000", "#FFFF00", "#00FFFF", "#FF00FF", "#0000FF")

nodeTypeValues <- c("inPathway", "outPathway")
nodeColors <- c("#FF0000", "#00FF00")

setNodeColorRule(cw, "moleculeType", nodeTypeValues, nodeColors, "lookup")
setEdgeColorRule(cw, "edgeType", edgeTypeValues, edgeColors, mode="lookup")

#g <- cw@graph # created above, in the section 'A minimal example'
#g <- initNodeAttribute(graph=g, attribute.name='moleculeType', attribute.type='char', default.value='undefined')
#g <- initNodeAttribute(graph=g, 'lfc', 'numeric', 0.0)
#
#nodeData(g, 'A', 'moleculeType') <- 'kinase'
#nodeData(g, 'B', 'moleculeType') <- 'TF'
#nodeData(g, 'C', 'moleculeType') <- 'cytokine'
#nodeData(g, 'A', 'lfc') <- -1.2
#nodeData(g, 'B', 'lfc') <- 1.8
#nodeData(g, 'C', 'lfc') <- 3.2
#
#cw <- setGraph(cw, g)

#nodeAttributeValues <- c('kinase', 'transcription factor', 'glycoprotein')
#lineWidths <- c(0, 8, 16)
#setNodeBorderWidthRule(cw, 'moleculeType', nodeAttributeValues, lineWidths)

displayGraph(cw)

layoutNetwork(cw, layout.name="force-directed")
redraw(cw)

fitContent(cw)

imageType <- "pdf" 
filename <- paste(getwd(), "allEdges.pdf", sep="/")
saveImage (cw, filename, "pdf", 2.0) 

