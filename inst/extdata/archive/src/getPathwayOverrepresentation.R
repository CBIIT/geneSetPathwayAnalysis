#' Run Pathway Analysis 
#' 
#' Conducts a pathway over-representation analysis using one of several
#' pathway gene sets.
#'  
#' @usage getPathwayOverrepresentation(geneClusters, pathwayList="all")
#'  
#' @param geneClusters: a list containing string vectors to be analyzed 
#' @param pathwayList: a string from the following list for the pathway set to be compared 
#'   "gsea", "mim", "pharmgkb", "lmp", "all"; where "all" is the combined lists
#' @param totalGenesNum: the size of the overall set of genes from the 
#'  geneClusters-associated gene sets are drawn.  If totalGeneNum is not set, it is 
#'  computed as the number of distinct genes over all elements of geneClusters.
#'
#' @details Adjusted p-values are calculated using the FDR correction. All results 
#'   are shown. 
#' 
#' @return getPathwayOverrepresentation produces a matrix of parameters  
#'   for a hypergeometric test, the resulting p-values, and FDR adjusted p-values
#'
#' @author Augustin Luna (augustin@mail.nih.gov)
#' @example 
#'#Example 1 
#'geneSet1 <- c("TP53", "TP63", "TP73")
#'geneSets <- list(geneSet1=geneSet1)
#'results <- getPathwayOverrepresentation(geneSets, pathwayList="gsea")
#'
#'#Example 2
#'load("~/Dropbox/drug_target_tmp/gene_set_pathway_analysis/clustermap.genes.Rdata")
#'geneClusters <- list()
#'for(i in 1:max(clustermap.genes)) {
#'	geneClusters[[i]] <- names(clustermap.genes[which(clustermap.genes == i)])
#'}
#'subGeneClusters <- geneClusters[100:138]
#'getPathwayOverrepresentation(subGeneClusters, pathwayList="lmp")
getPathwayOverrepresentation <- function(geneClusters, pathwayList="all", totalGenesNum=NULL, verbose=FALSE) {
	source(file.path(.lmp, "gene_set_pathway_analysis", "src", "readGmtPathwayList.R"))

	verbose <- FALSE 
	
	if(pathwayList == "mim") {
		input_file <- file.path(.lmp, "gene_set_pathway_analysis", "data", "mim_gene_sets.txt")
	} else if(pathwayList == "pharmgkb") {
		input_file <- file.path(.lmp, "gene_set_pathway_analysis", "data", "pharmgkb_gene_sets.txt")
	} else if(pathwayList == "lmp") {
		input_file <- file.path(.lmp, "gene_set_pathway_analysis", "data", "lmp_gene_sets.txt")
	} else if(pathwayList == "gsea") {
		input_file <- file.path(.lmp, "gene_set_pathway_analysis", "data", "c2.cp.v4.0.symbols.gmt")
	} else if(pathwayList == "all") {
		input_file <- file.path(.lmp, "gene_set_pathway_analysis", "data", "all_gene_sets.txt")
	}
	
	pathways_list <- readGmtPathwayList(input_file)
	
	for(pathway in 1:length(pathways_list)) {
		pathways_list[[pathway]] <- pathways_list[[pathway]][which(pathways_list[[pathway]] 
			!= "" & pathways_list[[pathway]] != "-")]
	}
	
	pvals <- NULL
	
	genes_in_pathway_num_vec <- NULL
	pathway_genes_num_vec <- NULL
	total_genes_num_vec <- NULL
	cluster_genes_num_vec <- NULL
	cluster_idx_vec <- NULL
	pathway_idx_vec <- NULL
	pathway_name_vec <- NULL
	genes_in_pathway_vec <- NULL
	pval_vec <- NULL

	# totalGenesNum is needed for phyper so it must be calculated prior to the 
	# iteration.
  if (is.null(totalGenesNum)){
    totalGenesNum <- length(unique(c(geneClusters, recursive=TRUE)))
  }
	
	for(cluster_idx in 1:length(geneClusters)) {
		for(pathway_idx in 1:length(pathways_list)) {
	#for(cluster_idx in 1:10) {
	#	for(pathway_idx in 1:10) {
		genes_in_pathway <- intersect(pathways_list[[pathway_idx]], 
				geneClusters[[cluster_idx]])
			genes_in_pathway_num <- length(genes_in_pathway)
		
			pathway_name <- names(pathways_list)[pathway_idx]	
				
			pathway_genes_num <- length(pathways_list[[pathway_idx]])
		
			cluster_genes_num <- length(geneClusters[[cluster_idx]])
			
			pval <- phyper(genes_in_pathway_num-1, 
				pathway_genes_num, 
				totalGenesNum-pathway_genes_num, 
				cluster_genes_num,
				lower.tail=FALSE)
    
      #-----[TEST]------------------------------------------------
# 		  require(testthat)
#       marginLabels <- list(c("in_pathway", "not_in_pathway"),
#                            c("in_cluster", "not_in_cluster"))
#       contTab <- matrix(0, nrow=2, ncol=2, dimnames=marginLabels)
#       contTab["in_pathway", "in_cluster"] <- genes_in_pathway_num
# 		  contTab["in_pathway", "not_in_cluster"] <- (pathway_genes_num - genes_in_pathway_num)
# 		  contTab["not_in_pathway", "in_cluster"] <- (cluster_genes_num - genes_in_pathway_num)
# 		  contTab["not_in_pathway", "not_in_cluster"] <- (totalGenesNum - cluster_genes_num - contTab["in_pathway", "not_in_cluster"])
#       stopifnot(sum(rowSums(contTab)) == totalGenesNum)
# 		  stopifnot(sum(rowSums(contTab)) == totalGenesNum)
# 		  expect_equal(fisher.test(contTab, alternative="greater")$p.value, pval)
      #-----------------------------------------------------------
				
			pvals <- c(pvals, pval)	
		
			if(verbose) {
				cat("GN: ", genes_in_pathway_num, 
					" PN: ", pathway_genes_num, 
					" TG: ", totalGenesNum, 
					" CG: ", cluster_genes_num,
					" P: ", pval, 
					" CI: ", cluster_idx,
					" PI: ", pathway_idx, 
					" P: ", pathway_name,
					" G: ", genes_in_pathway,
					"\n")
			}
							
			# Store results in individual vectors
			genes_in_pathway_num_vec <- c(genes_in_pathway_num_vec, genes_in_pathway_num)
			pathway_genes_num_vec <- c(pathway_genes_num_vec, pathway_genes_num) 
			total_genes_num_vec <- c(total_genes_num_vec, totalGenesNum)
			cluster_genes_num_vec <- c(cluster_genes_num_vec, cluster_genes_num)
			cluster_idx_vec <- c(cluster_idx_vec, cluster_idx)
			pathway_idx_vec <- c(pathway_idx_vec, pathway_idx)
			pathway_name_vec <- c(pathway_name_vec, pathway_name)
			genes_in_pathway_vec <- c(genes_in_pathway_vec, paste(genes_in_pathway, collapse=" "))
			pval_vec <- c(pval_vec, pval)
		}
	}
				
	pval_adjust <- p.adjust(pvals, "fdr")
	
	results <- data.frame(genes_in_pathway_num=genes_in_pathway_num_vec, 
				pathway_genes_num=pathway_genes_num_vec, 
				total_genes_num=total_genes_num_vec, 
				cluster_genes_num=cluster_genes_num_vec,
				cluster_idx=cluster_idx_vec,
				pathway_idx=pathway_idx_vec,
				pathway_name=pathway_name_vec,
				genes_in_pathway=genes_in_pathway_vec,
				pval=pval_vec, 
				pval_adjust=pval_adjust)

	#save(results, file=output_file)
	
	return(results)
}

