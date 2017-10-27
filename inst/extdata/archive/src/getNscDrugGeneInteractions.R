#' Get NSC Drug Gene Interactions from CTD, PharmGKB, and DrugBank
#'
#' @param nsc a string, an NSC ID
#' @param db a string, a database of drug-gene interactions (default: ctd)
#' @return a data.frame with the NSC, drug name, and gene symbol;
#'   for the "ctd" database a column of InteractionActions is provided to
#'   give additional details on the interaction.
#'
#' @examples
#' getNscDrugGeneInteractions("94600")
#'
#' @export
getNscDrugGeneInteractions <- function(nsc, db=c("ctd", "pharmgkb", "drugbank")) {
  require(RSQLite)

  db <- match.arg(db)

  con <- dbConnect(SQLite(), dbname=file.path(.db, "lmp_drug_target_big.sqlite"))

  results <- NULL

  returnValues <- paste0(" dtp_drugs.NSC, dtp_drugs.DrugName, ", db, "_genes.GeneSymbol")

  if(db == "ctd") {
    returnValues <- paste0(returnValues, ", ctd_drug_target.InteractionActions")
  }

  query <- paste0("SELECT", returnValues, "
    FROM dtp_drugs, nsc_cas, ", db, "_drugs, ", db, "_drug_target, ", db, "_genes
  	WHERE dtp_drugs.NSC='", nsc, "'
  	AND nsc_cas.NSC=dtp_drugs.NSC
  	AND nsc_cas.CasRN=", db, "_drugs.CasRN
  	AND ", db, "_drugs.ChemicalID=", db, "_drug_target.ChemicalID
  	AND ", db, "_drug_target.GeneID=", db, "_genes.GeneID")

  cat("NSC: ", nsc, "\n")
  tmp <- dbGetQuery(con, query)

  results <- rbind(results, tmp)
  results <- unique(results)

  return(results)
}
