#--------------------------------------------------------------------------------------------------
# PREPARE AND SAVE GENE SETS LIST
#--------------------------------------------------------------------------------------------------
library(rcmUtils)

maDnaRep <- sort(unique(read.table(file="inst/extdata/aladjem_dnarep.txt",
                                   header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                   comment.char="", quote="")$GENE_NAME))

geneSets <- list()
geneSets[["All Gene Sets"]] <- ""
geneSets[["ABC Transporters"]] <- rcmUtils::cancerGeneSets[["abc_transporters"]]
geneSets[["Apoptosis"]] <- rcmUtils::cancerGeneSets[["apoptosis"]]
geneSets[["Cell Signaling"]] <- rcmUtils::cancerGeneSets[["cell_signaling"]]

geneSets[["DDR"]] <- rcmUtils::cancerGeneSets[["ddr"]]
geneSets[["DDR (BER)"]]   <- rcmUtils::ddrGeneSets[["BER"]]
geneSets[["DDR (NER)"]]   <- rcmUtils::ddrGeneSets[["NER"]]
geneSets[["DDR (MMR)"]]   <- rcmUtils::ddrGeneSets[["MMR"]]
geneSets[["DDR (FA)"]]    <- rcmUtils::ddrGeneSets[["FA (Fanconi anemia pathway)"]]
geneSets[["DDR (HR)"]]    <- rcmUtils::ddrGeneSets[["HR (Homologous Recombination)"]]
geneSets[["DDR (NHEJ)"]]  <- rcmUtils::ddrGeneSets[["NHEJ"]]
geneSets[["DDR (Direct Repair)"]]   <- rcmUtils::ddrGeneSets[["Direct Repair"]]
geneSets[["DDR (G1-S checkpoint)"]]   <- rcmUtils::ddrGeneSets[["G1-S checkpoint"]]
geneSets[["DDR (DNA replication)"]]   <- sort(unique(c(rcmUtils::ddrGeneSets[["DNA replication"]],
                                                       maDnaRep)))
geneSets[["DDR (TLS)"]]   <- rcmUtils::ddrGeneSets[["TLS"]]
geneSets[["DDR (G2-M checkpoint)"]]   <- rcmUtils::ddrGeneSets[["G2-M checkpoint"]]
geneSets[["DDR (Chromatin)"]]   <- rcmUtils::ddrGeneSets[["Chromatin remodelling"]]

geneSets[["EMT (Epithelial)"]] <- sort(unique(c(
  "ESRP1", "C1orf172", "RAB25", "C1orf210", "MARVELD3", "CLDN7",
  "IRF6", "TMC4", "PRRG2", "PPP1R14D", "GRHL2", "ST14", "TMEM125",
  "ATP2C2", "CCDC64B", "S100A14", "CBLC", "SOWAHB", "LOC100288748",
  "SPINT1", "MAPK15", "CNKSR1", "CDC42BPG", "ELMO3", "BSPRY", "TJP3",
  "ILDR1", "PRSS8", "ESRP2", "EPHA1", "CLDN4")))

geneSets[["EMT (Mesenchymal)"]] <- sort(unique(c(
  "LIX1L", "VIM", "CCDC88A", "EMP3", "MSN", "CMTM3", "IKBIP", "QKI", "AP1M1",
  "FRMD8P1", "STARD9", "BICD2", "CHST10", "FAM126A", "SYDE1", "GNB4", "LRP12",
  "SACS", "DYRK3", "MAP7D1", "RECK", "NR3C1", "ST3GAL3", "SLC35B4", "SOAT1",
  "MAP7D3", "SPG20", "ELOVL5", "ETS1")))


geneSets[["Oncogenes"]] <- rcmUtils::cancerGeneSets[["oncogenes"]]
geneSets[["Protein Kinases"]] <- rcmUtils::cancerGeneSets[["protein_kinases"]]
geneSets[["Solute Carriers"]] <- rcmUtils::cancerGeneSets[["solute_carriers"]]
geneSets[["Transcription Factors"]] <- rcmUtils::cancerGeneSets[["transcription_factors"]]
geneSets[["Tumor Suppressors"]] <- rcmUtils::cancerGeneSets[["tumor_suppressors"]]

geneSets[["All Gene Sets"]] <- setdiff(sort(unique(c(geneSets, recursive = TRUE))), "")

geneSets <- lapply(geneSets, FUN = unique)

stopifnot(all(sort(vapply(geneSets, length, integer(1)))) > 0)

save(geneSets, file = "data/geneSets.RData")

#--------------------------------------------------------------------------------------------------
