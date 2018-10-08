#--------------------------------------------------------------------------------------------------
# PREPARE AND SAVE GENE SETS LIST
#--------------------------------------------------------------------------------------------------
load(file = "data/geneSets.RData")
nold=names(geneSets)
#names(geneSets)
# [1] "All Gene Sets"         "ABC Transporters"      "Apoptosis"             "Cell Signaling"
# [5] "DDR"                   "DDR (BER)"             "DDR (NER)"             "DDR (MMR)"
# [9] "DDR (FA)"              "DDR (HR)"              "DDR (NHEJ)"            "DDR (Direct Repair)"
# [13] "DDR (G1-S checkpoint)" "DDR (DNA replication)" "DDR (TLS)"             "DDR (G2-M checkpoint)"
# [17] "DDR (Chromatin)"       "EMT (Epithelial)"      "EMT (Mesenchymal)"     "Oncogenes"
# [21] "Protein Kinases"       "Solute Carriers"       "Transcription Factors" "Tumor Suppressors"
#
names(geneSets)[5]="DNA damage repair"
names(geneSets)[6]="DNA damage repair, break excision repair"
names(geneSets)[7]="DNA damage repair, nucleotide excision repair"
names(geneSets)[8]="DNA damage repair, mismatch repair"
names(geneSets)[9]="DNA damage repair, Fanconi anemia"
names(geneSets)[10]="DNA damage repair, homologous recombination"
names(geneSets)[11]="DNA damage repair, non-homologous end joining"
names(geneSets)[12]="DNA damage repair, direct repair"
names(geneSets)[13]="DNA damage repair, G1-S checkpoint"
names(geneSets)[14]="DNA damage repair, DNA replication"
names(geneSets)[15]="DNA damage repair, Translesion synthesis"
names(geneSets)[16]="DNA damage repair, G2-M checkpoint"
names(geneSets)[17]="DNA damage repair, chromatin"
names(geneSets)[18]="Epithelial"
names(geneSets)[19]="Mesenchymal"

#names(geneSets)

# [1] "All Gene Sets"                                 "ABC Transporters"
# [3] "Apoptosis"                                     "Cell Signaling"
# [5] "DNA damage repair"                             "DNA damage repair, break excision repair"
# [7] "DNA damage repair, nucleotide excision repair" "DNA damage repair, mismatch repair"
# [9] "DNA damage repair, Fanconi anemia"             "DNA damage repair, homologous recombination"
# [11] "DNA damage repair, non-homologous end joining" "DNA damage repair, direct repair"
# [13] "DNA damage repair, G1-S checkpoint"            "DNA damage repair, DNA replication"
# [15] "DNA damage repair, Translesion synthesis"      "DNA damage repair, G2-M checkpoint"
# [17] "DNA damage repair, chromatin"                  "EMT (Epithelial)"
# [19] "EMT (Mesenchymal)"                             "Oncogenes"
# [21] "Protein Kinases"                               "Solute Carriers"
# [23] "Transcription Factors"                         "Tumor Suppressors"

nnew=names(geneSets)
#write.csv(cbind(old=nold,new=nnew),"~/Documents/GeneSets_names.csv",row.names=F)
#
# not yet
save(geneSets, file = "data/geneSets.RData")

#--------------------------------------------------------------------------------------------------
