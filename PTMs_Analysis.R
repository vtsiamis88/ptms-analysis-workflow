#########################################################################
# A script for Protein and PTM quantitative analysis.

#1) Place the analysis and function scripts at a working directory
#2) Create three sub-folders: "1. Quality Control", "2. Protein Quantification" and "3. PTM Quantification"
#   respectively for PSM report quality control and merging, protein quantification and PTM quantification.
#3) Place the PSM reports at "1. Quality Control" folder and follow the analysis instructions in the present script.

#########################################################################
#1) Merge PD reports of fractions of multiple MS runs.

repo <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(repo)

#install.packages("protr")
library(protr)
#install.packages("stringi")
library(stringi)
source("meRgreeMS.R")
source("Functions.R")

setwd(paste0(repo, "/1. Quality Control"))

# PSM report files of the multi-run MS analysis.
files <- c("TMT1 nonmod psms.xlsx" ,"TMT2 nonmod psms.xlsx","TMT1 LA allmods psms.xlsx", "TMT2 LA allmods psms.xlsx", "TMT1 phos allmod psms.xlsx", "TMT2 phos allmod psms.xlsx")

# For merged.Dataset() function will merge the PD reports from different MS experiments.
# The function generates a QC report per PSM report on PSM and peptide level and for the merged dataset on peptide level.
# The QC reports can be found at: .../1. Quality Control/QC plots
# Input description for merged.Dataset() can be found at meRgreeMS.R file.
merged.Dataset <- datasets.Merge(files = files, 
                                 runIDs = c("TMT1", "TMT2"), 
                                 pepC = 5, #Column index in PD files
                                 ptmC = 6, #Column index in PD files
                                 intensitiesC = 30:39, #Column indices in PD files
                                 referenceC = 10, #
                                 log2 = F, 
                                 statC = 46, 
                                 statisticalThreshold = 0.01, 
                                 localC = 48, 
                                 localizationThreshold = 90,
                                 additionalC = 40, 
                                 label = "No Quan Labels", 
                                 mods = c("Acetyl", "Phospho"), 
                                 colnames = c("Peptide", "C_0_R_1", "C_0_R_2", "C_1_R_1", "C_1_R_2", "C_2_R_1", "C_2_R_2", "C_3_R_1", "C_4_R_1", "C_5_R_1",
                                              "C_0_R_3", "C_1_R_3", "C_2_R_3", "C_3_R_2", "C_3_R_3", "C_4_R_2", "C_4_R_3", "C_5_R_2", "C_5_R_3"),
                                 norm = "none")

# Export the merged peptide report.
write.csv(merged.Dataset, file = "Peptide_report.csv", row.names = FALSE, quote = FALSE)

#########################################################################

#########################################################################
#2) Protein Quantification.
# Use VIQoR for protein quantification: http://computproteomics.bmb.sdu.dk:8192/app_direct/VIQoR/
# Input: Peptide_report.csv and the protein sequence database file in .fasta  format
# Ouput: Protein_Expressions.csv and Modified_Peptides.csv

setwd(paste0(repo, "/2. Protein Quantification"))

# Load Protein expression file.
ProteinExprs <- read.csv("Protein_Expressions.csv", stringsAsFactors = F)

# Protein expression normalization and QC report at: .../2. Protein Quantification/QC plots.
ProteinExprs <- protein.Report.QC(data = ProteinExprs, normalization = "median")

# Export the normalized protein expressions.
write.csv(ProteinExprs, file = "Protein_Expressions_Norm.csv", row.names = FALSE)

# Use PolySTest for statistical analysis": http://computproteomics.bmb.sdu.dk:8192/app_direct/PolySTest/
# Input: Protein_Expression_Norm.csv
# Ouput: Protein_Statistics.csv

# Load Protein Statistics file.
ProteinStats <- read.csv("Protein_Statistics.csv", stringsAsFactors = F)

# Keep only Log2 ratios and PolySTest results.
ProteinStats <- ProteinStats[,5:14]

# The following function will keep only the significant protein features and createa barplot for the expressed features.
# The figure can be found at: .../2. Protein Quantification/Figures
# According to the analysis of PolySTest, the suggested log FC is ±0.0469393788945891 and FDR of 0.002.
ProteinStats <- regulated.Features(data = ProteinStats, FC.c = c(1,5), FC.thresholds = c(-0.0469393788945891,0.0469393788945891), Stats.c = c(6,10), Stats.threshold = c(0.002), label = "protein")

# Or it could be changed. For istance: q-value higher than 0.05 and half/double FC relatively to the reference condition. 
#ProteinStats <- regulated.Proteins(data = ProteinStats, FC.c = c(1,5), FC.thresholds = c(-1,1), Stats.c = c(6,10), Stats.threshold = c(0.05))

# Update the protein expression data frame and add the statistic columns.
ProteinExprs <- cbind(ProteinExprs[ProteinStats$Selected,], ProteinStats[ProteinStats$Selected, -11])
rownames(ProteinExprs) <- 1:nrow(ProteinExprs)

#Export the significant regulated protein group features for clustering analysis.
write.csv(ProteinExprs[,c(1,5:22)], file = "Protein_Expressions_Significant.csv", row.names = FALSE)

# Use VSClust for clustering the regulated protein features": http://computproteomics.bmb.sdu.dk:8192/app_direct/VSClust/
# Input: Protein_Expressions_Significant.csv
# Ouput: Protein_Clusters.csv

# Import the output of VSClust.
ProteinClusters <- read.csv("Protein_Clusters.csv", stringsAsFactors = F)
ProteinClusters <- ProteinClusters[order(ProteinClusters$X),]

# Add the cluster number to proteins and keep only the protein features that are members of the clusters.
ProteinExprs <- ProteinExprs[order(ProteinExprs$Protein),]
ProteinClusters <- ProteinClusters[order(ProteinClusters$X),]

ProteinExprs$Cluster <- ProteinClusters$cluster

ProteinExprs <- ProteinExprs[ProteinClusters$isClusterMember,]
rownames(ProteinExprs) <- 1:nrow(ProteinExprs)

# Export the final protein table.
write.csv(ProteinExprs, file = "Proteins_Final.csv", row.names = FALSE)

#...
#... Data ready for go-term and pathway enrichment analysis.
#...

#########################################################################

#########################################################################
#3) PTMs Quantification
setwd(paste0(repo, "/3. PTM Quantification"))

# Load the modified peptide data (output of VIQoR).
ModifiedPeptides <- read.csv(paste0(repo,"/2. Protein Quantification/Modified_Peptides.csv"), stringsAsFactors = F)

# Load the .fasta file.
Fasta <- protr::readFASTA(file = paste0(repo,"/1. Quality Control/uniprot-(taxonomy_9606)+AND+reviewed_yes.fasta"), legacy.mode = TRUE, seqonly = FALSE)

# Load the normalized protein expression data.
ProteinExprsNorm <- read.csv(paste0(repo, "/2. Protein Quantification/Protein_Expressions_Norm.csv"), stringsAsFactors = F)

# Normalize PTMs relatively to the protein expressions and create a small QC report.
# The QC report is located at: .../2. Protein Quantification/QC plots
ModifiedPeptides <- PTMs.Normalization(peptides = ModifiedPeptides, proteins = ProteinExprsNorm, fasta = Fasta, quantitative.columns = 18, normalization = "median")

# Export the normalized modified peptides to be used in PolySTest for statistical analysis.
write.csv(ModifiedPeptides$mapped, file = "Normalized_PTMs.csv", row.names = FALSE)

# Export the unmapped modified peptides.
write.csv(ModifiedPeptides$unmapped, file = "UnNormalized_PTMs.csv", row.names = FALSE)

# Use PolySTest for statistical analysis": http://computproteomics.bmb.sdu.dk:8192/app_direct/PolySTest/
# Input: Normalized_PTMs.csv
# Ouput: PTMs_Statistics.csv

# Load results of PolySTest on modified peptides.
PTMsStats <- read.csv("PTMs_Statistics.csv", stringsAsFactors = F)

# Keep only FC and PolySTest statistics.
PTMsStats <- PTMsStats[,4:13]

# The following function will keep only the significant modified peptide features and create a barplot for the expressed features.
# The figure can be found at: .../2. Protein Quantification/Figures
# According to the analysis of PolySTest, the suggested log FC is ±-0.185208164958878 and FDR of 0.002.
PTMsStats <- regulated.Features(data = PTMsStats, FC.c = c(1,5), FC.thresholds = c(-0.-0.185208164958878, 0.185208164958878), Stats.c = c(6,10), Stats.threshold = 0.002, label = "PTM")

PTMsSignificant <- cbind(ModifiedPeptides$mapped[PTMsStats$Selected,], PTMsStats[PTMsStats$Selected, -11])
rownames(PTMsSignificant) <- 1:nrow(PTMsSignificant)

# Separate only phosphorylated peptides.
PhosphoSignificant <- PTMsSignificant[stringi::stri_detect(PTMsSignificant$Peptide, regex = "^(?!.*Acetyl).*Phospho.*$"), ]
rownames(Phospho.Significant) <- 1:nrow(Phospho.Significant)

# Export significant phosphopeptide features for clustering analysis.
write.csv(PhosphoSignificant[ ,c(1,4:21)], file = "Phospho_Significant.csv", row.names = FALSE)

# Separate only acetylated peptides.
AcetylSignificant <- PTMsSignificant[stringi::stri_detect(PTMsSignificant$Peptide, regex = "^(?!.*Phospho).*Acetyl.*$"), ]
rownames(AcetylSignificant) <- 1:nrow(AcetylSignificant)

# Export the significant acetylated peptide features for clustering analysis.
write.csv(AcetylSignificant[,c(1,4:21)], file = "Acetyl_Significant.csv", row.names = FALSE)

# Use VSClust for clustering the PTM features: http://computproteomics.bmb.sdu.dk:8192/app_direct/VSClust/
# Input: Phospho_Significant.csv and Acetyl_Significant.csv
# Ouput: Phospho_Clusters.csv and Acetyl_Clusters.csv

# Import the clustered phosphorylated peptide features.
PhosphoClusters <- read.csv("Phospho_Clusters.csv", stringsAsFactors = F)
PhosphoClusters <- PhosphoClusters[order(PhosphoClusters$X),]

# Add the cluster number to the phosphopeptides and keep only the significant features that are members of the clusters.
PhosphoSignificant <- PhosphoSignificant[-sapply(setdiff(PhosphoSignificant$Peptide, PhosphoClusters$X), function(x) which(PhosphoSignificant$Peptide == x)), ]
PhosphoSignificant <- PhosphoSignificant[order(PhosphoSignificant$Peptide),]
PhosphoSignificant$Cluster <- PhosphoClusters$cluster
PhosphoSignificant <- PhosphoSignificant[PhosphoClusters$isClusterMember,]
rownames(PhosphoSignificant) <- 1:nrow(PhosphoSignificant)

# Export the final phosphopeptide table.
write.csv(PhosphoSignificant, file = "Phospho_Final.csv", row.names = FALSE)

# Import the clustered acetylated peptide features.
AcetylClusters <- read.csv("Acetyl_Clusters.csv", stringsAsFactors = F)
AcetylClusters <- AcetylClusters[order(AcetylClusters$X),]

# Add the cluster number to the acetylated peptides and keep only the significant features that are members of the clusters.
AcetylSignificant <- AcetylSignificant[order(AcetylSignificant$Peptide),]
AcetylSignificant$Cluster <- AcetylClusters$cluster
AcetylSignificant <- AcetylSignificant[AcetylClusters$isClusterMember,]
rownames(AcetylSignificant) <- 1:nrow(AcetylSignificant)

# Export the final table of the acetylated peptides.
write.csv(PhosphoSignificant, file = "Acetyl_Final.csv", row.names = FALSE)

#...
#... Data ready for go-term and pathway enrichment analysis.
#...

############################################



