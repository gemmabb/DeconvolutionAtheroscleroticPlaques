###################################################################################
#                                                                                 #
# CREATING DIFFERENT DATA INPUTS                                                  #
# 1) Our bulk data (data transformations: lin, log, sqrt and vst)                 #
# 2) Our sc data (original, quality control for cells, qc for cells and genes +   #
#         data transformations: lin, log, sqrt and vst)                           #
# 3) Carotid sc data from Huize Pan, 2020                                         #
#                                                                                 #
# Gemma Bel Bordes, Nov 20th 2022                                                 #
#                                                                                 #
###################################################################################

library(magrittr)
library(tidyverse)
library(Seurat) #4.1.0
library(Biobase)
library(rstudioapi)
#library(DESeq2) #1.32.0

#R version 4.1.2

setwd(dirname(getActiveDocumentContext()$path)) #set wd to the directory of this code file
source("helperFunctions/functions_CreatingDifferentDataInputs.R")

# 1. OUR BULK DATA ####
## a) First bulk dataset with 656 patients: ####
bulkDataNR <- read_delim("input/bulk_RNAseq_raw_counts.txt.minRib.txt.PC.txt") #only well annotated and non-ribosomal genes
bulkDataNR
#HGNC_ComplData <- read.delim("input/HGNC_database.txt") #downloaded on May 13th
HGNC_ComplData <- read.delim("input/hgnc_complete_set.txt") #downloaded on July 14th
head(HGNC_ComplData)
bulkDataNR_correct <- bulkDataNR %>% filter(symbol %in% HGNC_ComplData$symbol) %>% mutate(Updated=symbol) #already in the new database
bulkDataNR_noOriginalSymbol <- bulkDataNR %>% filter(is.na(symbol)) %>% select(!symbol) %>% 
  inner_join(HGNC_ComplData %>% select(c("ensembl_gene_id", "symbol")), by="ensembl_gene_id") %>% mutate(Updated=symbol) #check if there's a relation with the ensembl and the new symbols for those that had no HGNC symbol
bulkDataNR_toCorrect <- bulkDataNR %>% filter(!(symbol %in% HGNC_ComplData$symbol) & !is.na(symbol)) 
newsymbs <- c()
for(s in bulkDataNR_toCorrect$symbol){
  print(s)
  newsymbs <- c(newsymbs, UpdateSymbolList(s, verbose = T))
}
bulkDataNR_toCorrect$Updated <- newsymbs
bulkDataNR <- rbind(bulkDataNR_noOriginalSymbol, rbind(bulkDataNR_correct, bulkDataNR_toCorrect)) %>% select(-c("symbol")) %>% 
  rename(symbol = Updated) #remove old symbols column and replace it with the updated symbols... 35 genes less (no HGNC symbol)

### remove all non-coding genes that are left:
ncGenes <- read_delim("input/HGNC_NonCodingRNA.txt") #all non coding genes, retrieved from https://www.genenames.org/tools/search/ filtering by non coding RNA (July 14th)
bulkDataNR %>% filter(symbol %in% ncGenes$Symbol) #29 genes are non-coding!
bulkDataNR <- bulkDataNR %>% filter(!(symbol %in% ncGenes$Symbol)) #we are left with 20975 genes
### correct for HGNC unique symbols (remove lower count duplicates):
bulkDataNR_HGNCcorrect <- correctForHGNCUniqueSymbols(bulkDataNR) #1992 from duplicates removed
### correct for max read count and UMI saturation AND convert it to matrix counts:
bulkDataNR_HGNCcorrect_MaxUmiCorrect <- correctReadCounts_Max_UMIsaturation(bulkDataNR_HGNCcorrect, 4095) #after correction, max read count is finally 34070
# left with 18988 genes (original bulk dataset)

## b) Second bulk dataset with 32 more patients (overlapping with the sc data, watch out! 4 femoral plaques) ####
#### already corrected for max and umi count!!!
bulkData32 <- read_delim("input/bulkFromSC.txt")
bulkData32
bulkData32 %<>% rename_at(vars(starts_with('AE')), lst(~paste0("ae", str_remove_all(.,"\\D")))) #remove all non digit characters (from AE4432.sam.counts --> ae4432)
max(bulkData32 %>% select(starts_with("ae"))) #check if the counts were already corrected (yes, max read count is 31230)
bulkData28 <- bulkData32 %>% select(-c("ae4499", "ae4518", "ae4531", "ae4557")) #these are femoral plaques! we're left with 28 samples
#### to make it more similar to the original bulk samples, take only the genes that were present in the original samples (here only ensembl ids are annotated)
bulkData28 <- bulkData28 %>% inner_join(bulkDataNR_HGNCcorrect_MaxUmiCorrect %>% select(ensembl_gene_id, symbol), 
                                        by = c("gene"="ensembl_gene_id")) #left with 18563 genes, all unique symbols, no duplicates

# ## c) Integration of the 656+28 patients (intersection ensembl ids) -> 684 patients ####
# bulkDataAllOurPatients <- bulkData28 %>% inner_join(bulkDataNR_HGNCcorrect_MaxUmiCorrect, by = c("gene"="ensembl_gene_id")) 

## d) Transformation (linear, log, sqrt) and to eset and txt (for scaden and for cibersort) ####

### Original 656 samples:
mat_bulkDataOriginal <- as.matrix(bulkDataNR_HGNCcorrect_MaxUmiCorrect %>% select("symbol", starts_with("ae")) %>% column_to_rownames("symbol")) #656 patients
bulkTransformed <- list()
transformations <- c("linear", "log", "sqrt")
folderBulk <- "output/ourBulk/onlyOriginalSamples/"
bulkLinear <- readRDS("output/ourBulk/onlyOriginalSamples/linearBulk_eset.rds")
for(t in transformations){
  bulkTransformed[[t]] <- transformTo(mat_bulkDataOriginal, t)
  asExpressionSet_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_eset.rds"))
  asCibersortTxtInput_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_forCiber.txt"))
  asScadenTxtInput_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_forScaden.txt"))
}
#### divide the origial bulk samples for sex (only using linear data):
clinicalData <- read_delim("input/clinical_data_good.txt") #654 patients... 2 patients missing, consent papers
countsFemale <- mat_bulkDataOriginal[,clinicalData %>% filter(sex=="female") %>% pull(study_number)]
countsMale <- mat_bulkDataOriginal[,clinicalData %>% filter(sex=="male") %>% pull(study_number)]
##### 169 female plaques + 485 male plaques = 654 
asExpressionSet_bulk(countsFemale, folderBulk, "female_Bulk_eset.rds")
asExpressionSet_bulk(countsMale, folderBulk, "male_Bulk_eset.rds")
asScadenTxtInput_bulk(countsFemale, folderBulk, "female_Bulk_forScaden.txt")
asScadenTxtInput_bulk(countsMale, folderBulk, "male_Bulk_forScaden.txt")
asCibersortTxtInput_bulk(countsFemale, folderBulk, "female_Bulk_forCiber.txt")
asCibersortTxtInput_bulk(countsMale, folderBulk, "male_Bulk_forCiber.txt")

### SC-matched bulk 28 samples (only linear transformation):
mat_bulkDataNew <- as.matrix(bulkData28 %>% select("symbol", starts_with("ae")) %>% column_to_rownames("symbol")) #28 patients
bulkTransformed <- list()
transformations <- "linear"
folderBulk <- "output/ourBulk/scMatchingSamples/"
for(t in transformations){
  bulkTransformed[[t]] <- transformTo(mat_bulkDataNew, t)
  asExpressionSet_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_eset.rds"))
  asCibersortTxtInput_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_forCiber.txt"))
  asScadenTxtInput_bulk(bulkTransformed[[t]], folderBulk, paste0(t, "Bulk_forScaden.txt"))
}

# 2. OUR SINGLE CELL DATA ####
load("input/20210811.46.patients.Koen.rds") #defined as seuset
length(rownames(GetAssayData(object = seuset[["RNA"]], slot = "counts"))) #20111 genes
sum(rownames(GetAssayData(object = seuset[["RNA"]], slot = "counts")) %in% HGNC_ComplData$symbol) #19855 genes --> some should be updated
# but note that for the deconvolution, only the genes overlapping with the bulk data will be used:
sum(setdiff(rownames(GetAssayData(object = seuset[["RNA"]], slot = "counts")), HGNC_ComplData$symbol) %in% bulkDataNR_HGNCcorrect_MaxUmiCorrect$symbol) 
### since the genes that are not updated are NOT present in the bulk dataset, we do not have to care about them (changing an entire seurat object can be very tricky)

## a) With no corrections ####
scReadCounts <- as.matrix(round(GetAssayData(object = seuset[["RNA"]], slot = "counts"), 0)) #20111 genes, 4948 cells

#check relationship between #genes expressed and #counts
n_counts <- colSums(scReadCounts)
n_genes <- colSums(scReadCounts>0) 
plot(n_counts, n_genes); title("Original sc data")

### Pheno data with different cell types clustering: ####
UMAPPlot(seuset, label=T, label.size = 2, repel = TRUE)

### Remove cell types with too few cells sequenced ####
seuset$cellType <- Idents(seuset)
pheno_scReadCounts <- seuset@meta.data %>% mutate(SampleID = rownames(seuset@meta.data), .before=1) %>% 
  rename(SubjectName = Patient)
cellTypesDF <- as.data.frame(sort(table(pheno_scReadCounts$cellType), decreasing = T)) %>% rename(CellType=Var1, NumCells=Freq)
par(mar=c(14,4,4,2))
ggplot(cellTypesDF, aes(x=CellType, y=NumCells))+geom_bar(stat = "identity")+
  geom_hline(aes(yintercept=50, color="limit"))+scale_color_manual(name="",values=c("limit"="blue"))+
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
#change h line to see the effect of changing how many #cells we want for each cell type
scReadCounts_rmvCellTypes <- 
  afterQCSingleCell_rmvCellTypes(scReadCounts, pheno_scReadCounts, minNumCells=50, reads=T) #monocytes and plasma B cells
pheno_scReadCounts_rmvCellTypes <- 
  afterQCSingleCell_rmvCellTypes(scReadCounts, pheno_scReadCounts, minNumCells=50, pheno=T)

dim(scReadCounts_rmvCellTypes) #20111 genes, 4903 cells
dim(pheno_scReadCounts_rmvCellTypes) #4903 cells as well (from all 46 patients)
table(pheno_scReadCounts_rmvCellTypes$Sex) #1931 cells from females, 2972 cells from males

#### 18 cell types ####
pheno_scReadCounts_rmvCellTypes$cellType <- 
  factor(pheno_scReadCounts_rmvCellTypes$cellType, levels = unique(pheno_scReadCounts_rmvCellTypes$cellType))#cell types

#### 9 cell types ####
clustering <- list()
clustering["FOXP3+ T Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "FOXP3+"))
clustering["CD3+ T Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "CD3+.*T"))
clustering["NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "NK Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Class-switched"))
clustering["CD68+ (Foam Cells and inflam/resident macrophages)"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Inflam|Foam|Resident"))
clustering["CD68+ Dendritic"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Dendritic"))
clustering["Endothelial Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "CD34+"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "ACTA2+"))
clustering["CD68+ Mast Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)

clType <- as.character(pheno_scReadCounts_rmvCellTypes$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scReadCounts_rmvCellTypes$clusterType_9 <- factor(clType, levels = names(clustering))

#### 6 cell types ####
clustering <- list()
clustering["T and NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "(?:T|NK) Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Class-switched"))
clustering["CD68+ Cells (no mast)"] <- list(str_subset(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Mast", negate = T), "CD68+"))
clustering["Endothelial Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "CD34+"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "ACTA2+"))
clustering["Mast Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)

clType <- as.character(pheno_scReadCounts_rmvCellTypes$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scReadCounts_rmvCellTypes$clusterType_6 <- factor(clType, levels = names(clustering))

## b) Filtering out genes expressed in less than 5% of the cells ####
load("input/20210811.46.patients.Koen.rds") #defined as seuset
scReadCounts <- as.matrix(round(GetAssayData(object = seuset[["RNA"]], slot = "counts"), 0)) #20111 genes, 4948 cells
scReadCounts_afterGeneF <- scReadCounts %>% QCSingleCell_detectableGenes(minPercent = 5) #12264 genes removed

seuset$cellType <- Idents(seuset)
pheno_scReadCounts <- seuset@meta.data %>% mutate(SampleID = rownames(seuset@meta.data), .before=1) %>% 
  rename(SubjectName = Patient) #we didn't remove any cell, so it's the same pheno data

### Remove cell types with too few cells sequenced ####
scReadCounts_afterGeneF_rmvCellTypes <- 
  afterQCSingleCell_rmvCellTypes(scReadCounts_afterGeneF, pheno_scReadCounts, minNumCells=50, reads=T) #monocytes and plasma B cells
pheno_scReadCounts_afterGeneF_rmvCellTypes <- 
  afterQCSingleCell_rmvCellTypes(scReadCounts_afterGeneF, pheno_scReadCounts, minNumCells=50, pheno=T)

dim(scReadCounts_afterGeneF_rmvCellTypes) #7847 genes, 4903 cells
dim(pheno_scReadCounts_afterGeneF_rmvCellTypes) #4903 cells as well (from all 46 patients)

#### 18 cell types ####
pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType <- 
  factor(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType, levels = unique(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType))#cell types

#### 9 cell types ####
clustering <- list()
clustering["FOXP3+ T Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "FOXP3+"))
clustering["CD3+ T Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "CD3+.*T"))
clustering["NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "NK Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Class-switched"))
clustering["CD68+ (Foam Cells and inflam/resident macrophages)"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Inflam|Foam|Resident"))
clustering["CD68+ Dendritic"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Dendritic"))
clustering["Endothelial Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "CD34+"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "ACTA2+"))
clustering["CD68+ Mast Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)

clType <- as.character(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scReadCounts_afterGeneF_rmvCellTypes$clusterType_9 <- factor(clType, levels = names(clustering))

#### 6 cell types ####
clustering <- list()
clustering["T and NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "(?:T|NK) Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Class-switched"))
clustering["CD68+ Cells (no mast)"] <- list(str_subset(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Mast", negate = T), "CD68+"))
clustering["Endothelial Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "CD34+"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "ACTA2+"))
clustering["Mast Cells"] <- list(str_subset(levels(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)

clType <- as.character(pheno_scReadCounts_afterGeneF_rmvCellTypes$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scReadCounts_afterGeneF_rmvCellTypes$clusterType_6 <- factor(clType, levels = names(clustering))

## d) Transformation (linear, log, sqrt) and to eset and txt (for scaden and for cibersort) ####
scTransformed <- list()
transformations <- c("linear", "log", "sqrt")
folderSC <- "output/ourSC/"
phenoListed <- list(pheno_scReadCounts_rmvCellTypes, pheno_scReadCounts_afterGeneF_rmvCellTypes)
scCountsListed <- list(scReadCounts_rmvCellTypes, scReadCounts_afterGeneF_rmvCellTypes)
outputFilesSCListed <- list("scAllGenes", "scDetectGenes")
i <- 1
for(i in 1:2){ #iterate over the 2 possible datasets (all and filtered genes)
  print(outputFilesSCListed[[i]])
  for(t in transformations){ #iterate over the possible data transformations 
    print(t)
    scTransformed[[t]] <- transformTo(scCountsListed[[i]], t)
    asExpressionSet_sc(scTransformed[[t]], phenoListed[[i]], folderSC, paste0(t, "_", outputFilesSCListed[i], "sc_eset.rds"))
    if(t=="linear"){
      asCibersortTxtInput_sc(scTransformed[[t]], phenoListed[[i]], folderSC, paste0(t, "_", outputFilesSCListed[i], "sc_forCiber.txt"))
      asScadenTxtInput_sc(scTransformed[[t]], phenoListed[[i]], folderSC, paste0(t, "_", outputFilesSCListed[i], "sc_forScaden.txt"))
      asScadenTxtInput_scByPatient(scTransformed[[t]], phenoListed[[i]], paste0(folderSC, "byPatient_Scaden/", outputFilesSCListed[i], "/", ""))
    }
  }
  i <- i+1
}

# write individual files for all genes sc for scaden (only linear data):
# asScadenTxtInput_scByPatient(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "byPatient_Scaden18/"), "cellType")
# asScadenTxtInput_scByPatient(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "byPatient_Scaden9/"), "clusterType_9")
# asScadenTxtInput_scByPatient(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "byPatient_Scaden6/"), "clusterType_6")

# for the special cluster (only linear, only all genes), write sc files:
#### 8 cell types (results from bulk simulation!) ####
clustering <- list()
clustering["T and NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "(?:T|NK) Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Class-switched"))
clustering["Macrophages"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Inflam|Foam|Resident"))
clustering["Dendritic Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Dendritic"))
clustering["Endothelial Cells I"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Endothelial Cells I$"))
clustering["Endothelial Cells II"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Endothelial Cells II"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "ACTA2+"))
clustering["Mastocytes"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)
clustering
clType <- as.character(pheno_scReadCounts_rmvCellTypes$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scReadCounts_rmvCellTypes$specialCluster <- factor(clType, levels = names(clustering))

# #### 7 cell types (same as before but merging ECs)
# clustering <- list()
# clustering["T and NK Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "(?:T|NK) Cells"))
# clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Class-switched"))
# clustering["Macrophages"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Inflam|Foam|Resident"))
# clustering["Dendritic Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Dendritic"))
# clustering["Endothelial Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Endothelial Cells"))
# clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "ACTA2+"))
# clustering["Mastocytes"] <- list(str_subset(levels(pheno_scReadCounts_rmvCellTypes$cellType), "Mast"))
# 
# unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)
# clustering
# clType <- as.character(pheno_scReadCounts_rmvCellTypes$cellType)
# for(i in 1:length(clustering)){
#   clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
# }
# pheno_scReadCounts_rmvCellTypes$specialCluster_bothECs <- factor(clType, levels = names(clustering))

asExpressionSet_sc(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "specialClustering/specialClustering_linear_scAllGenessc_eset.rds"))
asCibersortTxtInput_sc(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "specialClustering/specialClustering_linear_scAllGenessc_forCiber.txt"), different3Clusters = F, "specialCluster")
asScadenTxtInput_scByPatient(scReadCounts_rmvCellTypes, pheno_scReadCounts_rmvCellTypes, paste0(folderSC, "byPatient_Scaden/scAllGenes/"), different3Clusters = F, "specialCluster")

## which samples are from females and from males?... needed for creating the 
paste0(paste0("ae", unique(pheno_scReadCounts %>% filter(Sex=="male") %>% pull(SubjectName)), "_counts.txt"), collapse = " ")
paste0(paste0("ae", unique(pheno_scReadCounts %>% filter(Sex=="male") %>% pull(SubjectName)), "_celltypes.txt"), collapse = " ")

paste0(paste0("ae", unique(pheno_scReadCounts %>% filter(Sex=="male") %>% pull(SubjectName)), ".h5ad"), collapse = ",")
paste0(paste0("ae", unique(pheno_scReadCounts %>% filter(Sex=="male") %>% pull(SubjectName)), "*"), collapse = " ")

# 3. OTHER CAROTID SINGLE CELL DATA (from Huize Pan et al, 2020) ####
load("input/Pan.seurat.RData")
seuset$SubjectName <- map_chr(colnames(seuset), ~{str_split(., pattern = "_")[[1]][2]}) #get subject names from the annotation of each cell id
counts_scHP <- as.data.frame(as.matrix(seuset@assays$RNA@counts)) #17818 genes x 8867 cells
unique(Idents(seuset))
### sanity checking (filtered dataset?) ####
#patient 1 and 2 are males, while patient 3 is female
#in the paper, they say they "genes expressed in less than 10 cells are excluded; then, they exclude the
#cells with +20000 UMI counts, with -200 or +4000 genes expressed and with +10% mapped reads of mitoch. genes
#--> which end up with 15796x2614 (subj1), 17397x3486 (subj2) and 15687x2767 (subj3)" ==> But we have a greater #genes in the seuset!!!
table(seuset$SubjectName) #how many cells per subject? --> same as what they said!
## check min number of cells in which genes are detected for each subject
min(rowSums(counts_scHP > 0)) # min of number of cells in which each gene is detected (it can be that after the filtering of cells, 
# we have a lower number of cells in which genes are detected OR that those additional genes we have in the seuset have 0s in all the cells, let's check that:
length(rownames(counts_scHP)[which(rowSums(counts_scHP > 0)==0)])
length(rownames(counts_scHP)[which(rowSums(counts_scHP > 0) < 10)]) 
## check max UMI count per cell:
max(colSums(counts_scHP)) #less than 20000, it seems to be after filtering
## check min and max number of expressed genes:
min(colSums(counts_scHP > 0)) #greater than 200, it seems to be after filtering
max(colSums(counts_scHP > 0)) #exactly 4000, it seems to be after filtering
## given that I have enough evidence for concluding it's the data after filtering, and that the gene IDs are not only HGNC, I skip the mitoch % checking
#subselect data for each subject and check number of genes if we avoid the ones full of 0 counts:
#s1
counts_scHP_s1 <- counts_scHP %>% select(contains("_1"))
dim(counts_scHP_s1) #17818 genes
sum(rowSums(counts_scHP_s1) == 0) #2022 genes have 0s in all cells
17818-2022 #15796 genes --> the ones they obtained after filtering! So it seems to be okay!
counts_scHP_s1 <- counts_scHP_s1[which(rowSums(counts_scHP_s1) != 0), ] #remove these null genes
#s2
counts_scHP_s2 <- counts_scHP %>% select(contains("_2"))
dim(counts_scHP_s2) #17818 genes
sum(rowSums(counts_scHP_s2) == 0) #422 genes have 0s in all cells
17818-422 #17396 genes --> the ones (-1) they obtained after filtering! So it seems to be okay!
counts_scHP_s2 <- counts_scHP_s2[which(rowSums(counts_scHP_s2) != 0), ] #remove these null genes
#s3
counts_scHP_s3 <- counts_scHP %>% select(contains("_3"))
dim(counts_scHP_s3) #17818 genes
sum(rowSums(counts_scHP_s3) == 0) #2131 genes have 0s in all cells
17818-2131 #15687 genes --> the ones they obtained after filtering! So it seems to be okay!
counts_scHP_s3 <- counts_scHP_s3[which(rowSums(counts_scHP_s3) != 0), ] #remove these null genes

## a) Properly arrange the gene annotations ####
sort(rownames(counts_scHP))
weirdGenes <- grep(".*\\.\\d*", sort(rownames(counts_scHP)), value = T) #I can't solve the problem removing them all because there are some HGNC genes hidden!
sort(rowSums(counts_scHP[weirdGenes,]), decreasing = T)[1:100] #these weird genes are expressed, I noticed they are long non coding genes BUT I ALSO HAVE SOME HGNC SYMBOLS

counts_scHP_tibble <- as_tibble(counts_scHP %>% rownames_to_column("gene")) #get back to tibble because we might get duplicates now (if we work with rownames, we can't)
counts_scHP_tibble %<>% mutate(geneWithoutVersion = map_chr(gene, ~{str_split(., pattern = "\\.")[[1]][1]})) 
sort(counts_scHP_tibble$geneWithoutVersion)

prev <- HGNC_ComplData$prev_symbol
nonHGNCupdatedGenes <- (counts_scHP_tibble %>% filter(!geneWithoutVersion %in% HGNC_ComplData$symbol))$geneWithoutVersion
length(nonHGNCupdatedGenes) #2420 genes 
#nonHGNCupdatedGenes contain both long non coding annotations AND non-updated HGNC genes
## weird long non coding (the ones not present in the old symbols):
allOldSymbols <- unlist(str_split(paste(prev, collapse = "|"), "\\|+")) #if +2 old symbols, they are annotated as gene1|gene2|... --> get array with all of them separatedly
setdiff(nonHGNCupdatedGenes, allOldSymbols) #it seems we only have the weird long coding genes here :)
notInHGNC <- setdiff(nonHGNCupdatedGenes, allOldSymbols) #only the weird long coding genes here :)
length(notInHGNC) #1707
counts_scHP_tibble %<>% filter(!(geneWithoutVersion %in% notInHGNC)) #remove weird annotated lnoncoding genes
## non-updated HGNC genes:
notUpdatedGenes <- intersect(nonHGNCupdatedGenes, allOldSymbols) #these are the ones I want to keep!

counts_scHP_tibbleCorrect <- counts_scHP_tibble %>% filter(!(geneWithoutVersion %in% notUpdatedGenes)) %>% mutate(Updated = geneWithoutVersion)
counts_scHP_tibbleToCorrect <- counts_scHP_tibble %>% filter(geneWithoutVersion %in% notUpdatedGenes) %>% 
  mutate(Updated = UpdateSymbolList(geneWithoutVersion)) #updated symbols found for 419 genes (and we had 421 to update!)
counts_scHP_tibbleToCorrect$Updated[which(!counts_scHP_tibbleToCorrect$Updated %in% HGNC_ComplData$symbol)] #the 2 that we're missing
#I searched them in the HGNC database and QARS and MUM1 now correspond to 2 and 3 new symbols... I also read they're not tissue specific, so I just delete them
counts_scHP_tibbleToCorrect %<>% filter(!(Updated %in% c("QARS", "MUM1")))
counts_scHP_tibbleUpdated <- rbind(counts_scHP_tibbleCorrect, counts_scHP_tibbleToCorrect) %>% mutate(symbol = Updated) %>% dplyr::select(-c(gene,Updated, geneWithoutVersion))
all(counts_scHP_tibbleUpdated$symbol %in% HGNC_ComplData$symbol) #we only have updated genes
sum(duplicated(counts_scHP_tibbleUpdated$symbol)) #we have duplicated gene symbols!

counts_scHP_tibbleUpdated <- correctForHGNCUniqueSymbols_onlyHGNC(counts_scHP_tibbleUpdated)
matrixcounts_scHP_Updated <- as.matrix(counts_scHP_tibbleUpdated %>% column_to_rownames("symbol"))
pheno_scHP <- seuset@meta.data

## b) Remove low populated and undefined cluster 13 ####
table(pheno_scHP$Author_Provided) #only 55 cells coming from this unknown cluster 13
index13 <- which(pheno_scHP$Author_Provided=="13") #these cells are the ones I want to remove
matrixcounts_scHP_Updated_without13 <- matrixcounts_scHP_Updated[,-index13]
all(rownames(pheno_scHP)[index13] == colnames(matrixcounts_scHP_Updated)[index13])
pheno_scHP_without13 <- pheno_scHP[-index13,]
pheno_scHP_without13$SampleID <- rownames(pheno_scHP_without13)

# to make it similar to our sc data, make one cluster for all macrophages
clustering <- list()
for(ctsame in c("SMC", "ICS", "Fibrochondrocyte", "Endothelial 1", "Endothelial 2", "T cell", "Fibroblast", "Mast cell", "Plasma cell")){
  clustering[ctsame] <- ctsame
}
clustering["Macrophages"] <- list(str_subset(levels(pheno_scHP_without13$Author_Provided), "Macrophage"))

unlist(unname(clustering)) #we still have 12 sub cell types (sanity check)

clType <- as.character(pheno_scHP_without13$Author_Provided)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
pheno_scHP_without13_oneMacro <- pheno_scHP_without13
pheno_scHP_without13_oneMacro$Author_Provided <- factor(clType, levels = names(clustering))
pheno_scHP_without13_oneMacro

## c) To eset and txt files (for cibersort and for scaden)
folderSCHP <- "output/huizePanSC/"
asExpressionSet_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "linear_sc_eset.rds")
asCibersortTxtInput_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "sc_forCiber.txt", different3Clusters = F, cluster = "Author_Provided")
asScadenTxtInput_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "sc_forScaden.txt", different3Clusters = F, cluster = "Author_Provided")
asScadenTxtInput_scByPatient(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, paste0(folderSCHP,"byPatient_Scaden/"), "Author_Provided")

folderSCHP <- "output/huizePanSC/"
asExpressionSet_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "oneMacro_linear_sc_eset.rds")
asCibersortTxtInput_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "oneMacro_sc_forCiber.txt", different3Clusters = F, cluster = "Author_Provided")
asScadenTxtInput_sc(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, folderSCHP, "oneMacro_sc_forScaden.txt", different3Clusters = F, cluster = "Author_Provided")
asScadenTxtInput_scByPatient(matrixcounts_scHP_Updated_without13, pheno_scHP_without13, paste0(folderSCHP,"byPatient_oneMacro/"), "Author_Provided")


