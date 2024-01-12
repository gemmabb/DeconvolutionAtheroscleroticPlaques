###################################################################################
#                                                                                 #
# AMSTERDAM SINGLE CELL DATASET                                                   #                                       
#                                                                                 #
# Gemma Bel Bordes                                                                #
#                                                                                 #
###################################################################################


library(rstudioapi)
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggcorrplot)
library(tableone)
library(forcats)
library(flextable)
#R version 4.1.2
setwd(dirname(getActiveDocumentContext()$path)) #set wd to the directory of this code file
source("helperFunctions/functions_CreatingDifferentDataInputs.R")

# 1. Single cell data to .txt files for Scaden ####
seuset <- readRDS("input/full.43p_10X.integrated.cleaned.seurat.RDS")
length(unique(seuset$Patient)) #49 patients
seuset$CellType <- Idents(seuset) #add cell identity to metadata
unique(seuset@meta.data %>% filter(Tissue=="plaque") %>% pull(CellType))
table(seuset@meta.data %>% filter(Tissue=="plaque") %>% pull(Patient)) #ae4675 has only 1 sequenced cell! Scaden needs +1 cell/patient --> remove it
table(seuset@meta.data %>% filter(Tissue=="plaque") %>% pull(CellType)) %>% as.data.frame() %>% 
  mutate(CellType = fct_reorder(Var1, desc(Freq))) %>% 
  ggplot(aes(x=Freq, y=CellType)) + geom_col() + geom_vline(xintercept = 50)
samples <- rownames(seuset@meta.data %>% filter(Tissue=="plaque") %>% 
                      filter(CellType!="CD79+ B Cells", CellType!="GYPA+ Erythroid Cells", CellType!="CD79+ Plasma B Cells") %>% 
                      filter(Patient != "4675"))
length(samples) #11715 cells
matrixCounts <- seuset@assays$RNA@counts[,samples] #all integers
metaData <- seuset@meta.data[samples,] %>% mutate(SubjectName = Patient) %>% 
  rownames_to_column("SampleID") 

clustering <- list()
clustering["T Cells"] <- list(str_subset(levels(metaData$CellType), "(?:T) Cells"))
clustering["NK Cells"] <- list(str_subset(levels(metaData$CellType), "(?:NK) Cells"))
clustering["SB Cells"] <- list(str_subset(levels(metaData$CellType), "B Cells"))
clustering["Inflammatory Macrophages"] <- list(str_subset(levels(metaData$CellType), "Inflammatory Monocyte")) #4 subtypes
clustering["Resident Macrophages"] <- list(str_subset(levels(metaData$CellType), "Resident")) #2 subtypes
clustering["LAM"] <- list(str_subset(levels(metaData$CellType), "Lipid Associated")) #1 subtype.. doublecheck!
clustering["Foamy Macrophages"] <- list(str_subset(levels(metaData$CellType), "Foamy")) #3 subtypes
clustering["Dendritic Cells"] <- list(str_subset(levels(metaData$CellType), "DC"))
clustering["Endothelial Cells"] <- list(str_subset(levels(metaData$CellType), "Endothelial Cells$"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(metaData$CellType), "ACTA2+"))
clustering["Mastocytes"] <- list(str_subset(levels(metaData$CellType), "Mast"))
unlist(unname(clustering)) #we still have 25 sub cell types (sanity check)
setdiff(unlist(unname(clustering)), unique(seuset@meta.data %>% filter(Tissue=="plaque") %>% pull(CellType)))

clustering
clType <- as.character(metaData$CellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
metaData$specialCluster <- factor(clType, levels = names(clustering))

rm(list=setdiff(ls(), c("metaData", "matrixCounts", "asScadenTxtInput_scByPatient"))) #free storage space
asScadenTxtInput_scByPatient(matrixCounts, metaData, "output/amsSC_4macros/byPatient_cellType/", different3Clusters = F, cluster = "specialCluster")
  ## median lib size = 2407
dim(matrixCounts)

library(Biobase)
metaData <- metaData %>% column_to_rownames("SampleID")
sc <- ExpressionSet(assayData = as.matrix(matrixCounts), 
              phenoData = AnnotatedDataFrame(data = metaData, varMetadata = data.frame(labelDescription = colnames(metaData), 
                                                                                    row.names = colnames(metaData))))
bulk <- readRDS("output/ourBulk/onlyOriginalSamples/linearBulk_eset.rds")
resMusic <- music_prop(bulk, sc, clusters = "specialCluster", samples = "Patient")
resMusic$Est.prop.weighted

### plot reduced umap:
seuset <- seuset[,samples]
ct <- Idents(seuset)
ct_reduced <- case_when(ct %in% c("CD1C+ cDC1", "CLEC9A+ cDC2")~"Dendritic Cells", 
                        ct %in% c("CD8+ T Cells", "FOXP3+ T Cells", "CD3+MKI67+ Proliferating T Cells", 
                                  "CD4+ T Cells", "CD56+KLRC1+ NK Cells", "CD56-CD16+ NK Cells")~"T and NK Cells", 
                        ct=="CD79+ Class-switched Memory B Cells"~"CD79+ Class-switched Memory B Cells", 
                        ct=="ACTA2+ Smooth Muscle Cells"~"ACTA2+ Smooth Muscle Cells", 
                        ct=="KIT+ Mast Cells"~"KIT+ Mast Cells", 
                        ct %in% c("CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages", 
                                  "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages",
                                  "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages",
                                  "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages")~"Inflammatory Macrophages", 
                        ct %in% c("CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", 
                                  "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages",
                                  "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages")~"Resident Macrophages", 
                        ct %in% c("CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages", 
                                  "CD14+TREM2-OLR1+ABCA+ Foamy Macrophages",
                                  "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages")~"Foamy Macrophages", 
                        ct == "CD34+ Endothelial Cells"~"CD34+ Endothelial Cells")
numCellType_reduced <- case_when(ct_reduced=="T and NK Cells"~1, 
                                 ct_reduced=="CD79+ Class-switched Memory B Cells"~2,
                                 ct_reduced=="Foamy Macrophages"~3,
                                 ct_reduced=="Inflammatory Macrophages"~4,
                                 ct_reduced=="Resident Macrophages"~5,
                                 ct_reduced=="Dendritic Cells"~6,
                                 ct_reduced=="ACTA2+ Smooth Muscle Cells"~7,
                                 ct_reduced=="CD34+ Endothelial Cells"~8,
                                 ct_reduced=="KIT+ Mast Cells"~9)
newLabel <- paste0(numCellType_reduced, ": ", ct_reduced)
Idents(seuset) <- numCellType_reduced
DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType_reduced), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), aspect.ratio = 1)


# 2. Analysis proportions:
boxplotAnalysis <- function(dataset, feature, cellType, pvalues=T, colour="white"){
  if(pvalues){
    require(rstatix)
    require(ggpubr)
    swr <- Vectorize(function(string, nwrap=20){paste(strwrap(string, width=nwrap), collapse="\n")}) #useful for compacting the names of the titles on the ggplot
    stat.test <- (dataset %>% filter(is.na(.data[[feature]])==FALSE)) %>% pairwise_wilcox_test(as.formula(paste0("`", cellType, "` ~ ", feature)))
    stat.test <- stat.test %>% add_xy_position(x = feature)
    return(ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
             geom_boxplot(fill=colour)+stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, step.increase = 0.02, size = 2)+
             theme_bw()+theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1.1))+ylab(paste0(swr(cellType), " proportion")))
  }
  swr <- Vectorize(function(string, nwrap=20){paste(strwrap(string, width=nwrap), collapse="\n")}) #useful for compacting the names of the titles on the ggplot
  return(ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
           geom_boxplot(fill=colour)+
           theme_bw()+theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+ylab(paste0(swr(cellType), " proportion")))
}

folderPredictFiles <- "results/scaden/amsSC/celltype/"
list.files(folderPredictFiles, "*_predictions")
propsRndSeeds_mean <- list()
for(f in list.files(folderPredictFiles, "*_predictions")){
  propsRndSeeds_mean[[f]] <- read_delim(paste0(folderPredictFiles, f)) %>% dplyr::rename(Patient = "...1") %>% 
    mutate(RndSeed=str_match(f, "seed(.*)_scaden")[2]) #str_match --> get random seed from file name
}
propsRndSeeds_mean<- bind_rows(propsRndSeeds_mean) %>% mutate(RndSeed=as.character(RndSeed)) #all the predictions within a single dataframe

for(cellType in sort(colnames(propsRndSeeds_mean)[2:10])){
  subdata <- propsRndSeeds_mean %>% select(Patient, cellType, RndSeed) %>% pivot_wider(names_from = RndSeed, values_from = cellType)  
  pearsons <- ggstatsplot::ggcorrmat(data = subdata, output = "dataframe", type = "parametric") #parametric for Pearson
  message(paste0("Mean correlation for ", cellType, ": ", round(mean(pearsons$estimate),3)))
}
propsRndSeeds_mean <- (propsRndSeeds_mean %>% group_by(Patient) %>% summarise_if(is.numeric, list(mean)))
propsRndSeeds_mean <- propsRndSeeds_mean %>% mutate_if(is.numeric, ~(.*100))

propsRndSeeds_mean %>% pivot_longer(-Patient, names_to = "CellType", values_to = "Percentage") %>% 
  ggplot(aes(x=CellType,y=Percentage,fill=CellType)) +
  geom_boxplot(outlier.size = 0.2) + ylim(0,100) + theme_light() +
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1), legend.position = "none") 
swr <- Vectorize(function(string, nwrap=20){paste(strwrap(string, width=nwrap), collapse="\n")}) #useful for compacting the names of the titles on the ggplot

propsRndSeeds_mean %>% select(Patient, ends_with("Macrophages")) %>% pivot_longer(-Patient, names_to = "CellType", values_to = "Percentage") %>% 
  mutate(CellType=swr(CellType)) %>% 
  ggplot(aes(x=CellType,y=Percentage,fill=CellType)) +
  geom_boxplot(outlier.size = 0.2) + theme_light() +
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "none") 

CreateTableOne(data=propsRndSeeds_mean %>% select(ends_with("Macrophages")))

histoClinical <- read_delim("input/clinical_data_good.txt") %>% 
  mutate(symptoms = case_when(symptoms_4g=="asymptomatic" ~ 0, symptoms_4g=="ocular" ~ 0, symptoms_4g=="TIA" ~ 1, symptoms_4g=="stroke"~1)) %>% 
  mutate(MACE = case_when(epmajor_3years=="Included" ~ 0, epmajor_3years=="Excluded" ~ 1)) #654 patients
props_cl <- merge(propsRndSeeds_mean, histoClinical, by.x = "Patient", by.y = "study_number")

props_Reduced <- propsRndSeeds_mean %>% mutate(Macro=propsRndSeeds_mean %>% select(ends_with("Macrophages")) %>% rowSums()) %>% 
  mutate(Dendritic=propsRndSeeds_mean %>% select(contains("cDC")) %>% rowSums()) %>% 
  mutate(`T and NK`=propsRndSeeds_mean %>% select(ends_with("T Cells"), ends_with("NK Cells")) %>% rowSums()) %>% 
  select(-ends_with("T Cells"), -ends_with("NK Cells"), -ends_with("Macrophages"), -contains("cDC")) 

props_Reduced %>% 
  pivot_longer(-Patient, names_to = "CellType", values_to = "Percentage") %>% 
  ggplot(aes(x=CellType,y=Percentage,fill=CellType)) +
  geom_boxplot(outlier.size = 0.2) + ylim(0,100) + theme_light() +
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.position = "none") 

## histological agreement:
props_Reduced_cl <- merge(props_Reduced, histoClinical, by.x = "Patient", by.y = "study_number")  
plot_grid(boxplotAnalysis(props_Reduced_cl %>% mutate(smc=factor(smc, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))), 
                          "smc", "Smooth Muscle Cells", pvalues = F, "#00A9FF"), boxplotAnalysis(props_Reduced_cl, "plaquephenotype", "Smooth Muscle Cells", F, "#00A9FF"))
plot_grid(boxplotAnalysis(props_Reduced_cl %>% mutate(macrophages=factor(macrophages, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))), 
                          "macrophages", "Macro", F, "#00BE67"), boxplotAnalysis(props_Reduced_cl, "plaquephenotype", "Macro", F, "#00BE67"))

for(subMacro in colnames(props_cl %>% select(ends_with("Macrophages", ignore.case = F)))){
  print(plot_grid(boxplotAnalysis(props_cl %>% mutate(macrophages=factor(macrophages, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))), 
                            "macrophages", subMacro, F, "#00BE67"), boxplotAnalysis(props_cl, "plaquephenotype", subMacro, F, "#00BE67")))
  
}

## association with MACE:
### reduced populations:
for(cell in colnames(props_cl)[2:9]){
  print(cell)
  print(paste0("p-value = ", summary(lm(as.formula(paste0("`", cell, "`~symptoms")), data = props_cl))$coefficients[2,"Pr(>|t|)"]))
  print(paste0("beta coeff (MACE) = : ", summary(lm(as.formula(paste0("`", cell, "`~symptoms")), data = props_cl))$coefficients[2,"Estimate"]))
  print("--------")
}

### look only macrophage subtypes:
table(props_cl$med_statin_derived)
for(dataset in list(props_cl %>% filter(med_statin_derived=="no"), props_cl %>% filter(med_statin_derived=="yes"))){
  dataset <- dataset %>% select(c("CD14+TREM2-OLR1+ABCA+ Foamy Macrophages", "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages", "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages", 
                                  "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages", "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages", "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages", "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages", 
                                  "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages", "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", "MACE", "symptoms"))
  univariate <- data.frame(row.names = NULL)
  for(cell in colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F)))){
    for(f in c("MACE", "symptoms")){
      coefficientsLm <- summary(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))$coefficients
      confInterv <- confint(lm(as.formula(paste("`", cell, "`", "~", f,  sep="")), dataset))
      for(categ in rownames(coefficientsLm)[-1]){
        univariate[categ, paste0(cell, ":beta [CI]")] <- paste0(round(coefficientsLm[categ, "Estimate"], 4), " [", paste(round(confInterv[categ,], 4), collapse = " to "), "]")
        univariate[categ, paste0(cell, ":pvalue")] <- round(coefficientsLm[categ, "Pr(>|t|)"], 4)
      }
    }
  }
  header_df = tibble(col_keys = c("Feature", names(univariate)), firstHead = c("Feature", rep(colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F))), each=2)), secondHead = c("Feature", rep(c("beta [CI]", "p-value"), length(colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F)))))))
  print(univariate %>% rownames_to_column("Feature") %>% mutate(Feature=c("MACE (yes vs no)", "Symptoms (yes vs no)"))%>% flextable() %>% set_header_df(header_df) %>%
    merge_v(part="header") %>% merge_h(part="header") %>% theme_box() %>% align(align = "center", part = "header") %>% autofit() %>% fontsize(size = 8, part = "all"))
  
}


# ## Lipid associated macrophages --> histology ?
# plot_grid(boxplotAnalysis(props_cl %>% mutate(macrophages=factor(macrophages, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))), 
#                           "macrophages", "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", F, "#00BE67"), 
#           boxplotAnalysis(props_cl, "plaquephenotype", "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", F, "#00BE67"), 
#           boxplotAnalysis(props_cl, "fat", "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", F, "#00BE67"))


# 3. Statines use or LDL association
univariate <- data.frame(row.names = NULL)
dataset <- props_cl %>% select(c("CD14+TREM2-OLR1+ABCA+ Foamy Macrophages", "CD14+TREM2-OLR1+NLRP3+ Inflammatory Foamy Macrophages", "CD14+TREM2-TIMP1+HSPA6+ Lipid-stress Activated Foamy Macrophages", 
                    "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages", "CD14+IL1B+SELL+MX1+ Interferon Activated Inflammatory Monocyte-derived Macrophages", "CD14+-IL1B+SELL+CD16+ Migrating Inflammatory Monocyte-derived Macrophages", "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages", 
                    "CD14+TNF+TREM2+FOLR2+ Inflammatory Resident-like Lipid Associated Macrophages", "CD14+IL1B-TREM2-FOLR2+ Resident-like Macrophages", "CD14+TREM2+FOLR2-ABCG+ Lipid Associated Macrophages", 
                    "med_statin_derived", "ldl_final"))
for(cell in colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F)))){
  for(f in c("med_statin_derived", "ldl_final")){
    coefficientsLm <- summary(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))$coefficients
    confInterv <- confint(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))
    for(categ in rownames(coefficientsLm)[-1]){
      univariate[categ, paste0(cell, ":beta [CI]")] <- paste0(round(coefficientsLm[categ, "Estimate"], 4), " [", paste(round(confInterv[categ,], 4), collapse = " to "), "]")
      univariate[categ, paste0(cell, ":pvalue")] <- round(coefficientsLm[categ, "Pr(>|t|)"], 4)
    }
  }
}
header_df = tibble(col_keys = c("Feature", names(univariate)), firstHead = c("Feature", rep(colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F))), each=2)), secondHead = c("Feature", rep(c("beta [CI]", "p-value"), length(colnames(dataset %>% select(ends_with("Macrophages", ignore.case = F)))))))
univariate %>% rownames_to_column("Feature") %>% mutate(Feature=c("statins use (yes vs no)", "LDL"))%>% flextable() %>% set_header_df(header_df) %>%
  merge_v(part="header") %>% merge_h(part="header") %>% theme_box() %>% align(align = "center", part = "header") %>% autofit() %>% fontsize(size = 8, part = "all")

## with Huize Pan data? ####
folderPredictFiles <- "results/scaden/ourBulk_hpSC/"
list.files(folderPredictFiles, "*_predictions")
propsRndSeeds_mean <- list()
for(f in list.files(folderPredictFiles, "*_predictions")){
  propsRndSeeds_mean[[f]] <- read_delim(paste0(folderPredictFiles, f)) %>% dplyr::rename(Patient = "...1") %>% 
    mutate(RndSeed=str_match(f, "seed(.*)_scaden")[2]) #str_match --> get random seed from file name
}
propsRndSeeds_mean<- bind_rows(propsRndSeeds_mean) %>% mutate(RndSeed=as.character(RndSeed)) #all the predictions within a single dataframe
propsRndSeeds_mean <- (propsRndSeeds_mean %>% group_by(Patient) %>% summarise_if(is.numeric, list(mean)))
propsRndSeeds_mean <- propsRndSeeds_mean %>% mutate_if(is.numeric, ~(.*100))

histoClinical <- read_delim("input/clinical_data_good.txt") %>% 
  mutate(symptoms = case_when(symptoms_4g=="asymptomatic" ~ 0, symptoms_4g=="ocular" ~ 0, symptoms_4g=="TIA" ~ 1, symptoms_4g=="stroke"~1)) %>% 
  mutate(MACE = case_when(epmajor_3years=="Included" ~ 0, epmajor_3years=="Excluded" ~ 1)) #654 patients
props_cl <- merge(propsRndSeeds_mean, histoClinical, by.x = "Patient", by.y = "study_number")

univariate <- data.frame(row.names = NULL)
dataset <- props_cl %>% select(c("Macrophage 1", "Macrophage 2", "Macrophage 3","med_statin_derived", "ldl_final"))
for(cell in colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F)))){
  for(f in c("med_statin_derived", "ldl_final")){
    coefficientsLm <- summary(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))$coefficients
    confInterv <- confint(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))
    for(categ in rownames(coefficientsLm)[-1]){
      univariate[categ, paste0(cell, ":beta [CI]")] <- paste0(round(coefficientsLm[categ, "Estimate"], 4), " [", paste(round(confInterv[categ,], 4), collapse = " to "), "]")
      univariate[categ, paste0(cell, ":pvalue")] <- round(coefficientsLm[categ, "Pr(>|t|)"], 4)
    }
  }
}
header_df = tibble(col_keys = c("Feature", names(univariate)), firstHead = c("Feature", rep(colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F))), each=2)), secondHead = c("Feature", rep(c("beta [CI]", "p-value"), length(colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F)))))))
univariate %>% rownames_to_column("Feature") %>% mutate(Feature=c("statines (yes vs no)", "LDL"))%>% flextable() %>% set_header_df(header_df) %>%
  merge_v(part="header") %>% merge_h(part="header") %>% theme_box() %>% align(align = "center", part = "header") %>% autofit() %>% fontsize(size = 8, part = "all")

univariate <- data.frame(row.names = NULL)
dataset <- props_cl %>% select(c("Macrophage 1", "Macrophage 2", "Macrophage 3","MACE", "symptoms"))
for(cell in colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F)))){
  for(f in c("MACE", "symptoms")){
    coefficientsLm <- summary(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))$coefficients
    confInterv <- confint(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), dataset))
    for(categ in rownames(coefficientsLm)[-1]){
      univariate[categ, paste0(cell, ":beta [CI]")] <- paste0(round(coefficientsLm[categ, "Estimate"], 4), " [", paste(round(confInterv[categ,], 4), collapse = " to "), "]")
      univariate[categ, paste0(cell, ":pvalue")] <- round(coefficientsLm[categ, "Pr(>|t|)"], 4)
    }
  }
}
header_df = tibble(col_keys = c("Feature", names(univariate)), firstHead = c("Feature", rep(colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F))), each=2)), secondHead = c("Feature", rep(c("beta [CI]", "p-value"), length(colnames(dataset %>% select(starts_with("Macrophage", ignore.case = F)))))))
univariate %>% rownames_to_column("Feature") %>% mutate(Feature=c("MACE (yes vs no)", "symptoms (yes vs no)"))%>% flextable() %>% set_header_df(header_df) %>%
  merge_v(part="header") %>% merge_h(part="header") %>% theme_box() %>% align(align = "center", part = "header") %>% autofit() %>% fontsize(size = 8, part = "all")


load("/Volumes/DATASHUR_P2/1stPart_redone/input/Pan.seurat.RData")
hpSC <- seuset
seuset@assays$RNA
#LAM markers:
VlnPlot(assay = "RNA", hpSC, features=c("TREM2", "CD9", "GPNMB"), idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3")) + theme(legend.position = "none")
#RES markers:
VlnPlot(assay = "RNA", hpSC, features=c("LYVE1", "MRC1", "FOLR2"), idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3")) + theme(legend.position = "none")
#FOAM markers:
VlnPlot(assay = "RNA", hpSC, features=c("OLR1", "ABCA1", "ABCG1"), idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3")) + theme(legend.position = "none")
#INFLAM markers:
VlnPlot(assay = "RNA", hpSC, features=c("IL1B", "TNF", "CASP1"), idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3")) + theme(legend.position = "none")
#MYELOID markers:
VlnPlot(assay = "RNA", hpSC, features=c("CD14", "FCGR1A", "FCGR3A"), idents = c("Macrophage 1", "Macrophage 2", "Macrophage 3")) + theme(legend.position = "none")



