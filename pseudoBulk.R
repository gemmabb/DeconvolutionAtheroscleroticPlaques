###########################################################################
#                                                                         #
# DECONVOLUTION ON ARTIFICIAL BULK DATA                                   #
#                                                                         #
# Gemma Bel Bordes, Nov 20th 2022                                         #
#                                                                         #
# Adapted from:                                                           #
# Avila Cobos F, Alquicira-Hernandez J, Powell JE, Mestdagh P and         #
# De Preter K. Benchmarking of cell type deconvolution pipelines for      #  
# transcriptomics data. (Nature Communications;                           #
# https://doi.org/10.1038/s41467-020-19015-1)                             #  
###########################################################################

library(Biobase)
library(tidyverse)
library(magrittr)
library(MuSiC)
library(BisqueRNA)
library(viridis)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("helperFunctions/functionsPseudoBulk.R")
source("helperFunctions/functions_CreatingDifferentDataInputs.R")

sc <- readRDS("output/ourSC/linear_scAllGenessc_eset.rds")
dim(exprs(sc)) #4903 cells
bulk <- readRDS("output/ourBulk/onlyOriginalSamples/linearBulk_eset.rds")

# 1 Pseudo bulk construction #### we'll use 50% of the sc dataset for the bulk simulation (training) and 50% for the posterior deconvolution (testing)
#scIntoTrainAndTest(sc, "pseudoBulkData/DataInput/") # 2329 cells and 46 patients for training -- 2574 cells and 46 patients for testing
scTrain <- readRDS("pseudoBulkData/DataInput/trainSC.eset")
dim(scTrain)

addTaskCallback(function(...) {set.seed(19599); TRUE}) #set seed for the entire R session... to make results reproducible
simulationAllSCGenes_18 <- simulateBulk(numMixtures = 700, numCells = 100, scEset = scTrain, cellCluster = "cellType")
saveEsetAndProps(simulationAllSCGenes_18, "pseudoBulkData/DataInput/")

## Loading the pseudo bulk data ####
bulk18_allGenes <- readRDS("pseudoBulkData/DataInput/pseudoBulk_18.eset") #20111 genes --> only using this one!

props18 <- read_delim("pseudoBulkData/DataInput/pseudoProps_18.txt")
props9 <- as_tibble(as.data.frame(t(as.data.frame(t(props18 %>% column_to_rownames("cellType")%>% select(-starts_with("numCells"), -starts_with("simNum")))) %>% 
                                      mutate(`CD3+ T Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")) %>% 
                                      #mutate(`FOXP3+ T Cells`=`FOXP3+ T Cells`) %>% select(-c("FOXP3+ T Cells")) %>% 
                                      mutate(`NK Cells`=select(., c("CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
                                      mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
                                      mutate(`CD68+ (Foam Cells and inflam/resident macrophages)`=select(., c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
                                      mutate(`CD68+ Dendritic`=`CD68+CD1C+ Dendritic Cells`) %>% select(-c("CD68+CD1C+ Dendritic Cells")) %>% 
                                      mutate(`Endothelial Cells`=select(., c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% rowSums()) %>% select(-c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% 
                                      mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
                                      mutate(`CD68+ Mast Cells`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`))) %>% rownames_to_column("cellType"))
props6 <- as_tibble(as.data.frame(t(as.data.frame(t(props18 %>% column_to_rownames("cellType")%>% select(-starts_with("numCells"), -starts_with("simNum")))) %>% 
  mutate(`T and NK Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
  mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
  mutate(`CD68+ Cells (no mast)`=select(., c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
  mutate(`Endothelial Cells`=select(., c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% rowSums()) %>% select(-c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% 
  mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
  mutate(`Mast Cells`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`))) %>% rownames_to_column("cellType"))

props18 %>% write_csv("dataPlots/pseudoProps18.csv")
props9 %>% write_csv("dataPlots/pseudoProps9.csv")
props6 %>% write_csv("dataPlots/pseudoProps6.csv")
#props9 <- read_delim("pseudoBulkData/AllSCGenes/pseudoProps_9.txt")
#props6 <- read_delim("pseudoBulkData/AllSCGenes/pseudoProps_6.txt")

### Reminder!!!! Deconvolutions must use the test sc dataset!
scTest <- readRDS("pseudoBulkData/DataInput/testSC.eset")
dim(scTest)

### Pseudo and reference data transformation ####
# linear (as it is):
bulk18_allGenes_linear <- bulk18_allGenes
scTest_linear <- scTest
rm(scTest, bulk18_allGenes) #remove from environment (+ storage)
# log transformed:
asExpressionSet_bulk(transformTo(exprs(bulk18_allGenes), "log"), "pseudoBulkData/DataInput/", "pseudoBulk_18_log.eset")
bulk18_allGenes_log <- readRDS("pseudoBulkData/DataInput/pseudoBulk_18_log.eset")
asExpressionSet_sc(transformTo(exprs(scTest), "log"), scTest@phenoData@data, "pseudoBulkData/DataInput", "testSC_log.eset")
scTest_log <- readRDS("pseudoBulkData/DataInput/testSC_log.eset")
# sqrt transformed:
asExpressionSet_bulk(transformTo(exprs(bulk18_allGenes), "sqrt"), "pseudoBulkData/DataInput/", "pseudoBulk_18_sqrt.eset")
bulk18_allGenes_sqrt <- readRDS("pseudoBulkData/DataInput/pseudoBulk_18_sqrt.eset")
asExpressionSet_sc(transformTo(exprs(scTest), "sqrt"), scTest@phenoData@data, "pseudoBulkData/DataInput/", "testSC_sqrt.eset")
scTest_sqrt <- readRDS("pseudoBulkData/DataInput/testSC_sqrt.eset")
# vstransformed:
asExpressionSet_bulk(transformTo(exprs(bulk18_allGenes), "vst"), "pseudoBulkData/DataInput/", "pseudoBulk_18_vst.eset")
bulk18_allGenes_vst <- readRDS("pseudoBulkData/DataInput/pseudoBulk_18_vst.eset")
asExpressionSet_sc(transformTo(exprs(scTest), "vst"), scTest@phenoData@data, "pseudoBulkData/DataInput/", "testSC_vst.eset")
scTest_vst <- readRDS("pseudoBulkData/DataInput/testSC_vst.eset")
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
## function: y = a/x + b, and a local regression fit was automatically substituted. Specify fitType='local' or 'mean' to avoid this message next time.
## VST gives us negative values for single cell --> we can't use that because it'll give us errors with the deconv.
hist(exprs(scTest_vst[,1])) #include in the supplements
hist(exprs(bulk18_allGenes_vst)[,1]) #compare it with bulk vstranformation

# 1. List of data and different parameters to analyze ####
listPseudoBulk_trans <- list(bulk18_allGenes_linear, bulk18_allGenes_log, bulk18_allGenes_sqrt) #to assess data transformation
rm(bulk18_allGenes_linear, bulk18_allGenes_log, bulk18_allGenes_sqrt, bulk18_allGenes_vst) #remove from environment (+ storage)
listSC_trans <- list(scTest_linear, scTest_log, scTest_sqrt)
rm(scTest_linear, scTest_log, scTest_sqrt, scTest_vst)
listTrans <- list("linear", "log", "sqrt")

listProps <- list(props18, props9, props6)
rm(props18, props9, props6) #remove from environment (+ storage)
listClustering <- list("cellType", "clusterType_9", "clusterType_6") #to assess cell clustering from the reference

# 2. MuSiC deconvolution ####
## Different data transformations:
for(transfIndex in 1:3){
  print(listTrans[[transfIndex]])
  for(clustIndex in 1:3){
    print(listClustering[[clustIndex]])
    musicDeconvolution(listPseudoBulk_trans[[transfIndex]], listSC_trans[[transfIndex]], listClustering[[clustIndex]]) %>% 
      write_csv(paste0("pseudoBulkData/DataTransformation_props/musicProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
  }
}

## Different clustering:
for(i in 1:3){
  musicDeconvolution(listPseudoBulk_trans[[1]], listSC_trans[[1]], listClustering[[i]]) %>% 
    write_csv(paste0("pseudoBulkData/Clustering_props/musicProps_", listClustering[[i]], ".csv"))
}
### to read:
musicMetrics_dataTrans <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                               data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(musicMetrics_dataTrans) <- c("metrics", "props")
musicMetrics_clust <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                               data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(musicMetrics_clust) <- c("metrics", "props")
for(clustIndex in 1:3){
  print(clustIndex)
  for(transfIndex in 1:3){
    musicResults <- read_csv(paste0("pseudoBulkData/DataTransformation_props/musicProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
    musicMetrics_dataTrans$metrics <- rbind(musicMetrics_dataTrans$metrics, musicResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
    musicMetrics_dataTrans$props <- rbind(musicMetrics_dataTrans$props, musicResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
  }
  musicResults <- read_csv(paste0("pseudoBulkData/Clustering_props/musicProps_", listClustering[[clustIndex]], ".csv"))
  musicMetrics_clust$metrics <- rbind(musicMetrics_clust$metrics, musicResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(Clustering=listClustering[[clustIndex]]))
  musicMetrics_clust$props <- rbind(musicMetrics_clust$props, musicResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
}
rm(musicResults)

# 2. NNLS deconvolution ####
## Different data transformations:
for(transfIndex in 1:3){
  print(listTrans[[transfIndex]])
  for(clustIndex in 1:3){
    print(listClustering[[clustIndex]])
    NNLSDeconvolution(listPseudoBulk_trans[[transfIndex]], listSC_trans[[transfIndex]], listClustering[[clustIndex]]) %>% 
      write_csv(paste0("pseudoBulkData/DataTransformation_props/nnlsProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
  }
}
## Different clustering:
for(i in 1:3){
  NNLSDeconvolution(listPseudoBulk_trans[[1]], listSC_trans[[1]], listClustering[[i]]) %>% 
    write_csv(paste0("pseudoBulkData/Clustering_props/nnlsProps_", listClustering[[i]], ".csv"))
}
### to read:
nnlsMetrics_dataTrans <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                               data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(nnlsMetrics_dataTrans) <- c("metrics", "props")
nnlsMetrics_clust <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                           data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(nnlsMetrics_clust) <- c("metrics", "props")
for(clustIndex in 1:3){
  print(clustIndex)
  for(transfIndex in 1:3){
    nnlsResults <- read_csv(paste0("pseudoBulkData/DataTransformation_props/nnlsProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
    nnlsMetrics_dataTrans$metrics <- rbind(nnlsMetrics_dataTrans$metrics, nnlsResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
    nnlsMetrics_dataTrans$props <- rbind(nnlsMetrics_dataTrans$props, nnlsResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
  }
  nnlsResults <- read_csv(paste0("pseudoBulkData/Clustering_props/nnlsProps_", listClustering[[clustIndex]], ".csv"))
  nnlsMetrics_clust$metrics <- rbind(nnlsMetrics_clust$metrics, nnlsResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(Clustering=listClustering[[clustIndex]]))
  nnlsMetrics_clust$props <- rbind(nnlsMetrics_clust$props, nnlsResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
}
rm(nnlsResults)

# 4. Bisque deconvolution ####
## Different data transformations:
for(transfIndex in 1:3){
  print(listTrans[[transfIndex]])
  for(clustIndex in 1:3){
    print(listClustering[[clustIndex]])
    bisqueDeconvolution(listPseudoBulk_trans[[transfIndex]], listSC_trans[[transfIndex]], listClustering[[clustIndex]]) %>% 
      write_csv(paste0("pseudoBulkData/DataTransformation_props/bisqueProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
  }
}
## Different clustering:
for(i in 1:3){
  bisqueDeconvolution(listPseudoBulk_trans[[1]], listSC_trans[[1]], listClustering[[i]]) %>% 
    write_csv(paste0("pseudoBulkData/Clustering_props/bisqueProps_", listClustering[[i]], ".csv"))
}
### to read:
bisqueMetrics_dataTrans <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                               data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(bisqueMetrics_dataTrans) <- c("metrics", "props")
bisqueMetrics_clust <- list(data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "dataTrans", "Clustering")))), 
                           data.frame(matrix(ncol=6, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "dataTrans", "Clustering")))))
names(bisqueMetrics_clust) <- c("metrics", "props")
for(clustIndex in 1:3){
  print(clustIndex)
  for(transfIndex in 1:3){
    bisqueResults <- read_csv(paste0("pseudoBulkData/DataTransformation_props/bisqueProps_", listTrans[[transfIndex]], "_", listClustering[[clustIndex]], ".csv"))
    bisqueMetrics_dataTrans$metrics <- rbind(bisqueMetrics_dataTrans$metrics, bisqueResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
    bisqueMetrics_dataTrans$props <- rbind(bisqueMetrics_dataTrans$props, bisqueResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(dataTrans=listTrans[[transfIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
  }
  bisqueResults <- read_csv(paste0("pseudoBulkData/Clustering_props/bisqueProps_", listClustering[[clustIndex]], ".csv"))
  bisqueMetrics_clust$metrics <- rbind(bisqueMetrics_clust$metrics, bisqueResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% checkMetrics() %>% mutate(Clustering=listClustering[[clustIndex]]))
  bisqueMetrics_clust$props <- rbind(bisqueMetrics_clust$props, bisqueResults %>% joinSimulProps(propsNumSimulation = listProps[[clustIndex]]) %>% mutate(Clustering=listClustering[[clustIndex]]))
}
rm(bisqueResults)

# 5. Cibersortx deconvolution ####
## bulk datasets to proper .txt file:
for(i in 1:3){
  as.data.frame(listBulkData_allGenes[[i]]@assayData$exprs) %>% rownames_to_column("GeneSymbol") %>% 
    write.table(file = paste0("pseudoBulkData/DataInput/CibersortInput/pseudobulkAllGenes_", listClustering[[i]], ".txt"), sep = "\t", row.names = T)
}
## single cell data: (instead of each cell ID we want the cell type) ####
for(clustering in tail(scTest@phenoData@varMetadata$labelDescription,3)){
  Counts <- as.data.frame(scTest@assayData$exprs)
  names(Counts) <- unlist(unname(scTest@phenoData@data[clustering]))
  Counts <- cbind(GeneSymbol = rownames(Counts), Counts)
  rownames(Counts) <- NULL
  write.table(Counts, file = paste0("pseudoBulkData/DataInput/CibersortInput/testSC_", clustering, ".txt"), sep = "\t", row.names = F)
}

## read and analyze results:
ciberMetrics_clust <- list(data.frame(matrix(ncol=4, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "Clustering")))), 
                            data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "Clustering")))))
names(ciberMetrics_clust) <- c("metrics", "props")
ciberFiles <- c("ciberProps_cellType.csv", "ciberProps_clusterType_9.csv", "ciberProps_clusterType_6.csv")
for(i in 1:3){
  ciberResults <- cibersortReadDeconvolution(paste0("pseudoBulkData/Clustering_props/", ciberFiles[[i]]))
  ciberMetrics_clust$metrics <- rbind(ciberMetrics_clust$metrics, ciberResults %>% joinSimulProps(propsNumSimulation = listProps[[i]]) %>% checkMetrics() %>% mutate(Clustering=str_match(ciberFiles[[i]], "Props_(.*).csv")[2]))
  ciberMetrics_clust$props <- rbind(ciberMetrics_clust$props, ciberResults %>% joinSimulProps(propsNumSimulation = listProps[[i]])%>% mutate(Clustering=str_match(ciberFiles[[i]], "Props_(.*).csv")[2]))
}
rm(ciberResults)

# 6. Scaden deconvolution ####
## bulk datasets to proper .txt file:
as.data.frame(bulk18_allGenes@assayData$exprs) %>%
    write.table(file = paste0("pseudoBulkData/DataInput/ScadenInput/pseudobulkAllGenes_18.txt"), sep = "\t", quote = F, col.names = NA)

## test sc dataset to proper .txt file:
### counts:
source("helperFunctions/functions_CreatingDifferentDataInputs.R")
asScadenTxtInput_scByPatient(scTest@assayData$exprs, scTest@phenoData@data, "pseudoBulkData/DataInput/ScadenInput/", different3Clusters = F, "cellType")
asScadenTxtInput_scByPatient(scTest@assayData$exprs, scTest@phenoData@data, "pseudoBulkData/DataInput/ScadenInput/", different3Clusters = F, "clusterType_9")
asScadenTxtInput_scByPatient(scTest@assayData$exprs, scTest@phenoData@data, "pseudoBulkData/DataInput/ScadenInput/", different3Clusters = F, "clusterType_6")
## median lib size: 2532.5
table(scTest@phenoData@data$SubjectName) #we have to remove from the directory those test samples with only 1 cell (because it can't simulate spare mixtures like with the others)
## --> remove ae4440 and ae4530 !!!

## read and analyze results:
scadenMetrics_clust <- list(data.frame(matrix(ncol=4, nrow=0, dimnames=list(NULL, c("resultsFor", "NRMSE", "Correlation", "Clustering")))), 
                           data.frame(matrix(ncol=5, nrow=0, dimnames=list(NULL, c("mixture", "cellType", "deconvProps", "simulProps", "Clustering")))))
names(scadenMetrics_clust) <- c("metrics", "props")
scadenFiles <- c("scadenProps_cellType.txt", "scadenProps_clusterType_9.txt", "scadenProps_clusterType_6.txt")
for(i in 1:3){
  scadenResults <- scadenReadDeconvolution(paste0("pseudoBulkData/Clustering_props/", scadenFiles[[i]]))
  scadenMetrics_clust$metrics <- rbind(scadenMetrics_clust$metrics, scadenResults %>% joinSimulProps(propsNumSimulation = listProps[[i]]) %>% checkMetrics() %>% mutate(Clustering=str_match(scadenFiles[[i]], "Props_(.*).txt")[2]))
  scadenMetrics_clust$props <- rbind(scadenMetrics_clust$props, scadenResults %>% joinSimulProps(propsNumSimulation = listProps[[i]])%>% mutate(Clustering=str_match(scadenFiles[[i]], "Props_(.*).txt")[2]))
}
rm(scadenResults)

# 7. Data transformation analysis:
dataTrans_results <- rbind(musicMetrics_dataTrans$metrics %>% mutate(Algorithm="MuSiC"), 
                           nnlsMetrics_dataTrans$metrics %>% mutate(Algorithm="NNLS"),
                           bisqueMetrics_dataTrans$metrics %>% mutate(Algorithm="Bisque"))
dataTrans_props <- rbind(musicMetrics_dataTrans$props %>% mutate(Algorithm="MuSiC"), 
                           nnlsMetrics_dataTrans$props %>% mutate(Algorithm="NNLS"),
                           bisqueMetrics_dataTrans$props %>% mutate(Algorithm="Bisque"))

dataTrans_results %>% filter(resultsFor!="all") %>% ggplot(aes(x=Algorithm, y=NRMSE, fill=dataTrans)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6)
dataTrans_results %>% filter(resultsFor=="all") %>% ggplot(aes(x=Algorithm, y=NRMSE, fill=dataTrans)) + geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.7) + theme_bw() #report

# 8. Clustering analysis: ####
## for 18 cell types
clust_results <- rbind(musicMetrics_clust$metrics %>% mutate(algorithm="MuSiC") %>% dplyr:::rename(cellType=resultsFor), 
                                      nnlsMetrics_clust$metrics %>% mutate(algorithm="NNLS") %>% dplyr:::rename(cellType=resultsFor),
                                      bisqueMetrics_clust$metrics %>% mutate(algorithm="Bisque") %>% dplyr:::rename(cellType=resultsFor), 
                                      ciberMetrics_clust$metrics %>% mutate(algorithm="CIBERSORTx") %>% dplyr:::rename(cellType=resultsFor), 
                                      scadenMetrics_clust$metrics %>% mutate(algorithm="Scaden") %>% dplyr:::rename(cellType=resultsFor))
clust_props <- rbind(musicMetrics_clust$props %>% mutate(algorithm="MuSiC"), 
                       nnlsMetrics_clust$props %>% mutate(algorithm="NNLS"),
                       bisqueMetrics_clust$props %>% mutate(algorithm="Bisque"), 
                       ciberMetrics_clust$props %>% mutate(algorithm="CIBERSORTx"), 
                       scadenMetrics_clust$props %>% mutate(algorithm="Scaden"))
clust_results %>% write_csv("dataPlots/pseudoBulk_Clust.csv")
## 18 cell types:
plot18 <- ggplot(clust_results %>% filter(Clustering=="cellType") %>% mutate(cellType=factor(cellType, levels=c("CD79A+ Class-switched Memory B Cells", "CD68+KIT+ Mast Cells", "CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I", "FOXP3+ T Cells", "CD3+ T Cells VI", "CD3+ T Cells V", "CD3+ T Cells IV", "CD3+ T Cells III", "CD3+ T Cells II", "CD3+ T Cells I", "CD34+ Endothelial Cells II", "CD34+ Endothelial Cells I", "ACTA2+ Smooth Muscle Cells", "all"))), 
       aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_fill_viridis_c(na.value = "white") + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## 9 cell types:
plot9 <- ggplot(clust_results %>% filter(Clustering=="clusterType_9") %>% mutate(cellType=factor(cellType, levels=c("Switched mem B Cells", "CD68+ Mast Cells", "CD68+ Dendritic", "null1", "CD68+ (Foam Cells and inflam/resident macrophages)", "null2", "null3", "NK Cells", "FOXP3+ T Cells", "null4", "null5", "null6", "CD3+ T Cells", "null7", "null8", "null9", "Endothelial Cells", "Smooth Muscle Cells", "all"))), 
       aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_y_discrete(drop=FALSE)+
  scale_fill_viridis_c(na.value = NA) + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## 6 cell types:
plot6 <- ggplot(clust_results %>% filter(Clustering=="clusterType_6") %>% mutate(cellType=factor(cellType, levels=c("Switched mem B Cells", "Mast Cells", "nullA", "nullB", "CD68+ Cells (no mast)", "nullC", "nullD", "nullE", "nullF", "nullG", "nullH", "nullI", "T and NK Cells", "nullJ", "nullK", "nullL", "Endothelial Cells", "Smooth Muscle Cells", "all"))), aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_y_discrete(drop=FALSE)+
  scale_fill_viridis_c(na.value = NA) + theme(axis.text.y = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## all plots combined... they must have the same 1/NRMSE legend!:
max(1/clust_results$NRMSE)
min(1/clust_results$NRMSE)
egg:::ggarrange(plot18 + ggtitle("   18 clusters") + theme(legend.position = "none", axis.ticks.y=element_blank(), plot.margin = margin(10,0,10,10))+geom_hline(yintercept = 18.5, color="black", size=1.5)+geom_hline(yintercept=seq(from = 1.5, to = 17.5, by = 1), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                   plot9 + ggtitle("   9 clusters") + theme(legend.position = "none", axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = margin(10,0,10,0))+geom_hline(yintercept = 18.5, color="black", size=1.5) + geom_hline(yintercept=c(1.5, 2.5, 3.5, 6.5, 8.5, 9.5, 15.5, 17.5), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                   plot6+ ggtitle("   6 clusters") + scale_y_discrete(labels=case_when(startsWith(levels(plot6$data$cellType), "null")~"", TRUE~levels(plot6$data$cellType)),position = "right", drop=F) + theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = margin(10,10,10,2.5))+geom_hline(yintercept = 18.5, color="black", size=1.5)+geom_hline(yintercept=c(1.5, 2.5, 6.5, 15.5, 17.5), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                   nrow = 1)

# 9. Boxplot deconvolution results vs true proportions ####
swr <- Vectorize(function(string, nwrap=20){paste(strwrap(string, width=nwrap), collapse="\n")}) #useful for compacting the names of the titles on the ggplot
for(clustering in listClustering){
  print(clustering)
  print(rbind(musicMetrics_clust$props %>% filter(Clustering==clustering) %>% mutate(algorithm="MuSiC"), 
        nnlsMetrics_clust$props %>% filter(Clustering==clustering)  %>% mutate(algorithm="NNLS"), 
        bisqueMetrics_clust$props %>% filter(Clustering==clustering)  %>% mutate(algorithm="Bisque"), 
        ciberMetrics_clust$props %>% filter(Clustering==clustering)  %>% mutate(algorithm="CIBERSORTx"), 
        scadenMetrics_clust$props %>% filter(Clustering==clustering)  %>% mutate(algorithm="Scaden")) %>%
    pivot_longer(c("deconvProps", "simulProps"), names_to="typeProps", values_to = "Props") %>%
    mutate(cellType=swr(cellType)) %>% 
    mutate(algorithm=case_when(typeProps=="simulProps"~"Simulation", TRUE~algorithm)) %>% 
    ggplot(aes(x=algorithm,y=Props)) + geom_boxplot(aes(fill=algorithm), outlier.size = 0.2) + theme_bw() +
    facet_wrap(~cellType, nrow=2) + theme(strip.text = element_text(size = 6), axis.title.x = element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank()) +
    ggtitle(paste0("Pseudo-bulk deconvolution, ", clustering)) + ylim(0,1)+ scale_fill_manual(values = c("#8DD3C7", "#FFFFB3","#BEBADA","#FDB462","#80B1D3","#C90606")))
}

# 9. Plot simul vs deconvoluted proportions ####
#### facet wrap cell types for music
musicMetrics_clust$props %>% filter(Clustering=="cellType") %>% mutate(cellType=swr(cellType)) %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free", ncol=3)+theme_bw()+theme(strip.text = element_text(size = 8), legend.position = "none") 

musicMetrics_clust$props %>% filter(Clustering=="clusterType_9") %>% mutate(cellType=swr(cellType)) %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free", ncol=3)+theme_bw()+theme(strip.text = element_text(size = 8), legend.position = "none") 

musicMetrics_clust$props %>% filter(Clustering=="clusterType_6") %>% mutate(cellType=swr(cellType)) %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free", ncol=3)+theme_bw()+theme(strip.text = element_text(size = 8), legend.position = "none") 

View(musicMetrics_clust$metrics %>% filter(Clustering=="cellType"))
View(musicMetrics_clust$metrics %>% filter(Clustering=="clusterType_9"))
View(musicMetrics_clust$metrics %>% filter(Clustering=="clusterType_6"))


# 10. Plot correlation of sc cell types ####
scCPM <- as.data.frame(t(apply(scTest@assayData$exprs, 2, function(x){x*10^6/sum(x)}))) %>% rownames_to_column("Sample")
scCPM_annotated<- as_tibble(merge(scCPM, scTest@phenoData@data %>% select(SampleID, cellType, clusterType_9, clusterType_6), by.x = "Sample", by.y = "SampleID"))

scCPM_cellTypes <- data.frame(row.names = colnames(scCPM %>% select(-Sample)))
for(ct in unique(scCPM_annotated$cellType)){
  print(ct)
  scCPM_cellTypes[[ct]] <- unname((scCPM_annotated %>% filter(cellType==ct))[,2:20112] %>% colMeans()) #remove columns like sample, celltype...
}

scCPM_cluster9 <- data.frame(row.names = colnames(scCPM %>% select(-Sample)))
for(ct in unique(scCPM_annotated$clusterType_9)){
  print(ct)
  scCPM_cluster9[[ct]] <- unname((scCPM_annotated %>% filter(clusterType_9==ct))[,2:20112] %>% colMeans())
}

scCPM_cluster6 <- data.frame(row.names = colnames(scCPM %>% select(-Sample)))
for(ct in unique(scCPM_annotated$clusterType_6)){
  print(ct)
  scCPM_cluster6[[ct]] <- unname((scCPM_annotated %>% filter(clusterType_6==ct))[,2:20112] %>% colMeans())
}

library(corrplot)
ggcorrplot(cor(scCPM_cellTypes), hc.order = TRUE,
           outline.col = "white", lab=TRUE, lab_size = 2) + 
  theme(axis.text.y = element_text(size = 7), 
        axis.text.x = element_text(size = 7)) + 
  scale_fill_gradient2(limit = c(0,1), low = "blue", high =  "#BD1313", mid = "white", midpoint = 0.5)+ 
  labs(fill = "Correlation")

# 11. Rereun MuSiC and Scaden deconvolution with the new clustering obtained from the results with pseudo bulk data ####
## Update the sc expression set:
clustering <- list()
clustering["T and NK Cells"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "(?:T|NK) Cells"))
clustering["Switched mem B Cells"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Class-switched"))
clustering["Macrophages"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Inflam|Foam|Resident"))
clustering["Dendritic Cells"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Dendritic"))
clustering["Endothelial Cells I"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Endothelial Cells I$"))
clustering["Endothelial Cells II"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Endothelial Cells II"))
clustering["Smooth Muscle Cells"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "ACTA2+"))
clustering["Mastocytes"] <- list(str_subset(levels(scTest@phenoData@data$cellType), "Mast"))

unlist(unname(clustering)) #we still have 18 sub cell types (sanity check)
clustering
clType <- as.character(scTest@phenoData@data$cellType)
for(i in 1:length(clustering)){
  clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
}
scTest@phenoData@data$specialCluster <- factor(clType, levels = names(clustering))

## Update the true proportions
propsSpecialClustering <- as_tibble(as.data.frame(t(as.data.frame(t(props18 %>% column_to_rownames("cellType")%>% select(-starts_with("numCells"), -starts_with("simNum")))) %>% 
                                      mutate(`T and NK Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
                                      mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
                                      mutate(`Macrophages`=select(., c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
                                      mutate(`Dendritic Cells`=`CD68+CD1C+ Dendritic Cells`) %>% select(-c("CD68+CD1C+ Dendritic Cells")) %>% 
                                      mutate(`Endothelial Cells I`=`CD34+ Endothelial Cells I`) %>% select(-c("CD34+ Endothelial Cells I")) %>% 
                                      mutate(`Endothelial Cells II`=`CD34+ Endothelial Cells II`) %>% select(-c("CD34+ Endothelial Cells II")) %>% 
                                      mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
                                      mutate(`Mastocytes`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`))) %>% rownames_to_column("cellType"))
propsSpecialClustering %>% select(-cellType) %>% colSums()

## Deconvolute with MuSiC:
musicResults <- musicDeconvolution(bulk18_allGenes, scTest, "specialCluster") %>% 
  write_csv(paste0("pseudoBulkData/AllSCGenes/musicProps_specialCluster.csv"))
musicMetrics_clust[["specialCluster"]] <- list(metrics=musicResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering) %>% checkMetrics(), 
                                                     props=musicResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering))
View(musicMetrics_clust$specialCluster$metrics)
ggplot(metrics, aes(x="music", y = resultsFor)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_fill_viridis_c(na.value = NA) + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Results for the deconvolution of the simulated mixtures"))

music$specialCluster$props %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ggtitle("MuSiC") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free")+theme(strip.text = element_text(size = 8))

## Deconvolute with Scaden:
asScadenTxtInput_scByPatient(scTest@assayData$exprs, scTest@phenoData@data, "pseudoBulkData/ScadenInput/", different3Clusters = F, cluster = "specialCluster")

scadenResults <- scadenReadDeconvolution("pseudoBulkData/AllSCGenes/scadenProps_specialCluster.txt")
resultsScaden_8 <- list(metrics=scadenResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering) %>% checkMetrics(), 
                                                      props=scadenResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering))
resultsScaden_8$props %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ggtitle("Scaden") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free")+theme(strip.text = element_text(size = 4))
## Deconvolute using NNLS:
nnlsResults <- NNLSDeconvolution(bulk18_allGenes, scTest, "specialCluster")
metrics <- nnlsResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering) %>% checkMetrics() 
props <- nnlsResults %>% joinSimulProps(propsNumSimulation = propsSpecialClustering)


props %>% ggplot(aes(x=simulProps, y=deconvProps)) + geom_point(aes(colour=cellType)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed") + ggtitle("NNLS") + ylab("Deconvoluted proportions") + xlab("Simulated proportions") +
  facet_wrap(vars(cellType), scales = "free")+theme(strip.text = element_text(size = 8))






