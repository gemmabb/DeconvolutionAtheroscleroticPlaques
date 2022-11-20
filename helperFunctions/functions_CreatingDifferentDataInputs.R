###########################################################################
#                                                                         #
# HELPER FUNCTIONS, DATA PROCESSING AND CREATING DECONV. INPUTS           #
#                                                                         #
# Gemma Bel Bordes, Nov 20th 2022                                         #
#                                                                         #
###########################################################################

### 1. DATA READING AND CORRECTIONS ####
correctForHGNCUniqueSymbols <- function(dataset){
  #delete entries with unknown symbol:
  numGenesBefore <- dim(dataset)[1]
  dataset %<>% filter(!is.na(symbol))
  message(paste0("Removing unkown genes: ", numGenesBefore-dim(dataset)[1], " genes were removed"))
  #delete duplicated genes with lower counts:
  numGenesBefore <- dim(dataset)[1]
  duplicates <- dataset$symbol[duplicated(dataset$symbol)]
  dataset_sumCounts <- dataset %>% mutate(Total = rowSums(dplyr::select(., starts_with("ae"))))
  for(s in unique(duplicates)){
    toDelete <- (dataset_sumCounts %>% filter(symbol==s) %>% arrange(desc(Total)))$ensembl_gene_id[-1] #select all but first gene (ordered from more to less counts)
    dataset_sumCounts <- dataset_sumCounts %>% filter(!(ensembl_gene_id %in% toDelete)) #delete rows without max total count 
  }
  message(paste0("Removing lower counts from duplicated genes: ", numGenesBefore-dim(dataset_sumCounts)[1], " genes were removed"))
  return(dataset_sumCounts %>% dplyr::select(-Total))
}

correctForHGNCUniqueSymbols_onlyHGNC <- function(dataset){ #for cases in which we don't have ensembl ids and HGNC symbols
  dataset %<>% mutate(geneNumber = 1:dim(dataset)[1])
  #delete duplicated genes with lower counts:
  numGenesBefore <- dim(dataset)[1]
  duplicates <- dataset$symbol[duplicated(dataset$symbol)]
  dataset_sumCounts <- dataset %>% mutate(Total = rowSums(dplyr::select(., -c("symbol", "geneNumber"))))
  for(s in unique(duplicates)){
    toDelete <- (dataset_sumCounts %>% filter(symbol==s) %>% arrange(desc(Total)))$geneNumber[-1] #select all but first gene (ordered from more to less counts)
    dataset_sumCounts <- dataset_sumCounts %>% filter(!(geneNumber %in% toDelete)) #delete rows without max total count 
  }
  message(paste0("Removing lower counts from duplicated genes: ", numGenesBefore-dim(dataset_sumCounts)[1], " genes were removed"))
  return(dataset_sumCounts %>% dplyr::select(-c("Total", "geneNumber")))
}

correctReadCounts_Max_UMIsaturation <- function(dataset, maxRC){
  library(tidyverse)
  library(magrittr)
  dataset %<>% mutate_at(vars(starts_with("ae")), ~map_dbl(., ~ifelse(.x>maxRC, maxRC, .x))) #max read count correction
  #correct for UMI saturation (there are maxRC+1 possible UMIs)
  message("Correcting for UMI saturation")
  dataset %<>% mutate_at(vars(starts_with("ae")), ~{round(-(maxRC+1)*log(1-(./(maxRC+1))))})
  message("After correction, max read count is finally ", max(dataset %>% select(starts_with("ae"))))
  return(dataset)
}

QCSingleCell_MitRibLibsize <- function(scMatrix, numMAD = 3, mitochondrialUniqueTreatment = F, mitochondrialMax = 0.05){
  library(ggplot2)
  #ribGenes <- suppressMessages(read_delim("/Volumes/DATASHUR_P2/preProcessing/input/ribosomalSymbols.txt")$`Approved symbol`) #all ribosomal genes, retrieved from https://www.genenames.org/tools/search/ filtering by ribosomal non coding (28/02/2022)
  #mitGenes <- suppressMessages(read_delim("/Volumes/DATASHUR_P2/preProcessing/input/mitochondrialSymbols.txt")$`Approved symbol`) #all mitochondrial genes, retrieved from https://www.genenames.org/tools/search/ filtering by mitochondrial encoded rnas(28/02/2022)
  #note that there are overlapping genes present in both mitGenes and ribGenes
  sc_df <- as.data.frame(scMatrix)
  libSize <- sc_df %>% colSums()
  #ribContentRatio <- (sc_df %>% rownames_to_column("symbol") %>% filter(symbol %in% ribGenes) %>% select(-1) %>% colSums())/libSize
  ribContentRatio <- (sc_df %>% rownames_to_column("symbol") %>% filter(grepl("^RPL|^RPS|_RPL|_RPS", symbol, ignore.case = T)) %>% select(-1) %>% colSums())/libSize
  mitContentRatio <- (sc_df %>% rownames_to_column("symbol") %>% filter(grepl("^MT-|_MT-", symbol, ignore.case = T)) %>% select(-1) %>% colSums())/libSize
  #mitContentRatio <- (sc_df %>% rownames_to_column("symbol") %>% filter(symbol %in% mitGenes) %>% select(-1) %>% colSums())/libSize
  numCellsBefore <- dim(sc_df)[2]
  listFeatures <- list(librarySize = libSize, ribosomalContentRatio = ribContentRatio, mitochondrialContentRatio = mitContentRatio)
  i <- 1
  for(feature in listFeatures){
    print(ggplot(data.frame(f = unname(feature)), aes(f)) + geom_histogram(bins = 30) + geom_vline(aes(xintercept=median(f), color="median")) + 
            geom_vline(aes(xintercept=median(f)+numMAD*mad(f), color="median ± 3 MAD")) + geom_vline(aes(xintercept=median(f)-numMAD*mad(f), color="median ± 3 MAD"))+
            scale_color_manual(name="statistics", values = c("median"="blue", "median ± 3 MAD"="red"))+
            ggtitle(paste0("Histogram of ", names(listFeatures)[i])) + xlab(names(listFeatures)[i]) + ylab("counts (cells)"))
    toRemove <- names(which(feature > median(feature) + numMAD*mad(feature) | 
                              feature < median(feature) - numMAD*mad(feature)))
    if(names(listFeatures)[i]=="mitochondrialContentRatio"){
      if(mitochondrialUniqueTreatment){
        print(ggplot(data.frame(f = unname(feature)), aes(f)) + geom_histogram(bins = 30) + geom_vline(aes(xintercept=mitochondrialMax, color="maxMitoContent")) + 
                scale_color_manual(name="statistics", values = c("maxMitoContent"="blue"))+
                ggtitle(paste0("Histogram of ", names(listFeatures)[i])) + xlab(names(listFeatures)[i]) + ylab("counts (cells)"))
        
        toRemove <- names(which(feature > mitochondrialMax))
      }
    }
    sc_df %<>% select(-(toRemove[toRemove %in% colnames(sc_df)])) #error if we try to remove a column that has been already deleted
    i <- i+1
  }
  if(mitochondrialUniqueTreatment){
    message(paste0("Removing cells with more than ", mitochondrialMax*100, "% counts coming from mitochondrial genes"))
    message(paste0("Removing cells with ribosomal content and library size further than ", 
                                                 numMAD, " mean absolute deviations... ", numCellsBefore-dim(sc_df)[2], " cells were deleted"))
  }else{
    message(paste0("Removing cells with ribosomal content, mitochondrial content or library size further than ", 
                   numMAD, " mean absolute deviations... ", numCellsBefore-dim(sc_df)[2], " cells were deleted"))
    }
  return(sc_df)
}

QCSingleCell_detectableGenes <- function(scMatrix, minPercent = 5){
  library(ggplot2)
  print(ggplot(data.frame(x=rowSums(scMatrix > 0)), aes(x)) + geom_histogram(bins = 30) + 
          geom_vline(aes(xintercept=round(minPercent*dim(scMatrix)[2]/100), color="limit")) +
          scale_color_manual(name="statistics", values = c("limit"="blue"))+
          ggtitle("Histogram of number of cells in which each gene is detected") + xlab("number of cells") + ylab("count (genes)"))
  numGenesBefore <- dim(scMatrix)[1]
  toKeep <- which(rowSums(scMatrix > 0) >= round(minPercent*dim(scMatrix)[2]/100))
  scMatrix <- scMatrix[toKeep,] #only those rows (genes) that are present in a % of cells
  message(paste0("Removing genes detected in less than ", minPercent, " % of the cells... ", 
                 numGenesBefore-dim(scMatrix)[1], " genes were removed"))
  return(scMatrix) 
}

afterQCSingleCell_rmvCellTypes <- function(scMatrix, phenoDataset, minNumCells = 50, reads=F, pheno=F){
  cellTypesToRemove <- names(which(table(phenoDataset[colnames(scMatrix),"cellType"]) < minNumCells))
  message(paste0("Removing cell types with less than ", minNumCells, " cells..."))
  message(cellTypesToRemove)
  cellsToKeep <- rownames(phenoDataset %>% filter(!cellType %in% cellTypesToRemove)) #get names of the cells we have to keep
  phenoDataset %<>% filter(!cellType %in% cellTypesToRemove) %>% mutate(cellType = factor(cellType)) #last step to change number of levels (from 20 to __)
  scMatrix <- scMatrix[,cellsToKeep]
  if(reads==T){
    return(scMatrix)
  }
  if(pheno==T){
    return(phenoDataset)
  }
}


### 2. DATA TRANSFORMATIONS ####

transformTo <- function(matrixCounts, typeTrans){
  if(typeTrans=="linear"){return(matrixCounts)}
  if(typeTrans=="log"){return(log1p(matrixCounts))}
  if(typeTrans=="sqrt"){return(sqrt(matrixCounts))}
  if(typeTrans=="vst"){
    suppressWarnings(suppressMessages(library(DESeq2)))
    return(varianceStabilizingTransformation(data.matrix(matrixCounts)))
    }
}

normalizeTo <- function(matrixCounts, normType, scPheno=NULL){
  if(normType=="column"){return(apply(matrixCounts, 2, function(col){col/sum(col)}))}
  if(normType=="globalMinMax"){return((matrixCounts-min(matrixCounts))/(max(matrixCounts)-min(matrixCounts)))}
  if(normType=="globalZscore"){return((matrixCounts-mean(data.matrix(matrixCounts)))/(sd(data.matrix(matrixCounts))))}
  if(normType=="mean"){return(apply(matrixCounts, 2, function(col){col-mean(col)}))}
  if(normType=="log"){
    suppressWarnings(suppressMessages(library(Seurat)))
    return(LogNormalize(matrixCounts, verbose = F))
    }
  if(normType=="TMM"){ #Counts per million from TMM-normalized lib size
    message("Reminder! Pheno data needed for single cell TMM normalization!")
    library(edgeR) #v4.1.1
    if(is.null(scPheno)){classes <- colnames(matrixCounts)} #bulk case
    else{classes <- scPheno$cellType} #sc case
    return(matrixCounts %>% DGEList(group = classes) %>% calcNormFactors(method="TMM") %>% cpm())
    }
  if(normType=="medianRatios"){
    message("Reminder! Pheno data needed for single cell TMM normalization!")
    library(DESeq2)
    if(is.null(scPheno)){classes <- colnames(matrixCounts)} #bulk case
    else{classes <- scPheno$cellType} #sc case
    metadata <- data.frame(classes=classes)
    return(matrixCounts %>% DESeq2::estimateSizeFactors %>% DESeq2::counts(normalized=T))
    }
  if(normType=="TPM"){
    require(SingleR)
    return(...)}
}

asExpressionSet_bulk <- function(matrixCounts, folder, outputFileName){
  saveRDS(ExpressionSet(assayData = as.matrix(matrixCounts), 
                        phenoData = AnnotatedDataFrame(data = data.frame(row.names = colnames(matrixCounts), 
                                                                         SampleID <- c(1:dim(matrixCounts)[2]),
                                                                         SampleNames <- colnames(matrixCounts)), 
                                                       varMetadata = data.frame(labelDescription = c("SampleID", "SampleNames"), 
                                                                                row.names = c("SampleID", "SampleNames")))),
                        file=paste0(folder, outputFileName))
  message(paste0("Bulk expression set saved as ", outputFileName))
}

asCibersortTxtInput_bulk <- function(matrixCounts, folder, outputFileName){
  as.data.frame(matrixCounts) %>% rownames_to_column("GeneSymbol") %>% 
    write.table(file = paste0(folder, outputFileName), sep = "\t", row.names = T)
  message(paste0("Bulk expression set saved as ", outputFileName))
}

asScadenTxtInput_bulk <- function(matrixCounts, folder, outputFileName){
  as.data.frame(matrixCounts) %>%
    write.table(file = paste0(folder, outputFileName), sep = "\t", quote = F, col.names = NA)
  message(paste0("Bulk expression set saved as ", outputFileName))
} 

asExpressionSet_sc <- function(matrixCounts, pheno, folder, outputFileName){
  saveRDS(ExpressionSet(assayData = as.matrix(matrixCounts), 
                        phenoData = AnnotatedDataFrame(data = pheno, varMetadata = data.frame(labelDescription = colnames(pheno), 
                                                                                              row.names = colnames(pheno)))),
                                                       file=paste0(folder, outputFileName))
  message(paste0("Single cell expression set saved as ", outputFileName))
}

asCibersortTxtInput_sc <- function(matrixCounts, pheno, folder, outputFileName, different3Clusters=T, cluster = NULL){
  if(different3Clusters){
    for(clustering in tail(colnames(pheno) ,3)){
      Counts <- as.data.frame(matrixCounts)
      names(Counts) <- unlist(unname(pheno[clustering])) #cell type instead of cell id
      Counts <- cbind(GeneSymbol = rownames(Counts), Counts)
      rownames(Counts) <- NULL
      write.table(Counts, file = paste0(folder, paste0(clustering, "_", outputFileName)), sep = "\t", row.names = F)
      message(paste0("Single cell expression sets saved as ", paste0(clustering, "_", outputFileName)))
    }
  }else{ #otherwise, get clustering column name
    Counts <- as.data.frame(matrixCounts)
    names(Counts) <- unlist(unname(pheno[cluster])) #cell type instead of cell id
    Counts <- cbind(GeneSymbol = rownames(Counts), Counts)
    rownames(Counts) <- NULL
    write.table(Counts, file = paste0(folder, paste0(outputFileName)), sep = "\t", row.names = F)
    message(paste0("Single cell expression sets saved as ", paste0(outputFileName)))
  }
}

asScadenTxtInput_sc <- function(matrixCounts, pheno, folder, outputFileName, different3Clusters=T, cluster = NULL){
  medianLibSize <- median(colSums(matrixCounts)) 
  message(paste0("Normalizing sc counts to median of lib size (", medianLibSize, ")"))
  sc_counts <- t(apply(matrixCounts, 2, function(x){x*medianLibSize/sum(x)})) #normalize by lib size 
  rownames(sc_counts) <- c(0:(dim(sc_counts)[1]-1))
  write.table(sc_counts, file = paste0(folder, outputFileName), sep = "\t", quote = F, col.names = NA)
  message(paste0("Single cell expression sets saved as ", paste0(folder, outputFileName)))
  ### cell types:
  if(different3Clusters){
    for(clustering in tail(colnames(pheno), 3)){
      write.table(data.frame(Celltype=unlist(unname(pheno[clustering]))), 
                  file = paste0(folder, paste0(clustering, "_", outputFileName)), sep = "\t", 
                  row.names = F, quote = F)
      message(paste0("Single cell annotations saved as ", paste0(folder, paste0(clustering, "_", outputFileName))))
    }
  }else{ #otherwise, get clustering column name
    write.table(data.frame(Celltype=unlist(unname(pheno[cluster]))), 
                file = paste0(folder, paste0("celltypes_", outputFileName)), sep = "\t", 
                row.names = F, quote = F)
    message(paste0("Single cell annotations saved as ", paste0("celltypes_", outputFileName)))
  }
  
} 

asScadenTxtInput_scByPatient <- function(matrixCounts, pheno, folder, different3Clusters=T, cluster = NULL){
  #inside the folder we should have subfolders with the different clusterings
  medianLibSize <- median(colSums(matrixCounts)) 
  message(paste0("Normalizing sc counts to median of lib size (", medianLibSize, ")"))
  sc_countsAll <- t(apply(matrixCounts, 2, function(x){x*medianLibSize/sum(x)})) #normalize by lib size
  phenoAll <- pheno
  
  if(different3Clusters){
    for(clustering in tail(colnames(phenoAll), 3)){
      print(clustering)
      message(paste0("saving in ", folder, clustering, "/"))
      for(patient in unique(phenoAll$SubjectName)){
        print(patient)
        sc_counts <- sc_countsAll[phenoAll %>% filter(SubjectName==patient) %>% pull(SampleID),,drop=F] #Select only those samples corresponding to the subject "Patient"
        pheno <- phenoAll %>% filter(SubjectName==patient)
        rownames(sc_counts) <- c(0:(dim(sc_counts)[1]-1))
        write.table(sc_counts, file = paste0(folder, clustering, "/", paste0("ae", patient, "_counts.txt")), sep = "\t", quote = F, col.names = NA)
        write.table(data.frame(Celltype=unlist(unname(pheno[clustering]))), 
                    file = paste0(folder, clustering, "/", paste0("ae", patient, "_celltypes.txt")), sep = "\t", 
                    row.names = F, quote = F)
      }
    }
  }else{ #otherwise, get clustering column name (cluster argument)
    print(cluster)
    message(paste0("saving in ", folder, cluster, "/"))
    for(patient in unique(phenoAll$SubjectName)){
      print(patient)
      sc_counts <- sc_countsAll[phenoAll %>% filter(SubjectName==patient) %>% pull(SampleID),,drop=F] #Select only those samples corresponding to the subject "Patient"
      pheno <- phenoAll %>% filter(SubjectName==patient)
      rownames(sc_counts) <- c(0:(dim(sc_counts)[1]-1))
      write.table(sc_counts, file = paste0(folder, cluster, "/", paste0("ae", patient, "_counts.txt")), sep = "\t", quote = F, col.names = NA)
      write.table(data.frame(Celltype=unlist(unname(pheno[cluster]))), 
                  file = paste0(folder, cluster, "/", paste0("ae", patient, "_celltypes.txt")), sep = "\t", 
                  row.names = F, quote = F)
    }
  }
}






