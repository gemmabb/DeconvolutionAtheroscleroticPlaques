###########################################################################
#                                                                         #
# HELPER FUNCTIONS, COMPARE DECONV RESULTS ANS SC-DERIVED PROPS           #
#                                                                         #
# Gemma Bel Bordes, Nov 20th 2022                                         #
#                                                                         #
###########################################################################

getTransformedPhenoOriginalSC <- function(scFile){
  scEset <- readRDS(scFile)
  pheno <- scEset@phenoData@data 
  pheno$SubjectName <- map_chr(pheno$SubjectName, ~{paste("ae", ., sep = "")})
  pheno %<>% filter(!cellType %in% c("CD3+ T Cells VI", "CD68+KIT+ Mast Cells", 
                                     "CD68+CD4+ Monocytes", "CD79+ Plasma B Cells")) %>% select(-c("clusterType_13", "clusterType_9"))
  pheno$cellType <- factor(pheno$cellType, levels = unique(pheno$cellType))
  #### 9 cell types ####
  clustering <- list()
  clustering["FOXP3+ T Cells"] <- list(str_subset(levels(pheno$cellType), "FOXP3+"))
  clustering["CD3+ T Cells I-IV"] <- list(str_subset(levels(pheno$cellType), "T Cells I"))
  clustering["CD3+ T Cells V"] <- list(str_subset(levels(pheno$cellType), "T Cells V$"))
  clustering["NK Cells"] <- list(str_subset(levels(pheno$cellType), "NK Cells"))
  clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno$cellType), "Class-switched"))
  clustering["CD68+ (Foam Cells and inflam/resident macrophages)"] <- list(str_subset(levels(pheno$cellType), "Inflam|Foam|Resident"))
  clustering["CD68+ Dendritic"] <- list(str_subset(levels(pheno$cellType), "Dendritic"))
  clustering["Endothelial Cells"] <- list(str_subset(levels(pheno$cellType), "CD34+"))
  clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno$cellType), "ACTA2+"))
  clType <- as.character(pheno$cellType)
  for(i in 1:length(clustering)){
    clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
  }
  pheno$clusterType_9 <- factor(clType, levels = names(clustering))
  #### 6 cell types 
  clustering <- list()
  clustering["T Cells"] <- list(str_subset(levels(pheno$cellType), "T Cells"))
  clustering["NK Cells"] <- list(str_subset(levels(pheno$cellType), "NK Cells"))
  clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno$cellType), "Class-switched"))
  clustering["CD68+ Cells"] <- list(str_subset(levels(pheno$cellType), "CD68+"))
  clustering["Endothelial Cells"] <- list(str_subset(levels(pheno$cellType), "CD34+"))
  clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno$cellType), "ACTA2+"))
  clType <- as.character(pheno$cellType)
  for(i in 1:length(clustering)){
    clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
  }
  pheno$clusterType_6 <- factor(clType, levels = names(clustering))
  return(pheno)
}


#' Get SC-derived proportions for the samples in the SC database
#'
#' @param scPheno Dataframe with the phenotypic data for each cell in the SC database
#' @param clustering String with the clustering that you want the proportions for (celltype, cluster6, cluster9)
#'
#' @return Dataframe with the proportions of each cell type/cluster
getScProportions <- function(scPheno, clustering){
  iter <- 1
  for(patient in unique(scPheno$SubjectName)){
    if(iter==1){
      sumCellTypesPerPatient <- data.frame(table(scPheno[[clustering]][scPheno$SubjectName==patient]), row.names = "Var1") #cell clusters as row names 
      colnames(sumCellTypesPerPatient) <- patient
    }
    else{
      df <- data.frame(table(scPheno[[clustering]][scPheno$SubjectName==patient]), row.names = "Var1") #cell clusters as row names 
      colnames(df) <- patient
      sumCellTypesPerPatient <- bind_cols(sumCellTypesPerPatient, df)
    }
    iter <- iter+1
  }
  print("Number of cells sequenced per patient:")
  print(colSums(sumCellTypesPerPatient))
  props <- as.data.frame(map_df(sumCellTypesPerPatient, ~{.x/sum(.)}))
  rownames(props) = rownames(sumCellTypesPerPatient)
  props <- props[order(rownames(props)), order(colnames(props))] #order names alphabetically
  return(props)
}

getMusicProportions <- function(musicResFile){
  props <- read_csv(musicResFile) 
  props <- tail(props, 32) %>% column_to_rownames("...1") #the last 32 are the new patients
  props <- as.data.frame(t(props))
  props <- props[order(rownames(props)), ] #order row names alphabetically
  return(props)
}

getBisqueProportions <- function(bisqueResFile){
  props <- read_csv(bisqueResFile) %>% column_to_rownames("...1") %>% 
    select(tail(names(.), 32)) #the last 32 are the new patients
  props <- props[order(rownames(props)), ] #order row names alphabetically
  return(props)
}

getScadenProportions <- function(scadenResFile){
  props <- read_delim(scadenResFile) 
  props <- tail(props, 32) %>% column_to_rownames("...1") #the last 32 are the new patients
  props <- as.data.frame(t(props))
  props <- props[order(rownames(props)), ] #order row names alphabetically
  return(props)
}

getCibersortProportions <- function(ciberResFile){
  props <- read_csv(ciberResFile) %>% select(-c("P-value", "Correlation", "RMSE")) %>% column_to_rownames("Mixture") %>% 
    tail(32) #the last 32 are the new patients
  props <- as.data.frame(t(props))
  return(props[order(rownames(props)), ])#order row names alphabetically
}

checkCorrelation <- function(props1, props2){
  if(dim(props1)[1] == dim(props2)[1]){
    for(celltype in rownames(props1)){
      message(paste0(celltype, " --> ", "Pearsons ",
                   cor(as.double((props1 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                            as.double((props2 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                            method = "pearson") %>% round(.,4), " / RMSE ", 
            ModelMetrics::rmse(as.double((props1 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                                                                                           as.double((props2 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]))%>% round(.,4)))
    }
  }
  else{print("Different amount of cell types between the datasets")}
}

compareProportions <- function(scProps, deconvProps, samplesToIgnore=NULL){
  samplesShared <- intersect(colnames(scProps), colnames(deconvProps))
  if(!is.null(samplesToIgnore)){
    scProps %<>% select(-samplesToIgnore)
    deconvProps %<>% select(-samplesToIgnore)
    samplesShared <- samplesShared[which(!samplesShared %in% samplesToIgnore)] #remove also from here
  }
  pheatmap::pheatmap(scProps %>% select(all_of(samplesShared)), main = "Proportions from SC data",  
                     cluster_cols = F, cluster_rows = F)
  pheatmap::pheatmap(deconvProps %>% select(all_of(samplesShared)), main="Proportions from deconvolution results", 
                     cluster_cols = F, cluster_rows = F)
  checkCorrelation(scProps, deconvProps) #Pearson's correlation coefficients for each cell type
  
  longDeconvProps <- deconvProps %>% select(all_of(samplesShared)) %>% rownames_to_column("CellType") %>% gather("Mixture", "DeconvPrediction", -1) #long dataframe to compute accuracy between props
  longSCProps <- scProps %>% select(all_of(samplesShared)) %>% rownames_to_column("CellType") %>% gather("Mixture", "scExpected", -1) #long dataframe to compute accuracy between props
  longAllProps <- longDeconvProps %>% right_join(longSCProps, by=c("CellType","Mixture"))
  longAllProps$DeconvPrediction <- round(longAllProps$DeconvPrediction, 3)
  longAllProps$scExpected <- round(longAllProps$scExpected, 3)
  
  print(longAllProps %>% summarise(RMSE = ModelMetrics::rmse(scExpected,DeconvPrediction) %>% round(.,4), 
                  Pearson=cor(scExpected,DeconvPrediction) %>% round(.,4)))
}

pearsonsAndRMSE_CellTypeAndDataInput <- function(phenoSC, clustering, folderResults, filesDeconvolution, annotDataUsed, deconvMethod){
  props1 <- getScProportions(phenoSC, clustering)
  CellTypeColumn <- c()
  DataInputColumn <- c()
  PearsonsColumn <- c()
  RmseColumn <- c()
  for(i in 1:length(filesDeconvolution)){
    if(deconvMethod=="music"){
      props2 <- getMusicProportions(paste0(folderResults, filesDeconvolution[i]))
    }else{
      if(deconvMethod=="cibersort"){
        props2 <- getCibersortProportions(paste0(folderResults, filesDeconvolution[i]))
      }else{
        if(deconvMethod=="bisque"){
          props2 <- getBisqueProportions(paste0(folderResults, filesDeconvolution[i]))
        }
        else{
          if(deconvMethod=="scaden"){
            props2 <- getScadenProportions(paste0(folderResults, filesDeconvolution[i]))
          }
        }
      }
    }
    DataInputColumn <- c(DataInputColumn, rep(annotDataUsed[i], dim(props1)[1]))
    for(celltype in rownames(props1)){
      CellTypeColumn <- c(CellTypeColumn, celltype)
      PearsonsColumn <- c(PearsonsColumn, cor(as.double((props1 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                                              as.double((props2 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                                              method = "pearson"))
      RmseColumn <- c(RmseColumn, ModelMetrics::rmse(as.double((props1 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ]), 
                                                     as.double((props2 %>% select(intersect(colnames(props1), colnames(props2))))[celltype, ])))
    }
  }
  results <- data.frame(CellType <- CellTypeColumn, DataInput <- DataInputColumn, Pearsons <- PearsonsColumn, RMSE <- RmseColumn)
  
  ggplot(results, aes(x = DataInput, y = CellType)) + geom_point(aes(size = 1/RMSE, fill = Pearsons), shape = 21) +
    scale_fill_viridis_c(na.value = NA) + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) +
    ggtitle(paste0(deconvMethod, " vs SC-dervied proportions")) 
}

## for when I have a df with proportions in the following way: columns are Patient, All Cell Types and Algorithm
pearsonsAndRMSE_allProps <- function(phenoSC, clustering, dataProps, dropPatient){
  props1 <- getScProportions(phenoSC, clustering)
  CellTypeColumn <- c()
  DataInputColumn <- c()
  PearsonsColumn <- c()
  RmseColumn <- c()
  for(algor in unique(dataProps$Algorithm)){
    DataInputColumn <- c(DataInputColumn, rep(algor, length(colnames(dataProps %>% select(-c("Patient", "Algorithm")))))) #rep as many times as cell types we have
    for(celltype in colnames(dataProps %>% select(-c("Patient", "Algorithm")))){
      CellTypeColumn <- c(CellTypeColumn, celltype)
      intersectPatients <- intersect(colnames(props1), dataProps$Patient)
      valuesSC <- as.double((props1 %>% select(intersectPatients[!intersectPatients %in% dropPatient]))[celltype, ])
      valuesDeconv <- pull((dataProps %>% filter(Patient %in% intersectPatients[!intersectPatients %in% dropPatient]) %>% 
                                   filter(Algorithm==algor) %>% arrange(Patient)), celltype)
      PearsonsColumn <- c(PearsonsColumn, cor(valuesSC, valuesDeconv, method = "pearson"))
      RmseColumn <- c(RmseColumn, ModelMetrics::rmse(valuesSC, valuesDeconv))
    }
  }
  results <- data.frame(CellType = CellTypeColumn, DataInput = DataInputColumn, Pearsons = PearsonsColumn, RMSE = RmseColumn)
  print(results[1:5,])
  print(ggplot(results, aes(x = DataInput, y = CellType)) + geom_point(aes(size = 1/RMSE, fill = Pearsons), shape = 21) +
    scale_fill_viridis_c(na.value = NA) + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Deconvolution vs SC-dervied proportions (", clustering, ")")))

  return(results)
}



