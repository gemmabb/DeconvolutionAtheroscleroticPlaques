###################################################################################
#                                                                                 #
# HELPER FUNCTIONS FOR THE INTEGRATION OF RESULTS                                 #                                       
#                                                                                 #
# Gemma Bel Bordes, Nov 20th 2022                                                 #
#                                                                                 #
###################################################################################


integrate_Music <- function(folder, clustering, differentInputs=F){
  if(differentInputs){ #files are written as "linear_scallGenes_cellType_PROPS.csv"
    i <- 1
    for(f in list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))){
      transf <- str_split(f, "_")[[1]][1]
      qc <- str_split(f, "_")[[1]][2]
      if(i==1){
        allRes <- read_csv(paste0(folder, f)) %>% mutate(Algorithm = "MuSiC", DataTrans = transf, ScQC = qc) %>% 
          rename(Patient = "...1") 
      }else{
        allRes <- rbind(allRes, read_csv(paste0(folder, f)) %>% mutate(Algorithm = "MuSiC", DataTrans = transf, ScQC = qc) %>% rename(Patient = "...1"))
      }
      i <- i+1
    }
    return(allRes)
  }
  else{
    f <- list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))
    return(read_csv(paste0(folder, f)) %>% mutate(Algorithm = "MuSiC") %>% rename(Patient = "...1"))
  }
}

integrate_NNLS <- function(folder, clustering, differentInputs=F){
  if(differentInputs){ #files are written as "linear_QC12_cellType_PROPS.csv"
    i <- 1
    for(f in list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))){
      transf <- str_split(f, "_")[[1]][1]
      qc <- str_split(f, "_")[[1]][2]
      if(i==1){
        allRes <- read_csv(paste0(folder, f)) %>% mutate(Algorithm = "NNLS", DataTrans = transf, ScQC = qc) %>% rename(Patient = "...1") 
      }else{
        allRes <- rbind(allRes, read_csv(paste0(folder, f)) %>% mutate(Algorithm = "NNLS", DataTrans = transf, ScQC = qc) %>% rename(Patient = "...1"))
      }
      i <- i+1
    }
    return(allRes)
  }
  else{
    f <- list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))
    return(read_csv(paste0(folder, f)) %>% mutate(Algorithm = "NNLS") %>% rename(Patient = "...1"))
  }
}

integrate_Bisque <- function(folder, clustering, differentInputs=F){
  if(differentInputs){
    i <- 1
    for(f in list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))){
      transf <- str_split(f, "_")[[1]][1]
      qc <- str_split(f, "_")[[1]][2]
      if(i==1){
        allRes <- as.data.frame(t(read_csv(paste0(folder, f)) %>% column_to_rownames("...1"))) %>% rownames_to_column("Patient") %>% 
          mutate(Algorithm = "Bisque", DataTrans = transf, ScQC = qc)
      }else{
        allRes <- rbind(allRes, as.data.frame(t(read_csv(paste0(folder, f)) %>% column_to_rownames("...1"))) %>% rownames_to_column("Patient") %>% 
                       mutate(Algorithm = "Bisque", DataTrans = transf, ScQC = qc))
      }
      i <- i+1
    }
    return(allRes)
  }
  else{
    f <- list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))
    return(as.data.frame(t(read_csv(paste0(folder, f)) %>% column_to_rownames("...1"))) %>% rownames_to_column("Patient") %>% 
             mutate(Algorithm = "Bisque"))
  }
}

integrate_Cibersort <- function(folder, clustering, differentInputs=F){
  #note that we have no data transformations here, data must be linear always
  if(differentInputs){
    i <- 1
    for(f in list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))){
      transf <- str_split(f, "_")[[1]][1]
      qc <- str_split(f, "_")[[1]][2]
      if(i==1){
        allRes <- read_csv(paste0(folder, f)) %>% rename(Patient=Mixture) %>% 
          select(-c("P-value", "Correlation", "RMSE")) %>% mutate(Algorithm = "CIBERSORTx", DataTrans = transf, ScQC = qc)
      }else{
        allRes <- rbind(allRes, read_csv(paste0(folder, f)) %>% rename(Patient=Mixture) %>% 
                       select(-c("P-value", "Correlation", "RMSE")) %>% mutate(Algorithm = "CIBERSORTx", DataTrans = transf, ScQC = qc))
      }
      i <- i+1
    }
    return(allRes)
  }else{
    f <- list.files(folder, pattern = paste0(clustering, "_PROPS.csv"))
    return(read_csv(paste0(folder, f)) %>% rename(Patient=Mixture) %>% 
             select(-c("P-value", "Correlation", "RMSE")) %>% mutate(Algorithm = "CIBERSORTx"))
  }
}

integrate_Scaden <- function(folder, clustering, differentInputs=F){
  #note that we have no data transformations here, data must be linear always
  if(differentInputs){
    i <- 1
    for(f in list.files(folder, pattern = paste0(clustering, "_PROPS.txt"))){
      transf <- str_split(f, "_")[[1]][1]
      qc <- str_split(f, "_")[[1]][2]
      if(i==1){
        allRes <- read_delim(paste0(folder, f)) %>% mutate(Algorithm = "Scaden", DataTrans = transf, ScQC = qc) %>% rename(Patient = "...1")
      }else{
        allRes <- rbind(allRes, read_delim(paste0(folder, f)) %>% mutate(Algorithm = "Scaden", DataTrans = transf, ScQC = qc) %>% rename(Patient = "...1"))
      }
      i <- i+1
    }
    return(allRes)
  }else{
    f <- list.files(folder, pattern = paste0(clustering, "_PROPS.txt"))
    return(read_delim(paste0(folder, f)) %>% mutate(Algorithm = "Scaden") %>% rename(Patient = "...1"))
  }
}

heatMapProportions <- function(datasetProps, Algor, Trans, QC){
  print(pheatmap::pheatmap(t(datasetProps %>% 
                               filter(Algorithm==Algor, DataTrans==Trans, ScQC==QC) %>% 
                               select(-c("Algorithm", "DataTrans", "ScQC")) %>% 
                               column_to_rownames("Patient")), 
                           show_colnames = F, main = paste0("Proportions using ", Algor, " ", Trans, " ", QC)))
}

boxPlotsProportions_differentQC <- function(datasetProps, trans){
  swr = function(string, nwrap=10) {
    paste(strwrap(string, width=nwrap), collapse="\n")
  }
  swr <- Vectorize(swr)
  df <- datasetProps  %>% filter(DataTrans==trans) %>% 
    pivot_longer(!c("Patient", "Algorithm", "DataTrans", "ScQC"), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=swr(CellType)) %>% rename(`Gene selection`=ScQC)
  print(df %>% 
          ggplot(aes(x=Algorithm,y=Proportion)) +
          geom_boxplot(aes(fill=`Gene selection`), outlier.size = 0.2, alpha=0.7) +
          facet_wrap(~CellType, ncol = ceiling(length(unique(df$CellType))/2)) +
          theme_bw() + theme(strip.text = element_text(size = 7.5), 
                axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) + scale_fill_manual(values=c("darkolivegreen1", "darkolivegreen4")) +
          ggtitle(paste0("Proportions using ", trans, " data")) + ylim(0,1))
}

boxPlotsProportions_differentTransf <- function(datasetProps){
  swr = function(string, nwrap=20) {
    paste(strwrap(string, width=nwrap), collapse="\n")
  }
  swr <- Vectorize(swr)
  df <- datasetProps %>% 
    pivot_longer(!c("Patient", "Algorithm", "DataTrans", "ScQC"), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=swr(CellType))
  print(df %>% 
          ggplot(aes(x=Algorithm,y=Proportion)) +
          geom_boxplot(aes(fill=DataTrans), outlier.size = 0.2, alpha=0.7) +
          facet_wrap(~CellType, ncol = ceiling(length(unique(df$CellType))/2)) +
          theme_bw() + scale_fill_manual(values=c("#fffa60", "#ffbe00", "#a77d00")) + 
          theme(strip.text = element_text(size = 6), 
                axis.text.x = element_text(angle = 90)) + ylim(0,1))
}

boxPlotsProportions_differentClustering <- function(datasetProps, trans){
  swr = function(string, nwrap=20) {
    paste(strwrap(string, width=nwrap), collapse="\n")
  }
  swr <- Vectorize(swr)
  df <- datasetProps  %>% filter(DataTrans==trans) %>% 
    pivot_longer(!c("Patient", "Algorithm", "DataTrans", "ClusteringDeconvolution"), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=swr(CellType))
  print(df %>% 
          ggplot(aes(x=Algorithm,y=Proportion)) +
          geom_boxplot(aes(fill=ClusteringDeconvolution), outlier.size = 0.2, alpha=0.7) +
          facet_wrap(~CellType, ncol = ceiling(length(unique(df$CellType))/2)) +
          theme_bw() + theme(strip.text = element_text(size = 6), 
                axis.text.x = element_text(angle = 90)) + scale_fill_manual(values=c("mistyrose1", "rosybrown3", "rosybrown4")) + 
          ggtitle(paste0("Proportions using ", trans, " data")) + ylim(0,1))
}

analyzeCorrelationProps_CellType <- function(dataset, cellType, arrayAlgorithms = NULL){
  library(ggcorrplot)
  library(pheatmap)
  library(ggstatsplot)
  library(cowplot)
  if(is.null(arrayAlgorithms)){
    subdata <- dataset %>% select(Patient, cellType, Algorithm) %>% pivot_wider(names_from = Algorithm, values_from = cellType)
  }
  else{
    subdata <- dataset %>% filter(Algorithm %in% arrayAlgorithms) %>% 
      select(Patient, cellType, Algorithm) %>% pivot_wider(names_from = Algorithm, values_from = cellType)
  }
  subdata_norm <- subdata %>% mutate_if(is.numeric, ~ (.-min(.))/(max(.)-min(.))) #range from 0 to 1 --> better for correlation detection
  #if max==min --> we will get infinite values!! 
  #try(print(pheatmap(t(subdata_norm %>% column_to_rownames("Patient")), main = paste0(cellType, " proportions (normalized 0-1)"), show_colnames = F)))
  pearsons <- ggcorrmat(data = subdata, output = "dataframe", type = "parametric") #parametric for Pearson
  #for the correlations it's exactly the same for normalized-range proportions or the original ones!
  print(ggcorrmat(data = subdata, pch=2, caption = "", output = "plot", type = "parametric", sig.level = 1, colors = c("#440154", "#21918c", "#fde725"), k=2, 
                      title = paste0(cellType), 
                      ggcorrplot.args = list(lab_size = 2.5)) + 
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 8), axis.text.y = ggplot2::element_text(size = 8), plot.margin = margin(10,0,10,10)
              ))
}

analyzeCategoricalFeatures <- function(dataset, listFeatures, cellType, titlePlot, plotOnlySignificant=F, pval=0.05, showpValue = T){
  library(rstatix)
  library(ggpubr)
  listPlots <- list()
  for(feature in listFeatures){
    print(feature)
    stat.test <- (dataset %>% filter(is.na(.data[[feature]])==FALSE)) %>% t_test(as.formula(paste0("`", cellType, "` ~ ", feature)))
    if(plotOnlySignificant){
      if(any(stat.test$p<pval)){
        print("--- something significant here")
        stat.test <- stat.test %>% add_xy_position(x = feature)
        listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
                geom_boxplot()+stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase = 0.02, size = 1.5)+
                ggtitle(paste0(titlePlot, paste0(" (n = ", sum(!is.na(dataset[[feature]])), ")")))+theme_bw()+
                theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)))
      }
    }
    else{
      if(showpValue){
        stat.test <- stat.test %>% add_xy_position(x = feature)
        listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
                                   geom_boxplot()+stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase = 0.02, size = 1.5)+
                                   ggtitle(paste0(titlePlot, paste0(" (n = ", sum(!is.na(dataset[[feature]])),")")))+theme_bw()+
                                   theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
                                   annotate("text", x=-Inf, y = Inf, vjust=1.5, hjust=22, label= paste0("n = ", sum(!is.na(dataset[[feature]])))))
        
      }else{
        listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
                                   geom_violin()+
                                   ggtitle(paste0(titlePlot, paste0(" (n = ", sum(!is.na(dataset[[feature]])),")")))+theme_bw()+
                                   theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
                                   annotate("text", x=-Inf, y = Inf, vjust=1.5, hjust=22, label= paste0("n = ", sum(!is.na(dataset[[feature]])))))
        
      }
    }
  }
  return(listPlots)
}

kendallResults <- function(dataset, ctProportions, histoFeature){
  print(paste0(ctProportions, histoFeature))
  corr <- round(cor.test(dataset[[ctProportions]], as.numeric(dataset[[histoFeature]]), method = "kendall")$estimate, 3)
  pval <- round(cor.test(dataset[[ctProportions]], as.numeric(dataset[[histoFeature]]), method = "kendall")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  return(paste0(corr, " (", psymbol, ")"))
}

pearsonResults <- function(dataset, ctProportions, histoFeature){
  print(paste0(ctProportions, histoFeature))
  corr <- round(cor.test(dataset[[ctProportions]], as.numeric(dataset[[histoFeature]]), method = "pearson")$estimate, 3)
  pval <- round(cor.test(dataset[[ctProportions]], as.numeric(dataset[[histoFeature]]), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  return(paste0(corr, " (", psymbol, ")"))
}

analyzeContinuousFeatures <- function(dataset, listFeatures, cellType, titlePlot){
  listPlots <- list()
  for(feature in listFeatures){
    print(feature)
    listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
            geom_point()+geom_smooth(method = "lm")+
            annotate("text", x=0.7*max(dataset[[feature]], na.rm = T), y=0.25, label= paste0("rho=", as.character(
              round(cor.test(dataset[[cellType]], as.numeric(dataset[[feature]]), method = "pearson")$estimate, 3))))+ggtitle(titlePlot))
  }
  return(listPlots)
}


getTransformedPhenoOriginalSC <- function(scFile){
  scEset <- readRDS(scFile)
  pheno <- scEset@phenoData@data 
  pheno$SubjectName <- map_chr(pheno$SubjectName, ~{paste("ae", ., sep = "")})
  pheno %<>% filter(!cellType %in% c("CD68+CD4+ Monocytes", "CD79+ Plasma B Cells")) %>% 
    select(-c("clusterType_9", "clusterType_6"))
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
  #### 8 cell types, special clustering
  clustering <- list()
  clustering["T and NK Cells"] <- list(str_subset(levels(pheno$cellType), "(?:T|NK) Cells"))
  clustering["Switched mem B Cells"] <- list(str_subset(levels(pheno$cellType), "Class-switched"))
  clustering["Macrophages"] <- list(str_subset(levels(pheno$cellType), "Inflam|Foam|Resident"))
  clustering["Dendritic Cells"] <- list(str_subset(levels(pheno$cellType), "Dendritic"))
  clustering["Endothelial Cells I"] <- list(str_subset(levels(pheno$cellType), "Endothelial Cells I$"))
  clustering["Endothelial Cells II"] <- list(str_subset(levels(pheno$cellType), "Endothelial Cells II"))
  clustering["Smooth Muscle Cells"] <- list(str_subset(levels(pheno$cellType), "ACTA2+"))
  clustering["Mastocytes"] <- list(str_subset(levels(pheno$cellType), "Mast"))
  clType <- as.character(pheno$cellType)
  for(i in 1:length(clustering)){
    clType[clType %in% unname(unlist(clustering[i]))] <- names(clustering[i])
  }
  pheno$specialCluster <- factor(clType, levels = names(clustering))
  return(pheno)
}

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
  print(data.frame(Patient = names(colSums(sumCellTypesPerPatient)), SequencedCells = unname(colSums(sumCellTypesPerPatient))) %>% 
    ggplot(aes_string(x="SequencedCells")) + geom_histogram(fill="gray65") + theme_bw())
  props <- as.data.frame(map_df(sumCellTypesPerPatient, ~{.x/sum(.)}))
  rownames(props) = rownames(sumCellTypesPerPatient)
  props <- props[order(rownames(props)), order(colnames(props))] #order names alphabetically
  return(as.data.frame(t(props)) %>% rownames_to_column("Patient") %>% mutate(Patient = paste0("ae", Patient)))
}

integrateSCandDeconvProps <- function(scProps, deconvProps_diffAlgorithms){
  deconvProps_diffAlgorithms_long <- deconvProps_diffAlgorithms %>% filter(Patient %in% intersect(scProps$Patient, deconvProps_diffAlgorithms$Patient)) %>% 
    pivot_longer(-c("Patient", "Algorithm"), names_to = "cellType", values_to = "deconvProps") 
  scProps_long <-scProps %>% filter(Patient %in% intersect(scProps$Patient, deconvProps_diffAlgorithms$Patient)) %>% 
    pivot_longer(-c("Patient"), names_to = "cellType", values_to = "scProps") 
  return(left_join(deconvProps_diffAlgorithms_long, scProps_long, by = c("Patient", "cellType")))
}

checkMetrics <- function(dataset){
  message("Computing metrics...")
  ct <- unique(dataset$cellType)
  metrics <- data.frame(resultsFor = c("all", ct))
  metrics$NRMSE <- c(with(dataset, {sqrt(mean((deconvProps-scProps)^2))/mean(scProps) %>% round(.,4)}), #for all the proportions
                    map_dbl(ct, function(ct){with(dataset %>% filter(cellType==ct),         #for each cell type
                                                  {sqrt(mean((deconvProps-scProps)^2))/mean(scProps) %>% round(.,4)})}))
  metrics$Correlation <- c(with(dataset, {cor(deconvProps, scProps) %>% round(.,4)}),           #for all the proportions
                           map_dbl(ct, function(ct){with(dataset %>% filter(cellType==ct),         #for each cell type
                                                         {cor(deconvProps, scProps) %>% round(.,4)})}))
  return(metrics)
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

boxPlot_histoVsProps <- function(dataset, score, cellType){
  return(ggplot(dataset %>% filter(!is.na(!!sym(score))), aes_string(x=score, y=paste0("`",cellType,"`")))+geom_boxplot()+
           theme_bw()+annotate("text", x=Inf, y = Inf, vjust=1, hjust=1, 
             label= paste0("n = ", sum(!is.na(dataset[[score]])))))
}

violinPlot_histoVsProps <- function(dataset, score, cellType){
  return(ggplot(dataset %>% filter(!is.na(!!sym(score))), aes_string(x=score, y=paste0("`",cellType,"`")))+geom_violin()+
           theme_bw()+annotate("text", x=Inf, y = Inf, vjust=1, hjust=1, 
                    label= paste0("n = ", sum(!is.na(dataset[[score]])))))
}

regrLine_markerExpression <- function(dataset, markerGene, cellType){
  library(ggpubr)
  p <- (dataset %>% ggplot(aes_string(markerGene, paste0("`",cellType,"`")))+geom_point()+geom_smooth(method = "lm", color="#DB6A13")+
          ylab(paste(cellType, " proportions"))+facet_wrap(~Algorithm, ncol = 5)+theme_bw()+theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)))
  resCor <- data.frame(Algorithm = sort(unique(all8_withExpr$Algorithm)))
  textCol <- c()
  for(i in 1:length(resCor$Algorithm)) {
    subdata <- dataset[dataset[, "Algorithm"] == resCor$Algorithm[i], ]
    corr <- round(cor.test(subdata[[cellType]], as.numeric(subdata[[markerGene]]), method = "pearson")$estimate, 3)
    pval <- round(cor.test(subdata[[cellType]], as.numeric(subdata[[markerGene]]), method = "pearson")$p.value, 4)
    textCol <- c(textCol, paste0("R=", corr, " , p=", pval))
  }
  resCor$text <- textCol
  p <- p + geom_label(data = resCor, 
                     aes(0.6 * max(dataset[[markerGene]], na.rm = TRUE),
                         0.99 * max(dataset[[cellType]], na.rm = TRUE),
                         label = text), size = 2.5)
  return(p)
}


