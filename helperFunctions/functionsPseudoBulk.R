###########################################################################
#                                                                         #
# HELPER FUNCTION FOR THE DECONVOLUTION ON ARTIFICIAL BULK DATA           #
#                                                                         #
# Gemma Bel Bordes, Nov 20th 2022                                         #
#                                                                         #
# Adapted from:                                                           #
# Avila Cobos F, Alquicira-Hernandez J, Powell JE, Mestdagh P and         #
# De Preter K. Benchmarking of cell type deconvolution pipelines for      #  
# transcriptomics data. (Nature Communications;                           #
# https://doi.org/10.1038/s41467-020-19015-1)                             #  
###########################################################################

scIntoTrainAndTest <- function(scEsetWhole, folder){
  cellsForTraining <- c() #cell names will be saved here
  tableCellTypeSubject <- as.data.frame.matrix(table(scEsetWhole$cellType, scEsetWhole$SubjectName))
  for(subject in colnames(tableCellTypeSubject)){
    for(ctype in rownames(tableCellTypeSubject)){
      cellsForSubjCelltype <- pull(scEsetWhole@phenoData@data %>% filter(SubjectName==subject, cellType==ctype), SampleID)
      if(length(cellsForSubjCelltype)<=1){ #if there's only 1 cell, randomly select if it's for train or test. If 0 cells, the cellsForTraining won't change
        if(sample(c(TRUE, FALSE), 1)){ #if random true, include cell into the test set
          cellsForTraining <- c(cellsForTraining, cellsForSubjCelltype)
        } #if false, do nothing
      }else{
        cellsForTraining <- c(cellsForTraining, 
                              sample(cellsForSubjCelltype, length(cellsForSubjCelltype)/2))
      }
    }
  }
  cellsForTesting <- pull(scEsetWhole@phenoData@data %>% filter(!SampleID %in% cellsForTraining), SampleID) #cell names
  message(paste0("Cells used for trainig: ", length(cellsForTraining), " / Patients in the subdataset: ", length(unique(sc@phenoData@data[cellsForTraining,]$SubjectName))))
  message(paste0("Cells used for testing: ", length(cellsForTesting), " / Patients in the subdataset: ", length(unique(sc@phenoData@data[cellsForTraining,]$SubjectName)))) 
  #possibility for the cells to not be exactly 50:50% because when getting numCells/2, we always get the number below --> less for training than for testing
  #create the 2 different esets:
  source("helperFunctions/functions_CreatingDifferentDataInputs.R")
  matrixCounts_training <- dim(exprs(scEsetWhole)[,cellsForTesting])
  asExpressionSet_sc(matrixCounts=exprs(scEsetWhole)[,cellsForTraining], pheno=scEsetWhole@phenoData@data[cellsForTraining,], 
                     folder=folder, outputFileName = "trainSC.eset")
  asExpressionSet_sc(matrixCounts=exprs(scEsetWhole)[,cellsForTesting], pheno=scEsetWhole@phenoData@data[cellsForTesting,], 
                     folder=folder, outputFileName = "testSC.eset")
}

simulateBulk <- function(numMixtures, numCells, scEset, cellCluster, selectedGenes=NULL){
  library(Biobase)
  pseudoBulk <- tibble(Gene = rownames(exprs(sc)))
  cellType_dataset = as_tibble(as.data.frame(table(phenoData(sc)[[cellCluster]]))) %>% rename(cellType=Var1, numCellsSC=Freq)
  cellTypes <- cellType_dataset$cellType
  for(i in 1:numMixtures){
    feasibleMixture = FALSE
    while(!feasibleMixture){
      cellType_dataset %<>% mutate(simProps=runif(length(cellTypes), 1, 99)) %>% mutate(simProps=simProps/sum(simProps)) 
      cellType_dataset %<>% mutate(simNumCells=round(simProps*numCells))
      if(all(cellType_dataset$numCellsSC >= cellType_dataset$simNumCells)){feasibleMixture=TRUE}
    }
    simNumCells <- cellType_dataset$simNumCells #specific for each iteration
    cellType_dataset %<>% rename(!!paste0("simProps", i):=simProps, !!paste0("simNumCells", i):=simNumCells) #save simulated proportions and number of cells
    #let's now create the simulated mixtures by adding up the counts of the specific number of cells selected above:
    cellsForSimulation <- str_c(map2_chr(cellTypes, simNumCells, function(cellTypes, simNumCells){
      str_c(sample(sc@phenoData@data %>% filter(!!as.symbol(cellCluster)==cellTypes) %>% pull(SampleID), simNumCells), collapse = ",")
    }), collapse = ",") %>% str_split(",+") %>% unlist() #combine all sample names in a string and then split 
    pseudoBulk %<>% mutate(!!paste0("mixture",i):=rowSums(as_tibble(exprs(sc)) %>% select(cellsForSimulation[which(!cellsForSimulation %in% "")]))) #imp to remove empty spaces from cellsFromSimulation
  }
  if(!is.null(selectedGenes)){
    pseudoBulk_genesSC %<>% filter(Gene %in% selectedGenes)
  }
  return(list(pseudoBulk=pseudoBulk, cellProps=cellType_dataset))
}

saveEsetAndProps <- function(listPseudoProps, folderOutput){
  source("helperFunctions/functions_CreatingDifferentDataInputs.R")
  asExpressionSet_bulk(as.matrix(listPseudoProps$pseudoBulk %>% column_to_rownames("Gene")), 
                       folderOutput, paste0("pseudoBulk_", nrow(listPseudoProps$cellProps), ".eset"))
  listPseudoProps$cellProps %>% write_delim(paste0(folderOutput, "pseudoProps_", nrow(listPseudoProps$cellProps), ".txt"))
  message(paste0("Simulated proportions saved as pseudoProps_",nrow(listPseudoProps$cellProps),".txt"))
}

musicDeconvolution <- function(bulkEset, scEset, clustering){
  library(MuSiC, quietly = TRUE)
  message("Deconvoluting...")
  estProp <- music_prop(bulk.eset = bulkEset, sc.eset = scEset, clusters = clustering, 
                        samples = "SubjectName", verbose = F)
  return(as.data.frame(estProp$Est.prop.weighted) %>% rownames_to_column("mixture") %>% 
           pivot_longer(!mixture,names_to = "cellType", values_to = "deconvProps") %>% arrange(cellType))
}

NNLSDeconvolution <- function(bulkEset, scEset, clustering){
  library(MuSiC, quietly = TRUE)
  message("Deconvoluting...")
  estProp <- music_prop(bulk.eset = bulkEset, sc.eset = scEset, clusters = clustering, 
                        samples = "SubjectName", verbose = F)
  return(as.data.frame(estProp$Est.prop.allgene) %>% rownames_to_column("mixture") %>% 
           pivot_longer(!mixture,names_to = "cellType", values_to = "deconvProps") %>% arrange(cellType))
}

bisqueDeconvolution <- function(bulkEset, scEset, clustering){
  library(BisqueRNA, quietly = TRUE)
  message("Deconvoluting...")
  estProp <- ReferenceBasedDecomposition(bulk.eset = bulkEset, sc.eset = scEset, markers = NULL, 
                                                           use.overlap = F, cell.types = clustering)
  return(as.data.frame(estProp$bulk.props) %>% rownames_to_column("cellType") %>% 
           pivot_longer(!cellType, names_to = "mixture", values_to = "deconvProps") %>% arrange(cellType))
}

cibersortReadDeconvolution <- function(fileProps){
  return(read_delim(fileProps) %>% dplyr:::rename(mixture=Mixture) %>% select(-`P-value`, -Correlation, -RMSE) %>% 
    pivot_longer(!mixture, names_to = "cellType", values_to = "deconvProps") %>% arrange(mixture, cellType))
}

scadenReadDeconvolution <- function(fileProps){
  return(read_delim(fileProps) %>% dplyr:::rename(mixture=...1) %>% 
           pivot_longer(!mixture, names_to = "cellType", values_to = "deconvProps") %>% 
           arrange(mixture, cellType))
}

joinSimulProps <- function(deconvProps, propsNumSimulation){
  return(deconvProps %>% right_join(propsNumSimulation %>% rename_at(vars(starts_with("simProps")), ~str_replace_all(., "simProps", "mixture")) %>% 
                                      select(cellType, starts_with("mixture")) %>% pivot_longer(!cellType, names_to = "mixture", values_to = "simulProps") 
                                    %>% arrange(cellType), by=c("cellType", "mixture")))
}

checkMetrics <- function(deconvAndSimProportions){
  message("Computing metrics...")
  ct <- unique(deconvAndSimProportions$cellType)
  metrics <- data.frame(resultsFor = c("all", ct))
  metrics$NRMSE <- c(with(deconvAndSimProportions, {sqrt(mean((deconvProps-simulProps)^2))/mean(simulProps) %>% round(.,4)}), #for all the proportions
                    map_dbl(ct, function(ct){with(deconvAndSimProportions %>% filter(cellType==ct),         #for each cell type
                                                  {sqrt(mean((deconvProps-simulProps)^2))/mean(simulProps) %>% round(.,4)})}))
  metrics$Correlation <- c(with(deconvAndSimProportions, {cor(deconvProps, simulProps) %>% round(.,4)}),           #for all the proportions
                           map_dbl(ct, function(ct){with(deconvAndSimProportions %>% filter(cellType==ct),         #for each cell type
                                                         {cor(deconvProps, simulProps) %>% round(.,4)})}))
  return(metrics)
}

reduceTo9 <- function(proportions){
  return(proportions %>% pivot_wider(names_from = cellType, values_from = deconvProps) %>% 
  mutate(`CD3+ T Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")) %>% 
  mutate(`NK Cells`=select(., c("CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
  mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
  mutate(`CD68+ (Foam Cells and inflam/resident macrophages)`=select(., c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
  mutate(`CD68+ Dendritic`=`CD68+CD1C+ Dendritic Cells`) %>% select(-c("CD68+CD1C+ Dendritic Cells")) %>% 
  mutate(`Endothelial Cells`=select(., c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% rowSums()) %>% select(-c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% 
  mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
  mutate(`CD68+ Mast Cells`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`) %>% 
  pivot_longer(-mixture, names_to = "cellType", values_to = "deconvProps"))
}

reduceTo6 <- function(proportions){
  return(proportions %>% pivot_wider(names_from = cellType, values_from = deconvProps) %>% 
    mutate(`T and NK Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
    mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
    mutate(`CD68+ Cells (no mast)`=select(., c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
    mutate(`Endothelial Cells`=select(., c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% rowSums()) %>% select(-c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% 
    mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
    mutate(`Mast Cells`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`) %>% 
    pivot_longer(-mixture, names_to = "cellType", values_to = "deconvProps"))
}

reduceTo8 <- function(proportions){
  return(proportions %>% pivot_wider(names_from = cellType, values_from = deconvProps) %>% 
    mutate(`T and NK Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
    mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
    mutate(`Macrophages`=select(., c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
    mutate(`Dendritic Cells`=`CD68+CD1C+ Dendritic Cells`) %>% select(-c("CD68+CD1C+ Dendritic Cells")) %>% 
    mutate(`Endothelial Cells I`=`CD34+ Endothelial Cells I`) %>% select(-c("CD34+ Endothelial Cells I")) %>% 
    mutate(`Endothelial Cells II`=`CD34+ Endothelial Cells II`) %>% select(-c("CD34+ Endothelial Cells II")) %>% 
    mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
    mutate(`Mastocytes`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`))
}
