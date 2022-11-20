###################################################################################
#                                                                                 #
# INTEGRATION OF RESULTS                                                          #                                       
#                                                                                 #
# Gemma Bel Bordes, Nov 20th 2022                                                 #
#                                                                                 #
###################################################################################

library(rstudioapi)
library(tidyverse)
library(magrittr)
library(cowplot)
library(Biobase)
library(grid)
library(gridExtra)
library(viridis)
#R version 4.1.2

setwd(dirname(getActiveDocumentContext()$path)) #set wd to the directory of this code file
source("helperFunctions/functions_IntegrationResults.R")
source("helperFunctions/functions_CreatingDifferentDataInputs.R")

# 1. Our bulk + Our SC ####
musicFolder <- "results/music/ourBulk_ourSC/onlyOriginalSamples/"
nnlsFolder <- "results/nnls/ourBulk_ourSC/onlyOriginalSamples/"
bisqueFolder <- "results/bisque/ourBulk_ourSC/onlyOriginalSamples/"
ciberFolder <- "results/cibersort/ourBulk_ourSC/onlyOriginalSamples/"
scadenFolder <- "results/scaden/ourBulk_ourSC/onlyOriginalSamples/"

## a) Integrate the results according to their clustering ####
clusters <- "cellType"
all18 <- rbind(integrate_Music(musicFolder, clusters, differentInputs=T), integrate_NNLS(nnlsFolder, clusters, differentInputs=T), integrate_Bisque(bisqueFolder, clusters, differentInputs=T), 
               integrate_Cibersort(ciberFolder, clusters, differentInputs=T), integrate_Scaden(scadenFolder, clusters, differentInputs=T)) %>% 
  rename(`CD68+IL18+ TLR4+TREM2+ Resident macrophages`=`CD68+IL18+TLR4+TREM2+ Resident macrophages`, `CD68+CASP1+ IL1B+SELL+ Inflammatory macrophages`=`CD68+CASP1+IL1B+SELL+ Inflammatory macrophages`, 
         `CD68+IL18+TLR4+ TREM2+ Resident macrophages`=`CD68+IL18+TLR4+TREM2+ Resident macrophages`, `CD68+ABCA1+ OLR1+TREM2+ Foam Cells`=`CD68+ABCA1+OLR1+TREM2+ Foam Cells`)
clusters <- "clusterType_9"
all9 <- rbind(integrate_Music(musicFolder, clusters, differentInputs=T), integrate_NNLS(nnlsFolder, clusters, differentInputs=T), integrate_Bisque(bisqueFolder, clusters, differentInputs=T), 
               integrate_Cibersort(ciberFolder, clusters, differentInputs=T), integrate_Scaden(scadenFolder, clusters, differentInputs=T))
clusters <- "clusterType_6"
all6 <- rbind(integrate_Music(musicFolder, clusters, differentInputs=T), integrate_NNLS(nnlsFolder, clusters, differentInputs=T), integrate_Bisque(bisqueFolder, clusters, differentInputs=T), 
               integrate_Cibersort(ciberFolder, clusters, differentInputs=T), integrate_Scaden(scadenFolder, clusters, differentInputs=T))

## b) Effect of qc for the sc data ####
boxPlotsProportions_differentQC(all18, "linear")
boxPlotsProportions_differentQC(all9, "linear")
boxPlotsProportions_differentQC(all6, "linear")

boxPlotsProportions_differentQC(all18, "log")
boxPlotsProportions_differentQC(all9, "log")
boxPlotsProportions_differentQC(all6, "log")

boxPlotsProportions_differentQC(all18, "sqrt")
boxPlotsProportions_differentQC(all9, "sqrt")
boxPlotsProportions_differentQC(all6, "sqrt")

## c) Effect of the data transformation ####
boxPlotsProportions_differentTransf(all18 %>% filter(ScQC=="AllGenes") %>% filter(Algorithm %in% c("NNLS", "MuSiC", "Bisque")))
boxPlotsProportions_differentTransf(all9 %>% filter(ScQC=="AllGenes") %>% filter(Algorithm %in% c("NNLS", "MuSiC", "Bisque")))
boxPlotsProportions_differentTransf(all6 %>% filter(ScQC=="AllGenes") %>% filter(Algorithm %in% c("NNLS", "MuSiC", "Bisque")))

## d) Effect of the cell type clustering ####
all18_reduced <- all18 %>% 
  mutate(`T and NK Cells`=select(., c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% rowSums()) %>% select(-c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI", "FOXP3+ T Cells", "CD3+CD56+ NK Cells II", "CD3+CD56+ NK Cells I")) %>% 
  mutate(`Switched mem B Cells`=`CD79A+ Class-switched Memory B Cells`) %>% select(-`CD79A+ Class-switched Memory B Cells`) %>% 
  mutate(`CD68+ Cells (no mast)`=select(., c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% rowSums()) %>% select(-c("CD68+CD1C+ Dendritic Cells", "CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")) %>% 
  mutate(`Endothelial Cells`=select(., c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% rowSums()) %>% select(-c("CD34+ Endothelial Cells I", "CD34+ Endothelial Cells II")) %>% 
  mutate(`Smooth Muscle Cells`=`ACTA2+ Smooth Muscle Cells`) %>% select(-`ACTA2+ Smooth Muscle Cells`) %>% 
  mutate(`Mast Cells`=`CD68+KIT+ Mast Cells`) %>% select(-`CD68+KIT+ Mast Cells`) %>% 
  filter(ScQC=="AllGenes") %>% select(-ScQC) %>% mutate(ClusteringDeconvolution="using 18 cell types")
  
all9_reduced <- all9 %>% 
  mutate(`T and NK Cells`=select(., c("CD3+ T Cells", "NK Cells", "FOXP3+ T Cells")) %>% rowSums()) %>% select(-c("CD3+ T Cells", "NK Cells", "FOXP3+ T Cells")) %>% 
  mutate(`CD68+ Cells (no mast)`=select(., c("CD68+ Dendritic", "CD68+ (Foam Cells and inflam/resident macrophages)")) %>% rowSums()) %>% select(-c("CD68+ Dendritic", "CD68+ (Foam Cells and inflam/resident macrophages)")) %>% 
  mutate(`Mast Cells`=`CD68+ Mast Cells`) %>% select(-`CD68+ Mast Cells`) %>%
  filter(ScQC=="AllGenes") %>% select(-ScQC) %>% mutate(ClusteringDeconvolution="using 9 clusters")

all6_reduced <- all6 %>% filter(ScQC=="AllGenes") %>% select(-ScQC) %>% mutate(ClusteringDeconvolution="using 6 clusters")

allReduced <- rbind(all18_reduced, all9_reduced, all6_reduced) %>% mutate(ClusteringDeconvolution=factor(ClusteringDeconvolution, levels=c("using 18 cell types", "using 9 clusters", "using 6 clusters")))
boxPlotsProportions_differentClustering(allReduced, "linear")

 #### for this analysis --> generation pseudo bulk data! (go to pseudoBulk.R)

## e) Visualize proportions for the chosen parameters so far (clustering: 8 populations -based on pseudoBulk analysis) ####
clusters <- "specialCluster" #only linear and all genes
all8 <- rbind(integrate_Music(musicFolder, clusters, differentInputs=T), integrate_NNLS(nnlsFolder, clusters, differentInputs=T), integrate_Bisque(bisqueFolder, clusters, differentInputs=T), 
              integrate_Cibersort(ciberFolder, clusters, differentInputs=T), integrate_Scaden(scadenFolder, clusters, differentInputs=T))
all8 %<>% filter(ScQC=="AllGenes", DataTrans=="linear") %>% select(-ScQC, -DataTrans)
all8 %>% write_csv("dataPlots/props8clustersAE.csv")

all8 %>% pivot_longer(!c("Patient", "Algorithm"), names_to = "CellType", values_to = "Proportion") %>% 
  ggplot(aes(x=Algorithm,y=Proportion)) + geom_boxplot(aes(fill=Algorithm), outlier.size = 0.2, alpha=0.7) +
  facet_wrap(~CellType, ncol = 4) + theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90)) +
  ggtitle(paste0("Proportions using 8 clusters")) + ylim(0,1) + scale_fill_brewer(palette = "Set3") + 
  theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("")

all8 %>% filter(Algorithm=="Scaden") %>% select(-Algorithm) %>% pivot_longer(-Patient, names_to = "CellType", values_to = "Proportion") %>% 
  ggplot(aes(x=CellType,y=Proportion)) + geom_boxplot(aes(fill=CellType), outlier.size = 0.2) + theme_bw() +
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")+ 
  ggtitle(paste0("Proportions for 8 clusters, scaden")) + ylim(0,1)

## e) Comparison with original (no qc) sc data-derived proportions... work with overlapping patients with both bulk and sc! ####
### results for only those patients with sc and bulk RNAseq data:
musicFolder <- "results/music/ourBulk_ourSC/scMatchingSamples/"
nnlsFolder <- "results/nnls/ourBulk_ourSC/scMatchingSamples/"
bisqueFolder <- "results/bisque/ourBulk_ourSC/scMatchingSamples/"
ciberFolder <- "results/cibersort/ourBulk_ourSC/scMatchingSamples/"
scadenFolder <- "results/scaden/ourBulk_ourSC/scMatchingSamples/"
clusters <- "specialCluster"
all8_MatchSamples <- rbind(integrate_Music(musicFolder, clusters), integrate_NNLS(nnlsFolder, clusters), integrate_Bisque(bisqueFolder, clusters), 
               integrate_Cibersort(ciberFolder, clusters), integrate_Scaden(scadenFolder, clusters)) 

### get proportions coming from the sc data:
scEset <- readRDS("output/ourSC/specialClustering/specialClustering_linear_scAllGenessc_eset.rds")
phenoSC <- pData(scEset)
scProps_8clust <- getScProportions(phenoSC, "specialCluster") #46 patients
length(intersect(scProps_8clust$Patient, unique(all8_MatchSamples$Patient))) #28 patients with both bulk and sc data 
scPropsAndDeconvProps_8 <- integrateSCandDeconvProps(scProps_8clust, all8_MatchSamples)
scPropsAndDeconvProps_8_allMetrics <- data.frame(resultsFor=character(), RMSE=double(), Correlation=double(), Algorithm=character())
for(algor in unique(scPropsAndDeconvProps_8$Algorithm)){
  scPropsAndDeconvProps_8_allMetrics <- rbind(scPropsAndDeconvProps_8_allMetrics, checkMetrics(scPropsAndDeconvProps_8 %>% filter(Algorithm==algor)) %>% 
                                                mutate(Algorithm=algor))
}
scPropsAndDeconvProps_8 %>% write_csv("dataPlots/propsScDec_28patients.csv")
scPropsAndDeconvProps_8_allMetrics %>% write_csv("dataPlots/metricsScDec_28patients.csv")

ggplot(scPropsAndDeconvProps_8_allMetrics %>% rename(cellType=resultsFor) %>%  mutate(cellType=factor(cellType, levels=c("Switched mem B Cells", "Mastocytes", "Macrophages", "Dendritic Cells", "T and NK Cells", "Endothelial Cells II", "Endothelial Cells I", "Smooth Muscle Cells", "all"))), 
       aes(x = Algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_fill_viridis_c(na.value = NA, limits=c(-1,1)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle(paste0("Deconvolution vs SC-dervied proportions")) + geom_hline(yintercept = 8.5, color="black", size=1.5)

scPropsAndDeconvProps_8 %>% pivot_longer(c("deconvProps", "scProps"), names_to="typeProps", values_to = "Props") %>%
  mutate(Algorithm=case_when(typeProps=="scProps"~"Simulation", TRUE~Algorithm)) %>% 
  ggplot(aes(x = Algorithm, y = Props)) + geom_boxplot(aes(fill=Algorithm), outlier.size = 0.2, alpha=0.7) +
  facet_wrap(~cellType, ncol = 4) + theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90)) +
  ggtitle("Deconvolution for overlapping patients vs sc proportions") + ylim(0,1) + scale_fill_manual(values = c("#8DD3C7", "#FFFFB3","#BEBADA","#FDB462","#80B1D3","#C90606"))

### are sc proportions reliable?
## Histological features for each sample:
table(phenoSC$SubjectName)
phenoSC %<>% filter(!SubjectName %in% c("ae4440", "ae4530", "ae4605", "ae4675")) #42 patients (as before)
phenoPlaque <- data.frame() #I want a df with histological info for each patient/sample
for (subject in unique(phenoSC$SubjectName)) {
  phenoPlaque <- rbind(phenoPlaque, phenoSC %>% filter(SubjectName==subject) %>% select(SubjectName,Phenotype,Sex,CD68.score,CD3.score,CD34.score, alpha.SMA.score, 
                                                                             SR.score,Glyc.c.score,Calcification) %>% head(1) %>% rename(Patient=SubjectName))
}
phenoPlaque <-  map_df(phenoPlaque, as.factor) %>% mutate(Patient=paste0("ae", Patient))
## relationship histological annotation - sc proportions (for 42 patients)
propsHisto_8clust <- as_tibble(merge(scProps_8clust, phenoPlaque, by="Patient"))
propsHisto_8clust %>% write_csv("dataPlots/scProps_histology.csv")
plot_grid(violinPlot_histoVsProps(propsHisto_8clust, "CD3.score", "T and NK Cells"),
          violinPlot_histoVsProps(propsHisto_8clust, "CD34.score", "Endothelial Cells I"),
          violinPlot_histoVsProps(propsHisto_8clust, "CD34.score", "Endothelial Cells II"),
          violinPlot_histoVsProps(propsHisto_8clust, "CD68.score", "Macrophages"),
          violinPlot_histoVsProps(propsHisto_8clust, "CD68.score", "Dendritic Cells"),
          violinPlot_histoVsProps(propsHisto_8clust, "alpha.SMA.score", "Smooth Muscle Cells"),
          ncol = 3)

## g) Relationship with histological data ####
histoClinical <- read_delim("input/clinical_data_good.txt")
length(intersect(histoClinical$study_number, unique(all8$Patient))) #654 patients with clinical data

featuresCategorical <- c("smc", "macrophages", "smc_macrophages_ratio",	"plaquephenotype", "fat",	"collagen", "media")
featuresContinuous <- c("vessel_density_averaged")

# histoClinical <- histoClinical %>% mutate(smc_macrophages_ratio=replace(smc_macrophages_ratio, smc_macrophages_ratio=="-888", NA)) %>% 
#   mutate(smc_macrophages_ratio=factor(smc_macrophages_ratio, levels = c("both not present", "macrophages dominant", "equal", "SMC dominant"))) %>% 
#   mutate(smc=factor(smc, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))) %>% 
#   mutate(macrophages=factor(macrophages, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))) %>% 
#   mutate(plaquephenotype=factor(plaquephenotype, levels = c("atheromatous", "fibroatheromatous", "fibrous"))) %>% 
#   mutate(fat=factor(fat, levels = c("no fat", "< 40% fat", "> 40% fat"))) %>% 
#   mutate(collagen=factor(collagen, levels = c("no staining", "minor staining", "moderate staining", "heavy staining"))) %>% 
#   mutate(media=factor(media, levels = c("no", "some remnants", "moderate presence media", "major parts media present"))) 

histoClinical %<>% mutate(ACTA2staining = factor(case_when(smc=="no staining"~0, smc=="minor staining"~1, smc=="moderate staining"~2, smc=="heavy staining"~3))) %>% 
  mutate(CD68staining = factor(case_when(macrophages=="no staining"~0, macrophages=="minor staining"~1, macrophages=="moderate staining"~2, macrophages=="heavy staining"~3))) %>% 
  mutate(PlaquePhenotype = factor(case_when(plaquephenotype=="atheromatous"~1, plaquephenotype=="fibroatheromatous"~2, plaquephenotype=="fibrous"~3))) %>% 
  mutate(ACTA2CD68ratio = factor(case_when(smc_macrophages_ratio=="macrophages dominant"~1, smc_macrophages_ratio=="equal"~2, smc_macrophages_ratio=="SMC dominant"~3))) %>% 
  mutate(Fat = factor(case_when(fat=="no fat"~0, fat=="< 40% fat"~1, fat=="> 40% fat"~2))) %>% 
  mutate(Collagen = factor(case_when(collagen=="no staining"~0, collagen=="minor staining"~1, collagen=="moderate staining"~2, collagen=="heavy staining"~3)))

colsTable <- sort(unique(all8$Algorithm))
correlationTable <- data.frame(matrix(nrow = 8, ncol = 5, dimnames=list(c("Macro_CD68staining", "Macro_fat", "Macro_ACTA2CD68ratio", "Macro_PlaquePhenotype", 
                                                                          "SMC_ACTA2staining", "SMC_collagen", "SMC_ACTA2CD68ratio", "SMC_PlaquePhenotype"), colsTable)))
for(algor in unique(all8$Algorithm)){
  propsHisto <- merge(all8 %>% filter(Algorithm==algor), histoClinical, by.x = "Patient", by.y = "study_number")
  resultsForAlgor <- c(kendallResults(propsHisto, "Macrophages", "CD68staining"), 
                       kendallResults(propsHisto, "Macrophages", "Fat"),
                       kendallResults(propsHisto, "Macrophages", "ACTA2CD68ratio"),
                       kendallResults(propsHisto, "Macrophages", "PlaquePhenotype"),
                       kendallResults(propsHisto, "Smooth Muscle Cells", "ACTA2staining"),
                       kendallResults(propsHisto, "Smooth Muscle Cells", "Collagen"),
                       kendallResults(propsHisto, "Smooth Muscle Cells", "ACTA2CD68ratio"),
                       kendallResults(propsHisto, "Smooth Muscle Cells", "PlaquePhenotype"))
  print(plot_grid(plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("CD68staining", "Fat", "ACTA2CD68ratio", "PlaquePhenotype"), "Macrophages", 
                                                                  algor, plotOnlySignificant = F, showpValue = F), nrow = 1), 
                  plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("ACTA2staining", "Collagen", "ACTA2CD68ratio", "PlaquePhenotype"), "Smooth Muscle Cells", 
                                                                  algor, plotOnlySignificant = F, showpValue = F), nrow = 1), 
                  nrow=2))
  correlationTable[[algor]] <- resultsForAlgor
}
View(correlationTable)

correlationTable2 <- data.frame(matrix(nrow = 2, ncol = 5, dimnames=list(c("EndoI_vesselDensity", "EndoII_vesselDensity"), colsTable)))
for(algor in unique(all8$Algorithm)){
  propsHisto <- merge(all8 %>% filter(Algorithm==algor), histoClinical, by.x = "Patient", by.y = "study_number")
  resultsForAlgor <- c(pearsonResults(propsHisto, "Endothelial Cells I", "vessel_density_averaged"), 
                       pearsonResults(propsHisto, "Endothelial Cells II", "vessel_density_averaged"))
  print(plot_grid(plot_grid(plotlist = analyzeContinuousFeatures(propsHisto, "vessel_density_averaged", "Endothelial Cells I",titlePlot = algor)),
                  plot_grid(plotlist = analyzeContinuousFeatures(propsHisto, "vessel_density_averaged", "Endothelial Cells II",titlePlot = algor)), ncol = 2))
  correlationTable2[[algor]] <- resultsForAlgor
}
View(correlationTable2)
rbind(correlationTable, correlationTable2)

## g) Analyze correlations of the proportions among different methods ####
clusters <- colnames(all8 %>% select(-c("Patient", "Algorithm", "DataTrans", "ScQC")))
listPlots <- list()
for(clust in clusters){
  if(clust != tail(clusters, 1)){
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all8, clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }else{
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all8, clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }
}
listPlots[["legend"]] <- ggplot(data.frame(Correlation=seq(-1, 1, 0.01)), aes(Correlation, Correlation, fill = Correlation))+
  geom_point() + scale_fill_viridis() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=12),
        panel.grid = element_blank(), panel.border = element_rect(colour = "white", fill='white', size=1), legend.key.size = unit(0.5, 'cm'))

plot_grid(plotlist = list(plot_grid(plotlist = listPlots[1:8], nrow=2), listPlots[[9]]), ncol = 2, rel_widths = c(8,1))


clusters <- colnames(all18 %>% filter(DataTrans=="linear", ScQC=="AllGenes") %>% select(-c("Patient", "Algorithm", "DataTrans", "ScQC")))
listPlots <- list()
for(clust in clusters){
  if(clust != tail(clusters, 1)){
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all18 %>% filter(DataTrans=="linear", ScQC=="AllGenes"), clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }else{
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all18 %>% filter(DataTrans=="linear", ScQC=="AllGenes"), clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }
}
listPlots[["legend"]] <- ggplot(data.frame(Correlation=seq(-1, 1, 0.01)), aes(Correlation, Correlation, fill = Correlation))+
  geom_point() + scale_fill_viridis() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=12),
        panel.grid = element_blank(), panel.border = element_rect(colour = "white", fill='white', size=1), legend.key.size = unit(0.5, 'cm'))

plot_grid(plotlist = list(plot_grid(plotlist = listPlots[1:18], nrow=3), listPlots[[19]]), ncol = 2, rel_widths = c(8,1))

## h) Relationship with marker genes expression ####
bulkExpr <- exprs(readRDS("output/ourBulk/onlyOriginalSamples/linearBulk_eset.rds")) #the original counts
normBulk <- normalizeTo(bulkExpr, "column")
all8_withExpr <- as_tibble(merge(all8, 
      as.data.frame(t(normBulk)) %>% rownames_to_column("Patient"),
      by="Patient"))
cowplot::plot_grid(regrLine_markerExpression(all8_withExpr, "CD68", "Macrophages"),
                   regrLine_markerExpression(all8_withExpr, "CD14", "Macrophages"))
cowplot::plot_grid(regrLine_markerExpression(all8_withExpr, "CD68", "Dendritic Cells"),
                   regrLine_markerExpression(all8_withExpr, "CD1C", "Dendritic Cells"))
cowplot::plot_grid(regrLine_markerExpression(all8_withExpr, "CD68", "Mastocytes"),
                   regrLine_markerExpression(all8_withExpr, "KIT", "Mastocytes"))
cowplot::plot_grid(regrLine_markerExpression(all8_withExpr, "CD3G", "T and NK Cells"),
                   regrLine_markerExpression(all8_withExpr, "CD3E", "T and NK Cells"),
                   regrLine_markerExpression(all8_withExpr, "CD3D", "T and NK Cells"),
                   regrLine_markerExpression(all8_withExpr, "FOXP3", "T and NK Cells"))
print(regrLine_markerExpression(all8_withExpr, "NCAM1", "T and NK Cells"))
print(regrLine_markerExpression(all8_withExpr, "ACTA2", "Smooth Muscle Cells"))
print(regrLine_markerExpression(all8_withExpr, "CD34", "Endothelial Cells I"))
print(regrLine_markerExpression(all8_withExpr, "CD34", "Endothelial Cells II"))
print(regrLine_markerExpression(all8_withExpr, "CD79A", "Switched mem B Cells"))


# 2. Our bulk + Huize Pan single cell ####
musicFolder <- "results/music/ourBulk_hpSC/"
bisqueFolder <- "results/bisque/ourBulk_hpSC/"
ciberFolder <- "results/cibersort/ourBulk_hpSC/"
scadenFolder <- "results/scaden/ourBulk_hpSC/groupedMacro/"
nnlsFolder <- "results/nnls/ourBulk_hpSC/"

## a) Integrate the results ####
all <- rbind(integrate_Music(musicFolder, clustering="", differentInputs=F), integrate_Bisque(bisqueFolder, clustering="", differentInputs=F), 
             integrate_NNLS(nnlsFolder, clustering="", differentInputs = F),
              integrate_Cibersort(ciberFolder, clustering="", differentInputs=F), integrate_Scaden(scadenFolder, clustering="", differentInputs=F)) 
## b) Visualize proportions for the methods ####
all %>% pivot_longer(!c("Patient", "Algorithm"), names_to = "CellType", values_to = "Proportion") %>% 
  ggplot(aes(x=Algorithm,y=Proportion)) + geom_boxplot(aes(fill=Algorithm), outlier.size = 0.2, alpha=0.7) +
  facet_wrap(~CellType, ncol = 4) + theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90)) +
  ggtitle(paste0("Proportions using 8 clusters")) + ylim(0,1) + scale_fill_brewer(palette = "Set3") + 
  theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("")

## c) Analyze correlations of the proportions among different methods ####
clusters <- colnames(all %>% select(-c("Patient", "Algorithm")))
listPlots <- list()
for(clust in clusters){
  if(clust != tail(clusters, 1)){
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all, clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }else{
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(all, clust) + theme(legend.position = "none", plot.margin = margin(0,0,0,0))
  }
}
listPlots[["legend"]] <- ggplot(data.frame(Correlation=seq(-1, 1, 0.01)), aes(Correlation, Correlation, fill = Correlation))+
  geom_point() + scale_fill_viridis() + 
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=12),
        panel.grid = element_blank(), panel.border = element_rect(colour = "white", fill='white', size=1), legend.key.size = unit(0.5, 'cm'))
plot_grid(plotlist = list(plot_grid(plotlist = listPlots[1:10], nrow=2), listPlots[[11]]), ncol = 2, rel_widths = c(8,1))

## d) Relationship with histological data ####
colsTable <- sort(unique(all$Algorithm))
correlationTable <- data.frame(matrix(nrow = 8, ncol = 5, dimnames=list(c("Macro_CD68staining", "Macro_fat", "Macro_ACTA2CD68ratio", "Macro_PlaquePhenotype", 
                                                                          "SMC_ACTA2staining", "SMC_collagen", "SMC_ACTA2CD68ratio", "SMC_PlaquePhenotype"), colsTable)))
for(algor in unique(all$Algorithm)){
  propsHisto <- merge(all %>% filter(Algorithm==algor), histoClinical, by.x = "Patient", by.y = "study_number")
  resultsForAlgor <- c(kendallResults(propsHisto, "Macrophages", "CD68staining"), 
                       kendallResults(propsHisto, "Macrophages", "Fat"),
                       kendallResults(propsHisto, "Macrophages", "ACTA2CD68ratio"),
                       kendallResults(propsHisto, "Macrophages", "PlaquePhenotype"),
                       kendallResults(propsHisto, "SMC", "ACTA2staining"),
                       kendallResults(propsHisto, "SMC", "Collagen"),
                       kendallResults(propsHisto, "SMC", "ACTA2CD68ratio"),
                       kendallResults(propsHisto, "SMC", "PlaquePhenotype"))
  print(plot_grid(plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("CD68staining", "Fat", "ACTA2CD68ratio", "PlaquePhenotype"), "Macrophages", 
                                                                  algor, plotOnlySignificant = F, showpValue = F), nrow = 1), 
                  plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("ACTA2staining", "Collagen", "ACTA2CD68ratio", "PlaquePhenotype"), "SMC", 
                                                                  algor, plotOnlySignificant = F, showpValue = F), nrow = 1), 
                  nrow=2))
  correlationTable[[algor]] <- resultsForAlgor
}
correlationTable2 <- data.frame(matrix(nrow = 2, ncol = 5, dimnames=list(c("EndoI_vesselDensity", "EndoII_vesselDensity"), colsTable)))
for(algor in unique(all$Algorithm)){
  propsHisto <- merge(all %>% filter(Algorithm==algor), histoClinical, by.x = "Patient", by.y = "study_number")
  resultsForAlgor <- c(pearsonResults(propsHisto, "Endothelial 1", "vessel_density_averaged"), 
                       pearsonResults(propsHisto, "Endothelial 2", "vessel_density_averaged"))
  print(plot_grid(plot_grid(plotlist = analyzeContinuousFeatures(propsHisto, "vessel_density_averaged", "Endothelial 1",titlePlot = algor)),
                  plot_grid(plotlist = analyzeContinuousFeatures(propsHisto, "vessel_density_averaged", "Endothelial 2",titlePlot = algor)), ncol = 2))
  correlationTable2[[algor]] <- resultsForAlgor
}
rbind(correlationTable, correlationTable2)
                       
## e) Relationship with marker genes expression ####
all_withExpr <- as_tibble(merge(all8, 
                                 as.data.frame(t(normBulk)) %>% rownames_to_column("Patient"),
                                 by="Patient"))
cowplot::plot_grid(regrLine_markerExpression(all_withExpr, "CD68", "Macrophages"),
                   regrLine_markerExpression(all_withExpr, "CD14", "Macrophages"))
cowplot::plot_grid(regrLine_markerExpression(all_withExpr, "CD68", "Dendritic Cells"),
                   regrLine_markerExpression(all_withExpr, "CD1C", "Dendritic Cells"))
cowplot::plot_grid(regrLine_markerExpression(all_withExpr, "CD68", "Mastocytes"),
                   regrLine_markerExpression(all_withExpr, "KIT", "Mastocytes"))
cowplot::plot_grid(regrLine_markerExpression(all_withExpr, "CD3G", "T and NK Cells"),
                   regrLine_markerExpression(all_withExpr, "CD3E", "T and NK Cells"),
                   regrLine_markerExpression(all_withExpr, "CD3D", "T and NK Cells"),
                   regrLine_markerExpression(all_withExpr, "FOXP3", "T and NK Cells"))
print(regrLine_markerExpression(all_withExpr, "NCAM1", "T and NK Cells"))
print(regrLine_markerExpression(all_withExpr, "ACTA2", "Smooth Muscle Cells"))
print(regrLine_markerExpression(all_withExpr, "CD34", "Endothelial Cells I"))
print(regrLine_markerExpression(all_withExpr, "CD34", "Endothelial Cells II"))
print(regrLine_markerExpression(all_withExpr, "CD79A", "Switched mem B Cells"))

## f) Correlation between these proportions and proportions in section (1) ####
colsTable <- sort(unique(all$Algorithm))
correlationTable <- data.frame(matrix(nrow = 6, ncol = 5, dimnames=list(c("Macrophages", "Smooth Muscle Cells", "Mastocytes", "T cells", 
                                                                          "Endothelial I", "Endothelial II"), colsTable)))

smcComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("Smooth Muscle Cells", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Smooth Muscle Cells"), 
                         all %>% filter(Algorithm==algor) %>% pull("SMC"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Smooth Muscle Cells"), 
                         all %>% filter(Algorithm==algor) %>% pull("SMC"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  smcComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

macroComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("Macrophages", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Macrophages"), 
                         all %>% filter(Algorithm==algor) %>% pull("Macrophages"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Macrophages"), 
                         all %>% filter(Algorithm==algor) %>% pull("Macrophages"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  macroComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

mastoComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("Mastocytes", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Mastocytes"), 
                         all %>% filter(Algorithm==algor) %>% pull("Mast cell"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Mastocytes"), 
                         all %>% filter(Algorithm==algor) %>% pull("Mast cell"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  mastoComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

tComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("T Cells", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("T and NK Cells"), 
                         all %>% filter(Algorithm==algor) %>% pull("T cell"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("T and NK Cells"), 
                         all %>% filter(Algorithm==algor) %>% pull("T cell"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  tComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

endoIComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("Endothelial Cells I", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Endothelial Cells I"), 
                         all %>% filter(Algorithm==algor) %>% pull("Endothelial 1"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Endothelial Cells I"), 
                         all %>% filter(Algorithm==algor) %>% pull("Endothelial 1"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  endoIComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

endoIIComparison <- data.frame(matrix(nrow = 1, ncol = 5, dimnames=list("Endothelial Cells II", colsTable)))
for(algor in unique(all$Algorithm)){
  corr <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Endothelial Cells II"), 
                         all %>% filter(Algorithm==algor) %>% pull("Endothelial 2"), method = "pearson")$estimate, 3)
  pval <- round(cor.test(all8 %>% filter(Algorithm==algor) %>% pull("Endothelial Cells II"), 
                         all %>% filter(Algorithm==algor) %>% pull("Endothelial 2"), method = "pearson")$p.value, 6)
  psymbol <-  case_when(pval>0.1~" ",  pval<0.01~"***", pval<0.05~"**", pval<0.1~"*")
  endoIIComparison[[algor]] <- (paste0(corr, " (", psymbol, ")"))
}

tableHP_AE <- rbind(tComparison, macroComparison, smcComparison, endoIComparison, endoIIComparison, mastoComparison)
tableHP_AE %>% rownames_to_column("Cell population")%>% write_csv("dataPlots/correlation_AEvsHPdeconv.csv", )

