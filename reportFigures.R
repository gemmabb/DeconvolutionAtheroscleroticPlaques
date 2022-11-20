###################################################################################
# REPORT FIGURES                                                                  #
#                                                                                 #
# Gemma Bel Bordes                                                                #
###################################################################################

library(rstudioapi)
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggcorrplot)
library(Seurat)
library(scales)
library(viridis)
library(pheatmap)
library(ggstatsplot)
library(rstatix)
library(survival)
library(ggfortify)
library(survminer)
#R version 4.1.2

setwd(dirname(getActiveDocumentContext()$path)) #set wd to the directory of this code file
swr = function(string, nwrap=10) {paste(strwrap(string, width=nwrap), collapse="\n")}
swr <- Vectorize(swr)

# FIG 1 ####
## A (UMAP, 18 cell types) ####
load("input/20210811.46.patients.Koen.rds")
seuset <- seuset[,!colnames(seuset) %in% colnames(seuset)[Idents(seuset) %in% c("CD79+ Plasma B Cells", "CD68+CD4+ Monocytes")]] #remove cell types that are not used for deconv.
ct <- Idents(seuset)
numCellType <- case_when(ct=="CD3+ T Cells I"~1,ct=="CD3+ T Cells II"~2,ct=="CD3+ T Cells III"~3,ct=="CD3+ T Cells IV"~4,
                         ct=="CD3+ T Cells V"~5,ct=="CD3+ T Cells VI"~6,ct=="FOXP3+ T Cells"~7,ct=="CD3+CD56+ NK Cells I"~8,ct=="CD3+CD56+ NK Cells II"~9,
                         ct=="CD79A+ Class-switched Memory B Cells"~10,ct=="CD68+ABCA1+OLR1+TREM2+ Foam Cells"~11,ct=="CD68+IL18+TLR4+TREM2+ Resident macrophages"~12,ct=="CD68+CASP1+IL1B+SELL+ Inflammatory macrophages"~13,
                         ct=="CD68+CD1C+ Dendritic Cells"~14,ct=="ACTA2+ Smooth Muscle Cells"~15,ct=="CD34+ Endothelial Cells I"~16,ct=="CD34+ Endothelial Cells II"~17,ct=="CD68+KIT+ Mast Cells"~18)
newLabel <- paste0(numCellType, ": ", ct)
Idents(seuset) <- numCellType
F1A <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), aspect.ratio = 1)

## B (% sequenced cells) ####
hue_pal()(18) #extract colours from the palette used in the umap
F1B <- as.data.frame(table(Idents(seuset))) %>% mutate(CellType=factor(Var1, levels = sort(as.numeric(Var1),decreasing = T))) %>% mutate(Proportion=100*Freq/sum(Freq)) %>% 
  ggplot(aes(x=Proportion, y=CellType))+geom_bar(stat="identity", fill = "gray")+theme_bw()+theme(legend.position = "none") + geom_vline(xintercept=5, linetype="dashed")+ylab("Cell type")+xlab("% sequenced cells")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,2), "cm"), axis.title.y =element_blank())+coord_cartesian(xlim = c(0, 15), ylim = c(1,18), clip = "off")+
  annotate("point", x = -2, y = 18, size = 3, color = "#F8766D")+  
  annotate("point", x = -2, y = 17, size = 3, color = "#E88526")+
  annotate("point", x = -2, y = 16, size = 3, color = "#D39200")+
  annotate("point", x = -2, y = 15, size = 3, color = "#B79F00")+
  annotate("point", x = -2, y = 14, size = 3, color = "#93AA00")+
  annotate("point", x = -2, y = 13, size = 3, color = "#5EB300")+
  annotate("point", x = -2, y = 12, size = 3, color = "#00BA38")+
  annotate("point", x = -2, y = 11, size = 3, color = "#00BF74")+
  annotate("point", x = -2, y = 10, size = 3, color = "#00C19F")+
  annotate("point", x = -2, y = 9, size = 3, color = "#00BFC4")+
  annotate("point", x = -2, y = 8, size = 3, color = "#00B9E3")+
  annotate("point", x = -2, y = 7, size = 3, color = "#00ADFA")+
  annotate("point", x = -2, y = 6, size = 3, color = "#619CFF")+
  annotate("point", x = -2, y = 5, size = 3, color = "#AE87FF")+
  annotate("point", x = -2, y = 4, size = 3, color = "#DB72FB")+
  annotate("point", x = -2, y = 3, size = 3, color = "#F564E3")+
  annotate("point", x = -2, y = 2, size = 3, color = "#FF61C3")+
  annotate("point", x = -2, y = 1, size = 3, color = "#FF699C")

## C (differences using all vs detect genes) ###
props18 <- read_csv("dataPlots/all18Props.csv")
df <- props18  %>% filter(DataTrans=="linear") %>% pivot_longer(!c("Patient", "Algorithm", "DataTrans", "ScQC"), names_to = "CellType", values_to = "Proportion") %>% 
  mutate(CellType=factor(swr(CellType), levels = c("CD3+ T\nCells I","CD3+ T\nCells II","CD3+ T\nCells III","CD3+ T\nCells IV","CD3+ T\nCells V","CD3+ T\nCells VI","FOXP3+ T\nCells","CD3+CD56+\nNK Cells\nI","CD3+CD56+\nNK Cells\nII",
                                              "CD79A+\nClass-switched\nMemory B\nCells","CD68+ABCA1+\nOLR1+TREM2+\nFoam\nCells","CD68+IL18+TLR4+\nTREM2+\nResident\nmacrophages","CD68+CASP1+\nIL1B+SELL+\nInflammatory\nmacrophages",
                                              "CD68+CD1C+\nDendritic\nCells","ACTA2+\nSmooth\nMuscle\nCells","CD34+\nEndothelial\nCells I","CD34+\nEndothelial\nCells II","CD68+KIT+\nMast\nCells"))) %>% 
  rename(`Gene selection`=ScQC)
F1C <- (df %>% 
        ggplot(aes(x=Algorithm,y=Proportion)) +
        geom_boxplot(aes(fill=`Gene selection`), outlier.size = 0.2, alpha=0.7) +
        facet_wrap(~CellType, ncol = ceiling(length(unique(df$CellType))/2)) +
        theme_bw() + theme(strip.text = element_text(size = 7.5), axis.title.x = element_blank(),
                           axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1), legend.position = "bottom") + scale_fill_manual(values=c("darkolivegreen1", "darkolivegreen4")) +
        ylim(0,1))

## Merge ####
png("Fig1.png", width = 900 , height = 1000)
plot_grid(plot_grid(F1A, F1B, nrow = 1, labels = c("A", "B"), rel_widths = c(1.5,1)),
          F1C, nrow = 2, labels = c("","D"), rel_heights = c(1,1.8))
dev.off()  

# FIG 2 ####
## A (differences AE data transformed)
df <- props18 %>% filter(Algorithm %in% c("NNLS","MuSiC","Bisque")) %>% 
    pivot_longer(!c("Patient", "Algorithm", "DataTrans", "ScQC"), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=factor(swr(CellType), levels = c("CD3+ T\nCells I","CD3+ T\nCells II","CD3+ T\nCells III","CD3+ T\nCells IV","CD3+ T\nCells V","CD3+ T\nCells VI","FOXP3+ T\nCells","CD3+CD56+\nNK Cells\nI","CD3+CD56+\nNK Cells\nII",
                                                   "CD79A+\nClass-switched\nMemory B\nCells","CD68+ABCA1+\nOLR1+TREM2+\nFoam\nCells","CD68+IL18+TLR4+\nTREM2+\nResident\nmacrophages","CD68+CASP1+\nIL1B+SELL+\nInflammatory\nmacrophages",
                                                   "CD68+CD1C+\nDendritic\nCells","ACTA2+\nSmooth\nMuscle\nCells","CD34+\nEndothelial\nCells I","CD34+\nEndothelial\nCells II","CD68+KIT+\nMast\nCells"))) 
F2A <- df %>% ggplot(aes(x=Algorithm,y=Proportion)) +
              geom_boxplot(aes(fill=DataTrans), outlier.size = 0.2, alpha=0.7) +
              facet_wrap(~CellType, ncol = ceiling(length(unique(df$CellType))/2)) +
              theme_bw() + scale_fill_manual(values=c("#fffa60", "#ffbe00", "#a77d00")) + 
              theme(strip.text = element_text(size = 6), 
                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.title.x = element_blank(), legend.position = "none") + ylim(0,1)

## B (NRMSE in pseudo bulk)
F2B <- read_csv("dataPlots/pseudoBulk_dataTrans.csv") %>% filter(resultsFor=="all") %>% 
  ggplot(aes(x=Algorithm, y=NRMSE, fill=dataTrans)) + geom_boxplot(alpha=0.8) + 
  scale_fill_manual(values=c("#fffa60", "#ffbe00", "#a77d00")) + labs(fill = "Data\ntransformation") + theme_bw() + theme(axis.title.x = element_blank()) #report

## Merge ####
tiff("Fig2.png", width = 900 , height = 700)
plot_grid(F2A, plot_grid(F2B, NULL, nrow = 1, labels = c("B", "")),
          nrow = 2, labels = c("A",""), rel_heights = c(1.8,1))
dev.off()  

# FIG 3 ####
## A ####
clust_results <- read_csv("dataPlots/pseudoBulk_Clust.csv")
clust_results %<>% filter(cellType!="all")
plot18 <- ggplot(clust_results %>% filter(Clustering=="cellType") %>% mutate(cellType=factor(cellType, 
                                                                                             levels=rev(c("CD3+ T Cells I","CD3+ T Cells II","CD3+ T Cells III","CD3+ T Cells IV",
                                                                                                                    "CD3+ T Cells V","CD3+ T Cells VI","FOXP3+ T Cells","CD3+CD56+ NK Cells I","CD3+CD56+ NK Cells II",
                                                                                                                    "CD79A+ Class-switched Memory B Cells","CD68+ABCA1+OLR1+TREM2+ Foam Cells","CD68+IL18+TLR4+TREM2+ Resident macrophages","CD68+CASP1+IL1B+SELL+ Inflammatory macrophages",
                                                                                                                    "CD68+CD1C+ Dendritic Cells","ACTA2+ Smooth Muscle Cells","CD34+ Endothelial Cells I","CD34+ Endothelial Cells II","CD68+KIT+ Mast Cells")))), 
                 aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_fill_viridis_c(na.value = "white") + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## 9 cell types:
plot9 <- ggplot(clust_results %>% filter(Clustering=="clusterType_9") %>% mutate(cellType=factor(cellType, levels=rev(c("null1","null2","CD3+ T Cells","null3",
                                                                                                                        "null4","null5","FOXP3+ T Cells","NK Cells","null6",
                                                                                                                        "Switched mem B Cells","null7","CD68+ (Foam Cells and inflam/resident macrophages)","null8", 
                                                                                                                        "CD68+ Dendritic","Smooth Muscle Cells","Endothelial Cells","null9","CD68+ Mast Cells")))), 
                aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_y_discrete(drop=FALSE)+
  scale_fill_viridis_c(na.value = NA) + theme(axis.text.x = element_text(angle = 20, vjust = 0.5, hjust=1)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## 6 cell types:
plot6 <- ggplot(clust_results %>% filter(Clustering=="clusterType_6") %>% mutate(cellType = case_when(cellType=="CD68+ Cells (no mast)"~"CD68+ Cells",T~cellType)) %>% mutate(cellType=factor(cellType, levels=rev(c("null1","null2","nullx","null3",
                                                                                                                        "T and NK Cells","null5","nully","nullz","nulla",
                                                                                                                        "Switched mem B Cells","null7","CD68+ Cells","null8", 
                                                                                                                        "nullw","Smooth Muscle Cells","Endothelial Cells","null9","Mast Cells")))), aes(x = algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_y_discrete(drop=FALSE)+
  scale_fill_viridis_c(na.value = NA) + theme(axis.text.y = element_text(angle = 20, vjust = 1, hjust=1)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), axis.title.y=element_blank(), axis.title.x=element_blank())
## all plots combined... they must have the same 1/NRMSE legend!:
max(1/clust_results$NRMSE)
min(1/clust_results$NRMSE)
F3A <- egg:::ggarrange(plot18 + ggtitle("18 clusters") + theme(plot.title = element_text(size=12),legend.position = "none", axis.ticks.y=element_blank(), plot.margin = margin(10,0,10,10))+geom_hline(yintercept=seq(from = 1.5, to = 17.5, by = 1), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                plot9 + ggtitle(" 9 clusters") + theme(plot.title = element_text(size=12),legend.position = "none", axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = margin(10,0,10,0)) + geom_hline(yintercept=c(1.5, 3.5, 4.5, 5.5, 8.5, 9.5, 11.5,12.5), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                plot6+ ggtitle("  6 clusters") + scale_y_discrete(labels=case_when(startsWith(levels(plot6$data$cellType), "null")~"", TRUE~levels(plot6$data$cellType)),position = "right", drop=F) + theme(plot.title = element_text(size=12),axis.title.y=element_blank(), axis.ticks.y=element_blank(), plot.margin = margin(10,10,10,2.5))+geom_hline(yintercept=c(1.5,3.5,4.5, 8.5,9.5), alpha=0.5)+scale_size_continuous(limits=c(0.4, 9.8)), 
                nrow = 1)

## B ####
ct_8 <- case_when(ct %in% c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI","CD3+CD56+ NK Cells I", "CD3+CD56+ NK Cells II","FOXP3+ T Cells")~"T and NK Cells", 
                  ct=="CD79A+ Class-switched Memory B Cells"~"Switched-memory B Cells", 
                  ct=="CD68+CD1C+ Dendritic Cells"~"Dendritic Cells", 
                  ct=="CD34+ Endothelial Cells I"~"Endothelial Cells I",
                  ct=="CD34+ Endothelial Cells II"~"Endothelial Cells II",
                  ct=="ACTA2+ Smooth Muscle Cells"~"Smooth Muscle Cells", ct=="CD68+KIT+ Mast Cells"~"Mastocytes", ct %in% c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")~"Macrophages")
numCellType_8 <- case_when(ct_8=="T and NK Cells"~1, ct_8=="Switched-memory B Cells"~2, ct_8=="Macrophages"~3, ct_8=="Dendritic Cells"~4,
                           ct_8=="Smooth Muscle Cells"~5, ct_8=="Endothelial Cells I"~6, ct_8=="Endothelial Cells II"~7,  ct_8=="Mastocytes"~8)
newLabel <- paste0(numCellType_8, ": ", ct_8)
Idents(seuset) <- numCellType_8
F3B <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "bottom", aspect.ratio = 1)+guides(colour=guide_legend(ncol=2, override.aes = list(size=3)))

png("Fig3.png", width = 600 , height = 700)
plot_grid(F3A, F3B, nrow = 2, labels = c("A", "B"), rel_heights = c(3,1))
dev.off()  

# FIG 4 (sc proportions) ####
scPropsAndDeconvProps_8_allMetrics <- read_csv("dataPlots/metricsScDec_28patients.csv") 
F4 <- ggplot(scPropsAndDeconvProps_8_allMetrics %>% rename(cellType=resultsFor) %>% 
               mutate(cellType=factor(case_when(cellType=="Switched mem B Cells"~"Switched-memory B Cells", TRUE~cellType), levels=rev(c("all", "T and NK Cells", "Switched-memory B Cells", "Macrophages", "Dendritic Cells", "Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes")))), 
       aes(x = Algorithm, y = cellType)) + geom_point(aes(size = 1/NRMSE, fill = Correlation), shape = 21) +
  scale_fill_viridis_c(na.value = NA, limits=c(-1,1)) + theme_bw() + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = 8.5, color="black", size=1.5)
png("Fig4.png", width = 400 , height = 300)
F4
dev.off() 

# FIG 5 (8 cluster-proportions vs histology) ####
props8 <- read_csv("dataPlots/props8clustersAE.csv") %>% rename(`Switched-memory\nB Cells`=`Switched mem B Cells`) %>% rename(`Smooth Muscle\nCells`=`Smooth Muscle Cells`)
# A ####
F5A <- props8 %>% pivot_longer(!c("Patient", "Algorithm"), names_to = "cellType", values_to = "Proportion") %>% 
  mutate(cellType=factor(cellType, levels=c("T and NK Cells", "Switched-memory\nB Cells", "Macrophages", "Dendritic Cells", "Smooth Muscle\nCells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes"))) %>% 
  ggplot(aes(x=Algorithm,y=Proportion)) + geom_boxplot(aes(fill=Algorithm), outlier.size = 0.2, alpha=0.7) +
  facet_wrap(~cellType, ncol = 8) + theme_bw() + theme(legend.position = "none", strip.text = element_text(size = 8), axis.text.x = element_text(angle = 90), axis.title.x = element_blank()) +
  ylim(0,1) + scale_fill_brewer(palette = "Set3") 


# B ####
props8 <- read_csv("dataPlots/props8clustersAE.csv") %>%  rename(`Smooth Muscle\nCells`=`Smooth Muscle Cells`, `Switched-memory\nB Cells`=`Switched mem B Cells`, `T and NK Cells\n`=`T and NK Cells`, `Macrophages\n`=`Macrophages`,
                   `Dendritic Cells\n`=`Dendritic Cells`, `Endothelial Cells I\n`=`Endothelial Cells I`, `Endothelial Cells II\n`=`Endothelial Cells II`, `Mastocytes\n`=`Mastocytes`)
analyzeCorrelationProps_CellType <- function(dataset, cellType){
  subdata <- dataset %>% select(Patient, cellType, Algorithm) %>% pivot_wider(names_from = Algorithm, values_from = cellType)
  pearsons <- ggcorrmat(data = subdata, output = "dataframe", type = "parametric") #parametric for Pearson
  return(ggcorrmat(data = subdata, pch=7, size=2,output = "plot", type = "parametric", colors = c("#440154", "#21918c", "#fde725"), k=2, 
                  ggcorrplot.args = list(lab_size = 2.5)) + ggtitle(cellType) + theme(axis.text.x = element_text(size = 8), 
                                                                  axis.text.y = element_text(size = 8), plot.margin = margin(1,1,1,1),
                                                                  legend.position = "none", plot.title = element_text(size = 10)))
}
listPlots <- list()
clusters <- c("T and NK Cells\n", "Switched-memory\nB Cells", "Macrophages\n", "Dendritic Cells\n", "Smooth Muscle\nCells", "Endothelial Cells I\n", "Endothelial Cells II\n", "Mastocytes\n")
for(clust in clusters){
    listPlots[[clust]] <- analyzeCorrelationProps_CellType(props8, clust)
}
listPlots[["legend1"]] <- ggplot(data.frame(Correlation=seq(-1, 1, 0.01)), aes(Correlation, Correlation, fill = Correlation))+
  geom_point() + scale_fill_viridis() + theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                                              legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=10), 
                                              panel.border = element_rect(colour = "white", fill='white', size=1),
                                              legend.key.size = unit(0.7, 'cm'), plot.margin = margin(0,0,0,0), legend.box = "horizontal", legend.direction = "horizontal") 
listPlots[["legend2"]] <- ggplot(data.frame(Correlation=seq(-1, 1, 0.01)), aes(Correlation, Correlation))+
  geom_point(color="white") + coord_cartesian(xlim = c(-1, 1), ylim = c(-1,1), clip = "off") + theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), aspect.ratio = 0.5) +
  annotate("point", x=-2, y=0, shape=7, size=7) + annotate("text", x=-0.02, y=0, label="= correlation non-significant (p < 0.05; Holm adjustment)", size=3.5)

F5B <- plot_grid(plotlist = list(plot_grid(plotlist = listPlots[1:8], nrow=2), 
                                 plot_grid(plotlist = listPlots[9:10], ncol=2)), ncol = 1, rel_heights = c(4,1))
# merge ####
png("Fig5.png", width = 650 , height = 700)
plot_grid(F5A, F5B, nrow = 2, labels = c("A", "B"), rel_heights = c(0.5,1))
dev.off()  

# FIG 6 (Sex differences) ####
scaden8props <- read_csv("dataPlots/props8scaden_10runs.csv") %>% rename(`Switched-memory B Cells`=`Switched mem B Cells`)
boxPlot_AllCellType_2groups <- function(dataset, cells, feature, removeNonSign = T, order=unique(cells)){
  data <- dataset %>% select(cells, "Patient", feature) %>% pivot_longer(-c("Patient", feature), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=factor(CellType, levels=order)) 
  stat_pvalues <- data %>% group_by(CellType) %>% pairwise_wilcox_test(as.formula(paste0("Proportion~", feature))) %>% 
    add_xy_position("CellType") %>% mutate(p=paste0("p = ", format.pval(p, eps=.001, digits=3, nsmall=3)))
  if(removeNonSign) {stat_pvalues %<>% filter(p.adj<0.1)}
  data %>% ggplot(aes(x=CellType,y=Proportion)) +
    geom_boxplot(aes_string(fill=feature), outlier.size = 0.2)  + theme_bw() + scale_fill_manual(values = c("#5D3A9B", "#E66001")) +
    theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          legend.position="right") + 
    xlab("") + stat_pvalue_manual(stat_pvalues, label = "p", size = 3.4) 
}
## A ####
F6A <- boxPlot_AllCellType_2groups(scaden8props, colnames(scaden8props)[2:9], "Sex", removeNonSign = F, 
                            order = c("T and NK Cells", "Switched-memory B Cells", "Macrophages", "Dendritic Cells", "Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes")) + ylim(0,1)
## B ####
scaden8props %<>% mutate(`SMC:Macro`=`Smooth Muscle Cells`/Macrophages) %>% 
  mutate(totalStruct= scaden8props %>% select(c("Endothelial Cells I", "Smooth Muscle Cells", "Endothelial Cells II")) %>% rowSums()) %>% 
  mutate(totalImmune= scaden8props %>% select(c("Macrophages", "Dendritic Cells", "Mastocytes", "Switched-memory B Cells", "T and NK Cells")) %>% rowSums()) %>% 
  mutate(`Struct:Immune` = totalStruct/totalImmune) 
F6B <- boxPlot_AllCellType_2groups(scaden8props, c("Struct:Immune"), "Sex", removeNonSign = F) + ylim(0,3) +
  ylab("Structural/Immune ratio") + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12))
## C ####
F6C <- boxPlot_AllCellType_2groups(scaden8props, c("SMC:Macro"), "Sex", removeNonSign = F) + ylim(0,8.2) +
  ylab("SMC/Macrophages ratio") + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12))
## merge ####
png("Fig6.png", width = 650, height = 600)
plot_grid(F6A, plot_grid(F6B, F6C, ncol=2, labels = c("B", "C")), nrow = 2, labels = c("A", ""), rel_heights = c(1,0.5))
dev.off()  

# FIG 7 (Cox hazard models) ####
macrosMace <- read_csv("dataPlots/macrosMace.csv") %>% mutate(Macrophages_tertile=factor(Macrophages_tertile, levels = c("low", "intermediate", "high")),
                                                              foam_tertile=factor(foam_tertile, levels = c("low", "intermediate", "high")),
                                                              res_tertile=factor(res_tertile, levels = c("low", "intermediate", "high")),
                                                              inf_tertile=factor(inf_tertile, levels = c("low", "intermediate", "high")))
  
## A ####
cox <- coxph(Surv(ep_major_t_3years, MACE)~Macrophages_tertile, data=macrosMace) 
F7A <- autoplot(survfit(Surv(ep_major_t_3years, MACE)~Macrophages_tertile, data=macrosMace)) +
  labs(x="Time (years)", y="Event-free survival", color = "Tertile", fill = "Tertile") + 
  theme_classic() + theme(legend.position = "bottom")  + ggtitle("Macrophage (single cluster)")+
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["Macrophages_tertilehigh","Pr(>|z|)"], 3)))
## B ####
cox <- coxph(Surv(ep_major_t_3years, MACE)~foam_tertile, data=macrosMace) 
F7B <- autoplot(survfit(Surv(ep_major_t_3years, MACE)~foam_tertile, data=macrosMace)) +
  labs(x="Time (years)", y="Event-free survival", color = "Foam cells tertile", fill = "Foam cells tertile") + 
  theme_classic() + theme(legend.position = "none")  +ggtitle("Foam Cells")+
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["foam_tertilehigh","Pr(>|z|)"], 3)))
## C ####
cox <- coxph(Surv(ep_major_t_3years, MACE)~res_tertile, data=macrosMace) 
F7C <- autoplot(survfit(Surv(ep_major_t_3years, MACE)~res_tertile, data=macrosMace)) +
  labs(x="Time (years)", y="Event-free survival", color = "Resident macrophages tertile", fill = "Resident macrophages tertile") + 
  theme_classic() + theme(legend.position = "none")  + ggtitle("Resident Macrophages")+
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["res_tertilehigh","Pr(>|z|)"], 3)))
## D ####
cox <- coxph(Surv(ep_major_t_3years, MACE)~inf_tertile, data=macrosMace) 
F7D <- autoplot(survfit(Surv(ep_major_t_3years, MACE)~inf_tertile, data=macrosMace)) +
  labs(x="Time (years)", y="Event-free survival", color = "Inflammatory macrophages tertile", fill = "Inflammatory macrophages tertile") + 
  theme_classic() + theme(legend.position = "none")  +ggtitle("Inflammatory Macrophages")+
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["inf_tertilehigh","Pr(>|z|)"], 3)))

## merge ####
png("Fig7.png", width = 650, height = 500)
plot_grid(plot_grid(NULL, F7A, NULL, ncol=3, rel_widths = c(0.2,1,0.2)), plot_grid(F7B, F7C, F7D, ncol=3, labels = c("B", "C", "D")), nrow = 2, labels = c("A", ""), rel_heights = c(1,0.5))
dev.off()  


# S. FIG 1 (possible gene marker)####
Idents(seuset) <- factor(numCellType, levels = sort(unique(numCellType)))
FS1 <- VlnPlot(seuset, features = "CD79A")+theme_bw()+theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"), axis.title.x =element_blank())+ylim(-4,5.3)+
  coord_cartesian(xlim = c(1,18), ylim = c(0,5), clip = "off") +
  annotate("point", y = -0.7, x = 18, size = 3, color = "#F8766D")+  
  annotate("point", y = -0.7, x = 17, size = 3, color = "#E88526")+
  annotate("point", y = -0.7, x = 16, size = 3, color = "#D39200")+
  annotate("point", y = -0.7, x = 15, size = 3, color = "#B79F00")+
  annotate("point", y = -0.7, x = 14, size = 3, color = "#93AA00")+
  annotate("point", y = -0.7, x = 13, size = 3, color = "#5EB300")+
  annotate("point", y = -0.7, x = 12, size = 3, color = "#00BA38")+
  annotate("point", y = -0.7, x = 11, size = 3, color = "#00BF74")+
  annotate("point", y = -0.7, x = 10, size = 3, color = "#00C19F")+
  annotate("point", y = -0.7, x = 9, size = 3, color = "#00BFC4")+
  annotate("point", y = -0.7, x = 8, size = 3, color = "#00B9E3")+
  annotate("point", y = -0.7, x = 7, size = 3, color = "#00ADFA")+
  annotate("point", y = -0.7, x = 6, size = 3, color = "#619CFF")+
  annotate("point", y = -0.7, x = 5, size = 3, color = "#AE87FF")+
  annotate("point", y = -0.7, x = 4, size = 3, color = "#DB72FB")+
  annotate("point", y = -0.7, x = 3, size = 3, color = "#F564E3")+
  annotate("point", y = -0.7, x = 2, size = 3, color = "#FF61C3")+
  annotate("point", y = -0.7, x = 1, size = 3, color = "#FF699C")

# take legend from UMAP with 18 cell types (F1A)
Legend_FS1 <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=12),
        panel.grid = element_blank(), panel.border = element_rect(colour = "white", fill='white', size=1), legend.key.size = unit(0.5, 'cm'))

### merge ####
png("SFig1.png", width = 600 , height = 300)
plot_grid(FS1, Legend_FS1)
dev.off()  

# S. FIG 2 (cell type similarity) ####
scCPM <- as.data.frame(t(apply(seuset@assays$RNA@counts, 2, function(x){x/sum(x)}))) %>% rownames_to_column("Sample")
seuset@meta.data$CellType <- ct
scCPM_annotated<- as_tibble(merge(scCPM, seuset@meta.data %>% rownames_to_column("Sample") %>% select(Sample, CellType), by = "Sample"))
scCPM_cellTypes <- data.frame(row.names = colnames(scCPM %>% select(-Sample)))
for(celltype in unique(scCPM_annotated$CellType)){
  print(celltype)
  scCPM_cellTypes[[celltype]] <- unname((scCPM_annotated %>% filter(CellType==celltype))[,2:20112] %>% colMeans()) #remove columns like sample, celltype...
}
FS2 <- ggcorrplot(cor(scCPM_cellTypes), hc.order = T, ggtheme = theme_bw(), 
           outline.col = "white", lab=TRUE, lab_size = 2) + 
  theme(axis.text.y = element_text(size = 7), 
        axis.text.x = element_text(size = 7)) + 
  scale_fill_gradient2(limit = c(0,1), low = "#440154", high =  "#fde725", mid = "#21918c", midpoint = 0.5)+ 
  labs(fill = "Correlation")

png("SFig2.png", width = 600 , height = 600)
FS2
dev.off() 

# S. FIG 3 (pseudo-bulk proportions) ####
## A ####
simulProps18_pseudobulk <- read_csv("dataPlots/pseudoProps18.csv")
FS3A  <- simulProps18_pseudobulk %>% select(cellType, starts_with("simProps")) %>% 
  mutate(cellType=factor(cellType, levels = c("CD3+ T Cells I","CD3+ T Cells II","CD3+ T Cells III","CD3+ T Cells IV","CD3+ T Cells V","CD3+ T Cells VI","FOXP3+ T Cells","CD3+CD56+ NK Cells I","CD3+CD56+ NK Cells II",
                                              "CD79A+ Class-switched Memory B Cells","CD68+ABCA1+OLR1+TREM2+ Foam Cells","CD68+IL18+TLR4+TREM2+ Resident macrophages","CD68+CASP1+IL1B+SELL+ Inflammatory macrophages",
                                              "CD68+CD1C+ Dendritic Cells","ACTA2+ Smooth Muscle Cells","CD34+ Endothelial Cells I","CD34+ Endothelial Cells II","CD68+KIT+ Mast Cells"))) %>% pivot_longer(-cellType, values_to = "Proportion", names_to = "Mixture") %>% 
  ggplot(aes(cellType, Proportion)) + geom_boxplot() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=40, hjust = 1))
## B ####
props18_pseudobulk <- read_csv("dataPlots/pseudoBulk_18Linear_props.csv") %>% mutate(cellType=case_when(cellType=="CD68+ABCA1+OLR1+TREM2+ Foam Cells"~"CD68+ABCA1+ OLR1+TREM2+ Foam Cells", 
                                                                                                        cellType=="CD68+CASP1+IL1B+SELL+ Inflammatory macrophages"~"CD68+CASP1+ IL1B+SELL+ Inflammatory macrophages",
                                                                                                        cellType=="CD68+IL18+TLR4+TREM2+ Resident macrophages"~"CD68+IL18+TLR4+ TREM2+ Resident macrophages",
                                                                                                        TRUE~cellType))
FS3B <- props18_pseudobulk %>% pivot_longer(c("deconvProps", "simulProps"), names_to="typeProps", values_to = "Props") %>%
  mutate(algorithm=case_when(typeProps=="simulProps"~"...simulated...", TRUE~algorithm)) %>% 
  mutate(cellType=factor(swr(cellType), levels = c("CD3+ T\nCells I","CD3+ T\nCells II","CD3+ T\nCells III","CD3+ T\nCells IV","CD3+ T\nCells V","CD3+ T\nCells VI","FOXP3+ T\nCells","CD3+CD56+\nNK Cells\nI","CD3+CD56+\nNK Cells\nII",
                                                    "CD79A+\nClass-switched\nMemory B\nCells","CD68+ABCA1+\nOLR1+TREM2+\nFoam\nCells","CD68+IL18+TLR4+\nTREM2+\nResident\nmacrophages","CD68+CASP1+\nIL1B+SELL+\nInflammatory\nmacrophages",
                                                   "CD68+CD1C+\nDendritic\nCells","ACTA2+\nSmooth\nMuscle\nCells","CD34+\nEndothelial\nCells I","CD34+\nEndothelial\nCells II","CD68+KIT+\nMast\nCells"))) %>% 
  ggplot(aes(x = algorithm, y = Props)) + geom_boxplot(aes(fill=algorithm), outlier.size = 0.2, alpha=0.7) + ylab("Proportion") +
  facet_wrap(~cellType, ncol = 9) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90)) +
  ylim(0,0.52) + scale_fill_manual(values = c("#C90606", "#8DD3C7", "#FFFFB3","#BEBADA","#FDB462","#80B1D3"))

## merge ####
png("SFig3.png", width = 800 , height = 900)
plot_grid(plot_grid(FS3A, NULL, ncol=2, rel_widths = c(2,1)), FS3B, nrow = 2, rel_heights = c(1,1.3), labels = c("A", "B"))
dev.off() 

# S. FIG 4 (from 18 cell types --> 9 clusters --> 6 clusters) ####
ct #scells IDs <-> cell type
ct_9 <- case_when(ct=="CD68+CD1C+ Dendritic Cells"~"Dendritic Cells", ct %in% c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")~"CD3+ T Cells", 
                  ct %in% c("CD3+CD56+ NK Cells I", "CD3+CD56+ NK Cells II")~"NK Cells", ct=="FOXP3+ T Cells"~"FOXP3+ T Cells", ct=="CD79A+ Class-switched Memory B Cells"~"Switched mem B Cells", 
                  ct=="ACTA2+ Smooth Muscle Cells"~"Smooth Muscle Cells", ct=="CD68+KIT+ Mast Cells"~"Mastocytes", ct %in% c("CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")~"Macrophages", 
                  ct %in% c("CD34+ Endothelial Cells I","CD34+ Endothelial Cells II")~"Endothelial Cells")
numCellType_9 <- case_when(ct_9=="CD3+ T Cells"~1,ct_9=="FOXP3+ T Cells"~2,ct_9=="NK Cells"~3,
                         ct_9=="Switched mem B Cells"~4, ct_9=="Macrophages"~5, ct_9=="Dendritic Cells"~6,
                         ct_9=="Smooth Muscle Cells"~7, ct_9=="Endothelial Cells"~8, ct_9=="Mastocytes"~9)
newLabel <- paste0(numCellType_9, ": ", ct_9)
Idents(seuset) <- numCellType_9
FS4_9 <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "right")
ct_6 <- case_when(ct %in% c("CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI","CD3+CD56+ NK Cells I", "CD3+CD56+ NK Cells II","FOXP3+ T Cells")~"T and NK Cells", 
                  ct=="CD79A+ Class-switched Memory B Cells"~"Switched mem B Cells", 
                  ct=="ACTA2+ Smooth Muscle Cells"~"Smooth Muscle Cells", ct=="CD68+KIT+ Mast Cells"~"Mastocytes", ct %in% c("CD68+CD1C+ Dendritic Cells","CD68+IL18+TLR4+TREM2+ Resident macrophages", "CD68+CASP1+IL1B+SELL+ Inflammatory macrophages", "CD68+ABCA1+OLR1+TREM2+ Foam Cells")~"CD68+ Cells", 
                  ct %in% c("CD34+ Endothelial Cells I","CD34+ Endothelial Cells II")~"Endothelial Cells")
numCellType_6 <- case_when(ct_6=="T and NK Cells"~1, ct_6=="Switched mem B Cells"~2, ct_6=="CD68+ Cells"~3,
                           ct_6=="Smooth Muscle Cells"~4, ct_6=="Endothelial Cells"~5, ct_6=="Mastocytes"~6)
newLabel <- paste0(numCellType_6, ": ", ct_6)
Idents(seuset) <- numCellType_6
FS4_6 <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "right")
Legend_F1A <- F1A+theme(legend.position = "bottom")+theme(axis.title = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                                                          legend.position = c(0.5, 0.5), legend.key = element_rect(fill='NA'), legend.title = element_text(size=12),
                                                      panel.grid = element_blank(), panel.border = element_rect(colour = "white", fill='white', size=1), legend.key.size = unit(0.5, 'cm'))+guides(colour=guide_legend(ncol=1, override.aes = list(size=3)))

png("SFig4.png", width = 800 , height = 600)
plot_grid(F1A+theme(legend.position = "none",aspect.ratio = 1), FS4_9, 
          Legend_F1A, FS4_6, ncol=2, labels = c("A","B","","C"))
dev.off()  

# S. FIG 5 (sc proportions) ####
## A ####
scPropsAndDeconvProps_8 <- read_csv("dataPlots/propsScDec_28patients.csv")
FS5A <- scPropsAndDeconvProps_8 %>% pivot_longer(c("deconvProps", "scProps"), names_to="typeProps", values_to = "Props") %>%
  mutate(Algorithm=case_when(typeProps=="scProps"~"...sc-props...", TRUE~Algorithm)) %>% 
  mutate(cellType=factor(case_when(cellType=="Switched mem B Cells"~"Switched-memory B Cells", TRUE~cellType), levels=c("T and NK Cells", "Switched-memory B Cells", "Macrophages", "Dendritic Cells", "Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes"))) %>% 
  ggplot(aes(x = Algorithm, y = Props)) + geom_boxplot(aes(fill=Algorithm), outlier.size = 0.2, alpha=0.7) + ylab("Proportion") +
  facet_wrap(~cellType, ncol = 4) + theme_bw() + theme(legend.position = "none", axis.title.x = element_blank(), strip.text = element_text(size = 6), axis.text.x = element_text(angle = 90)) +
  ylim(0,1) + scale_fill_manual(values = c("#C90606", "#8DD3C7", "#FFFFB3","#BEBADA","#FDB462","#80B1D3"))
## B ####
propsHisto_8clust <-  read_csv("dataPlots/scProps_histology.csv") %>% mutate_at(vars(10:18), as.factor)
boxPlot_histoVsProps <- function(dataset, score, cellType){
  return(ggplot(dataset %>% filter(!is.na(!!sym(score))), aes_string(x=score, y=paste0("`",cellType,"`")))+geom_boxplot()+
           theme_bw())
}
FS5B <- plot_grid(boxPlot_histoVsProps(propsHisto_8clust, "CD3.score", "T and NK Cells"),
          boxPlot_histoVsProps(propsHisto_8clust, "CD34.score", "Endothelial Cells I"),
          boxPlot_histoVsProps(propsHisto_8clust, "CD34.score", "Endothelial Cells II"),
          boxPlot_histoVsProps(propsHisto_8clust, "CD68.score", "Macrophages"),
          boxPlot_histoVsProps(propsHisto_8clust, "CD68.score", "Dendritic Cells"),
          boxPlot_histoVsProps(propsHisto_8clust, "alpha.SMA.score", "Smooth Muscle Cells"),
          ncol = 3)
## merge ####
png("SFig5.png", width = 800 , height = 900)
plot_grid(FS5A, FS5B, ncol= 1, labels = c("A","B"), rel_heights = c(1.5,1))
dev.off()

# S. FIG 6 (histology) ####
histoClinical %<>% mutate(ACTA2 = factor(case_when(smc=="no staining"~0, smc=="minor staining"~1, smc=="moderate staining"~2, smc=="heavy staining"~3))) %>% 
  mutate(CD68 = factor(case_when(macrophages=="no staining"~0, macrophages=="minor staining"~1, macrophages=="moderate staining"~2, macrophages=="heavy staining"~3))) %>% 
  mutate(Phenotype = factor(case_when(plaquephenotype=="atheromatous"~1, plaquephenotype=="fibroatheromatous"~2, plaquephenotype=="fibrous"~3))) %>% 
  mutate(`ACTA2.CD68` = factor(case_when(smc_macrophages_ratio=="macrophages dominant"~1, smc_macrophages_ratio=="equal"~2, smc_macrophages_ratio=="SMC dominant"~3))) %>% 
  mutate(Fat = factor(case_when(fat=="no fat"~0, fat=="< 40% fat"~1, fat=="> 40% fat"~2))) %>% 
  mutate(Collagen = factor(case_when(collagen=="no staining"~0, collagen=="minor staining"~1, collagen=="moderate staining"~2, collagen=="heavy staining"~3)))
analyzeCategoricalFeatures <- function(dataset, listFeatures, cellType){
  listPlots <- list()
  for(feature in listFeatures){
    if(feature==listFeatures[1]){ #print y axis (cell type)
      listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
                                 geom_boxplot())+theme_bw()+
        theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), axis.title.y = element_text(size=11), axis.title.x = element_text(size=8))
    }else{ #remove y axis 
      listPlots[[feature]] <- (ggplot(dataset %>% filter(is.na(.data[[feature]])==FALSE), aes_string(x=feature, y=paste0("`", cellType, "`")))+
                                 geom_boxplot())+theme_bw()+
        theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), axis.title.x = element_text(size=11), axis.title.y = element_blank())
      
    }
  }
  return(listPlots)
}
propsHisto <- 
plotsAlgor <- list()
for(algor in sort(unique(props8$Algorithm))){
  propsHisto <- merge(props8 %>% filter(Algorithm==algor), histoClinical, by.x = "Patient", by.y = "study_number")
  plotsAlgor[[algor]] <- (plot_grid(plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("CD68", "Fat", "ACTA2.CD68", "Phenotype"), "Macrophages"), nrow = 1), 
                                    plot_grid(plotlist = analyzeCategoricalFeatures(propsHisto, c("ACTA2", "Collagen", "ACTA2.CD68", "Phenotype"), "Smooth Muscle Cells"), nrow = 1), 
                                    nrow=1))
}
png("SFig6.png", width = 800 , height = 1000)
plot_grid(plotlist = plotsAlgor, nrow = 5, labels=c("A", "B", "C", "D", "E"))
dev.off()

# S. FIG 7 ####
## A (UMAP, macro degrouped) ####
ct <- Idents(seuset)
ct_macroDegrouped <- case_when(ct=="CD68+CD1C+ Dendritic Cells"~"Dendritic Cells", ct %in% c("FOXP3+ T Cells","CD3+CD56+ NK Cells I", "CD3+CD56+ NK Cells II","CD3+ T Cells I", "CD3+ T Cells II", "CD3+ T Cells III", "CD3+ T Cells IV", "CD3+ T Cells V", "CD3+ T Cells VI")~"T and NK Cells", 
                  ct=="CD79A+ Class-switched Memory B Cells"~"Switched-memory B Cells", 
                  ct=="ACTA2+ Smooth Muscle Cells"~"Smooth Muscle Cells", ct=="CD68+KIT+ Mast Cells"~"Mastocytes", ct=="CD68+IL18+TLR4+TREM2+ Resident macrophages"~"Resident Macrophages", ct=="CD68+CASP1+IL1B+SELL+ Inflammatory macrophages"~"Inflammatory Macrophages", 
                  ct=="CD68+ABCA1+OLR1+TREM2+ Foam Cells"~"Foam Cells", 
                  ct=="CD34+ Endothelial Cells I"~"Endothelial Cells I", ct=="CD34+ Endothelial Cells II"~"Endothelial Cells II")
numCellType_macroDegrouped <- case_when(ct_macroDegrouped=="T and NK Cells"~1,
                           ct_macroDegrouped=="Switched-memory B Cells"~2, ct_macroDegrouped=="Foam Cells"~3, ct_macroDegrouped=="Resident Macrophages"~4, 
                           ct_macroDegrouped=="Inflammatory Macrophages"~5, ct_macroDegrouped=="Dendritic Cells"~6,
                           ct_macroDegrouped=="Smooth Muscle Cells"~7, ct_macroDegrouped=="Endothelial Cells I"~8, ct_macroDegrouped=="Endothelial Cells II"~9, ct_macroDegrouped=="Mastocytes"~10)
newLabel <- paste0(numCellType_macroDegrouped, ": ", ct_macroDegrouped)
Idents(seuset) <- numCellType_macroDegrouped
FS7A <- DimPlot(seuset, reduction = "umap", label = T, label.size = 4, label.box = F, repel = T, order = sort(unique(numCellType_macroDegrouped), decreasing = T))+theme_bw()+
  scale_color_discrete(labels=unique(newLabel)[order(as.numeric(gsub(":.*", "", unique(newLabel))))])+
  theme(panel.grid.minor = element_line(colour = "black", size = 0), 
        panel.grid.major = element_line(colour = "black", size = 0), 
        axis.title.x = element_text(hjust = 0), axis.title.y = element_text(angle = -90, hjust = 1), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), aspect.ratio = 1)
## B (proportions) ####
degrProps <- read_csv("dataPlots/propsMacroUngrouped.csv") %>% mutate_if(is.numeric, ~(./100)) %>% rename(`Endothelial Cells I`=`CD34+ Endothelial Cells I`, `Foam Cells`=`Foam Macrophages`, `Switched-memory B Cells`=`Switched mem B Cells`, 
                                                                                                          `Resident Macrophages`=`Res Macrophages`, `Inflammatory Macrophages`=`Inflam Macrophages`)
FS7B <- degrProps %>% pivot_longer(-Patient, names_to = "CellType", values_to = "Proportion") %>% 
  mutate(CellType=factor(CellType, levels = c("T and NK Cells", "Switched-memory B Cells", "Foam Cells", "Resident Macrophages", "Inflammatory Macrophages", "Dendritic Cells", "Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes"))) %>% 
  ggplot(aes(x=CellType,y=Proportion)) + geom_boxplot(outlier.size = 0.2) + theme_bw() +
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "none")+ 
  ylim(0,1) + xlab("")
## merge ####
png("SFig7.png", width = 500 , height = 500)
plot_grid(FS7A, FS7B, nrow = 2, labels=c("A", "B"))
dev.off()







