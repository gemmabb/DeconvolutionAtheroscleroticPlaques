###################################################################################
#                                                                                 #
# CLINICAL ASSOCIATIONS                                                           #                                       
#                                                                                 #
# Gemma Bel Bordes, Nov 20th 2022                                                 #
#                                                                                 #
###################################################################################

library(rstudioapi)
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggcorrplot)
library(tableone)
library(survival)
library(ggfortify)
library(survminer)
library(questionr)
#R version 4.1.2

setwd(dirname(getActiveDocumentContext()$path)) #set wd to the directory of this code file

# DATASET ####
### Scaden results (deconvolution performed 10 times using 10 different fixed rnd seeds [final predictions: mean]):
folderPredictFiles <- "output/ourSC/experimentSimulation/bothSex/"
list.files(folderPredictFiles, "both.*_predictions.txt")
propsRndSeeds <- list()
for(f in list.files(folderPredictFiles, "bothSex.*_predictions")){
  propsRndSeeds[[f]] <- read_delim(paste0(folderPredictFiles, f)) %>% rename(Patient = "...1") %>% 
    mutate(RndSeed=str_match(f, "bothSex(.*)_scaden")[2]) #str_match --> get seed from file name
}
propsRndSeeds <- bind_rows(propsRndSeeds) %>% mutate(RndSeed=as.character(RndSeed)) #all the predictions within a single dataframe
propsRndSeeds_mean <- propsRndSeeds %>% group_by(Patient) %>% summarise_if(is.numeric, list(mean)) %>% mutate_if(is.numeric, ~(.*100)) #get mean for the 10 rnd seeds
propsRndSeeds_mean %>% write_csv("dataPlots/props8scaden_10runs.csv")

### Clinical and histological annotation:
histoClinical <- read_delim("input/clinical_data_good.txt") %>% 
  mutate(symptoms_cat = case_when(symptoms_4g=="asymptomatic" ~ "asymptomatic", symptoms_4g=="ocular" ~ "asymptomatic", symptoms_4g=="TIA" ~ "symptomatic", symptoms_4g=="stroke"~"symptomatic")) %>% 
  mutate(symptoms = case_when(symptoms_4g=="asymptomatic" ~ 0, symptoms_4g=="ocular" ~ 0, symptoms_4g=="TIA" ~ 1, symptoms_4g=="stroke"~1)) %>% 
  mutate(MACE_cat = case_when(epmajor_3years=="Included" ~ "No MACE", epmajor_3years=="Excluded" ~ "MACE")) %>% 
  mutate(MACE = case_when(epmajor_3years=="Included" ~ 0, epmajor_3years=="Excluded" ~ 1)) #654 patients

### Combination of clinical/histological annotation with deconvolved proportions (from 656 to 654 patients!):
propsRndSeeds_mean <- merge(propsRndSeeds_mean, histoClinical, by.x = "Patient", by.y = "study_number") %>% rename(Sex=sex) #this is the dataset to be used for the plots
# Univariate proportions - baseline ####
propsRndSeeds_mean %<>% mutate(`SMC:Macro`=`Smooth Muscle Cells`/Macrophages) %>% 
  mutate(totalStruct= propsRndSeeds_mean %>% select(c("Endothelial Cells I", "Smooth Muscle Cells", "Endothelial Cells II")) %>% rowSums()) %>% 
  mutate(totalImmune= propsRndSeeds_mean %>% select(c("Macrophages", "Dendritic Cells", "Mastocytes", "Switched mem B Cells", "T and NK Cells")) %>% rowSums()) %>% 
  mutate(`Struct:Immune` = totalStruct/totalImmune) 
propsRndSeeds_mean %<>% mutate(Sex=factor(Sex, levels=c("male", "female")))
baseline <- c("Sex", "smokercurrent", "dm_composite", "hypertension_composite", "symptoms_cat", "age", 
              "bmi", "tc_final", "tg_final", "ldl_final", "hdl_final", "cad_history", "MACE_cat")
cells <- c(colnames(propsRndSeeds_mean)[2:9], "SMC:Macro", "Struct:Immune")
univariate <- data.frame(row.names = cells)
for(cell in cells){
  for(f in baseline){
    coefficientsLm <- summary(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), propsRndSeeds_mean))$coefficients
    confInterv <- confint(lm(as.formula(paste("`", cell, "`", "~", f, sep="")), propsRndSeeds_mean))
    for(categ in rownames(coefficientsLm)[-1]){
      univariate[cell, paste0(categ, ":pvalue")] <- round(coefficientsLm[categ, "Pr(>|t|)"], 4)
      univariate[cell, paste0(categ, ":beta [95% CI]")] <- paste0(round(coefficientsLm[categ, "Estimate"], 4), "[", paste(round(confInterv[categ,], 4), collapse = ","),"]")
    }
  }
}
View(t(univariate))

# Sex differences ####
baseline_bySex <- CreateTableOne(data=propsRndSeeds_mean, vars=baseline, strata="Sex")
write.csv(print(baseline_bySex, formatOptions = list(big.mark = ","), nonnormal = cells, contDigits = 3), 
          file = "dataPlots/baselineBySex.csv")
## T-test sex differences
boxPlot_AllCellType_2groups <- function(dataset, cells, feature, removeNonSign = T, order=unique(cells)){
  require(rstatix)
  require(ggpubr)
  data <- dataset %>% select(cells, "Patient", feature) %>% pivot_longer(-c("Patient", feature), names_to = "CellType", values_to = "Proportion") %>% 
    mutate(CellType=factor(CellType, levels=order)) 
  stat_pvalues <- data %>% group_by(CellType) %>% pairwise_wilcox_test(as.formula(paste0("Proportion~", feature))) %>% 
    add_xy_position("CellType") %>% mutate(p=paste0("p = ", format.pval(p, eps=.001, digits=3, nsmall=3)))
  if(removeNonSign) {stat_pvalues %<>% filter(p.adj<0.1)}
  data %>% ggplot(aes(x=CellType,y=Proportion)) +
    geom_boxplot(aes_string(fill=feature), outlier.size = 0.2)  + theme_bw() + scale_fill_manual(values = c("#5D3A9B", "#E66001")) +
    theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.justification=c(1,1),
          legend.position="bottom") + 
    xlab("") + stat_pvalue_manual(stat_pvalues, label = "p", size = 3.4) 
}
boxPlot_AllCellType_2groups(propsRndSeeds_mean %>% rename(`Switched-memory\nB Cells`=`Switched mem B Cells`), colnames(propsRndSeeds_mean)[2:9], "Sex", removeNonSign = F, 
                            order = c("T and NK Cells", "Switched-memory\nB Cells", "Macrophages", "Dendritic Cells", "Smooth Muscle Cells", "Endothelial Cells I", "Endothelial Cells II", "Mastocytes")) + ylim(0,1)

propsRndSeeds_mean %<>% mutate(`SMC:Macro`=`Smooth Muscle Cells`/Macrophages) %>% 
  mutate(totalStruct= propsRndSeeds_mean %>% select(c("Endothelial Cells I", "Smooth Muscle Cells", "Endothelial Cells II")) %>% rowSums()) %>% 
  mutate(totalImmune= propsRndSeeds_mean %>% select(c("Macrophages", "Dendritic Cells", "Mastocytes", "Switched mem B Cells", "T and NK Cells")) %>% rowSums()) %>% 
  mutate(`Struct:Immune` = totalStruct/totalImmune) 
boxPlot_AllCellType_2groups(propsRndSeeds_mean, c("Struct:Immune"), "Sex", removeNonSign = F) + ylim(0,3) +
  ylab("Structural/Immune ratio") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12))
boxPlot_AllCellType_2groups(propsRndSeeds_mean, c("SMC:Macro"), "Sex", removeNonSign = F) + ylim(0,8.2) +
  ylab("SMC/Macrophages ratio") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y = element_text(size = 12))

# Symptoms ####
baseline_bySymptoms <- CreateTableOne(data=propsRndSeeds_mean, vars=baseline, strata="symptoms")
write.csv(print(baseline_bySymptoms, formatOptions = list(big.mark = ","), nonnormal = cells, contDigits = 3), 
          file = "dataPlots/baselineBySymptoms.csv")
odds.ratio((glm(symptoms~`Macrophages`, propsRndSeeds_mean, family = "binomial")), level = 0.95)
odds.ratio(glm(symptoms~Macrophages+age+bmi, propsRndSeeds_mean, family = "binomial"), 0.95)

propsRndSeeds_mean %>% filter(!is.na(symptoms), !is.na(MACE)) %>% select(Macrophages, symptoms, MACE) %>% 
  mutate(symptoms=case_when(symptoms==1~"present", symptoms==0~"absent"), MACE=case_when(MACE==1~"present", MACE==0~"absent")) %>% 
  pivot_longer(-Macrophages, names_to = "Feature", values_to = "Presence") %>% 
  ggplot(aes(Feature, Macrophages, fill=Presence)) + geom_boxplot() + stat_compare_means(aes(group=Presence), 
                                                                                         label = "p.signif") +
  theme_classic() + ylab("Macrophage content") +
  theme(legend.title = element_blank(), axis.title.x = element_blank())

# MACE ####
baseline_byMACE <- CreateTableOne(data=propsRndSeeds_mean, vars=baseline, strata="MACE")
write.csv(print(baseline_byMACE, formatOptions = list(big.mark = ","), nonnormal = cells, contDigits = 3), 
          file = "dataPlots/baselineByMACE.csv")
### Macrophages as continuous variable ####
#uni:
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophages`, data=propsRndSeeds_mean)) 
#multi:
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophages`+age+hypertension_composite+dm_composite+hdl_final, data=propsRndSeeds_mean))

### Macrophages into tertiles ####
propsRndSeeds_mean %<>% mutate(Macrophages_tertile=as.factor(ntile(`Macrophages`, 3))) %>% 
  mutate(Macrophages_tertile = factor(case_when(Macrophages_tertile==1~"low", Macrophages_tertile==2~"intermediate", Macrophages_tertile==3~"high"), levels = c("low", "intermediate", "high")))
quantile(propsRndSeeds_mean$Macrophages, probs = seq(0, 1, 1/3))
table(propsRndSeeds_mean$Macrophages_tertile)
#uni:
cox <- coxph(Surv(ep_major_t_3years, MACE)~Macrophages_tertile, data=propsRndSeeds_mean)
summary(cox)
autoplot(survfit(Surv(ep_major_t_3years, MACE)~Macrophages_tertile, data=propsRndSeeds_mean)) +
  labs(x="Time (years)", y="Event-free survival", color = "Tertile", fill = "Tertile") + 
  theme_classic() + theme(legend.position = c(0.18,0.2))  +
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["Macrophages_tertilehigh","Pr(>|z|)"], 3)))
#multi:
cox <- coxph(Surv(ep_major_t_3years, MACE)~Macrophages_tertile+age+hypertension_composite+dm_composite+hdl_final, data=propsRndSeeds_mean)
summary(cox)
autoplot(survfit(Surv(ep_major_t_3years, MACE)~Macrophages_tertile, data=propsRndSeeds_mean)) +
  labs(x="Time (years)", y="Event-free survival", color = "Tertile", fill = "Tertile") + 
  theme_classic() + theme(legend.position = c(0.18,0.2)) +
  annotate("text", x = 1.5, y = 1, label = paste0("p = ", round(summary(cox)$coefficients["Macrophages_tertilehigh","Pr(>|z|)"], 3)))

# UNGROUPED MACRO REFERENCE ####
folderPredictFiles <- "output/ourSC/byPatient_Scaden/scAllGenes/cellType/"
list.files(folderPredictFiles, "*_predictions")
ungroupedProps <- list()
for(f in list.files(folderPredictFiles, "*_predictions")){
  ungroupedProps[[f]] <- read_delim(paste0(folderPredictFiles, f)) %>% dplyr::rename(Patient = "...1") %>% 
    mutate(RndSeed=str_match(f, "seed(.*)_scaden")[2]) #str_match --> get random seed from file name
}
ungroupedProps <- bind_rows(ungroupedProps) %>% mutate(RndSeed=as.character(RndSeed)) #all the predictions within a single dataframe
ungroupedProps <- (ungroupedProps %>% group_by(Patient) %>% summarise_if(is.numeric, list(mean))) %>% mutate_if(is.numeric, ~(.*100)) %>% 
  rename(`Foam Macrophages`=`CD68+ABCA1+OLR1+TREM2+ Foam Cells`, `Res Macrophages`=`CD68+IL18+TLR4+TREM2+ Resident macrophages`, `Inflam Macrophages`=`CD68+CASP1+IL1B+SELL+ Inflammatory macrophages`)
ungroupedProps %>% write_csv("dataPlots/propsMacroUngrouped.csv")
ungroupedProps <- merge(ungroupedProps, histoClinical, by.x = "Patient", by.y = "study_number") %>% rename(Sex=sex) #this is the dataset to be used for the plots
## symptoms ####
odds.ratio(glm(symptoms~`Foam Macrophages`, ungroupedProps, family = "binomial"), 0.95)
odds.ratio(glm(symptoms~`Res Macrophages`, ungroupedProps, family = "binomial"), 0.95)
odds.ratio(glm(symptoms~`Inflam Macrophages`, ungroupedProps, family = "binomial"), 0.95)

odds.ratio(glm(symptoms~`Foam Macrophages`+age+bmi, ungroupedProps, family = "binomial"), 0.95)
odds.ratio(glm(symptoms~`Res Macrophages`+age+bmi, ungroupedProps, family = "binomial"), 0.95)
odds.ratio(glm(symptoms~`Inflam Macrophages`+age+bmi, ungroupedProps, family = "binomial"), 0.95)

## MACE  ####
ungroupedProps %<>% mutate(foam_tertile=as.factor(ntile(`Foam Macrophages`, 3))) %>% 
  mutate(foam_tertile = factor(case_when(foam_tertile==1~"low", foam_tertile==2~"intermediate", foam_tertile==3~"high"), levels = c("low", "intermediate", "high"))) %>% 
  mutate(res_tertile=as.factor(ntile(`Res Macrophages`, 3))) %>% 
  mutate(res_tertile = factor(case_when(res_tertile==1~"low", res_tertile==2~"intermediate", res_tertile==3~"high"), levels = c("low", "intermediate", "high"))) %>% 
  mutate(inf_tertile=as.factor(ntile(`Inflam Macrophages`, 3))) %>% 
  mutate(inf_tertile = factor(case_when(inf_tertile==1~"low", inf_tertile==2~"intermediate", inf_tertile==3~"high"), levels = c("low", "intermediate", "high"))) 
quantile(ungroupedProps$`Foam Macrophages`, probs = seq(0, 1, 1/3))
quantile(ungroupedProps$`Res Macrophages`, probs = seq(0, 1, 1/3))
quantile(ungroupedProps$`Inflam Macrophages`, probs = seq(0, 1, 1/3))

summary(coxph(Surv(ep_major_t_3years, MACE)~`Foam Macrophages`, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Res Macrophages`, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Inflam Macrophages`, data=ungroupedProps)) 

summary(coxph(Surv(ep_major_t_3years, MACE)~foam_tertile, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~res_tertile, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~inf_tertile, data=ungroupedProps)) 

summary(coxph(Surv(ep_major_t_3years, MACE)~`Foam Macrophages`+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Res Macrophages`+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Inflam Macrophages`+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 

summary(coxph(Surv(ep_major_t_3years, MACE)~foam_tertile+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~res_tertile+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~inf_tertile+age+hypertension_composite+dm_composite+hdl_final, data=ungroupedProps)) 

# HUIZE PAN REFERENCE ####
folderPredictFiles <- "results/scaden/ourBulk_hpSC/allCelltypes/"
list.files(folderPredictFiles, "*_predictions")
hpProps <- list()
for(f in list.files(folderPredictFiles, "*_predictions")){
  hpProps[[f]] <- read_delim(paste0(folderPredictFiles, f)) %>% dplyr::rename(Patient = "...1") %>% 
    mutate(RndSeed=str_match(f, "seed(.*)_scaden")[2]) #str_match --> get random seed from file name
}
hpProps <- bind_rows(hpProps) %>% mutate(RndSeed=as.character(RndSeed)) #all the predictions within a single dataframe
hpProps <- (hpProps %>% group_by(Patient) %>% summarise_if(is.numeric, list(mean))) %>% mutate_if(is.numeric, ~(.*100))
hpProps <- merge(hpProps, histoClinical, by.x = "Patient", by.y = "study_number") %>% rename(Sex=sex) #this is the dataset to be used for the plots
## symptoms ####
odds.ratio(glm(symptoms~`Macrophage 1`, hpProps, family = "binomial"), .95)
odds.ratio(glm(symptoms~`Macrophage 2`, hpProps, family = "binomial"), .95)
odds.ratio(glm(symptoms~`Macrophage 3`, hpProps, family = "binomial"), .95)

odds.ratio(glm(symptoms~`Macrophage 1`+age+bmi, hpProps, family = "binomial"), .95)
odds.ratio(glm(symptoms~`Macrophage 2`+age+bmi, hpProps, family = "binomial"), .95)
odds.ratio(glm(symptoms~`Macrophage 3`+age+bmi, hpProps, family = "binomial"), .95)

## MACE  ####
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 1`, data=hpProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 2`, data=hpProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 3`, data=hpProps)) 

summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 1`+age+hypertension_composite+dm_composite+hdl_final, data=hpProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 2`+age+hypertension_composite+dm_composite+hdl_final, data=hpProps)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~`Macrophage 3`+age+hypertension_composite+dm_composite+hdl_final, data=hpProps)) 

## for plots in the report:
merge(propsRndSeeds_mean %>% select(Patient, Macrophages_tertile, MACE, ep_major_t_3years), 
      ungroupedProps %>% select(Patient, foam_tertile, res_tertile, inf_tertile), 
      by = "Patient") %>% write_csv("dataPlots/macrosMace.csv")

# Added value in deconvolution? ####
bulk <- exprs(readRDS("output/ourBulk/onlyOriginalSamples/linearBulk_eset.rds"))
norm_bulkMarkers <- as.data.frame(t(apply(bulk, 2,  function(x) 1000000*x/sum(x)))) %>% 
  rownames_to_column("study_number") %>% merge(histoClinical, "study_number")

#In general for macrophages:
odds.ratio(glm(symptoms~CD68, norm_bulkMarkers, family = "binomial"), .95) 
odds.ratio(glm(symptoms~CD14, norm_bulkMarkers, family = "binomial"), .95) 

summary(coxph(Surv(ep_major_t_3years, MACE)~CD68, data=norm_bulkMarkers)) 
summary(coxph(Surv(ep_major_t_3years, MACE)~CD14, data=norm_bulkMarkers))

#For foam cells:
odds.ratio(glm(symptoms~ABCA1, norm_bulkMarkers, family = "binomial"), .95) 
odds.ratio(glm(symptoms~ABCA1+age+bmi, norm_bulkMarkers, family = "binomial"))
summary(coxph(Surv(ep_major_t_3years, MACE)~ABCA1, data=norm_bulkMarkers))
summary(coxph(Surv(ep_major_t_3years, MACE)~ABCA1+age+hypertension_composite+dm_composite+hdl_final, data=norm_bulkMarkers)) 

#For inflam macrophages:
odds.ratio(glm(symptoms~IL1B, norm_bulkMarkers, family = "binomial"), .95) #higher in symptomatic, p=0.2
odds.ratio(glm(symptoms~IL1B+age+bmi, norm_bulkMarkers, family = "binomial")) #higher in symptomatic, p=0.2
summary(coxph(Surv(ep_major_t_3years, MACE)~IL1B, data=norm_bulkMarkers))
summary(coxph(Surv(ep_major_t_3years, MACE)~IL1B+age+hypertension_composite+dm_composite+hdl_final, data=norm_bulkMarkers)) 

#For resident macrophages:
odds.ratio(glm(symptoms~LYVE1, norm_bulkMarkers, family = "binomial"), .95) #higher in symptomatic, p=0.2
odds.ratio(glm(symptoms~LYVE1+age+bmi, norm_bulkMarkers, family = "binomial")) #higher in symptomatic, p=0.2
summary(coxph(Surv(ep_major_t_3years, MACE)~LYVE1, data=norm_bulkMarkers))
summary(coxph(Surv(ep_major_t_3years, MACE)~LYVE1+age+hypertension_composite+dm_composite+hdl_final, data=norm_bulkMarkers)) 

# associations with sex
summary(lm(ACTA2~sex, norm_bulkMarkers %>% mutate(sex=factor(sex, levels = c("male", "female"))))) 
confint(lm(ACTA2~sex, norm_bulkMarkers %>% mutate(sex=factor(sex, levels = c("male", "female")))))

summary(lm(CD79A~sex, norm_bulkMarkers %>% mutate(sex=factor(sex, levels = c("male", "female"))))) 
confint(lm(CD79A~sex, norm_bulkMarkers %>% mutate(sex=factor(sex, levels = c("male", "female")))))

hist(norm_bulkMarkers$ACTA2)
odds.ratio(glm(symptoms~macrophages, norm_bulkMarkers, family = "binomial"), .95) 
summary(coxph(Surv(ep_major_t_3years, MACE)~macrophages+age+hypertension_composite+dm_composite+hdl_final, data=norm_bulkMarkers)) 

