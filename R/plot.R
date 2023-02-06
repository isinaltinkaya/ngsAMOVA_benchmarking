################################################################################
# PLOT
################################################################################


################################################################################
# LOAD LIBRARIES
################################################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(latex2exp)


################################################################################
# 230205



################################################################################
# MODEL: Model 1

################################################################################
## Phi statistic

# raw simulated genotypes
true_phi_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_rawTruth_amova_phi.csv"
true_phi<-read.csv(true_phi_fn,header=T)
true_phi %>% filter(Model=="model1") -> true_phi
true_phi %>% filter(Contig %in% c(1,100)) -> true_phi
true_phi$Measure <- as.factor(true_phi$Measure)
true_phi$Level <- as.factor(true_phi$Level)
true_phi$Model <- as.factor(true_phi$Model)
true_phi$Contig <- as.factor(true_phi$Contig)
true_phi$Rep <- as.factor(true_phi$Rep)




# genotype call
gc_phi_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_calledGt_amova_phi.csv"
gc_phi<-read.csv(gc_phi_fn,header=T)
gc_phi %>% filter(Model=="model1") -> gc_phi
gc_phi %>% filter(Contig %in% c(1,100)) -> gc_phi
gc_phi$Measure <- as.factor(gc_phi$Measure)
gc_phi$Level <- as.factor(gc_phi$Level)
gc_phi$Model <- as.factor(gc_phi$Model)
gc_phi$Contig <- as.factor(gc_phi$Contig)
gc_phi$Rep <- as.factor(gc_phi$Rep)
gc_phi$Depth <- as.factor(gc_phi$Depth)







# check if there is any missing data
# anti_join(true_phi,gc_phi,by=c("Measure","Level","Model","Contig","Rep"))

gle_phi_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_gle_tole5_amova_phi.csv"
gle_phi<-read.csv(gle_phi_fn,header=T)
gle_phi %>% filter(Model=="model1") -> gle_phi
gle_phi$Measure <- as.factor(gle_phi$Measure)
gle_phi$Level <- as.factor(gle_phi$Level)
gle_phi$Model <- as.factor(gle_phi$Model)
gle_phi$Contig <- as.factor(gle_phi$Contig)
gle_phi$Rep <- as.factor(gle_phi$Rep)
gle_phi$Depth <- as.factor(gle_phi$Depth)










phi_true_gc<-full_join(true_phi,gc_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gc"))
phi_true_gc$diff2_true_gc <- (phi_true_gc$Value_true - phi_true_gc$Value_gc)^2
phi_true_gc %>% group_by(Measure,Level,Model,Contig, Depth) %>% summarize(RMSE=sqrt(mean(diff2_true_gc))) -> phi_true_gc_rmse


phi_true_gle<-full_join(true_phi,gle_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gle"))
phi_true_gle$diff2_true_gle <- (phi_true_gle$Value_true - phi_true_gle$Value_gle)^2
phi_true_gle %>% group_by(Measure,Level,Model,Contig,Depth) %>% summarize(RMSE=sqrt(mean(diff2_true_gle))) -> phi_true_gle_rmse



gle_gc_rmse<-rbind(phi_true_gc_rmse,phi_true_gle_rmse)


library(ggplot2)



















################################################################################
## Variance components

# raw simulated genotypes
raw_vc_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_rawTruth_amova_varianceComponent.csv"
raw_vc<-read.csv(raw_vc_fn,header=T)
raw_vc %>% filter(Model=="model1") -> raw_vc

# genotype call
gc_vc_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_calledGt_amova_varianceComponent.csv"
gc_vc<-read.csv(gc_vc_fn,header=T)
gc_vc %>% filter(Model=="model1") -> gc_vc

# genotype likelihood estimation


gle_vc_fn<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_gle_tole5_amova_varianceComponent.csv"
gle_vc<-read.csv(gle_vc_fn,header=T)
gle_vc %>% filter(Model=="model1") -> gle_vc
gle_vc$Measure <- as.factor(gle_vc$Measure)
gle_vc$Level <- as.factor(gle_vc$Level)
gle_vc$Model <- as.factor(gle_vc$Model)
gle_vc$Contig <- as.factor(gle_vc$Contig)
gle_vc$Rep <- as.factor(gle_vc$Rep)
gle_vc$Depth <- as.factor(gle_vc$Depth)

# HEATMAP

# Coefficicient of variation (CV)
