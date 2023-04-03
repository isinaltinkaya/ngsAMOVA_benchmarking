############################################################################################################
# read_data.R
# Functions to read in data from genotype calling, GLE, and true phi
#
# isinaltinkaya@gmail.com
# 230307
#
############################################################################################################
# LOAD LIBRARIES
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(latex2exp)
library(scales)
library(viridis)
############################################################################################################
# FUNCTIONS

## genotype calling
read_gc_phi <- function(gc_phi_fn){
  gc_phi <- read.csv(gc_phi_fn, header=TRUE)
  gc_phi$Measure <- as.factor(gc_phi$Measure)
  gc_phi$Level <- as.factor(gc_phi$Level)
  gc_phi$Model <- as.factor(gc_phi$Model)
  gc_phi$Contig <- as.factor(gc_phi$Contig)
  gc_phi$Rep <- as.factor(gc_phi$Rep)
  gc_phi$Depth <- as.factor(gc_phi$Depth)
  gc_phi$AnalysisType <- as.factor(gc_phi$AnalysisType)
  gc_phi$doPost <- as.factor(gc_phi$doPost)
  gc_phi$postCutoff <- as.factor(gc_phi$postCutoff)
  gc_phi$Method <- as.factor(paste0("GTC_", gc_phi$doPost, "_", gc_phi$postCutoff))
  gc_phi[gc_phi$Value=="NoResult",]$Value<-NA
  gc_phi$Value<-as.double(gc_phi$Value)
  gc_phi$Level <- str_replace(gc_phi$Level, "Region_in_Total", "Phi_CT")
  gc_phi$Level <- str_replace(gc_phi$Level, "Population_in_Region", "Phi_SC")
  gc_phi$Level <- str_replace(gc_phi$Level, "Population_in_Total", "Phi_ST")
  gc_phi$Level<-as.factor(gc_phi$Level)
  return(gc_phi)
}


## true phi (from msprime)
read_true_phi <- function(true_phi_fn){
  true_phi <- read.csv(true_phi_fn, header=TRUE)
  true_phi$Measure <- as.factor(true_phi$Measure)
  true_phi$Level <- as.factor(true_phi$Level)
  true_phi$Model <- as.factor(true_phi$Model)
  true_phi$Contig <- as.factor(true_phi$Contig)
  true_phi$Rep <- as.factor(true_phi$Rep)
  true_phi$AnalysisType <-  as.factor(true_phi$AnalysisType)
  true_phi$Level <- str_replace(true_phi$Level, "Region_in_Total", "Phi_CT")
  true_phi$Level <- str_replace(true_phi$Level, "Population_in_Region", "Phi_SC")
  true_phi$Level <- str_replace(true_phi$Level, "Population_in_Total", "Phi_ST")
  true_phi$Level<-as.factor(true_phi$Level)
  return(true_phi)
}

## genotype likelihood estimation
read_gle_phi <- function(gle_phi_fn){
  gle_phi <- read.csv(gle_phi_fn, header=TRUE)
  gle_phi$Measure <- as.factor(gle_phi$Measure)
  gle_phi$Level <- as.factor(gle_phi$Level)
  gle_phi$Model <- as.factor(gle_phi$Model)
  gle_phi$Contig <- as.factor(gle_phi$Contig)
  gle_phi$Rep <- as.factor(gle_phi$Rep)
  gle_phi$Depth <- as.factor(gle_phi$Depth)
  gle_phi$AnalysisType <- as.factor(gle_phi$AnalysisType)
  gle_phi$Method <- as.factor(paste0("GLE_", gle_phi$Tole))
  gle_phi$Tole<-as.factor(gle_phi$Tole)
  gle_phi$Level <- str_replace(gle_phi$Level, "Region_in_Total", "Phi_CT")
  gle_phi$Level <- str_replace(gle_phi$Level, "Population_in_Region", "Phi_SC")
  gle_phi$Level <- str_replace(gle_phi$Level, "Population_in_Total", "Phi_ST")
  gle_phi$Level<-as.factor(gle_phi$Level)
  return(gle_phi)
}


join_true_gc <- function(true_phi, gc_phi){
  phi_true_gc<-full_join(true_phi,gc_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gc"), multiple="all")
  phi_true_gc$diff2_true <- (phi_true_gc$Value_true - phi_true_gc$Value_gc)^2
  phi_true_gc_rmse<-aggregate(diff2_true ~ Measure + Level + Model + Contig + Depth + Method, phi_true_gc, function(x) sqrt(mean(x)))
  return(phi_true_gc_rmse)
}

join_true_gle <- function(true_phi, gle_phi){
  phi_true_gle<-full_join(true_phi,gle_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gle"), multiple="all")
  phi_true_gle$diff2_true <- (phi_true_gle$Value_true - phi_true_gle$Value_gle)^2
  phi_true_gle_rmse<-aggregate(diff2_true ~ Measure + Level + Model + Contig + Depth + Method, phi_true_gle, function(x) sqrt(mean(x)))
  return(phi_true_gle_rmse)
}

get_dt <- function(phi_true_gc_rmse, phi_true_gle_rmse){
  dt<-rbind(phi_true_gc_rmse,phi_true_gle_rmse)
  dt<-dt%>%rename("RMSE"="diff2_true")
  return(dt)
}

# number of sites used in analysis
get_nSites_used <- function(nsu_gc_fn){
  nsu<-read.csv(nsu_gc_fn,header=TRUE)
  dns<-nsu %>% group_by(Model,Contig,Depth) %>% summarize(meanNSites=floor(mean(mean_nSites_used)))
  return(dns)
}

get_df <- function(true_phi_fn, gc_phi_fn, gle_phi_fn, nsu_gc_fn){
  true_phi<-read_true_phi(true_phi_fn)
  gle_phi<-read_gle_phi(gle_phi_fn)
  gc_phi<-read_gc_phi(gc_phi_fn)
  phi_true_gc_rmse<-join_true_gc(true_phi, gc_phi)
  phi_true_gle_rmse<-join_true_gle(true_phi, gle_phi)
  dt<-get_dt(phi_true_gc_rmse, phi_true_gle_rmse)
  dns<-get_nSites_used(nsu_gc_fn)
  df<-merge(dt, dns, by=c("Model","Contig","Depth"))
  ### convert names to phi notations
  #reg in tot =phi_ct
  #pop in reg =phi_sc
  #pop in tot =phi_st
  df$Level <- str_replace(df$Level, "Region_in_Total", "Phi_CT")
  df$Level <- str_replace(df$Level, "Population_in_Region", "Phi_SC")
  df$Level <- str_replace(df$Level, "Population_in_Total", "Phi_ST")
  df$Level<-as.factor(df$Level)
  return(df)
}

get_df2 <- function(true_phi, gc_phi, gle_phi, nsites_used){
  phi_true_gc_rmse<-join_true_gc(true_phi, gc_phi)
  phi_true_gle_rmse<-join_true_gle(true_phi, gle_phi)
  dt<-get_dt(phi_true_gc_rmse, phi_true_gle_rmse)
  df<-merge(dt, nsites_used, by=c("Model","Contig","Depth"))
  ### convert names to phi notations
  #reg in tot =phi_ct
  #pop in reg =phi_sc
  #pop in tot =phi_st
  df$Level <- str_replace(df$Level, "Region_in_Total", "Phi_CT")
  df$Level <- str_replace(df$Level, "Population_in_Region", "Phi_SC")
  df$Level <- str_replace(df$Level, "Population_in_Total", "Phi_ST")
  df$Level<-as.factor(df$Level)
  return(df)
}
