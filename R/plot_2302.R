

source("/home/isin/Projects/AMOVA/AMOVA_paper_analyses/R/shared.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
#library(reshape2)
library(latex2exp)
library(akima)
library(fields)



read_process_files <- function(true_phi_fn, gc_phi_fn, gle_phi_fn) {
  # read in true phi data
  true_phi <- read.csv(true_phi_fn, header=TRUE)
  true_phi$Measure <- as.factor(true_phi$Measure)
  true_phi$Level <- as.factor(true_phi$Level)
  true_phi$Model <- as.factor(true_phi$Model)
  true_phi$Contig <- as.factor(true_phi$Contig)
  true_phi$Rep <- as.factor(true_phi$Rep)

  # read in genotype call phi data
  gc_phi <- read.csv(gc_phi_fn, header=TRUE)
  gc_phi$Measure <- as.factor(gc_phi$Measure)
  gc_phi$Level <- as.factor(gc_phi$Level)
  gc_phi$Model <- as.factor(gc_phi$Model)
  gc_phi$Contig <- as.factor(gc_phi$Contig)
  gc_phi$Rep <- as.factor(gc_phi$Rep)
  gc_phi$Depth <- as.factor(gc_phi$Depth)

  # read in GLE phi data
  gle_phi <- read.csv(gle_phi_fn, header=TRUE)
  gle_phi$Measure <- as.factor(gle_phi$Measure)
  gle_phi$Level <- as.factor(gle_phi$Level)
  gle_phi$Model <- as.factor(gle_phi$Model)
  gle_phi$Contig <- as.factor(gle_phi$Contig)
  gle_phi$Rep <- as.factor(gle_phi$Rep)
  gle_phi$Depth <- as.factor(gle_phi$Depth)

  return(list(true_phi=true_phi, gc_phi=gc_phi, gle_phi=gle_phi))

}

file_list <- read_process_files(true_phi_fn="/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_rawTruth_amova_phi.csv",
                                gc_phi_fn="/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_calledGt_amova_phi.csv",
                                gle_phi_fn="/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova1_gle_tole5_amova_phi.csv")

#sanity check
sum(file_list$true_phi$Model=="model1")
sum(file_list$true_phi$Model=="model2")
sum(file_list$gc_phi$Model=="model1")
sum(file_list$gc_phi$Model=="model2")
sum(file_list$gle_phi$Model=="model1")
sum(file_list$gle_phi$Model=="model2")







true_phi=file_list$true_phi
gc_phi=file_list$gc_phi
gle_phi=file_list$gle_phi

anti_join(true_phi,gc_phi,by=c("Measure","Level","Model","Contig","Rep"))
phi_true_gc<-full_join(true_phi,gc_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gc"), multiple="all")
check_na_col(phi_true_gc)
phi_true_gc$diff2_true_gc <- (phi_true_gc$Value_true - phi_true_gc$Value_gc)^2
phi_true_gc_rmse<-aggregate(diff2_true_gc ~ Measure + Level + Model + Contig + Depth, phi_true_gc, function(x) sqrt(mean(x)))
phi_true_gc_rmse$Method<-"GenotypeCalling"
colnames(phi_true_gc_rmse)<-c("Measure","Level","Model","Contig","Depth","RMSE","Method")
phi_true_gc_rmse$Method<-as.factor(phi_true_gc_rmse$Method)

# check if any discordance
anti_join(true_phi,gle_phi,by=c("Measure","Level","Model","Contig","Rep"))

phi_true_gle<-full_join(true_phi,gle_phi,by=c("Measure","Level","Model","Contig","Rep"),suffix=c("_true","_gle"),multiple="all")
phi_true_gle$diff2_true_gle <- (phi_true_gle$Value_true - phi_true_gle$Value_gle)^2
phi_true_gle_rmse <- aggregate(diff2_true_gle ~ Measure + Level + Model + Contig + Depth, phi_true_gle, function(x) sqrt(mean(x)))

phi_true_gle_rmse$Method<-"GenotypeLikelihood"
phi_true_gle_rmse$Method<-as.factor(phi_true_gle_rmse$Method)
colnames(phi_true_gle_rmse)<-c("Measure","Level","Model","Contig","Depth","RMSE","Method")

check_na_col(phi_true_gle_rmse)
check_na_col(phi_true_gc_rmse)

gle_gc_rmse<-rbind(phi_true_gc_rmse,phi_true_gle_rmse)


ninds<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_calledGt_avgNIndPerSite.txt"
nis<-read.csv(ninds,header=T)
nis$Model<-as.factor(nis$Model)
nis$Contig<-as.factor(nis$Contig)
nis$Rep<-as.factor(nis$Rep)
nis$Depth<-as.factor(nis$Depth)

nsites<-"/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_feb23/sim_demes_v2_doAmova2_calledGt_nSites.txt"
nsts<- read.csv(nsites,header=T)
nsts$Model<-as.factor(nsts$Model)
nsts$Contig<-as.factor(nsts$Contig)
nsts$Rep<-as.factor(nsts$Rep)
nsts$Depth<-as.factor(nsts$Depth)


# mean nsites among replicates
# this is the total number of sites included in analysis
# i.e. if one site exists only in two individuals (==1 individual pair) it counts it
# so there is not much difference between 0.1 and 0.2 since it is across 40 individuals
# dns : mean number of (across reps) sites used in analysis
#nsts %>% filter(Model=="model2") %>% group_by(Model,Contig,Depth) %>% summarize(meanNSites=floor(mean(nSites))) -> dns
nsts %>% group_by(Model,Contig,Depth) %>% summarize(meanNSites=floor(mean(nSites))) -> dns

check_na_col(dns)


# dnips: mean number of (across reps) average number of individuals per site
# nis %>% filter(Model=="model2") %>% group_by(Model,Contig,Depth) %>% summarize(meanNInds=floor(mean(avgNIndPerSite))) -> dnips
nis %>% group_by(Model,Contig,Depth) %>% summarize(meanNInds=floor(mean(avgNIndPerSite))) -> dnips


df<-merge(gle_gc_rmse, dns, by=c("Model","Contig","Depth"))
df<-merge(df, dnips, by=c("Model","Contig","Depth"))
check_na_col(df)




################################################################################



ggplot(df,aes(x=Depth, y=meanNSites,color=RMSE))+
  geom_point()

ndf<-data.frame(Depth=as.numeric(as.character(df$Depth)), meanNSites=as.numeric(as.character(df$meanNSites)),RMSE=df$RMSE )

ddf<-ndf[
  order( ndf[,"Depth"], ndf[,"meanNSites"],ndf[,"RMSE"] ),
]


x=scale(ddf$Depth)
y=scale(ddf$meanNSites)
z=ddf$RMSE

intp<-interp(ddf$Depth,y,z,duplicate = "mean")


## interpolate data
#fld <- with(ddf, interp(x = Depth, y = y, z = RMSE,duplicate = "mean"))

## prepare data in long format
library(reshape2)
#meld <- melt(fld$z, na.rm = TRUE)
#names(meld) <- c("x", "y", "Dewpoint")

# g <- ggplot(data = meld, aes(x = x, y = y, z =RMSE))  +
#   labs(x = "depth", y = "nsites",
#        color = "RMSE")
# g
# g + stat_contour(aes(color = ..level.., fill = RMSE))



length(y)
length(ndf$meanNSites)
length(ndf$Depth)


z=dfgl$RMSE

dfgl_n$Depth
dfgl_n$meanNSites

t<-interp(x,y,z,duplicate ="mean")
#t <- interp(x,y,z)

tt<-data.frame(t)


contour(dfgl_n$Depth,dfgl_n$meanNSites,t(dfgl$RMSE),add=T)
contour(dfgl_n$Depth,dfgl_n$meanNSites,dfgl$RMSE,add=T)

contour(intp$x,intp$y,intp$z)

contour(intp,nlevels = 10)

filled.contour(t,nlevels = 5,zlim = c(0,0.01))

filled.contour(t,nlevels = 5)
filled.contour(t,levels=c(min(unlist(t$z),na.rm = T),5e-4,6e-4,1e-3,5e-3,1e-2,2e-2,0.1,max(unlist(t$z),na.rm = T)))



filled.contour(t,nlevels = 5,zlim = c(0,0.01))
t2<-as.data.frame(interp2xyz(t))
ggplot(t2, aes(x=x, y=y, z=z)) + geom_contour()
ggplot(t2, aes(x=x, y=y, z=z)) + geom_contour_filled(binwidth = 0.005)





unscale <- function(scaled_data){
  scale(scale(scaled_data,center=FALSE,scale=1/attr(scaled_data,'scaled:scale')),
        center=-attr(scaled_data,'scaled:center'),scale=FALSE)
}


