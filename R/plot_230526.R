#230526
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(latex2exp)
library(viridis)
library(ggpubr)
library(gridExtra)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
#install.packages(c("ggplot2","dplyr","tidyr","stringr","reshape2","latex2exp","viridis","ggpubr"))

resdir<-"/home/isin/Mount/AMOVA/simulations/msprime/"

fn_gli<-"simulations/sim_demes_v2/results/collect_results_genotype_likelihood_inferredMajorMinor.csv"
fn_glf<-"simulations/sim_demes_v2/results/collect_results_genotype_likelihood_fixedMajorMinor.csv"
fn_truth<-"simulations/sim_demes_v2/results/collect_results_truth.csv"
fn_gci<-"simulations/sim_demes_v2/results/collect_results_genotype_calling_inferredMajorMinor.csv"
fn_gcf<-"simulations/sim_demes_v2/results/collect_results_genotype_calling_fixedMajorMinor.csv"

gli<-read.csv(paste0(resdir,fn_gli),header=T,sep=",")
glf<-read.csv(paste0(resdir,fn_glf),header=T,sep=",")
truth<-read.csv(paste0(resdir,fn_truth),header=T,sep=",")
truth<-truth[c("AnalysisType","Model","Contig","Rep","Measure","Level","Value")]
gci<-read.csv(paste0(resdir,fn_gci),header=T,sep=",")
gcf<-read.csv(paste0(resdir,fn_gcf),header=T,sep=",")

gl<-rbind(gli,glf)
gl[is.infinite(gl$Value),]$Value<-NA
rm(gli,glf)
gc()

gc<-rbind(gci,gcf)
gc[is.infinite(gc$Value),]$Value<-NA
rm(gci,gcf)
gc()
############################################################################################################



# > colnames(gl)
# [1] "AnalysisType" "Model"        "Contig"       "Rep"          "Depth"       
# [6] "Measure"      "Level"        "Value"     
gld<-merge(gl,truth, by=c("Model","Contig","Rep","Measure","Level"),suffixes=c("_gl","_truth"))
gld<-gld[!is.na(gld$Value_gl),]
if(0 != sum(is.na(gld$Value_truth))) stop("Truth has NA values")


gld$diff2<-(gld$Value_gl-gld$Value_truth)^2


gcd<-merge(gc,truth, by=c("Model","Contig","Rep","Measure","Level"),suffixes=c("_gc","_truth"))
gcd<-gcd[!is.na(gcd$Value_gc),]
if(0 != sum(is.na(gcd$Value_truth))) stop("Truth has NA values")

gcd$diff2<-(gcd$Value_gc-gcd$Value_truth)^2


############################################################################################################



dt<-merge(gld,gcd,by=c("Model","Contig","Rep","Measure","Level","Depth"),suffixes=c("_gl","_gc"))
dt$Method<-NA
dt$Method_gl<-NA
dt$Method_gc<-NA
dt[dt$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5",]$Method_gl<-"GL_5_500"
dt[dt$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5_doMajorMinor1",]$Method_gl<-"GL_5_500_inferMajorMinor"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1",]$Method_gc<-"GC_095_1"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2",]$Method_gc<-"GC_095_2"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1",]$Method_gc<-"GC_03_1"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2",]$Method_gc<-"GC_03_2"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1_doMajorMinor1",]$Method_gc<-"GC_095_1_inferMajorMinor"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2_doMajorMinor1",]$Method_gc<-"GC_095_2_inferMajorMinor"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1_doMajorMinor1",]$Method_gc<-"GC_03_1_inferMajorMinor"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2_doMajorMinor1",]$Method_gc<-"GC_03_2_inferMajorMinor"

dt$Method<-paste0(dt$Method_gl,"_",dt$Method_gc)


sum(is.na(dt$diff2_gl))
sum(is.na(dt$diff2_gt))
sum(is.infinite(dt$diff2_gl))
sum(is.infinite(dt$diff2_gt))

dt$diff2<-dt$diff2_gl+dt$diff2_gc
dt$Depth<-as.factor(dt$Depth)
dt%>%filter(Depth!=0.1 & Model == "model1")%>%filter(Depth %in% c(10,20,50)) %>%
ggplot(aes(x=Depth, y=diff2, color=Contig))+
 geom_boxplot()+
    theme_bw()+theme(legend.position="bottom")+
    facet_wrap(~Method, scales="free")+
    labs(x="Depth",y="Difference of squared distance from truth (GL-GC)")
ggsave("diffOfSquaredDistanceFromTruth_betweenMethods.png",width=12,height=10,units="in",dpi=300)

rm(dt)
gc()

############################################################################################################


sum(is.infinite(gcd$diff2))
sum(is.infinite(gld$diff2))
sum(is.na(gcd$diff2))
sum(is.na(gld$diff2))


gcd$Method<-NA
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1",]$Method<-"GC_095_1"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2",]$Method<-"GC_095_2"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1",]$Method<-"GC_03_1"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2",]$Method<-"GC_03_2"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1_doMajorMinor1",]$Method<-"GC_095_1_inferMajorMinor"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2_doMajorMinor1",]$Method<-"GC_095_2_inferMajorMinor"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1_doMajorMinor1",]$Method<-"GC_03_1_inferMajorMinor"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2_doMajorMinor1",]$Method<-"GC_03_2_inferMajorMinor"

gld$Method<-NA
gld[gld$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5",]$Method<-"GL_5_500"
gld[gld$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5_doMajorMinor1",]$Method<-"GL_5_500_inferMajorMinor"

gcd$Depth<-as.factor(gcd$Depth)
gld$Depth<-as.factor(gld$Depth)

# calculate RMSE among replicates (N is replicates)
gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE=sqrt(sum(diff2)/N))%>%filter(Depth!=0.1 & Model == "model1")%>%filter(Depth %in% c(10,20,50))%>%ggplot(aes(x=Depth, y=RMSE, color=Contig))+
 geom_boxplot()+
    theme_bw()+theme(legend.position="bottom")+
    facet_wrap(~Method, scales="free")+
    labs(x="Depth",y="RMSE")

ggsave("RMSE_betweenMethods.png",width=12,height=10,units="in",dpi=300)

d_rmse_gcd<-gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n())
# here we lose some of the replicates altogether due to not having enough data
# all these are GC_095_1 so it is due to using postCutoff 0.95 and doPost 1
d_rmse_gcd[d_rmse_gcd$N!=20,]

d_rmse_gld<-gld%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n())
# OK; all are equal to 20
d_rmse_gld[d_rmse_gld$N!=20,]

rm(d_rmse_gcd,d_rmse_gld)
gc()


############################################################################################################

rmse_gcd<-gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE = sqrt(sum(diff2)/N))
rmse_gld<-gld%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE = sqrt(sum(diff2)/N))

df<-rbind(rmse_gcd,rmse_gld)
df$Contig <- as.factor(df$Contig)
df$Method <- as.factor(df$Method)

df %>% filter(Depth!=0.1 & Model == "model1") %>% filter(Method %in% c("GC_095_1","GC_095_2","GC_03_1","GC_03_2", "GL_5_500")) %>%
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method, group=interaction(Contig,Method)), size=2)+
    geom_line(aes(x=Depth, y=RMSE, linetype=Contig, group=interaction(Contig,Method),color=Method), size=1)+
    theme_bw()+theme(legend.position="bottom")+
    scale_linetype_manual(values = c("dotted","longdash","solid"))+
    facet_wrap(~Level, scales="free")


df %>% filter(Model == "model1") %>% filter(Method %in% c("GC_095_1","GC_095_2","GC_03_1","GC_03_2", "GL_5_500")) %>%
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method, group=interaction(Contig,Method)), size=2)+
    geom_line(aes(x=Depth, y=RMSE, linetype=Contig, group=interaction(Contig,Method),color=Method), size=1)+
    theme_bw()+theme(legend.position="bottom")+
    scale_linetype_manual(values = c("dotted","longdash","solid"))+
    facet_wrap(~Level, scales="free")

#   scale_y_continuous(n.breaks =20)+

############################################################################################################
############################################################################################################
# do not use inferMajorMinor results for now
############################################################################################################
############################################################################################################
rm(list=ls())
gc()


resdir<-"/home/isin/Mount/AMOVA/simulations/msprime/"
fn_glf<-"simulations/sim_demes_v2/results/collect_results_genotype_likelihood_fixedMajorMinor.csv"
fn_truth<-"simulations/sim_demes_v2/results/collect_results_truth.csv"
# fn_gcf<-"simulations/sim_demes_v2/results/collect_results_genotype_calling_fixedMajorMinor.csv"
fn_gcf<-"simulations/sim_demes_v2/results/collect_results_genotype_calling_fixedMajorMinor_run2.csv"

gl<-read.csv(paste0(resdir,fn_glf),header=T,sep=",")
truth<-read.csv(paste0(resdir,fn_truth),header=T,sep=",")
truth<-truth[c("AnalysisType","Model","Contig","Rep","Measure","Level","Value")]
gc<-read.csv(paste0(resdir,fn_gcf),header=T,sep=",")

gl[is.infinite(gl$Value),]$Value<-NA
gl$Contig<-as.factor(gl$Contig)
gl$Rep<-as.factor(gl$Rep)

gc[is.infinite(gc$Value),]$Value<-NA
gc$Contig<-as.factor(gc$Contig)
gc$Rep<-as.factor(gc$Rep)

## 
modelNames<-c("Model 1","Model 2")
modelList<-c("model1","model2")
names(modelNames)<-modelList
phiStatList<-c("Phi_ST","Phi_CT","Phi_SC")
contig_lens_lut<-c(1,2,10,20,50,100)
names(contig_lens_lut)<-c("1e6","2e6","10e6","20e6","50e6","100e6")
method_labels <- c("Genotype Likelihood", "Genotype Calling", "Genotype Calling", "Genotype Calling", "Genotype Calling")
methodList<-c("GL_5_500", "GC_095_1", "GC_095_2", "GC_03_1", "GC_03_2")
names(method_labels) <- methodList


color_map_gl5_gc<-c("Genotype Likelihood" = "blue",
  "Genotype Calling" = "#FC4E07",
  "True value" = "black")


phiStatLevelLut<-c("Phi_ST"="Population_in_Total", "Phi_CT"="Region_in_Total", "Phi_SC"="Population_in_Region")

# maketex label generator for latex formatting
maketex <- function(string){
#   string <- gsub("Phi_CT", expression("$\\\u03A6_{CT}$"), string)
#   string <- gsub("Phi_SC", expression("$\\\u03A6_{SC}$"), string)
#   string <- gsub("Phi_ST", expression("$\\\u03A6_{ST}$"), string)
  string <- gsub(phiStatLevelLut["Phi_CT"][[1]], expression("$\\\u03A6_{CT}$"), string)
  string <- gsub(phiStatLevelLut["Phi_SC"][[1]], expression("$\\\u03A6_{SC}$"), string)
  string <- gsub(phiStatLevelLut["Phi_ST"][[1]], expression("$\\\u03A6_{ST}$"), string)
  TeX(string)
}

makepretex <- function(string){
#   string <- gsub("Phi_CT", expression("$\\\u03A6_{CT}$"), string)
#   string <- gsub("Phi_SC", expression("$\\\u03A6_{SC}$"), string)
#   string <- gsub("Phi_ST", expression("$\\\u03A6_{ST}$"), string)
    string <- gsub(phiStatLevelLut["Phi_CT"][[1]], expression("$\\\u03A6_{CT}$"), string)
    string <- gsub(phiStatLevelLut["Phi_SC"][[1]], expression("$\\\u03A6_{SC}$"), string)
    string <- gsub(phiStatLevelLut["Phi_ST"][[1]], expression("$\\\u03A6_{ST}$"), string)
    string
}

ContigsToUse<-c(1,10,100)
ContigsToUseId<-"contigs1-10-100"

colormap<-c("#FC4E07","blue","black")


##

gld<-merge(gl,truth, by=c("Model","Contig","Rep","Measure","Level"),suffixes=c("_gl","_truth"))
gld<-gld[!is.na(gld$Value_gl),]
if(0 != sum(is.na(gld$Value_truth))) stop("Truth has NA values")


gld$diff2<-(gld$Value_gl-gld$Value_truth)^2


gcd<-merge(gc,truth, by=c("Model","Contig","Rep","Measure","Level"),suffixes=c("_gc","_truth"))
gcd<-gcd[!is.na(gcd$Value_gc),]
if(0 != sum(is.na(gcd$Value_truth))) stop("Truth has NA values")

gcd$diff2<-(gcd$Value_gc-gcd$Value_truth)^2

dt<-merge(gld,gcd,by=c("Model","Contig","Rep","Measure","Level","Depth"),suffixes=c("_gl","_gc"))
dt$Method<-NA
dt$Method_gl<-NA
dt$Method_gc<-NA
dt[dt$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5",]$Method_gl<-"GL_5_500"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1",]$Method_gc<-"GC_095_1"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2",]$Method_gc<-"GC_095_2"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1",]$Method_gc<-"GC_03_1"
dt[dt$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2",]$Method_gc<-"GC_03_2"

dt$Method<-paste0(dt$Method_gl,"_",dt$Method_gc)


sum(is.na(dt$diff2_gl))
sum(is.na(dt$diff2_gt))
sum(is.infinite(dt$diff2_gl))
sum(is.infinite(dt$diff2_gt))

dt$diff2<-dt$diff2_gl+dt$diff2_gc
dt$Depth<-as.factor(dt$Depth)
dt%>%filter(Depth!=0.1 & Model == "model1")%>%filter(Depth %in% c(10,20,50)) %>%
ggplot(aes(x=Depth, y=diff2, color=Contig))+
 geom_boxplot()+
    theme_bw()+theme(legend.position="bottom")+
    facet_wrap(~Method, scales="free")+
    labs(x="Depth",y="Difference of squared distance from truth (GL-GC)")

ggsave("diffOfSquaredDistanceFromTruth_betweenMethods_fixedMajorMinor.png",width=10,height=10,units="in",dpi=300)


############################################################################################################


sum(is.infinite(gcd$diff2))
sum(is.infinite(gld$diff2))
sum(is.na(gcd$diff2))
sum(is.na(gld$diff2))

gcd$Method<-NA
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost1",]$Method<-"GC_095_1"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff095_doPost2",]$Method<-"GC_095_2"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost1",]$Method<-"GC_03_1"
gcd[gcd$AnalysisType_gc=="genotype_calling_postCutoff03_doPost2",]$Method<-"GC_03_2"

gld$Method<-NA
gld[gld$AnalysisType_gl=="genotype_likelihood_maxIter500_tole5",]$Method<-"GL_5_500"

gcd$Depth<-as.factor(gcd$Depth)
gld$Depth<-as.factor(gld$Depth)

# calculate RMSE among replicates (N is replicates)
gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE=sqrt(sum(diff2)/N))%>%filter(Depth!=0.1 & Model == "model1")%>%filter(Depth %in% c(10,20,50))%>%ggplot(aes(x=Depth, y=RMSE, color=Contig))+
 geom_boxplot()+
    theme_bw()+theme(legend.position="bottom")+
    facet_wrap(~Method, scales="free")+
    labs(x="Depth",y="RMSE")

ggsave("RMSE_betweenMethods.png",width=12,height=10,units="in",dpi=300)

d_rmse_gcd<-gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n())
# here we lose some of the replicates altogether due to not having enough data
# OK; all these are postCutoff 0.95 with either doPost 1 or 2 at low depths as expected
d_rmse_gcd[d_rmse_gcd$N!=20,]

d_rmse_gld<-gld%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n())
# OK; all are equal to 20
d_rmse_gld[d_rmse_gld$N!=20,]

rm(d_rmse_gcd,d_rmse_gld)
gc()


############################################################################################################

rmse_gcd<-gcd%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE = sqrt(sum(diff2)/N))
rmse_gld<-gld%>%group_by(Model,Contig,Depth,Method,Measure,Level)%>%summarize(N=n(),RMSE = sqrt(sum(diff2)/N))

df<-rbind(rmse_gcd,rmse_gld)
df$Contig <- as.factor(df$Contig)
df$Method <- as.factor(df$Method)


# methodList[-1] to exclude GL, which is already plotted in all
plotOutDir="results_230526"
for (method_i in methodList[-1]){
    for (model_i in modelList){

        df %>% 
            filter(Method %in% c(method_i,"GL_5_500")) %>% 
            filter(Model == model_i) %>%
            ggplot()+
            geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
            geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
            theme_bw()+
            labs(x="Depth",y="RMSE",
                color="Method")+
            facet_wrap(~Level,scale="free",labeller= as_labeller(maketex, default=label_parsed))+
            scale_y_continuous(n.breaks = 20)+
            theme(legend.position="bottom")+
            scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
            scale_color_manual(name="Method", labels=method_labels, values = colormap)
        
        ggsave(paste0(plotOutDir,"/RMSE_",method_i,"_",model_i, ".png"),width=10,height=10,units="in",dpi=300)

        df %>% 
            filter(Method %in% c(method_i,"GL_5_500")) %>% 
            filter(Model == model_i) %>%
            ggplot()+
            geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
            geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
            theme_bw()+
            labs(x="Depth",y="RMSE (log10)",
                color="Method")+
            facet_wrap(~Level,scale="free",labeller= as_labeller(maketex, default=label_parsed))+
            scale_y_continuous(n.breaks = 20)+
            theme(legend.position="bottom")+
            scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
            scale_color_manual(name="Method", labels=method_labels, values = colormap)

        ggsave(paste0(plotOutDir,"/RMSE_",method_i,"_",model_i, "_log10.png"),width=10,height=10,units="in",dpi=300)
    }
}




method_i<-methodList[-1][4]
model_i<-modelList[1]
colormap<-c("#FC4E07","blue","black")
# choose the method that seems to favor the other method the most
p_legend<-
  df %>% 
  filter(Model=="model1") %>%
  filter(Method %in% c(method_i,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="bottom")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="A")+
  theme(plot.tag = element_text(size=20, face="bold"), 
    legend.box="vertical", legend.margin=margin()
  )


legend<- get_legend(p_legend)

p1 <-
  df %>% 
  filter(Model=="model1") %>%
  filter(Method %in% c(method_i,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="A")+
  theme(plot.tag = element_text(size=20, face="bold"))
p1


p2<- 
  df %>% 
  filter(Model=="model2") %>%
  filter(Method %in% c(method_i,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="B")+
  theme(plot.tag = element_text(size=20, face="bold"))
p2  

# ggarrange(p1,p2,ncol=1, nrow=2, common.legend = TRUE, legend = "bottom")
grid<-gridExtra::grid.arrange(p1,p2, legend, layout_matrix=rbind(1,2,3), heights=c(1,1,0.2))
grid
ggsave(paste0(plotOutDir,"/","grid_RMSE_",method_i,"_",model_i,"_log10.png"),width=10,height=10, plot=grid)





## 

## multiple methods panel


rm(method_i)
method_i1<-methodList[-1][1]

p1<-df %>% 
  filter(Method %in% c(method_i1,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="RMSE (log10)",
      color="Method")+
  facet_wrap(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="A")+
  theme(plot.tag = element_text(size=20, face="bold"))
p1

method_i2<-methodList[-1][2]

p2<-df %>% 
  filter(Method %in% c(method_i2,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="",
      color="Method")+
  facet_wrap(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="B")+
  theme(plot.tag = element_text(size=20, face="bold"))
p2  



method_i3<-methodList[-1][3]

p3<-df %>% 
  filter(Method %in% c(method_i3,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_wrap(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="C")+
  theme(plot.tag = element_text(size=20, face="bold"))
p3


method_i4<-methodList[-1][4]

p4<-df %>% 
  filter(Method %in% c(method_i4,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="",
      color="Method")+
  facet_wrap(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="D")+
  theme(plot.tag = element_text(size=20, face="bold"))
p4

gridExtra::grid.arrange(p1,p2,p3,p4)
panel



panel<-gridExtra::grid.arrange(p1,p2,p3,p4)

ggsave(paste0(plotOutDir,"/","panel_RMSE_4methods_",model_i,"_log10.png"),width=10,height=10, plot=panel)



#### # panel with facet_grid

rm(method_i)
method_i1<-methodList[-1][1]

p1<-df %>% 
  filter(Method %in% c(method_i1,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="RMSE (log10)",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="A")+
  theme(plot.tag = element_text(size=20, face="bold"))
p1

method_i2<-methodList[-1][2]

p2<-df %>% 
  filter(Method %in% c(method_i2,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="B")+
  theme(plot.tag = element_text(size=20, face="bold"))
p2  



method_i3<-methodList[-1][3]

p3<-df %>% 
  filter(Method %in% c(method_i3,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="C")+
  theme(plot.tag = element_text(size=20, face="bold"))
p3


method_i4<-methodList[-1][4]

p4<-df %>% 
  filter(Method %in% c(method_i4,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="D")+
  theme(plot.tag = element_text(size=20, face="bold"))
p4

gridExtra::grid.arrange(p1,p2,p3,p4)



panel<-gridExtra::grid.arrange(p1,p2,p3,p4)

ggsave(paste0(plotOutDir,"/","v2_panel_RMSE_4methods_",model_i,"_log10.png"),width=10,height=10, plot=panel)



### version 3 

method_i1<-methodList[-1][1]

p1<-df %>% 
  filter(Method %in% c(method_i1,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="RMSE (log10)",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  theme(strip.text.y=element_blank())+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="A")+
  theme(plot.tag = element_text(size=20, face="bold"))
p1

method_i2<-methodList[-1][2]

p2<-df %>% 
  filter(Method %in% c(method_i2,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="",y="",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="B")+
  theme(plot.tag = element_text(size=20, face="bold"))
p2  



method_i3<-methodList[-1][3]

p3<-df %>% 
  filter(Method %in% c(method_i3,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10)",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  theme(strip.text.x=element_blank())+
  theme(strip.text.y=element_blank())+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="C")+
  theme(plot.tag = element_text(size=20, face="bold"))
p3


method_i4<-methodList[-1][4]

p4<-df %>% 
  filter(Method %in% c(method_i4,"GL_5_500")) %>% 
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="",
      color="Method")+
  facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  theme(strip.text.x=element_blank())+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) +
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="D")+
  theme(plot.tag = element_text(size=20, face="bold"))
p4

gridExtra::grid.arrange(p1,p2,p3,p4, legend, ncol=2, layout_matrix = rbind(c(1,2),c(3,4),5), heights=c(10,10,3))
panel<-gridExtra::grid.arrange(p1,p2,p3,p4, legend, ncol=2, layout_matrix = rbind(c(1,2),c(3,4),5), heights=c(10,10,3))


# panel<-gridExtra::grid.arrange(p1,p2,p3,p4)

ggsave(paste0(plotOutDir,"/","v3_panel_RMSE_4methods_",model_i,"_log10.png"),width=10,height=10, plot=panel)


## try doing it in one plot


ggplot()+
  geom_point(data=df[df$Method%in%c(method_i1,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(data=df[df$Method%in%c(method_i1,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  geom_point(data=df[df$Method%in%c(method_i2,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(data=df[df$Method%in%c(method_i2,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
geom_point(data=df[df$Method%in%c(method_i3,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
geom_line(data=df[df$Method%in%c(method_i3,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
    geom_point(data=df[df$Method%in%c(method_i4,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
    geom_line(data=df[df$Method%in%c(method_i4,"GL_5_500"),],aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="",
      color="Method")+
  facet_grid(Method~Level+Model,labeller= as_labeller(maketex, default=label_parsed))+
#   facet_grid(Model~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  scale_y_continuous(n.breaks =10)+
  theme(legend.position="none")+
  theme(strip.text.x=element_blank())+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) 
  scale_color_manual(name="Method", labels=method_labels, values = colormap)+
  labs(tag="D")+
  theme(plot.tag = element_text(size=20, face="bold"))

