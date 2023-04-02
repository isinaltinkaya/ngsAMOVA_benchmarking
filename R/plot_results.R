############################################################################################################
# plot_results.R
# 
# isinaltinkaya@gmail.com
# 230307
#
############################################################################################################
library(latex2exp)

source("R/read_data.R")


## --------------------------------------------------------------------------------------------
## MODEL 1

model_id="model1"


resdir="/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_march23/"

gc_phi_fn<-paste0(resdir,"sim_demes_v2_genotype_calling_phi_",model_id,".csv")
true_phi_fn<-paste0(resdir,"sim_demes_v2_truth_phi_",model_id,".csv")
gle_phi_fn<-paste0(resdir,"sim_demes_v2_genotype_likelihood_phi_",model_id,".csv")
nsu_gc_fn<-paste0(resdir,"sim_demes_v2_genotype_calling_mean_nSites_used_",model_id,".csv")
df <- get_df(true_phi_fn,gc_phi_fn, gle_phi_fn, nsu_gc_fn)

contig_lens<-as.integer(levels(df$Contig))*1e6

linetype3=scale_linetype_manual(values = c("dotted","dashed","solid"))
lut <- c("Genotype Likelihood", "Genotype Calling 1 (0.95)", "Genotype Calling 2 (0.95)", "Genotype Calling 1 (0.3)", "Genotype Calling 2 (0.3)")
# names(lut) <-c("GLE_5", "GTC_1_95", "GTC_2_95", "GTC_1_3", "GTC_2_3")
names(lut) <-c("GLE_5", "GTC_1_0.95", "GTC_2_0.95", "GTC_1_0.3", "GTC_2_0.3")


gc_phi <- read_gc_phi(gc_phi_fn)
gle_phi <- read_gle_phi(gle_phi_fn)
true_phi <- read_true_phi(true_phi_fn)
gc_103<-gc_phi%>%filter(Method == "GTC_1_0.3")
gc_195<-gc_phi%>%filter(Method == "GTC_1_0.95")
gc_203<-gc_phi%>%filter(Method == "GTC_2_0.3")
gc_295<-gc_phi%>%filter(Method == "GTC_2_0.95")

rm_gc_103<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_103, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))
rm_gc_195<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_195, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))
rm_gc_203<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_203, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))
rm_gc_295<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_295, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))
rm_gle_phi<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gle_phi, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))
rm_real_phi<-aggregate(Value ~ Measure + Level + Model + Contig, true_phi, function(x) mean(x)) %>%filter(Contig %in% c(1,10,100))

# method_labels <- c("Genotype Likelihood", "Genotype Calling 1 (0.95)", "Genotype Calling 2 (0.95)", "Genotype Calling 1 (0.3)", "Genotype Calling 2 (0.3)")
# method_labels <- c("Genotype Likelihood", "Genotype Calling -doPost 1 -postCutoff 0.95", "Genotype Calling -doPost 2 -postCutoff 0.95", "Genotype Calling -doPost 1 -postCutoff 0.3", "Genotype Calling -doPost 2 -postCutoff 0.3")
method_labels <- c("Genotype Likelihood", "Genotype Calling", "Genotype Calling", "Genotype Calling", "Genotype Calling")
names(method_labels) <-c("GLE_5", "GTC_1_0.95", "GTC_2_0.95", "GTC_1_0.3", "GTC_2_0.3")

# labeller for latex formatting
phi.labs <- c("Phi_CT","Phi_SC", "Phi_ST")
names()

maketex <- function(string){
  # convert Phi_CT to $\Phi_{CT}$
  string <- gsub("Phi_CT", "$\\\\Phi_{CT}$", string)
  string <- gsub("Phi_SC", "$\\\\Phi_{SC}$", string)
  string <- gsub("Phi_ST", "$\\\\Phi_{ST}$", string)
  TeX(string)
}


modelNames<-c("Model 1","Model 2")
names(modelNames)<-c("model1","model2")
#maketex(df$Level)

color_map_gl5_gt203_rmse<-c("GLE_5" = "blue",
  "GTC_2_0.3" = "#FC4E07")

phiStatList<-c("Phi_ST","Phi_CT","Phi_SC")
modelList<-c("model1","model2")

color_map_gl5_gc<-c("Genotype Likelihood" = "blue",
  "Genotype Calling" = "#FC4E07",
  "True value" = "black")

for(model in modelList){
  # scale free
  df %>% filter(Method %in% c("GTC_2_0.3","GLE_5")) %>% filter(Contig%in%c(1,10,100)) %>%
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE",
        color="Method")+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    linetype3 + 
    scale_color_manual(name="Method", labels=method_labels, values = color_map_gl5_gt203_rmse)

  ggsave(paste0("figures/phi_rmse_gtc203_gle5_scalefree_contig110100_",model_id,".png"))

  # scale free, y log
  df %>% filter(Method %in% c("GTC_2_0.3","GLE_5")) %>% filter(Contig%in%c(1,10,100)) %>%
    ggplot()+
    geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE (log10 transformed)",
        color="Method")+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    linetype3 + 
    #change the names of the methods in the legend
    scale_color_manual(name = "Method", labels = method_labels, values = color_map_gl5_gt203_rmse)
  ggsave(paste0("figures/phi_rmse_gtc203_gle5_scalefree_contig110100_",model_id,"_log.png"))
}
 

gc_data_list<-c("rm_gc_295", "rm_gc_195", "rm_gc_203", "rm_gc_103")
names(gc_data_list)<-c("gc295","gc195", "gc203", "gc103")
get_gc_data<-function(str){
  return(get(gc_data_list[str]))
}

for(gc_data_i in names(gc_data_list)){
  ggplot()+
    geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood"),shape=1)+
    geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
    color="Genotype Likelihood",linetype=Contig))+
    geom_point(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,color="Genotype Calling"),shape=1)+
    geom_line(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,group=interaction(Contig,Level),
    color="Genotype Calling",linetype=Contig))+
    geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value"))+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    theme_bw()+
    scale_y_continuous(n.breaks =20)+
    scale_color_manual(values=color_map_gl5_gc)+
    theme(legend.position="bottom")+
    labs(colour="",
        y="Phi Estimate")+
    linetype3

  ggsave(paste0("figures/phi_estimates_",gc_data_i,"_gle5_scalefree_contig110100_",model_id,".png"))
}

for(gc_data_i in names(gc_data_list)){

  ggplot()+
    geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood",shape=Contig))+
    geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
    color="Genotype Likelihood",linetype=Contig))+
    geom_point(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,color="Genotype Calling",shape=Contig))+
    geom_line(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,group=interaction(Contig,Level),
    color="Genotype Calling",linetype=Contig))+
    geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value",shape=Contig))+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    theme_bw()+
    scale_y_continuous(n.breaks =20)+
    scale_color_manual(values=color_map_gl5_gc)+
    theme(legend.position="bottom")+
    labs(colour="",
        y="Phi Estimate")+
    linetype3
  # ggplot()+
  #   geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood"),shape=1)+
  #   geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  #   color="Genotype Likelihood",linetype=Contig))+
  #   geom_point(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,color="Genotype Calling"),shape=1)+
  #   geom_line(data=get_gc_data(gc_data_i),aes(x=Depth,y=Value,group=interaction(Contig,Level),
  #   color="Genotype Calling",linetype=Contig))+
  #   geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value"))+
  #   facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  #   theme_bw()+
  #   scale_y_continuous(n.breaks =20)+
  #   scale_color_manual(values=color_map_gl5_gc)+
  #   theme(legend.position="bottom")+
  #   labs(colour="",
  #       y="Phi Estimate")+
  #   linetype3

  ggsave(paste0("figures/phi_estimates_",gc_data_i,"_gle5_scalefree_contig110100_",model_id,".png"))
}



color_map<-c("Genotype Likelihood" = "blue",
  "Genotype Calling (1,0.95)" = "green",
  "Genotype Calling (2,0.95)" = "red",
  "Genotype Calling (1,0.3)" = "blue",
  "Genotype Calling (2,0.3)" = "orange",
  "True value" = "black")



ggplot()+
  geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood"),shape=1)+
  geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Likelihood",linetype=Contig))+
  geom_point(data=rm_gc_195,aes(x=Depth,y=Value,color="Genotype Calling (1,0.95)"),shape=1)+
  geom_line(data=rm_gc_195,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling (1,0.95)",linetype=Contig))+
  geom_point(data=rm_gc_295,aes(x=Depth,y=Value,color="Genotype Calling (2,0.95)"),shape=1)+
  geom_line(data=rm_gc_295,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling (2,0.95)",linetype=Contig))+
  geom_point(data=rm_gc_103,aes(x=Depth,y=Value,color="Genotype Calling (1,0.3)"),shape=1)+
  geom_line(data=rm_gc_103,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling (1,0.3)",linetype=Contig))+
  geom_point(data=rm_gc_203,aes(x=Depth,y=Value,color="Genotype Calling (2,0.3)"),shape=1)+
  geom_line(data=rm_gc_203,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling (2,0.3)",linetype=Contig))+
  geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value"))+
  facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  theme_bw()+
  scale_y_continuous(n.breaks =20)+
  scale_color_manual(values=color_map)+
  theme(legend.position="bottom")+
  labs(colour="",
      y="Phi Estimate")+
      #  title="Mean Phi Estimates Across 20 Replicates",
      #  subtitle=model_id)+
  linetype3






ggplot()+
  geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="GenotypeLikelihood"),shape=1)+
  geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig),color="GenotypeLikelihood",linetype=Contig))+
  geom_point(data=rm_gc_203,aes(x=Depth,y=Value,color="GenotypeCalling"),shape=1)+
  geom_line(data=rm_gc_203,aes(x=Depth,y=Value,group=interaction(Contig),color="GenotypeCalling",linetype=Contig))+
  geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="Truth"))+
  facet_wrap(~Level,scale="free")+
  theme_bw()+
  scale_color_manual(values=c("red","blue","black"))+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="bottom")+
  labs(colour="",
       y="Estimated value of Phi statistic",
       title="Mean estimated Phi statistic values across 20 replicates",
       subtitle=model_id)+
  linetype3



##


phi_stat<-"Phi_ST"
method_id<-"GLE_5"

phi_stat_lut<-c("$\\Phi_{ST}$","$\\Phi_{CT}$","$\\Phi_{SC}$")
names(phi_stat_lut)<-c("Phi_ST","Phi_CT","Phi_SC")

# n.b. nSites is mean number of sites
df%>%filter(Method==method_id)%>% filter(Level==phi_stat) %>%  group_by(Contig,Depth) %>% summarize(mmeanNSites=floor(mean(meanNSites)),RMSE=RMSE) %>%
  ggplot(aes(x=Depth,y=Contig,fill=RMSE))+
  geom_tile(color="white")+
  scale_fill_viridis()+
  coord_equal()+
  theme_bw()+
  geom_text(aes(x=Depth,y=Contig,label=paste0("nSites: ",format(mmeanNSites,big.mark=","),"\nRMSE: ",format(RMSE))),check_overlap=TRUE,color="white",fontface="bold")+
  labs( y="Number of sites simulated")+
  # ggtitle(TeX(paste0("RMSE of ",lut[method_id][[1]]," for ",phi_stat_lut[phi_stat]))) +
  scale_y_discrete(labels=contig_lens)

ggsave(paste0("figures/rmse_heatmap_",method_id,"_",phi_stat,".png"))


rmse_heatmap_plotter <- function(df,method_id,phi_stat){
  df%>%filter(Method==method_id)%>% filter(Level==phi_stat) %>%  group_by(Contig,Depth) %>% summarize(mmeanNSites=floor(mean(meanNSites)),RMSE=RMSE) %>%
    ggplot(aes(x=Depth,y=Contig,fill=RMSE))+
    geom_tile(color="white")+
    scale_fill_viridis()+
    coord_equal()+
    theme_bw()+
    geom_text(aes(x=Depth,y=Contig,label=paste0("nSites: ",format(mmeanNSites,big.mark=","),"\nRMSE: ",format(RMSE))),check_overlap=TRUE,color="white",fontface="bold")+
    labs( y="Number of sites simulated")+
    scale_y_discrete(labels=contig_lens)
    ggsave(paste0("figures/rmse_heatmap_",method_id,"_",phi_stat,"_",model_id,".png"))
}

for(model in modelList){
  for (phiStat in phiStatList){
    for (method in c("GLE_5","GTC_1_0.3","GTC_2_0.3","GTC_1_0.95","GTC_2_0.95")){
      rmse_heatmap_plotter(df,method,phiStat)
    }
  }
}


  

###################################################
df%>% filter(Contig%in%c(1,10,100))%>%
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10 transformed)",
       title="RMSE of Phi Estimates",
       color="Method",
       subtitle=model_id)+
  facet_wrap(~Level)+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="bottom")+
  linetype3

# scale free, y log
df%>% filter(Contig%in%c(1,10,100))%>%
  ggplot()+
  geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE (log10 transformed)",
       title="RMSE of Phi Estimates",
       color="Method",
       subtitle=model_id)+
  facet_wrap(~Level,scale="free")+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="bottom")+
  linetype3
###################################################



