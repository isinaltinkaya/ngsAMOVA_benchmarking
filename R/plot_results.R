############################################################################################################
# plot_results.R
# 
# isinaltinkaya@gmail.com
# 230307
#
############################################################################################################

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
names(lut) <-c("GLE_5", "GTC_1_95", "GTC_2_95", "GTC_1_3", "GTC_2_3")


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



# scale free
df%>% filter(Contig%in%c(1,10,100)) %>%
  ggplot()+
  geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
  geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
  theme_bw()+
  labs(x="Depth",y="RMSE",
       title="RMSE of Phi Estimates",
       color="Method",
       subtitle=model_id)+
  facet_wrap(~Level,scale="free")+
  scale_y_continuous(n.breaks =20)+
  theme(legend.position="bottom")+
  linetype3


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
# 

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




color_map<-c("Genotype Likelihood" = "blue",
  "Genotype Calling 1 0.95" = "red",
  "Genotype Calling 2 0.95" = "green",
  "Genotype Calling 1 0.3" = "orange",
  "Genotype Calling 2 0.3" = "purple",
  "True" = "black")


ggplot()+
  geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood"),shape=1)+
  geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Likelihood",linetype=Contig))+
  # geom_point(data=rm_gc_195,aes(x=Depth,y=Value,color="Genotype Calling 1"),shape=1)+
  # geom_line(data=rm_gc_195,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling 1 0.95",linetype=Contig))+
  # geom_point(data=rm_gc_295,aes(x=Depth,y=Value,color="Genotype Calling 2"),shape=1)+
  # geom_line(data=rm_gc_295,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling 2 0.95",linetype=Contig))+
  geom_point(data=rm_gc_103,aes(x=Depth,y=Value,color="Genotype Calling 2"),shape=1)+
  geom_line(data=rm_gc_103,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling 1 0.3",linetype=Contig))+
  geom_point(data=rm_gc_203,aes(x=Depth,y=Value,color="Genotype Calling 2"),shape=1)+
  geom_line(data=rm_gc_203,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Calling 2 0.3",linetype=Contig))+
  geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True"))+
  facet_wrap(~Level,scale="free")+
  theme_bw()+
  scale_y_continuous(n.breaks =20)+
  scale_color_manual(values=color_map)+
  theme(legend.position="bottom")+
  labs(colour="",
       y="Estimated value of Phi statistic",
       title="Mean estimated Phi statistic values across 20 replicates",
       subtitle=model_id)+
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


df%>%filter(Method==method_id)%>% filter(Level==phi_stat) %>%  group_by(Contig,Depth) %>% summarize(mmeanNSites=floor(mean(meanNSites)),RMSE=RMSE) %>%
  ggplot(aes(x=Depth,y=Contig,fill=RMSE))+
  geom_tile()+
  scale_fill_viridis()+
  coord_fixed()+
  theme_bw()+
  geom_text(aes(x=Depth,y=Contig,label=paste0(format(mmeanNSites,big.mark=","),"\n",format(RMSE))),check_overlap=TRUE,color="white",fontface="bold")+
  labs(title=paste0("RMSE of ",lut[method_id][[1]]," for ",phi_stat," using ", model_id),
  y="Number of sites simulated")+
  scale_y_discrete(labels=contig_lens)




############################################################################################################
## MODEL 2

model_id="model2"

