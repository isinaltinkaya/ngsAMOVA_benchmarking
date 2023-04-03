############################################################################################################
# plot_results.R
# 
# isinaltinkaya@gmail.com
# 230307
#
############################################################################################################

source("R/read_data.R")

############################################################################################################

modelNames<-c("Model 1","Model 2")
names(modelNames)<-c("model1","model2")
phiStatList<-c("Phi_ST","Phi_CT","Phi_SC")
modelList<-c("model1","model2")
contig_lens_lut<-c(1,2,10,20,50,100)
names(contig_lens_lut)<-c("1e6","2e6","10e6","20e6","50e6","100e6")
method_labels <- c("Genotype Likelihood", "Genotype Calling", "Genotype Calling", "Genotype Calling", "Genotype Calling")
names(method_labels) <-c("GLE_5", "GTC_1_0.95", "GTC_2_0.95", "GTC_1_0.3", "GTC_2_0.3")


color_map_gl5_gc<-c("Genotype Likelihood" = "blue",
  "Genotype Calling" = "#FC4E07",
  "True value" = "black")

# maketex label generator for latex formatting
maketex <- function(string){
  string <- gsub("Phi_CT", expression("$\\\u03A6_{CT}$"), string)
  string <- gsub("Phi_SC", expression("$\\\u03A6_{SC}$"), string)
  string <- gsub("Phi_ST", expression("$\\\u03A6_{ST}$"), string)
  TeX(string)
}

makepretex <- function(string){
  string <- gsub("Phi_CT", expression("$\\\u03A6_{CT}$"), string)
  string <- gsub("Phi_SC", expression("$\\\u03A6_{SC}$"), string)
  string <- gsub("Phi_ST", expression("$\\\u03A6_{ST}$"), string)
  string
}


############################################################################################################


resdir="/home/isin/Mount/AMOVA/simulations/msprime/simulations/sim_demes_v2/collected_results/collect_march23/"


ContigsToUse<-c(1,10,100)
ContigsToUseId<-"contigs1-10-100"

# contig_lens<-as.integer(levels(df$Contig))*1e6

for(model_id in modelList){
  gc_phi_fn<-paste0(resdir,"sim_demes_v2_genotype_calling_phi_",model_id,".csv")
  true_phi_fn<-paste0(resdir,"sim_demes_v2_truth_phi_",model_id,".csv")
  gle_phi_fn<-paste0(resdir,"sim_demes_v2_genotype_likelihood_phi_",model_id,".csv")
  nsu_gc_fn<-paste0(resdir,"sim_demes_v2_genotype_calling_mean_nSites_used_",model_id,".csv")

  gc_phi <- read_gc_phi(gc_phi_fn) %>% filter(Contig %in% ContigsToUse, Model == model_id)
  gle_phi <- read_gle_phi(gle_phi_fn) %>% filter(Contig %in% ContigsToUse, Model == model_id)
  true_phi <- read_true_phi(true_phi_fn) %>% filter(Contig %in% ContigsToUse, Model == model_id)
  nsu<-get_nSites_used(nsu_gc_fn)

  df <- get_df2(true_phi = true_phi, gle_phi = gle_phi, gc_phi = gc_phi, nsites_used = nsu)

  gc_103<-gc_phi%>%filter(Method == "GTC_1_0.3")
  gc_195<-gc_phi%>%filter(Method == "GTC_1_0.95")
  gc_203<-gc_phi%>%filter(Method == "GTC_2_0.3")
  gc_295<-gc_phi%>%filter(Method == "GTC_2_0.95")

  rm_gc_103<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_103, function(x) mean(x))
  rm_gc_195<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_195, function(x) mean(x))
  rm_gc_203<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_203, function(x) mean(x)) 
  rm_gc_295<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gc_295, function(x) mean(x))
  rm_gle_phi<-aggregate(Value ~ Measure + Level + Model + Contig + Depth, gle_phi, function(x) mean(x))
  rm_real_phi<-aggregate(Value ~ Measure + Level + Model + Contig, true_phi, function(x) mean(x)) 




for(method_i in c("GTC_1_0.3","GTC_1_0.95","GTC_2_0.3","GTC_2_0.95")){
  
  colormap<-c("#FC4E07","blue","black")
  df %>% 
    filter(Method %in% c(method_i,"GLE_5")) %>% 
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE",
        color="Method")+
    facet_wrap(~Level,labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
    scale_color_manual(name="Method", labels=method_labels, values = colormap)
  ggsave(paste0("figures/phi_rmse_",method_i,"_gle5_",ContigsToUseId,"_",model_id,".png"), width=8, height=6, units="in", dpi=300)

  # y log
  df %>% 
    filter(Method %in% c(method_i,"GLE_5")) %>% 
    ggplot()+
    geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE (log10 transformed)",
        color="Method")+
    facet_wrap(~Level,labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
    scale_color_manual(name="Method", labels=method_labels, values = colormap)
  ggsave(paste0("figures/phi_rmse_",method_i,"_gle5_",ContigsToUseId,"_",model_id,"_log.png"), width=8, height=6, units="in", dpi=300)


  # scale free
  df %>% 
    filter(Method %in% c(method_i,"GLE_5")) %>% 
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE",
        color="Method")+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
    scale_color_manual(name="Method", labels=method_labels, values = colormap)

  ggsave(paste0("figures/phi_rmse_",method_i,"_gle5_scalefree_",ContigsToUseId,"_",model_id,".png"), width=8, height=6, units="in", dpi=300)



  # scale free, y log
  df %>% 
  filter(Method %in% c(method_i,"GLE_5")) %>% 
    ggplot()+
    geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE (log10 transformed)",
        color="Method")+
    facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) + 
    scale_color_manual(name="Method", labels=method_labels, values = colormap)
  ggsave(paste0("figures/phi_rmse_",method_i,"_gle5_scalefree_",ContigsToUseId,"_",model_id,"_log.png"), width=8, height=6, units="in", dpi=300)


  gc_data_list<-c("rm_gc_295", "rm_gc_195", "rm_gc_203", "rm_gc_103")
  names(gc_data_list)<-c("gc295","gc195", "gc203", "gc103")

  for(gc_data_i in names(gc_data_list)){
    ggplot()+
      geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,color="Genotype Likelihood"),shape=1)+
      geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),color="Genotype Likelihood",linetype=Contig))+
      geom_point(data=get(gc_data_list[gc_data_i]),aes(x=Depth,y=Value,color="Genotype Calling"),shape=1)+
      geom_line(data=get(gc_data_list[gc_data_i]),aes(x=Depth,y=Value,group=interaction(Contig,Level),color="Genotype Calling",linetype=Contig))+
      geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value"))+
      facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
      theme_bw()+
      scale_y_continuous(n.breaks =20)+
      scale_color_manual(values=color_map_gl5_gc)+
      theme(legend.position="bottom")+
      labs(colour="",
          y="Phi Estimate")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8"))
    ggsave(paste0("figures/phi_estimates_",gc_data_i,"_gle5_scalefree_",ContigsToUseId,"_",model_id,".png"), width=8, height=6, units="in", dpi=300)

  }
}



# color_map<-c("Genotype Likelihood" = "blue",
#   "Genotype Calling (1,0.95)" = "green",
#   "Genotype Calling (2,0.95)" = "red",
#   "Genotype Calling (1,0.3)" = "blue",
#   "Genotype Calling (2,0.3)" = "orange",
#   "True value" = "black")



ggplot()+
  geom_point(data=rm_gle_phi,aes(x=Depth,y=Value,
  color="Genotype Likelihood"),shape=1)+
  geom_line(data=rm_gle_phi,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  color="Genotype Likelihood",linetype=Contig))+
  geom_point(data=rm_gc_195,aes(x=Depth,y=Value,
  # color="Genotype Calling (1,0.95)"),shape=1)+
  color="GTC_1_0.95"),shape=1)+
  geom_line(data=rm_gc_195,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling (1,0.95)",linetype=Contig))+
  color="GTC_1_0.95",linetype=Contig))+
  geom_point(data=rm_gc_295,aes(x=Depth,y=Value,
  # color="Genotype Calling (2,0.95)"),shape=1)+
  color="GTC_2_0.95"),shape=1)+
  geom_line(data=rm_gc_295,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling (2,0.95)",linetype=Contig))+
  color="GTC_2_0.95",linetype=Contig))+
  geom_point(data=rm_gc_103,aes(x=Depth,y=Value,
  # color="Genotype Calling (1,0.3)"),shape=1)+
  color="GTC_1_0.3"),shape=1)+
  geom_line(data=rm_gc_103,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling (1,0.3)",linetype=Contig))+
  color="GTC_1_0.3",linetype=Contig))+
  geom_point(data=rm_gc_203,aes(x=Depth,y=Value,
  # color="Genotype Calling (2,0.3)"),shape=1)+
  color="GTC_2_0.3"),shape=1)+
  geom_line(data=rm_gc_203,aes(x=Depth,y=Value,group=interaction(Contig,Level),
  # color="Genotype Calling (2,0.3)",linetype=Contig))+
  color="GTC_2_0.3",linetype=Contig))+
  geom_point(data=rm_real_phi,aes(x=7.3,y=Value,color="True value"))+
  facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
  theme_bw()+
  scale_y_continuous(n.breaks =20)+
  # scale_color_manual(values=color_map)+
  theme(legend.position="bottom")+
  labs(colour="",
      y="Phi Estimate")+
      #  title="Mean Phi Estimates Across 20 Replicates",
      #  subtitle=model_id)+
  scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8"))
  ggsave(paste0("figures/phi_estimates_all_gle5_scalefree_",ContigsToUseId,"_",model_id,".png"), width=10, height=6, units="in", dpi=300)





phi_stat_lut<-c("$\\Phi_{ST}$","$\\Phi_{CT}$","$\\Phi_{SC}$")
names(phi_stat_lut)<-c("Phi_ST","Phi_CT","Phi_SC")

  methods_lut <- c("Genotype Likelihood", "Genotype Calling 1 (0.95)", "Genotype Calling 2 (0.95)", "Genotype Calling 1 (0.3)", "Genotype Calling 2 (0.3)")
  names(methods_lut) <-c("GLE_5", "GTC_1_0.95", "GTC_2_0.95", "GTC_1_0.3", "GTC_2_0.3")

  for (phi_stat in phiStatList){
    for (method_id in c("GLE_5","GTC_1_0.3","GTC_2_0.3","GTC_1_0.95","GTC_2_0.95")){
      df%>%
        filter(Method==method_id)%>% 
        filter(Level==phi_stat) %>%
        group_by(Contig,Depth) %>% 
        summarize(mmeanNSites=floor(mean(meanNSites)),RMSE=RMSE) %>%
        ggplot(aes(x=Depth,y=Contig,fill=RMSE))+
        geom_tile(color="white")+
        scale_fill_viridis()+
        coord_equal()+
        theme_bw()+
        geom_text(aes(x=Depth,y=Contig,label=paste0("nSites: ",format(mmeanNSites,big.mark=","),"\nRMSE: ",format(RMSE))),check_overlap=TRUE,color="white",fontface="bold")+
        labs( y="Contig size",x="Depth",fill="RMSE")+
        scale_y_discrete(labels=c("1"="1e6", "10"="1e7", "100"="1e8"))
        # ggtitle(paste0("RMSE of ",methods_lut[method_id][[1]]," Estimations"),subtitle=maketex(phi_stat))

        ggsave(paste0("figures/rmse_heatmap_",method_id,"_",phi_stat,"_",model_id,".png"), width=20, height=10, units="in", dpi=300)
    }
  }


  # rmse with all methods
  df %>% 
    ggplot()+
    geom_point(aes(x=Depth, y=RMSE, color=Method,group=interaction(Method,Contig)))+
    geom_line(aes(x=Depth, y=RMSE, color=Method, group=interaction(Method,Contig),linetype=Contig))+
    theme_bw()+
    labs(x="Depth",y="RMSE",
        color="Method")+
    facet_wrap(~Level,labeller= as_labeller(maketex, default=label_parsed))+
    scale_y_continuous(n.breaks =20)+
    theme(legend.position="bottom")+
    scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8"))

  ggsave(paste0("figures/phi_rmse_all_gle5_",ContigsToUseId,"_",model_id,".png"), width=10, height=6, units="in", dpi=300)

    df %>% 
      ggplot()+
      geom_point(aes(x=Depth, y=log10(RMSE), color=Method,group=interaction(Method,Contig)))+
      geom_line(aes(x=Depth, y=log10(RMSE), color=Method, group=interaction(Method,Contig),linetype=Contig))+
      theme_bw()+
      labs(x="Depth",y="RMSE (log10 transformed)",
          color="Method")+
      facet_wrap(~Level,scale="free", labeller= as_labeller(maketex, default=label_parsed))+
      scale_y_continuous(n.breaks =20)+
      theme(legend.position="bottom")+
      scale_linetype_manual(name="Contig size",values=c("1"="dotted","10"="dashed","100"="solid"), labels=c("1"="1e6", "10"="1e7", "100"="1e8")) 

    ggsave(paste0("figures/phi_rmse_all_gle5_",ContigsToUseId,"_",model_id,"_log.png"), width=10, height=6, units="in", dpi=300)




}
