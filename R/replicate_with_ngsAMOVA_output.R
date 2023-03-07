library("readr")

dmf<-"/home/isin/Projects/AMOVA/feb23/ngsAMOVA/amovaput.distance_matrix.csv"
dm<-as.data.frame(t(read.csv(dmf,header=F,sep = ',',colClasses = "double")))
colnames(dm)<-"distance"
dm$PairIdx<-seq(1,length(dm$distance))

jgcdf<-"/home/isin/Projects/AMOVA/feb23/ngsAMOVA/amovaput.joint_geno_count_dist.csv"

d<-read.csv(jgcdf,header = F)

get_dij<-function(d){
  1 - (((d$cA/d$nSites) + (d$cI/d$nSites)) + (((d$cB/d$nSites) + (d$cD/d$nSites) + (d$cE/d$nSites) + (d$cF/d$nSites) + (d$cH/d$nSites))/2))
}

colnames(d)<-c("PairIdx","cA","cB","cC","cD","cE","cF","cG","cH","cI","nSites")
#colnames(d)<-c("PairIdx","cA","cD","cG","cB","cE","cH","cC","cF","cI","nSites")

d$dij2<-get_dij(d)^2
all.equal(d$dij2,dm$distance)

mtdf<-"/home/isin/Projects/AMOVA/feb23/ngsAMOVA/wd/metadata_2lvl_with_header.tsv"
mtd<-read.csv(mtdf,header = T,sep="\t")

dcol="dij2"
nInd=40
m<-matrix(NA,nInd,nInd)
diag(m)<-0


m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]
dd.d<-as.dist(m)
dd.pops<-factor(mtd$Population)
dd.regs<-factor(mtd$Region)

amv<-pegas::amova(dd.d ~ dd.regs/dd.pops,is.squared=TRUE,nperm = 0)

