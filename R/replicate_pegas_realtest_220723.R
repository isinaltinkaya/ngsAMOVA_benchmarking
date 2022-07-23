d <- read_csv("sim2-OutOfAfrica_3G09-chr22-rep0.sfs.csv")

d<-set_pops(d)              
ncols<-c("Ind1","Ind2")
icols<-c("I_g1","I_g2")
gcols<-c("Pop1","Pop2")
dcol<-"Fij"

d[icols]<-"sameI"


nInd<-length(unique(c(unique(d[ncols][1][[1]]),unique(d[ncols][2][[1]]))))

m<-matrix(NA,nInd,nInd)
diag(m)<-0

m[lower.tri(m,diag=FALSE)]<- d[dcol][[1]]

#print(m,upper=F,diag=F)

#validate correct positions of pairwise distances
m[40,][2]
d%>%filter(Ind1=="pop1_ind2",Ind2=="pop1_ind40")%>%select(all_of(dcol))



dd.d<-as.dist(m)

dd.pops<-factor(c(rep("pop1",50),rep("pop2",50),rep("pop3",50)))



pegas::amova(dd.d ~ dd.pops,is.squared=FALSE,nperm=0)


dd.I<-factor(c(rep("sameI",nInd)))
pegas::amova(dd.d ~ dd.I/dd.pops,is.squared=FALSE,nperm=0)



df_ag<-get_df_ag(d=d, G_cols=gcols)
ssd_ag<-get_ssd_ag(d=d,d_col=dcol,G_cols=gcols,I_cols=icols,N_cols=ncols)
msd_ag<-get_msd(ssd = ssd_ag,df=df_ag)

df_wp<-get_df_wp(d=d, N_cols=ncols,G_cols=gcols, I_cols=icols)
ssd_wp<-get_ssd_wp(d=d,d_col=dcol,N_cols=ncols,G_cols=gcols,I_cols=icols)
msd_wp<-get_msd(ssd = ssd_wp,df=df_wp)


df_total<-get_df_total(d=d, N_cols=ncols)
ssd_total<-get_ssd_total(d=d,d_col=dcol,N_cols=ncols)
msd_total<-get_msd(ssd = ssd_total,df=df_total)

n0<-get_n_0(d=d,G_cols=gcols,I_cols=icols)
n1<-get_n_1(d=d,G_cols=gcols,I_cols=icols)
n2<-get_n_2(d=d,N_cols=ncols,G_cols=gcols,I_cols=icols)

#(50+50+50)-(((50^2)+(50^2)+(50^2))/(50+50+50))

sigmasq_a<-get_sigmasq_a(msd_ag=msd_ag,msd_wp=msd_wp,n2=n2)

phi_stats<-get_phi_stats(sigmasq_a=sigmasq_a, msd_wp=msd_wp)


#row1 (dd.pops)
#ssd
ssd_ag
msd_ag
df_ag

#row2 (Error)
ssd_wp
msd_wp
df_wp

#row3 (Total)
ssd_total
msd_total
df_total

#Variance components:
#dd.pops
sigmasq_a
#Error
msd_wp

#Phi-statistics
#dd.pops.in.GLOBAL
phi_stats


#Variance coefficients
#a
n2
