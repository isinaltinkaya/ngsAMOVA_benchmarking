args = commandArgs(trailingOnly=TRUE)




df<-read.csv(args[1],header = F)
df2<-read.csv(args[2],header = F)

model<-args[3]
contig<-args[4]
rep<-args[5]
depth<-args[6]




df[df==".."]<-NA
df2[df2==".."]<-NA

if(nrow(df)==0){
    error("df is empty")
}
if(nrow(df2)==0){
    error("df2 is empty")
}


# print(nrow(df))
# print(nrow(df2))
# if(nrow(df)!=nrow(df2)){
#     stop("different number of rows")
# }

nr<-nrow(df)
# print(nr)
# out<-data.frame(matrix(NA,nr,7))
nt<-0
nf<-0
nn<-0

if(nrow(df)>nrow(df2)){
    stop("df1 has more rows than df2")
}

nr_true<-nrow(df2)

for (i in 1:nr_true){
    if(sum(df$V1==df2[i,]$V1)>0){
        if(sum(df$V1==df2[i,]$V1)>0){
            if(length(df2[i,][-1])!=length(df[df$V1==df2[i,]$V1,][-1])){
                stop(paste0("different number of columns in row ",i))
            }
            cmp<-df2[i,][-1]==df[df$V1==df2[i,]$V1,][-1]
            nt<-nt+length(cmp[cmp==TRUE & !is.na(cmp)])
            nf<-nf+length(cmp[cmp==FALSE& !is.na(cmp)])
            nn<-nn+length(cmp[is.na(cmp)])
        }
    }
}

rt<-nt/(nt+nf+nn)
rf<-nf/(nt+nf+nn)
rn<-nn/(nt+nf+nn)
cat(paste0(rt,",",rf,",",rn,",",nt,",",nf,",",nn,",",nt+nf+nn,",",model,",",contig,",",rep,",",depth))
