#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Loading required packages for the phylogenetic analysis

# library(seqinr)
library(ape)
# library(msa) # this package is available through the Bioconductor packages. If you do not have the BiocManager package installed, you could install this using the following lines of code... then load it using library(msa)
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
	# install.packages("BiocManager")
#
# # BiocManager::install("adagenet")
# # BiocManager::install("phangorn")
#
# # library(adagenet)
# library(phangorn)
#
# d<- read.phyDat(faFile,format="fasta")

# library(msa)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Supply first arg input", call.=FALSE)
} else if (length(args)==1) {
  stop("Supply second arg output prefix", call.=FALSE)
}

faFile<-args[1]
outprefix<-args[2]
outnewick<-paste0(outprefix,".newick")
outplot<-paste0(outprefix,".png")

dna <- read.dna(file = faFile ,format = "fasta")
D<-dist.dna(dna)

# D <- dist.dna(dna, model = "TN93")
# temp <- as.data.frame(as.matrix(D))
# table.paint(temp, cleg = 0, clabel.row = 0.5, clabel.col = 0.5)


# tre <- bionj(D)
tre<-nj(D)
# tre<-ladderize(tre)

write.tree(tre,file=outnewick)

png(outplot)
plot(tre)
dev.off()






