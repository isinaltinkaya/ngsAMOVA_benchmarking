###
library('ade4')

#################
# ADE4
#################

# EXAMPLE DATA
# ------------
#
# Source: https://rdrr.io/cran/ade4/man/humDNAm.html
#
# Details from ade4
# humDNAm: human mitochondrial DNA restriction data
#
# * Frequencies of haplotypes of mitochondrial DNA restriction data in ten populations all over the world.
# * Distances among the haplotypes.
# 
# distances
# is an object of class dist with 56 haplotypes. These distances are computed by counting the number of differences in restriction sites between two haplotypes.
# 
# samples
# is a data frame with 56 haplotypes, 10 abundance variables (populations). These variables give the haplotype abundance in a given population.
# 
# structures
# is a data frame with 10 populations, 1 variable (classification). This variable gives the name of the continent in which a given population is located.
#
#
# More details about the data from Arlequin:
#
# This example is about an AMOVA analysis based on an
# distance matrix given in an external file. The
# haplotypic composition of the haplotypes 
# is thus only used through the distance matrix 
#
# Data source: Excoffier, L., Smouse, P., and Quattro, J., 1992, 
# Analysis of molecular variance inferred from metric distances 
# among DNA haplotypes: Application to human mitochondrial DNA restriction data, 
# Genetics 131:479-491.

library(ade4)
data(humDNAm)

amovahum <- ade4::amova(humDNAm$samples, sqrt(humDNAm$distances), humDNAm$structures)

print(amovahum)
