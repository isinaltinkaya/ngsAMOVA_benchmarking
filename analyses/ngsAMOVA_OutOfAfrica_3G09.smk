import msprime
import tszip
import tskit
import stdpopsim

from stdpopsim import models

import math
from itertools import product, combinations
import numpy as np
import pandas as pd
import sys
import os
import subprocess

SIMULATION_ID="sim2"

###################################################
# Version for running / not active develop dir
ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA/ngsAMOVA"


# Average per site depth
# DEPTH=[100,20,10,5,2,1,0.5,0.2,0.1,0.01]
DEPTH=[20,10,5,2,1]
# DEPTH=[20,10,5,2,1,0.5,0.2,0.1]
# DEPTH=DEPTH[int(config["ND"])]

# Number of replicates
n_reps=200

REP=[*range(n_reps)]
# REP=REP[int(config["NN"])]
REP=REP[:20]


# Set seed
np.random.seed(42)
SEED=np.random.randint(1,2**32-1,n_reps)
# print(SEED)


species=stdpopsim.get_species("HomSap")


# exclude_chr_list=['chrY']
exclude_chr_list=['chrX','chrY']
# CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]
# CONTIGID="chr21"
CONTIGID="chr22"
# CONTIGID=["chr20","chr21","chr22"]
# CONTIGID=CONTIGID[int(config["NC"])]

MODEL=["OutOfAfrica_3G09"]


n_samples_popi=50
ploidy=2
n_pops=3

samples_per_pop=[100]*3

# TOLE=[9,11]
TOLE=10


# Include sites if exists for both in the pair for each pair
PAIRWISE=[2]

def get_nInd(samples_per_pop,ploidy,minInd_p,minInd,i):
	minInd[i]=int((sum(samples_per_pop)/ploidy)*(minInd_p[i]/100))


# MININD_P=[50,80,90]
# MININD_P=[90,100]
MININD_P=[100]
MININD=[0]*len(MININD_P)
list(map(lambda i:get_nInd(samples_per_pop,ploidy,MININD_P,MININD,i),range(0,len(MININD))))
# print(MININD)
if 0 in MININD:
	print("MININD contains 0")

###################################################



rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
				simid=SIMULATION_ID,
				model_id=MODEL,
				contig=CONTIGID,
				tole=TOLE,
				minInd=[2],
				rep=REP,
				depth=DEPTH),



rule run_ngsAMOVA_sfs_var_minInd2:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	params:
		ngsAMOVA=ngsAMOVA,
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} \
				-minInd {wildcards.minInd} \
				-out {params.outprefix} -doAMOVA 3 2> {log}
		"""



rule run_ngsAMOVA_sfs_var_minInd2_2122:
	input:
		"simulations/{simid}/model_{model_id}/contig_chr2122/masked_vcfgl_var_concat/{simid}-{model_id}-chr2122-rep{rep}-d{depth}.bcf",
	output:
		"simulations/{simid}/model_{model_id}/contig_chr2122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr2122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	params:
		ngsAMOVA=ngsAMOVA,
		outprefix="simulations/{simid}/model_{model_id}/contig_chr2122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr2122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_chr2122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr2122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} \
				-minInd {wildcards.minInd} \
				-out {params.outprefix} -doAMOVA 3 2> {log}
		"""






rule run_ngsAMOVA_sfs_var_minInd2_202122:
	input:
		"simulations/{simid}/model_{model_id}/contig_chr202122/masked_vcfgl_var_concat/{simid}-{model_id}-chr202122-rep{rep}-d{depth}.bcf",
	output:
		"simulations/{simid}/model_{model_id}/contig_chr202122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr202122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	params:
		ngsAMOVA=ngsAMOVA,
		outprefix="simulations/{simid}/model_{model_id}/contig_chr202122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr202122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_chr202122/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-chr202122-rep{rep}-d{depth}-tole{tole}-minInd{minInd}.sfs.csv",
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} \
				-minInd {wildcards.minInd} \
				-out {params.outprefix} -doAMOVA 3 2> {log}
		"""





#
# rule run_ngsAMOVA_sfs_var_minInd_2122_trueSFS:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_chr2122/masked_vcfgl_var_concat/{simid}-{model_id}-chr2122-rep{rep}-d100.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_chr2122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr2122-rep{rep}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_chr2122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr2122-rep{rep}",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_chr2122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr2122-rep{rep}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 \
				# -minInd 150 \
				# -doTest 0 -out {params.outprefix} -doAMOVA 2 2> {log}
		# """
#


#
# rule run_ngsAMOVA_sfs_var_minInd_202122_trueSFS:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_chr202122/masked_vcfgl_var_concat/{simid}-{model_id}-chr202122-rep{rep}-d100.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_chr202122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr202122-rep{rep}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_chr202122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr202122-rep{rep}",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_chr202122/ngsAMOVA_trueSFS_masked_var/{simid}-{model_id}-chr202122-rep{rep}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 \
				# -minInd 150 \
				# -doTest 0 -out {params.outprefix} -doAMOVA 2 2> {log}
		# """
#
#
#
#
# #true GT SFS including all sites
# rule run_ngsAMOVA_sfs_masked_var_trueSFS:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_trueSFSgt_masked_var/{simid}-{model_id}-{contig}-rep{rep}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_trueSFSgt_masked_var/{simid}-{model_id}-{contig}-rep{rep}",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsAMOVA_trueSFSgt_masked_var/{simid}-{model_id}-{contig}-rep{rep}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -doAMOVA 2 -out {params.outprefix} -doTest 0 2> {log}
		# """




