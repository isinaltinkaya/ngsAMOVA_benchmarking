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

# SIMULATION_ID="sim2"
# SIMULATION_ID="sim3"
SIMULATION_ID="simwf"

###################################################
# Version for running / not active develop dir
ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA/ngsAMOVA"


# Average per site depth
# DEPTH=[20,10,5,2,1,0.5,0.2,0.1,0.01]
DEPTH=[20,10,5,2,1,0.5,0.2,0.1]
# DEPTH=[20,10,5,2,1]
# DEPTH=[0.5,0.2,0.1]
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
		expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/maxIter1e3/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
				simid=SIMULATION_ID,
				model_id=MODEL,
				contig=CONTIGID,
				tole=TOLE,
				minInd=[2],
				rep=REP,
				dist=[1],
				depth=DEPTH),
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# minInd=[2],
				# rep=REP,
				# dist=[1],
				# depth=DEPTH),
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# minInd=[2],
				# rep=REP,
				# dist=[1],
				# depth=DEPTH),
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/rawSFS/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-dist{dist}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# rep=REP,
				# dist=[2]),

#
# rule run_ngsAMOVA_sfs_var_minInd2_100K_maxIter1000:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/masked_vcfgl_var/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}",
		# metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} -minInd {wildcards.minInd} -P 5 -maxIter 1000 -out {params.outprefix} -doAMOVA 3 -doDist {wildcards.dist}  -m {params.metadata} 2> {log}
		# """
#

rule run_ngsAMOVA_sfs_var_minInd2_100K_doDist1_sqDist0_dij_maxIter1000:
	input:
		"simulations/{simid}/model_{model_id}/contig_all/subsample_100K/masked_vcfgl_var/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}.bcf",
	output:
		"simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/maxIter1e3/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
	params:
		ngsAMOVA=ngsAMOVA,
		outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/maxIter1e3/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}",
		metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_sfs_masked_var/dist{dist}/maxIter1e3/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-tole{tole}-minInd{minInd}-dist{dist}.sfs.csv",
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} -minInd {wildcards.minInd} -P 5 -maxIter 1000 -out {params.outprefix} -doAMOVA 3 -doDist {wildcards.dist}  -sqDist 0 -m {params.metadata} 2> {log}
		"""



# #true GT SFS including all sites
# rule run_ngsAMOVA_sfs_masked_var_trueSFS:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/{simid}-{model_id}-contig_all_100K-rep{rep}.vcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/rawSFS/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-dist{dist}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/rawSFS/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-dist{dist}",
		# metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/rawSFS/dist{dist}/{simid}-{model_id}-contig_all_100K-rep{rep}-dist{dist}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -minInd 0 -out {params.outprefix} -doAMOVA 2 -doDist {wildcards.dist}  -m {params.metadata} 2> {log}
		# """
