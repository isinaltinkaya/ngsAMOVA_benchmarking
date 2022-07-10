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
DEPTH=[100,20,10,5,2,1,0.5,0.2,0.1,0.01]
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
CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]
# CONTIGID="chr21"
# CONTIGID="chr22"
# CONTIGID=CONTIGID[int(config["NC"])]

MODEL=["OutOfAfrica_3G09"]


n_samples_popi=50
ploidy=2
n_pops=3

samples_per_pop=[100]*3

# TOLE=[7,8,9,10,11]
TOLE=10



def get_nInd(samples_per_pop,ploidy,minInd_p,minInd,i):
	minInd[i]=(sum(samples_per_pop)/ploidy)*(minInd_p[i]/100)
	

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
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}_sfs.csv",
				simid=SIMULATION_ID,
				model_id=MODEL,
				minInd=MININD,
				contig=CONTIGID,
				tole=TOLE,
				rep=REP,
				depth=DEPTH),


# #shared_sfs= only use sites that contains data for all individuals
# #220702 update: now this requires -onlyShared 1
# rule run_ngsAMOVA_shared_sfs_var:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var_shared_sfs/tole{tole}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/masked_vcfgl_var_shared_sfs/tole{tole}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -onlyShared 1 -tole 1e-{wildcards.tole} > {output} 2> {log}
		# """



#
# # evaluate the inclusion of sites per ind pair basis
# # if a site is nonmissing for a specific individual pair, it is included
# # whereas with shared_sfs it should be non-missing for all individuals
# rule run_ngsAMOVA_sfs_var:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var_sfs/tole{tole}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/masked_vcfgl_var_sfs/tole{tole}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -onlyShared 0 -tole 1e-{wildcards.tole} > {output} 2> {log}
		# """




# set minind
# # test em tole vals
rule run_ngsAMOVA_sfs_var_minInd:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}_sfs.csv",
	params:
		ngsAMOVA=ngsAMOVA,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsAMOVA_sfs_masked_var/tole{tole}/minInd{minInd}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-tole{tole}-minInd{minInd}_sfs.csv",
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -tole 1e-{wildcards.tole} -minInd {wildcards.minInd} > {output} 2> {log}
		"""


