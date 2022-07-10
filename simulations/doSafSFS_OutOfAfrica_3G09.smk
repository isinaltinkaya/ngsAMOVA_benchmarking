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


## Disable warnings
# There were bunch of warnings due to msprime version issues
import warnings
warnings.filterwarnings('ignore')



# SIMULATION_ID="sim1"
SIMULATION_ID="sim2"



angsd="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/angsd"
splitgl="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/misc/splitgl"
realSFS="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/misc/realSFS"



###################################################
# MSTOGLF 
###################################################

# average per site depth
DEPTH=[100,20,10,5,2,1,0.5,0.1]
DEPTH=DEPTH[int(config["ND"])]

###################################################



DEF_POPS= {
		"1" : "1 50",
		"2" : "51 100",
		"3" : "101 150"
		}

# Number of replicates
n_reps=200


REP=[*range(n_reps)]
REP=REP[int(config["NN"])]
# REP=REP[:20]


# Set seed
np.random.seed(42)
SEED=np.random.randint(1,2**32-1,n_reps)
# print(SEED)


species=stdpopsim.get_species("HomSap")


recombination_map_id="HapMapII_GRCh37"


# exclude_chr_list=['chrY']
# CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]
CONTIGID="chr22"
# CONTIGID=CONTIGID[int(config["NC"])]

MODEL=["OutOfAfrica_3G09"]



n_samples_popi=50
ploidy=2
n_pops=3

samples_per_pop=[100]*3
# samples_per_pop=[n_samples_popi*ploidy]*3


POP_IND_ID=list( "_".join(comb) for comb in product(list(map('pop{}'.format,range(1,n_pops+1))),list(map('ind{}'.format,range(1,n_samples_popi+1)))))
# POP_IND_ID=list( "_".join(comb) for comb in product(list(map('pop{}'.format,range(1,3))),list(map('ind{}'.format,range(1,2)))))
IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(POP_IND_ID,2))
# print(POP_IND_ID)



rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				pop_id=DEF_POPS.keys(),
				model_id=MODEL),
#
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.glf.gz",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				pop_id=DEF_POPS.keys(),
				ind_id=range(1,n_samples_popi+1),
				model_id=MODEL),
#
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.saf.idx",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				pop_id=DEF_POPS.keys(),
				ind_id=range(1,n_samples_popi+1),
				model_id=MODEL),

		expand("simulations/{simid}/model_{model_id}/contig_{contig}/sfs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_indpair_{pop_ind_pair}.sfs",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL,
				pop_ind_pair=IND_PAIRS),


rule split_glf_to_pops:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
	params:
		pop_size=30,
		pop_range=lambda wildcards: DEF_POPS[wildcards.pop_id],
		splitgl="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd_doSaf_fixed/angsd/misc/splitgl"
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/glf/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
	shell:
		"""
		({params.splitgl} {input} {params.pop_size} {params.pop_range} > {output} ) 2> {log}
		"""


rule split_glf_to_pops_var:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
	params:
		tot_size=n_samples_popi*n_pops,
		pop_range=lambda wildcards: DEF_POPS[wildcards.pop_id],
		splitgl=splitgl,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/glf/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
	shell:
		"""
		({params.splitgl} {input} {params.tot_size} {params.pop_range} > {output} ) 2> {log}
		"""


rule split_glf_to_inds:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/pops/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}.glf.gz",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.glf.gz",
	params:
		pop_size=n_samples_popi,
		ind_range="{ind_id} {ind_id}",
		splitgl=splitgl,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/glf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.glf.gz",
	shell:
		"""
		({params.splitgl} {input} {params.pop_size} {params.ind_range} > {output}) 2> {log}
		"""

rule doSaf_perInd_var:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.glf.gz",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.saf.idx",
	params:
		outbase="simulations/{simid}/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}",
		angsd=angsd,
		ref="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/ref/ref.fa",
		fai="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/ref/ref.fa.fai",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.saf.idx",
	shell:
		"""
		({params.angsd} \
				-glf {input} \
				-doSaf 1 \
				-anc {params.ref} \
				-isSim 1 \
				-out {params.outbase}) 2> {log}

		"""

rule realSFS_2dsfs_ind2ind:
	input:
		ind1=lambda wildcards: "simulations/{simid}/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-"+ wildcards.pop_ind_pair.split("-")[0] + ".saf.idx",
		ind2=lambda wildcards: "simulations/{simid}/model_{model_id}/contig_{contig}/saf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-"+ wildcards.pop_ind_pair.split("-")[1] + ".saf.idx",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/sfs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_indpair_{pop_ind_pair}.sfs",
	params:
		realSFS=realSFS,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/sfs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_indpair_{pop_ind_pair}.sfs",
	shell:
		"""
		({params.realSFS} {input.ind1} {input.ind2} > {output} )2>{log}
		"""




rule doSaf_perInd_var_rmap:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.glf.gz",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/saf_var_rmap/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.saf.idx",
	params:
		outbase="simulations/{simid}/model_{model_id}/contig_{contig}/saf_var_rmap/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}",
		angsd=angsd,
		ref="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/ref/ref.fa",
		fai="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/ref/ref.fa.fai",
		rf="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/resources/genetic_map_mask/inclusive_chr22.bed",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/saf_var_rmap/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-pop{pop_id}_ind{ind_id}.saf.idx",
	shell:
		"""
		({params.angsd} \
				-vcf-gl {input} \
				-doSaf 1 \
				-anc {params.ref} \
				-isSim 1 \
				-rf {params.rf} \
				-out {params.outbase}) 2> {log}

		"""


rule realSFS_2dsfs_ind2ind_rmap:
	input:
		ind1=lambda wildcards: "simulations/{simid}/model_{model_id}/contig_{contig}/saf_var_rmap/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-"+ wildcards.pop_ind_pair.split("-")[0] + ".saf.idx",
		ind2=lambda wildcards: "simulations/{simid}/model_{model_id}/contig_{contig}/saf_var_rmap/inds/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-"+ wildcards.pop_ind_pair.split("-")[1] + ".saf.idx",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/sfs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_indpair_{pop_ind_pair}.sfs",
	params:
		realSFS=realSFS,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/sfs/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_indpair_{pop_ind_pair}.sfs",
	shell:
		"""
		({params.realSFS} {input.ind1} {input.ind2} > {output} )2>{log}
		"""


