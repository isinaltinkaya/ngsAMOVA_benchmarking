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


###################################################
# FUNCTIONS

######
# from stdpopsim_dev

def read_bed(mask_fpath, chrom):
	lines = np.loadtxt(
			mask_fpath,
			dtype={"names": ("chrom", "left", "right"), "formats": (object, int, int)},
			delimiter="\t",
			usecols=(0, 1, 2),
			)
	in_chrom = lines["chrom"] == f"{chrom}"
	lines = lines.compress(in_chrom)
	intervals = np.array([lines["left"], lines["right"]]).T
	return intervals


def mask_tree_sequence(ts, mask_intervals):
	ts = ts.delete_intervals(mask_intervals)
	return ts


######


###################################################
# CONFIG


vcfgl="/maps/projects/lundbeck/scratch/pfs488/vcf-gl/vcfgl"


SIMULATION_ID="sim2"

###################################################

# average per site depth
DEPTH=[100,20,10,5,2,1,0.5,0.2,0.1,0.01]
# DEPTH=DEPTH[int(config["NN"])]

###################################################


# Number of replicates
n_reps=200
REP=[*range(n_reps)]
REP=REP[:10]

# Set seed
ULTIMATE_ANSWER=42
np.random.seed(ULTIMATE_ANSWER)
SEED=np.random.randint(1,2**32-1,n_reps)


###################################################
# stdpopsim config
#


species=stdpopsim.get_species("HomSap")
recombination_map_id="HapMapII_GRCh37"

exclude_chr_list=['chrY']
CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]
# CONTIGID="chr22"
CONTIGID=CONTIGID[int(config['NC'])]

MODEL=["AncientEurasia_9K19"]



###################################################
# define populations
#
ploidy=2

# define population with Nploid individuals
# <pop_id> : <n_individuals_in_pop>
DEF_POPS={
		1 : 50,
		2 : 50,
		4 : 50,
		5 : 50,
		6 : 50,
		7 : 50,
		8 : 50,
		}


#note that msprime>1 and stdpopsim_dev version
#assumes ploidy=2 by default
#so this may be subject to change in future 

haplo_list=[]
indv_names=[]

for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value*ploidy)
	for ind in range(value):

		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")


# END CONFIG
###################################################



rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL),

		expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL),

		expand("simulations/{simid}/model_{model_id}/contig_{contig}/bcf_gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL),



		r_map_downloaded_flag= "."+recombination_map_id+"_downloaded"
# Function source:
# https://github.com/popsim-consortium/analysis/blob/86e3a884309c35f3962971a0ca8abbd04edda626/two_population_analysis/Snakefile
rule download_genetic_map:
	output: r_map_downloaded_flag
	run:
		if recombination_map_id is not None:
			r_map = species.get_genetic_map(recombination_map_id)
			if not r_map.is_cached():
				r_map.download()
			with open(output[0], "w") as f:
				print("File to indicate genetic map has been downloaded", file=f)


rule download_mask_bed:
	output:
		"simulations/resources/genetic_map_mask/HapmapII_GRCh37.mask.bed"
	shell:
		"""
		wget https://raw.githubusercontent.com/popsim-consortium/analysis/master/n_t/masks/HapmapII_GRCh37.mask.bed -O {output}
		"""


rule simulation:
	input:
		r_map_downloaded_flag,
	output: 
		"simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	run:
		model=species.get_demographic_model(wildcards.model_id)
		contig = species.get_contig(wildcards.contig, genetic_map=recombination_map_id)
		samples = model.get_samples(*haplo_list)
		engine = stdpopsim.get_engine("msprime")
		ts = engine.simulate(model,contig, samples, seed=SEED[int(wildcards.rep)])
		with open(output[0],"w") as tsout:
			ts.dump(tsout)


# using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
	input:
		rules.simulation.output,
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	run:
		ts=tskit.load(input[0])
		with open(output[0],"w") as vcfout:
			ts.write_vcf(output=vcfout,
					ploidy=ploidy,
					contig_id=str(wildcards.contig),
					position_transform="legacy",
					individual_names=indv_names)


			# ## masking tree for tree-based statistics
# rule tree_to_masked:
	# input:
		# rules.simulation.output,
	# output:
		# ""
	# run:
		# ts=tskit.load(input[0])
		# masked_ts=mask_tree_sequence(ts,read_bed(str(input["mask"]),chrom=str(wildcards.contig)))
			# masked_ts.dump(tsout)


rule vcf_to_gl:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/bcf_gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	params:
		vcfgl=vcfgl,
		prefix=lambda wildcards, output: output[0].rsplit('.',1)[0],
		seed=lambda wildcards: SEED[int(wildcards.rep)],
		# bcf output
		mode="b",
		error_rate=0.002,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/bcf_gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	shell:
		"""
		({params.vcfgl} -in {input} -out {params.prefix} \
				-depth {wildcards.depth} \
				-err {params.error_rate}  \
				-isSim 1 \
				-mode {params.mode} \
				-seed {params.seed}) 2> {log}
		"""
