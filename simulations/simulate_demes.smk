
import math
from itertools import product, combinations
import numpy as np
import pandas as pd
import sys
import os
import subprocess

import tskit,msprime
import demes, demesdraw


## Disable warnings
# There were bunch of warnings due to msprime version issues
import warnings
warnings.filterwarnings('ignore')


###################################################
# CONFIG


vcfgl="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"

# SIMULATION_ID="sim3"
SIMULATION_ID="simwf"

###################################################
# DEPTH 

# average per site depth
# DEPTH=[100,20,10,5,2,1,0.5,0.2,0.1,0.01]
DEPTH=[20,10,5,2,1,0.5,0.2,0.1,0.01]

###################################################


# Number of replicates
n_reps=200
REP=[*range(n_reps)]
REP=REP[:20]



###################################################
# define populations
#
ploidy=2

# define population with diploid individuals
# <pop_id> : <n_individuals_in_pop>
DEF_POPS={"B1C1":10, 
		"B1C2":10, 
		"B2C1":10,
		"B2C2":10}


#note that msprime>1 and stdpopsim_dev version
#assumes ploidy=2 by default

haplo_list=[]
indv_names=[]


for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value*ploidy)
	for ind in range(value):

		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")

IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(indv_names,2))




mutation_rate=1.29e-08
recombination_rate=1.14856e-08

# END CONFIG
###################################################






rule all:
	input:
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				# simid="sim_models",
				# model_id=["model1","model2"],
				# contig=[1],
				# rep=REP),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
				# simid="sim_models",
				# model_id=["model1","model2"],
				# contig=[1],
				# rep=REP),
#
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid="sim_models",
				# model_id=["model1","model2"],
				# contig=[1],
				# depth=DEPTH,
				# rep=REP),
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				simid="sim_models",
				model_id=["model1","model2"],
				contig=[10],
				rep=REP),
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
				simid="sim_models",
				model_id=["model1","model2"],
				contig=[10],
				rep=REP),

		expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				simid="sim_models",
				model_id=["model1","model2"],
				contig=[10],
				depth=DEPTH,
				rep=REP),







#############################################
### Simulation


rule simulation:
	output: 
		"simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	params:
		seedfile="simulations/{simid}/model_{model_id}/contig_{contig}/trees/.seed.{simid}-{model_id}-{contig}-rep{rep}.trees",
	run:
		model=demes.load(str(wildcards.model_id)+".yaml")
		demography=msprime.Demography.from_demes(model)
		samples=DEF_POPS
		seedval=np.random.randint(2**20) + (100*int(wildcards.rep))  + int(wildcards.contig)
		with open(params.seedfile,"w") as seedout:
			print("Contig :"+str(wildcards.contig)+", rep:"+str(wildcards.rep),file=seedout)
			print("Seed value: "+str(seedval), file=seedout)
		sequence_length=int(wildcards.contig)*100000
		ts = msprime.sim_ancestry([msprime.SampleSet(n, population=p) for p,n in samples.items()],
				demography=demography,
				sequence_length=sequence_length,
				random_seed=seedval,
				recombination_rate=recombination_rate)
		mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seedval)
		with open(output[0],"w") as tsout:
			mts.dump(tsout)



## Convert tree sequence to VCF file
# Using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
	input:
		rules.simulation.output,
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	run:
		ts=tskit.load(input[0])
		with open(output[0],"w") as vcfout:
			ts.write_vcf(output=vcfout,
					contig_id=str(wildcards.contig),
					position_transform="legacy",
					individual_names=indv_names)
					# ploidy=ploidy,




## Using only variable sites
rule vcf_to_vcfgl_var:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	params:
		vcfgl=vcfgl,
		prefix="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		mode="b",
		error_rate=0.002,
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	shell:
		"""
		RNG_SEED=$(( {wildcards.rep} + 1 ))

		({params.vcfgl} -in {input} -out {params.prefix} \
				-depth {wildcards.depth} \
				-err {params.error_rate}  \
				-mode {params.mode} \
				-explode 0 \
				-seed ${{RNG_SEED}} ) 2> {log}
		"""


