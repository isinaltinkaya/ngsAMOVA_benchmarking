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


# POP_IND_ID=list( "_".join(comb) for comb in product(list(map('pop{}'.format,range(1,4))),list(map('ind{}'.format,range(1,11)))))
# IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(POP_IND_ID,2))



SIMULATION_ID="sim2"
# SIMULATION_ID="sim3"

#sim2: 50 ind per pop
#sim1: 10 ind per pop

###################################################
# MSTOGLF 
###################################################

# average per site depth
DEPTH=[100,20,10,5,2,1,0.5,0.1]

###################################################


conda_env="/maps/projects/lundbeck/scratch/pfs488/AMOVA/env/simulation_env.yml"

DEF_POPS= {
		"1" : "1 50",
		"2" : "51 100",
		"3" : "101 150"
		}



# Number of replicates
n_reps=200
REP=[*range(n_reps)]

# Set seed
np.random.seed(42)
SEED=np.random.randint(1,2**32-1,n_reps)


species=stdpopsim.get_species("HomSap")
recombination_map_id="HapMapII_GRCh37"

exclude_chr_list=['chrY']
CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]


MODEL=["OutOfAfrica_3G09"]



# samples_per_pop=[20,20,20]
samples_per_pop=[100]*3



rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL),
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/ms/{simid}-{model_id}-{contig}-rep{rep}.ms",
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/glf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
				# depth=DEPTH,
				# simid=SIMULATION_ID,
				# contig=CONTIGID,
				# rep=REP,
				# model_id=MODEL),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
				# depth=DEPTH,
				# simid=SIMULATION_ID,
				# contig=CONTIGID,
				# rep=REP,
				# model_id=MODEL),
#
rule generate_reference:
	output:
		"simulations/ref/ref.fa"
	shell:
		"""
		Rscript -e 'cat(">ref\n",paste(rep("A",2e9),sep="", collapse=""),"\n",sep="")' > {output}
		samtools faidx {output}
		"""


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


rule simulation:
	input:
		r_map_downloaded_flag
	output: 
		"simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees"
	run:
		model=species.get_demographic_model(wildcards.model_id)
		contig = species.get_contig(wildcards.contig, genetic_map=recombination_map_id)
		samples = model.get_samples(*samples_per_pop)
		engine = stdpopsim.get_engine("msprime")
		ts = engine.simulate(model,contig, samples, seed=SEED[int(wildcards.rep)])
		ts.dump(output[0])


rule msprime_to_ms:
	input:
		rules.simulation.output
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/ms/{simid}-{model_id}-{contig}-rep{rep}.ms",
	run:
		ts=tskit.load(input[0])
		ms=open(output[0],"w")
		tskit.write_ms(ts,ms)


# rule tree2svg:
	# input:
		# rules.simulation.output
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/tree/svg/{simid}-{model_id}-{contig}-rep{rep}.svg",
#
# rule tree_to_sequence_length:
	# input:
		# rules.simulation.output
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/trees/sequence_length/{simid}-{model_id}-{contig}-rep{rep}.sequence_length",
	# run:
		# ts=tskit.load(input[0])
		# of=open(output[0],"w")
		# print(int(ts.sequence_length),file=of)
#
# rule tree_to_num_sites:
	# input:
		# rules.simulation.output
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/trees/num_sites/{simid}-{model_id}-{contig}-rep{rep}.num_sites",
	# run:
		# ts=tskit.load(input[0])
		# of=open(output[0],"w")
		# print(ts.num_sites,file=of)

rule tree_to_:
	input:
		rules.simulation.output
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/testing/seg_sites/{simid}-{model_id}-{contig}-rep{rep}",
	run:
		ts=tskit.load(input[0])
		print(ts.sites())

# rule ms_to_var_glf:
	# input:
		# ms="simulations/{simid}/model_{model_id}/contig_{contig}/ms/{simid}-{model_id}-{contig}-rep{rep}.ms",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	# params:
		# msToGlf="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/misc/msToGlf",
		# outbase="simulations/{simid}/model_{model_id}/contig_{contig}/glf_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		# seed=42,
		# error_rate=0.002,
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/glf_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	# shell:
		# """
		# ({params.msToGlf} -in {input.ms} -out {params.outbase} \
				# -singleOut 1 \
				# -depth {wildcards.depth} \
				# -err {params.error_rate}  \
				# -regLen 0 \
				# -seed {params.seed}) 2> {log}
		# """


# rule ms_to_glf:
	# input:
		# ms="simulations/{simid}/model_{model_id}/contig_{contig}/ms/{simid}-{model_id}-{contig}-rep{rep}.ms",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/glf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	# params:
		# msToGlf="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/angsd/misc/msToGlf",
		# outbase="simulations/{simid}/model_{model_id}/contig_{contig}/glf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		# seed=42,
		# error_rate=0.002,
		# sequence_length=lambda wildcards: species.genome.get_chromosome(wildcards.contig).length,
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/glf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.glf.gz",
	# shell:
		# """
		# ({params.msToGlf} -in {input.ms} -out {params.outbase} \
				# -singleOut 1 \
				# -depth {wildcards.depth} \
				# -err {params.error_rate}  \
				# -regLen {params.sequence_length} \
				# -seed {params.seed}) 2> {log}
		# """
