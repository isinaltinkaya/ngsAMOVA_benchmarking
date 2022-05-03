import msprime
import tszip
import stdpopsim

from stdpopsim import models

import math
from itertools import product, combinations
import numpy as np
import pandas as pd

import sys
import os
import subprocess


POP_IND_ID=list( "_".join(comb) for comb in product(list(map('pop{}'.format,range(1,4))),list(map('ind{}'.format,range(1,11)))))
IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(POP_IND_ID,2))


# REP=range(1,10,1)
REP=1


SIMULATION_ID="sim1"
#SIMULATION_ID="sim2"+"_L"+config['region_length']

###################################################
# MSPRIME
###################################################

# REGION_LENGTH=range(5,70,5)
REGION_LENGTH=1



###################################################
# MSTOGLF 
###################################################

# average per site depth
# DEPTH=[100,20,10,5,2,1,0.5,0.1]
DEPTH=[1]

#0.1 to 5
#DEPTH=[100,20,10,0.01]+ [ x / pow(.1, -1) for x in range(1, 50 + 1) ]


###################################################


conda_env="/maps/projects/lundbeck/scratch/pfs488/AMOVA/env/simulation_env.yml"

# DEF_POPS= {
		# "1" : "1 10",
		# "2" : "11 20",
		# "3" : "21 30"
		# }



#
# {
    # "seed" : 12345,
    # "num_samples_per_population" : [20, 20, 0],
    # "replicates" : 1,
    # "species" : "HomSap",
    # "model" : "OutOfAfrica_3G09",
    # "genetic_map" : "HapmapII_GRCh37",
    # "chrm_list" : "chr22",
    # "mask_file" : "masks/HapmapII_GRCh37.mask.bed"
# }

# Number of replicates
n_reps=3


REP=[*range(n_reps)]

# Set seed
np.random.seed(42)
# SEED=np.random.randn(n_reps)
SEED=np.random.randint(1,2**32-1,n_reps)
# print(SEED)


species=stdpopsim.get_species("HomSap")


recombination_map_id="HapMapII_GRCh37"

#
# All chromosomes
# CONTIG=[c for c in species.genome.chromosomes]
# CONTIG="chr22"

# CONTIGID="all"
CONTIGID="chr22"

if CONTIGID != "all":
	CONTIG=CONTIGID
else:
	CONTIG=[c for c in species.genome.chromosomes]


# MODEL=["OutOfAfrica_3G09","OutOfAfrica_2T12","AncientEurasia_9K19","AshkSub_7G19"]
MODEL=["OutOfAfrica_3G09"]


# print(SEED[0])

samples_per_pop=[20,20,20]


		# seed=SEED

		# expand("simulations/{simid}-model_{model_id}-contig_{contig}-rep_{rep}.trees",
rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL)
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/msprime/trees/{sid}_{srlen}MB_rep{rep}.trees",
				# rep=REP,
				# sid=SIMULATION_ID,
				# srlen=REGION_LENGTH),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/msprime/ms/{sid}_{srlen}MB_rep{rep}.ms",
				# rep=REP,
				# sid=SIMULATION_ID,
				# srlen=REGION_LENGTH),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz",
				# rep=REP,
				# srlen=REGION_LENGTH,
				# sid=SIMULATION_ID,
				# depth=DEPTH),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz",
				# rep=REP,
				# srlen=REGION_LENGTH,
				# sid=SIMULATION_ID,
				# depth=DEPTH,
				# pop_id=DEF_POPS.keys()),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz",
				# rep=REP,
				# srlen=REGION_LENGTH,
				# sid=SIMULATION_ID,
				# depth=DEPTH,
				# pop_id=DEF_POPS.keys(),
				# ind_id=range(1,11)),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.saf.idx",
				# rep=REP,
				# srlen=REGION_LENGTH,
				# sid=SIMULATION_ID,
				# depth=DEPTH,
				# pop_id=DEF_POPS.keys(),
				# ind_id=range(1,11)),
		# expand("simulations/{sid}_{srlen}MB_rep{rep}/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs",
				# rep=REP,
				# srlen=REGION_LENGTH,
				# sid=SIMULATION_ID,
				# depth=DEPTH,
				# pop_ind_pair=IND_PAIRS),

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
		contig = species.get_contig(CONTIG, genetic_map=recombination_map_id)
		samples = model.get_samples(*samples_per_pop)
		engine = stdpopsim.get_engine("msprime")
		ts = engine.simulate(model,contig, samples, seed=SEED[int(wildcards.rep)])
		ts.dump(output[0])


	# log:
		# "simulations/logs/{simid}-{model}-{contig}-{seed}.trees",
	# conda:
		# conda_env
rule msprime2ms:
	input:
		rules.simulation.output
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/ms/{simid}-{model_id}-{contig}-rep{rep}.ms",
	run:
		ts=tskit.load(input[0])
		ms=open(output[0],"w")
		tskit.write_ms(ts,ms)
		
	# conda:
		# conda_env
	# shell:
		# """
		# python3.10 scripts/write_ms.py {input} {output}
		# """
#

#
#
#
#
#
# rule msprime_simulation:
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/msprime/trees/{sid}_{srlen}MB_rep{rep}.trees"
	# params:
		# rlen=lambda wildcards: int(float(wildcards.srlen)*1e6),
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/msprime/{sid}_{srlen}MB_rep{rep}.ms"
	# conda:
		# conda_env
	# shell:
		# """
		# ( python3.10 scripts/sim_demography.py {output} 10 10 10 {params.rlen} )2> {log}
		# """
#
# rule msprime2ms:
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/msprime/ms/{sid}_{srlen}MB_rep{rep}.ms"
	# input:
		# rules.msprime_simulation.output
	# conda:
		# conda_env
	# shell:
		# """
		# python3.10 scripts/write_ms.py {input} {output}
		# """
#
#
#
#
# rule ms_to_glf:
	# input:
		# "simulations/{sid}_{srlen}MB_rep{rep}/msprime/ms/{sid}_{srlen}MB_rep{rep}.ms"
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz"
	# params:
		# msToGlf="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/msToGlf",
		# outbase="simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}",
		# seed=42,
		# error_rate=0.002,
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz.log"
	# shell:
		# """
		# ({params.msToGlf} -in {input} -out {params.outbase} \
				# -singleOut 1 \
				# -depth {wildcards.depth} \
				# -err {params.error_rate}  \
				# -nSites 0 \
				# -seed {params.seed}) 2> {log}
		# """
#
# rule split_glf_to_pops:
	# input:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz"
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz",
	# params:
		# pop_size=30,
		# pop_range=lambda wildcards: DEF_POPS[wildcards.pop_id],
		# sim_glf="simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz",
		# outbase="simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}",
		# splitgl="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/splitgl",
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz.log",
	# shell:
		# """
		# ({params.splitgl} {params.sim_glf} {params.pop_size} {params.pop_range} > {params.outbase}.glf.gz ) 2> {log}
		# """
#
#
# rule split_glf_to_inds:
	# input:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz"
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz"
	# params:
		# pop_size=10,
		# ind_range="{ind_id} {ind_id}",
		# splitgl="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/splitgl",
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz.log"
	# shell:
		# """
		# ({params.splitgl} {input} {params.pop_size} {params.ind_range} > {output}) 2> {log}
		# """
#
#
#
# rule doSaf_perInd:
	# input:
		# "simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz"
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.saf.idx"
	# params:
		# outbase="simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}",
		# angsd="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/angsd/angsd",
		# ref="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/ms/ref.fa",
		# fai="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/ms/ref.fa.fai",
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz.log"
	# shell:
		# """
		# ({params.angsd} \
				# -glf {input} \
				# -doSaf 1 \
				# -nInd 1 \
				# -ref {params.ref} \
				# -fai {params.fai} \
				# -isSim 1 \
				# -out {params.outbase}) 2> {log}
#
		# """
#
#
# rule realSFS_2dsfs_ind2ind:
	# input:
		# ind1=lambda wildcards: "simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_"+ wildcards.pop_ind_pair.split("-")[0] + ".saf.idx",
		# ind2=lambda wildcards: "simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_"+ wildcards.pop_ind_pair.split("-")[1] + ".saf.idx",
	# output:
		# "simulations/{sid}_{srlen}MB_rep{rep}/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs",
	# params:
		# realSFS="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/angsd/misc/realSFS",
	# log:
		# "simulations/{sid}_{srlen}MB_rep{rep}/logs/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs.log",
	# shell:
		# """
		# ({params.realSFS} {input.ind1} {input.ind2} > {output} )2>{log}
		# """
