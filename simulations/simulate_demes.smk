
import math
from itertools import product, combinations
import numpy as np
import pandas as pd
import sys
import os
import subprocess

import tskit,msprime
import demes, demesdraw

import csv


## Disable warnings
# There were bunch of warnings due to msprime version issues
import warnings
warnings.filterwarnings('ignore')


###############################################################################
# BEGIN CONFIG
###############################################################################

vcfgl="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"
ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA"

SIMULATION_ID="sim_demes_v2"

###############################################################################
# DEPTH 
# average per site depth
# DEPTH=[20,10,5,2,1,0.5,0.2,0.1,0.01]
# DEPTH=[20,10,5,2,1,0.5,0.2,0.1]
DEPTH=[10,5,2,1,0.5,0.2,0.1]


###############################################################################
# Number of replicates
n_reps=200
REP=[*range(n_reps)]
REP=REP[:20]



###############################################################################
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

# max number of iterations for EM algorithm
MAXIT=[500]

AMOVASET={"noAMOVA":-1}

# CONTIGS=[1,2,10,50,100,200,1000]
CONTIGS=[1,2,10,50,100]


MODELS=["model1","model2"]


###############################################################################
# END CONFIG
###############################################################################


rule all:
	input:
		"simulations/"+SIMULATION_ID+"/model_all/distanceMatrixList_amovaResultsList.txt",
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/distanceMatrices/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}_emIter{it}.distance_matrix.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP,
				# it=[i for i in range(1,501,1)],
				# tole=[5,10],
				# allow_missing=True),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv_emIter1.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP,
				# tole=[5,10]),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP,
				# tole=[5,10]),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
		#221213 with per iter log
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=10,
				# rep=REP),
		#221214 tole10 + tole5
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=[10,5],
				# rep=REP),
		#221214 benchmark time
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=5,
				# rep=REP),
# TODO CHECK THIS
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/amova_results/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=10,
				# rep=REP),
# COALESCENCE SIMULATION ##########################
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP),
# MISC ############################################
		# expand("simulations/{simid}/resources/{contig}.fa",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
		# expand("simulations/{simid}/resources/{contig}.fa.fai",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
		# expand("simulations/{simid}/resources/{contig}.sites",
				# simid=SIMULATION_ID,
				# contig=CONTIGS)
		# expand("simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
				# simid=SIMULATION_ID,
				# model_id=MODELS),
# TODO CHECK BELOW ##################################
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# it=MAXIT,
				# depth=DEPTH,
				# doamv=["noAMOVA"],
				# nthr=THRSET,
				# rep=REP),
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/iterLike/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# it=MAXIT,
				# depth=DEPTH,
				# doamv=["noAMOVA"],
				# nthr=THRSET,
				# rep=REP),





# ###############################################################################
# ### Simulation
#
#
# # rule generate_reference:
	# # output:
		# # "simulations/{simid}/resources/{contig}.fa"
	# # params:
		# # contig_length=lambda wildcards: int(int(wildcards.contig)*1e6)
	# # shell:
		# # """
		# # dd if=/dev/zero bs=1 count={params.contig_length} | tr '\\0' 'A' > {output}
		# # sed -i -e '$a\\' {output}
		# # sed -i -e '1s/^/>{wildcards.contig}\\n/' {output}
		# # """
# #
# # rule index_reference:
	# # input:
		# # "simulations/{simid}/resources/{contig}.fa"
	# # output:
		# # "simulations/{simid}/resources/{contig}.fa.fai"
	# # shell:
		# # """
		# # samtools faidx {input}
		# # """
# #
# # rule generate_sites:
	# # input:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# # output:
	# # shell:
		# # """
		# # bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n'  sim_demes_v2-model2-1-rep0-d1.bcf > sites.txt
		# # """
# #
# # rule generate_sites_contig:
	# # output:
		# # "simulations/{simid}/resources/{contig}.sites"
	# # params:
		# # contig_length=lambda wildcards: int(int(wildcards.contig)*1e6)
	# # shell:
		# # """
		# # awk -v X={wildcards.contig} -v XN={params.contig_length} 'BEGIN{{for(c=1;c<XN;c++) printf "%s\\t%d\\tA\\tC\\n",X,c}}' > {output}
		# # {ANGSD} sites index {output}
		# # """
# #
# #
# # ###############################################################################
# # # ./angsd
# # # -doMajorMinor 3 # Use pre-specified major and minor
# # # -sites sites.txt
# # # -doMaf 1  # Calculate frequencies using fixed major and minor
# # # -doGeno 31 # 1 + 2 + 4 + 8 + 16 [below]
# # # -doPost 1 # Estimate the posterior genotype probability based on allele frequency as prior
# # # -vcf-gl {input.vcf}
# # ###############################################################################
# # ## doGeno
# # # 1: print out major minor
# # # 2: print the called genotype as -1,0,1,2 (count of minor)
# # # 4: print the called genotype as AA, AC, AG, ...
# # # 8: print all 3 posts (major,major),(major,minor),(minor,minor)
# # # 16: print the posterior of the called genotype
# # ###############################################################################
# # rule angsd_call_genotypes:
	# # input:
		# # vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# # sites="simulations/{simid}/resources/{contig}.sites"
	# # output:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/simid}-{model_id}-{contig}-rep{rep}-d{depth}.geno.gz"
	# # params:
		# # prefix="simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	# # shell:
		# # """
		# # {ANGSD} -doMajorMinor 3 -doMaf 1 -doGeno 31 -doPost 1 -vcf-gl {input.vcf} -sites {input.sites} -o {params.prefix}
		# # """
# #
# #
# #
# # rule simulation:
	# # output:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	# # params:
		# # seedfile="simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}_simulation_params.txt",
	# # run:
		# # model=demes.load("demes_conf/"+str(wildcards.model_id)+".yaml")
		# # demography=msprime.Demography.from_demes(model)
		# # samples=DEF_POPS
		# # sequence_length=int(wildcards.contig)*1e6
		# # seedval=np.random.randint(2**20) + (100*int(wildcards.rep))  + int(wildcards.contig)
		# # with open(params.seedfile,"w") as seedout:
			# # print("Simulation_ID:"+str(wildcards.simid)+
					# # "\nReplicate:"+str(wildcards.rep)+
					# # "\nSeed:"+str(seedval)+
					# # "\nContig:"+str(wildcards.contig)+
					# # "\nSequence_length:"+str(sequence_length)+
					# # "\nMutation_rate:"+str(mutation_rate)+
					# # "\nRecombination_rate:"+str(recombination_rate),
					# # file=seedout)
			# # ts = msprime.sim_ancestry([msprime.SampleSet(n, population=p) for p,n in samples.items()],
					# # demography=demography,
					# # sequence_length=sequence_length,
					# # random_seed=seedval,
					# # recombination_rate=recombination_rate)
			# # mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seedval, model="binary")
		# # with open(output[0],"w") as tsout:
			# # mts.dump(tsout)
# #
# #
# #
# # ## Convert tree sequence to VCF file
# # # Using legacy format to avoid multiple instances of sites
# # rule tree_to_vcf:
	# # input:
		# # rules.simulation.output,
	# # output:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	# # run:
		# # ts=tskit.load(input[0])
		# # with open(output[0],"w") as vcfout:
			# # ts.write_vcf(output=vcfout,
					# # contig_id=str(wildcards.contig),
					# # position_transform="legacy",
					# # individual_names=indv_names)
# #
# #
# #
# #
			# # ## Using only variable sites
# # rule vcf_to_vcfgl_var:
	# # input:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
	# # output:
		# # "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	# # params:
		# # vcfgl=vcfgl,
		# # prefix="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		# # mode="b",
		# # error_rate=0.002,
	# # log:
		# # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	# # shell:
		# # """
		# # RNG_SEED=$(( {wildcards.rep} + 1 ))
# #
		# # ({params.vcfgl} -in {input} -out {params.prefix} \
				# # -depth {wildcards.depth} \
				# # -err {params.error_rate}  \
				# # -mode {params.mode} \
				# # -explode 0 \
				# # -seed ${{RNG_SEED}} ) 2> {log}
		# # """
# #
# #
# #
# #
# #
# # #
# # # rule benchmark_run_ngsAMOVA_sfs_var_maxIterSet:
	# # # input:
		# # # "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# # # output:
		# # # "simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
	# # # params:
		# # # outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}",
		# # # iterout="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}_perIterD.csv",
		# # # doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
	# # # log:
		# # # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/doAmova/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
	# # # shell:
		# # # """
		# # # ( {ngsAMOVA} -in {input} -isSim 1 -P {wildcards.nthr} -mEmIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1 > {params.iterout} )  2> {log}
		# # # """
# # #
# # #
# # # rule benchmark_run_ngsAMOVA_sfs_var_maxIterSet_printLike:
	# # # input:
		# # # "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# # # output:
		# # # "simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/iterLike/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
	# # # params:
		# # # outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/iterLike/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}",
		# # # iterout="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova/iterLike/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}_perIterLike.csv",
		# # # doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
	# # # log:
		# # # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/doAmova/iterLike/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-i{it}-t{nthr}-{doamv}.sfs.csv",
	# # # shell:
		# # # """
		# # # ( {ngsAMOVA} -in {input} -isSim 1 -P {wildcards.nthr} -mEmIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1 > {params.iterout} )  2> {log}
		# # # """
# #
# #
# # #################################
# # # individual names are in form of
# # # full_ind_id=popid_indid
# # #
# # # metadata should be in form of
# # # full_ind_id\tpopid
# # #################################
# # rule prepare_metadata:
	# # input:
		# # expand("simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# # simid=SIMULATION_ID,
				# # model_id=MODELS[0],
				# # contig=CONTIGS[0],
				# # depth=DEPTH[0],
				# # rep=REP[0]),
	# # output:
		# # "simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# # shell:
		# # """
		# # paste <(bcftools query -l {input}) <(bcftools query -l {input}|cut -d_ -f1) > {output}
		# # """
# #
# #
# #
# # # rule run_ngsAMOVA_sfs_var_gtgl_forAmovaResults:
	# # # input:
		# # # bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# # # metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# # # output:
		# # # amv="simulations/{simid}/model_{model_id}/contig_{contig}/amova_results/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
		# # # iterlog="simulations/{simid}/model_{model_id}/contig_{contig}/amova_results/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.csv",
	# # # params:
		# # # outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/amova_results/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}",
	# # # threads:
		# # # 1
	# # # log:
		# # # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/amova_results/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
	# # # shell:
		# # # """
		# # # (
		# # # {ngsAMOVA} -in {input.bcf} -isSim 1 -P {threads} -mEmIter 500 -out {params.outprefix} -tole 1e-{wildcards.tole} -doAMOVA 3 -doDist 1 -m {input.metadata} > {output.iterlog}
		# # # )  2> {log}
		# # # """
#
# #
# # rule run_ngsAMOVA_sfs_var_gtgl_forAmovaResults_v2:
	# # input:
		# # bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# # metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# # output:
		# # amv="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
	# # params:
		# # outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}",
	# # threads:
		# # 1
	# # log:
		# # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
	# # shell:
		# # """
		# # (
		# # {ngsAMOVA} -in {input.bcf} -isSim 1 -P {threads} -mEmIter 500 -out {params.outprefix} -tole 1e-{wildcards.tole} -doAMOVA 3 -doDist 1 -m {input.metadata} -printMatrix 1  -sqDist 1
		# # )  2> {log}
		# # """
#
#
# # # rule run_ngsAMOVA_sfs_var_gtgl_forAmovaResults_v3:
	# # # input:
		# # # bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# # # metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# # # output:
		# # # amv="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
		# # # iterlog="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.csv",
	# # # params:
		# # # outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}",
	# # # threads:
		# # # 1
	# # # log:
		# # # "simulations/{simid}/logs/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
	# # # shell:
		# # # """
		# # # (
		# # # {ngsAMOVA} -in {input.bcf} -isSim 1 -P {threads} -mEmIter 500 -out {params.outprefix} -tole 1e-{wildcards.tole} -doAMOVA 3 -doDist 1 -m {input.metadata} -printMatrix 0  -sqDist 1 > {output.iterlog}
		# # # )  2> {log}
		# # # """
#
#





#
# rule collect_results_1:
	# input:
		# expand("simulations/{simid}/benchmark/model_{model_id}/contig_{contig}/doAmova1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=5,
				# rep=REP)
	# output:
		# "simulations/sim_demes_v2/collected_results/benchmark/sim_demes_v2.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header=['s', 'h_m_s', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load', 'cpu_time', 'fid', 'simid', 'model', 'contig', 'rep', 'depth', 'tole']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile,delimiter='\t')
					# next(reader)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# depth = filename.split('_')[2].split('-')[4].split('d')[1]
					# tole = filename.split('_')[3].split('-')[0].split('tole')[1].split('.')[0]
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# row.append(tole)
						# writer.writerow(row)
#
#
#
#
#
#
# rule collect_results_2:
	# input:
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/iterLog/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=10,
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_doAmova3_gl_gt_iterlog.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header=['pair','iter','d','A','D','G','B','E','H','C','F','I', 'fid', 'simid', 'model', 'contig', 'rep', 'depth', 'tole']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# depth = filename.split('_')[2].split('-')[4].split('d')[1]
					# tole = filename.split('_')[3].split('-')[0].split('tole')[1].split('.')[0]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# row.append(tole)
						# writer.writerow(row)
#
#
#
# rule collect_results_3:
	# input:
		# #221214 tole10 + tole5
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=[10,5],
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_doAmova3_gl_gt_amova.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header = ['method','df_AG', 'ssd_AG', 'msd_AG', 'df_AIWG', 'ssd_AIWG', 'msd_AIWG', 'df_TOTAL', 'ssd_TOTAL', 'msd_TOTAL', 'coef_n', 'sigmasq_a', 'sigmasq_b', 'phi_a','fid', 'simid', 'model', 'contig', 'rep', 'depth', 'tole']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# depth = filename.split('_')[2].split('-')[4].split('d')[1]
					# tole = filename.split('_')[3].split('-')[0].split('tole')[1].split('.')[0]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# row.append(tole)
						# writer.writerow(row)
#
#
#
# #221217
# rule run_ngsAMOVA_sfs_var_calledgt:
	# input:
		# bcf="simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# params:
		# outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	# threads:
		# 1
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# shell:
		# """
		# (
		# {ngsAMOVA} -in {input.bcf} -isSim 0 -P {threads} -out {params.outprefix} -doAMOVA 2 -doDist 1 -m {input.metadata} -printMatrix 0  -sqDist 1
		# )  2> {log}
		# """


# rule run_ngsAMOVA_sfs_var_rawgt:
	# input:
		# vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
		# metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}.sfs.csv",
	# params:
		# outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}",
	# threads:
		# 1
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	# shell:
		# """
		# (
		# {ngsAMOVA} -in {input.vcf} -isSim 0 -P {threads} -out {params.outprefix} -doAMOVA 2 -doDist 1 -m {input.metadata} -printMatrix 0  -sqDist 1
		# )  2> {log}
		# """



# rule collect_results_4:
	# input:
		# #221214 tole10 + tole5
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_doAmova2_called_gt_amova.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header = ['method','df_AG', 'ssd_AG', 'msd_AG', 'df_AIWG', 'ssd_AIWG', 'msd_AIWG', 'df_TOTAL', 'ssd_TOTAL', 'msd_TOTAL', 'coef_n', 'sigmasq_a', 'sigmasq_b', 'phi_a','fid', 'simid', 'model', 'contig', 'rep', 'depth']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# depth = filename.split('.amova.csv')[0].split('-')[4].split('d')[1]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# writer.writerow(row)


# rule collect_results_5:
	# input:
		# #221214 tole10 + tole5
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_truth_raw_all_amova.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header = ['method','df_AG', 'ssd_AG', 'msd_AG', 'df_AIWG', 'ssd_AIWG', 'msd_AIWG', 'df_TOTAL', 'ssd_TOTAL', 'msd_TOTAL', 'coef_n', 'sigmasq_a', 'sigmasq_b', 'phi_a','fid', 'simid', 'model', 'contig', 'rep']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('.')[0].split('-')[3].split('rep')[1]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# writer.writerow(row)


# rule collect_results_6:
	# input:
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova2/called_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_doAmova2_called_gt_sfs.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header=[ 'Method',    'Ind1',    'Ind2',    'A',    'D',    'G',    'B',    'E',    'H',    'C',    'F',    'I',    'n_em_iter',    'shared_nSites',    'Delta',    'Tole',    'Sij',    'Fij',    'Fij2',    'IBS0',    'IBS1',    'IBS2',    'R0',    'R1',    'Kin', 'fid', 'simid', 'model', 'contig', 'rep', 'depth']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# next(reader)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# # depth = filename.split('-')[4].split('d')[1]
					# depth = filename.split('.sfs.csv')[0].split('-')[4].split('d')[1]
					# # print(depth)
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# writer.writerow(row)


# rule collect_results_7:
	# input:
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/truth_raw_all/{simid}-{model_id}-{contig}-rep{rep}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_truth_raw_all_sfs.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header=[ 'Method',    'Ind1',    'Ind2',    'A',    'D',    'G',    'B',    'E',    'H',    'C',    'F',    'I',    'n_em_iter',    'shared_nSites',    'Delta',    'Tole',    'Sij',    'Fij',    'Fij2',    'IBS0',    'IBS1',    'IBS2',    'R0',    'R1',    'Kin', 'fid', 'simid', 'model', 'contig', 'rep']
#
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# next(reader)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('.')[0].split('-')[3].split('rep')[1]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# writer.writerow(row)




# rule collect_results_8:
	# input:
		# expand("simulations/{simid}/model_{model_id}/contig_{contig}/doAmova3/gl_gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# depth=DEPTH,
				# tole=[10,5],
				# rep=REP),
	# output:
		# "simulations/sim_demes_v2/collected_results/sim_demes_v2_doAmova3_gl_gt_sfs.csv"
	# run:
		# rows = []
		# with open(output[0], "w") as outfile:
			# header=[ 'Method',    'Ind1',    'Ind2',    'A',    'D',    'G',    'B',    'E',    'H',    'C',    'F',    'I',    'n_em_iter',    'shared_nSites',    'Delta',    'Tole',    'Sij',    'Fij',    'Fij2',    'IBS0',    'IBS1',    'IBS2',    'R0',    'R1',    'Kin', 'fid', 'simid', 'model', 'contig', 'rep', 'depth', 'tole']
			# writer = csv.writer(outfile)
			# writer.writerow(header)
			# for fi in input:
#
				# with open(fi, "r") as infile:
					# reader = csv.reader(infile)
					# next(reader)
					# filename = fi.split('/')[-1]
					# simid = filename.split('-')[0]
					# model_id = filename.split('-')[1]
					# contig = filename.split('-')[2]
					# rep = filename.split('-')[3].split('rep')[1]
					# depth = filename.split('_')[2].split('-')[4].split('d')[1]
					# tole = filename.split('_')[3].split('-')[0].split('tole')[1].split('.')[0]
#
					# for row in reader:
						# row.append(fi)
						# row.append(simid)
						# row.append(model_id)
						# row.append(contig)
						# row.append(rep)
						# row.append(depth)
						# row.append(tole)
						# writer.writerow(row)






###############################################################################
# PER ITERATION DISTANCE MATRIX
###############################################################################


###############################################################################
# print per iteration distance matrices
# input: bcf file, metadata file
# output: per iteration distance matrix
# 	for all pairs of individuals and all iterations
# 
# 	columns in output:
# 	pair_index,n_em_iter,sq_dij,log10_d
rule run_printPerIterDistanceMatrix_devmode_stdout:
	input:
		bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		metadata="simulations/{simid}/model_{model_id}/{simid}_{model_id}_metadata.tsv",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv"
	params:
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv"
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv"
	shell:
		"""
		( {ngsAMOVA} -in {input.bcf} -doEM 1 -isSim 1 -P 1 -maxIter 500  -out {params.outprefix} -doAMOVA 1 -doDist 1 -printMatrix 0 -tole 1e-{wildcards.tole} -dev 1 -m {input.metadata} > {output} ) 2> {log}
		"""

###############################################################################
# collect per iteration distance matrix
# input: per iteration distance matrix
# output: distance matrix for each iteration
rule collect_run_printPerIterDistanceMatrix_devmode_stdout:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv"
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/perIterationDist/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_tole{tole}.perIterDistances.csv_emIter1.csv"
	params:
		emptyFiles="simulations/{simid}/logs/perIterationDist_emptyFiles.txt"
	shell:
		# get column 2, split file into multiple files based on column 2
		"""
		# if file is not empth
		if [[ -s {input} ]]; then 
			awk -F"," '{{print > FILENAME"_emIter"$2".csv"}}' {input}
		else
		# if file is empty
			echo "{input}" >> {params.emptyFiles}
			echo "empty" > {output}
		fi

		"""

################################################################################
# transform per iteration distance matrix to distance matrix input for ngsAMOVA
# input: distance matrix for each iteration
# script output: distance matrix input to be used by ngsAMOVA
# snakemake output: list of distance matrices
rule transform_data_to_distanceMatrix:
	output:
		"simulations/"+SIMULATION_ID+"/model_all/distanceMatrixList.txt"
	shell:
		"""
		bash transform_data_to_distanceMatrix.sh > {output}
		"""
	
###############################################################################
# use distance matrices with ngsAMOVA to get AMOVA results for each iteration
run_ngsAMOVA_perIterationDistList="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/run_ngsAMOVA_perIterationDistList.sh"
rule run_ngsAMOVA_with_distanceMatrix:
	input:
		"simulations/"+SIMULATION_ID+"/model_all/distanceMatrixList.txt",
	output:
		"simulations/"+SIMULATION_ID+"/model_all/distanceMatrixList_amovaResultsList.txt",
	log:
		"simulations/"+SIMULATION_ID+"/logs/model_all/distanceMatrixList_amovaResultsList.txt",
	shell:
		"""
		( bash {run_ngsAMOVA_perIterationDistList} {input} {output} ) 2> {log}
		"""
