import math
from itertools import product, combinations
import numpy as np
import pandas as pd
import sys
import os
import subprocess
import tskit, msprime
import demes, demesdraw
import csv

## Disable warnings (due to msprime version issues)
import warnings

warnings.filterwarnings("ignore")

###############################################################################
# BEGIN CONFIG
###############################################################################
SIMULATION_ID = "sim_demes_v3"

###############################################################################
# paths to programs
ANGSD = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_angsd/angsd/angsd"
vcfgl = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"
ngsAMOVA = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA"
bcftools = "/opt/software/bcftools/1.16/bin/bcftools"

###############################################################################
# simulation models
MODELS = ["model1", "model2"]

metadata_file = "sim/sim_demes_v3/metadata_2lvl_with_header.tsv"

# average per site depth
DEPTH = [50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1]
# DEPTH = [0.1, 1, 2, 5, 10, 20, 50]

# contig id x corresponds to contig length x * 1e6
# CONTIGS = [1, 2, 10, 50, 100]
# CONTIGS = [1, 10]
CONTIGS = [1, 10, 100, 500]


# Number of replicates
n_reps = 20
REP = [*range(n_reps)]

###############################################################################
# define populations
#
ploidy = 2

# define population with diploid individuals
# <pop_id> : <n_individuals_in_pop>
DEF_POPS = {"B1C1": 10, "B1C2": 10, "B2C1": 10, "B2C2": 10}


# note that msprime>1 and stdpopsim_dev version
# assumes ploidy=2 by default
haplo_list = []
indv_names = []
for i, (key, value) in enumerate(DEF_POPS.items()):
	haplo_list.append(value * ploidy)
	for ind in range(value):
		# Using PLINK-like format: <Family-ID>_<Individual-ID>
		# to store <Population-ID>_<Individual-ID>
		indv_names.append(f"pop{key}_ind{str(ind+1)}")
IND_PAIRS = list("-".join(map(str, comb)) for comb in combinations(indv_names, 2))


mutation_rate = 1.29e-08
recombination_rate = 1.14856e-08


###############################################################################
# AMOVA parameters


# END CONFIG
###############################################################################

###############################################################################
# BEGIN RULES

postCutoffDict = {"095": " -postCutoff 0.95 ", "03": " "}


rule all:
	input:
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var_baseCounts/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.baseCounts.tsv",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/data/gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/data/gt_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# postCutoff=postCutoffDict.keys(),
				# doPost=["1", "2"]
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/data/gt_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# postCutoff=postCutoffDict.keys(),
				# doPost=["1", "2"]
				# ),
		# expand("sim/{simid}/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.fasta",
				# simid=SIMULATION_ID,
				# model_id=MODELS,
				# contig=CONTIGS,
				# rep=REP,
				# depth=DEPTH,
				# ),
		expand("sim/{simid}/model_{model_id}/contig_{contig}/afa_nj_tree/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.newick",
				simid=SIMULATION_ID,
				model_id=MODELS,
				contig=CONTIGS,
				rep=REP,
				depth=DEPTH,
				),

		# expand(
			# simid=SIMULATION_ID,
			# model_id=MODELS,
			# contig=CONTIGS,
			# rep=REP,
			# depth=DEPTH,
			# postCutoff=postCutoffDict.keys(),
			# doPost=["1", "2"],
		# ),


###############################################################################
###############################################################################
# SIMULATION



rule simulation:
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
	params:
		seedfile="sim/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}_simulation_params.txt",
	run:
		model = demes.load("demes_conf/" + str(wildcards.model_id) + ".yaml")
		demography = msprime.Demography.from_demes(model)
		samples = DEF_POPS
		sequence_length = int(wildcards.contig) * 1e6
		seedval = (
				np.random.randint(2**20)
				+ (100 * int(wildcards.rep))
				+ int(wildcards.contig)
				)
		with open(params.seedfile, "w") as seedout:
			print(
					"Simulation_ID:"
					+ str(wildcards.simid)
					+ "\nReplicate:"
					+ str(wildcards.rep)
					+ "\nSeed:"
					+ str(seedval)
					+ "\nContig:"
					+ str(wildcards.contig)
					+ "\nSequence_length:"
					+ str(sequence_length)
					+ "\nMutation_rate:"
					+ str(mutation_rate)
					+ "\nRecombination_rate:"
					+ str(recombination_rate),
					file=seedout,
					)
			ts = msprime.sim_ancestry(
					[msprime.SampleSet(n, population=p) for p, n in samples.items()],
					demography=demography,
					sequence_length=sequence_length,
					random_seed=seedval,
					recombination_rate=recombination_rate,
					)
			mts = msprime.sim_mutations(
					ts, rate=mutation_rate, random_seed=seedval, model="binary"
					)
			with open(output[0], "w") as tsout:
				mts.dump(tsout)

###############################################################################
# GET VCF FROM TREE

## Convert tree sequence to VCF file
# Using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
	input:
		rules.simulation.output,
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	run:
		ts = tskit.load(input[0])
		with open(output[0], "w") as vcfout:
			ts.write_vcf(
					output=vcfout,
					contig_id=str(wildcards.contig),
					position_transform="legacy",
					individual_names=indv_names,
					)


			###############################################################################
# SIMULATE GENOTYPE LIKELIHOODS


## Using only variable sites
rule vcf_to_vcfgl_var:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	params:
		vcfgl=vcfgl,
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		mode="b",
		error_rate=0.002,
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	shell:
		"""
		RNG_SEED=$(( {wildcards.rep} + 1 ))

		({params.vcfgl} -in {input} -out {params.outprefix} \
				-depth {wildcards.depth} \
				-err {params.error_rate}  \
				-mode {params.mode} \
				-explode 0 \
				-seed ${{RNG_SEED}} ) 2> {log}
		"""


# use the same seeds to avoid rerunning all analyses
# print basecounts in the vcfgl simulations
rule vcf_to_vcfgl_var_printBaseCounts:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var_baseCounts/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.baseCounts.tsv",
	params:
		vcfgl=vcfgl,
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var_baseCounts/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
		mode="b",
		error_rate=0.002,
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var_baseCounts/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.baseCounts.tsv",
	shell:
		"""
		RNG_SEED=$(( {wildcards.rep} + 1 ))

		({params.vcfgl} -in {input} -out {params.outprefix} \
				-depth {wildcards.depth} \
				-err {params.error_rate}  \
				-mode {params.mode} \
				-explode 0 \
				-printBaseCounts 1 \
				-seed ${{RNG_SEED}} ) 2> {log}
		"""


# call fasta for phylogenetic tree construction using genotypes
# calls most common; if tie randomly sample
GETFASTA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/simulations/getFasta/getFasta"


rule baseCounts_to_fasta:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var_baseCounts/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.baseCounts.tsv",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.fasta",
	params:
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.fasta",
	shell:
		"""
		RNG_SEED=$(( {wildcards.rep} + 1 ))
		({GETFASTA} -i {input} -o {params.outprefix} -s ${{RNG_SEED}}) 2> {log}
		"""


#
# MUSCLE="./muscle"
# /maps/projects/lundbeck/scratch/pfs488/Programs/muscle
# rule align_pseudoHaplo_fasta_sequences:
	# input:
		# "sim/{simid}/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.fasta",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/aligned_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.afa",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/aligned_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.afa",
	# shell:
		# """
		# {MUSCLE} -align {input} -output {output}
		# """


getTree="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/R/get_njTree_from_afa.R"
#
# rule get_nj_tree_plot_from_pseudohaplo_afa:
	# input:
		# "sim/{simid}/model_{model_id}/contig_{contig}/aligned_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.afa",
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/afa_nj_tree/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.newick",
	# params:
		# outprefix="sim/{simid}/model_{model_id}/contig_{contig}/afa_nj_tree/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	# shell:
		# """
		# Rscript {getTree} {input} {params.outprefix}
		# """
#
#

rule get_nj_tree_plot_from_fastas:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/pseudo_haplo_fasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.fasta",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/afa_nj_tree/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.newick",
	params:
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/afa_nj_tree/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	shell:
		"""
		/opt/software/R/4.3.1/bin/Rscript {getTree} {input} {params.outprefix} 
		"""


###############################################################################
# PREPARE FILES

# remove simulated true genotypes from gl vcfs
# to prepare vcf with 10 gls ready for analyses
rule prepare_glvcfs_remove_gts:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	shell:
		"""
		{bcftools} annotate -x FORMAT/GT {input} -O b -o {output}
		"""


# get gl vcfs with 3 gls where majorminor is set to true majorminor allelic states
rule get_vcf_gl3_trueMajorMinor:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	params:
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/data/gl_3_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	shell:
		"""
		({ANGSD} -doMaf 1 -doPost 1 -doBcf 1 -vcf-gl {input} -out {params.outprefix} )2>{log}
		"""


###############################################################################
# ESTIMATE MAJOR MINOR

# estimate majorminor with angsd
# get gl vcfs with 3 gls where majorminor is set to estimated majorminor
rule get_vcf_gl3_estimMajorMinor:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	params:
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/data/gl_3_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	shell:
		"""
		({ANGSD} -doMajorMinor 1 -doBcf 1 -vcf-gl {input} -out {params.outprefix} )2>{log}
		"""

###############################################################################
# GENOTYPE CALLING

rule angsd_genotype_calling_postCutoffvar_doPostvar_estimMajorMinor:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gt_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/data/gt_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf"
	params:
		cutoffarg=lambda wildcards: str(postCutoffDict[str(wildcards.postCutoff)]),
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/data/gt_estimMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}"
	shell:
		"""
		({ANGSD} -doMaf 1 -doGeno 31 -doBcf 1 -vcf-gl {input} -out {params.outprefix} -doPost {wildcards.doPost} {params.cutoffarg} )2>{log}
		bcftools annotate -x "FORMAT/GL,FORMAT/GP" {params.outprefix}.bcf -O b -o {params.outprefix}.tmp.bcf
		mv {params.outprefix}.tmp.bcf {params.outprefix}.bcf
		"""

rule angsd_genotype_calling_postCutoffvar_doPostvar_trueMajorMinor:
	input:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gl_3_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
		"sim/{simid}/model_{model_id}/contig_{contig}/data/gt_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf"
	log:
		"sim/{simid}/logs/model_{model_id}/contig_{contig}/data/gt_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf"
	params:
		cutoffarg=lambda wildcards: str(postCutoffDict[str(wildcards.postCutoff)]),
		outprefix="sim/{simid}/model_{model_id}/contig_{contig}/data/gt_trueMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}"
	shell:
		"""
		({ANGSD} -doMaf 1 -doGeno 31 -doBcf 1 -vcf-gl {input} -out {params.outprefix} -doPost {wildcards.doPost} {params.cutoffarg} )2>{log}
		bcftools annotate -x "FORMAT/GL,FORMAT/GP" {params.outprefix}.bcf -O b -o {params.outprefix}.tmp.bcf
		mv {params.outprefix}.tmp.bcf {params.outprefix}.bcf
		"""




### SAFE ABOVE

#
# rule run_ngsAMOVA_genotype_likelihood_estimMajorMinor:
#
# rule run_ngsAMOVA_genotype_likelihood_trueMajorMinor:
#
# rule run_ngsAMOVA_genotype_calling_trueMajorMinor:
#
#
# # run with python to handle exit errors from ngsAMOVA
# # exit errors may occur if there isn't enough data
#
# rule run_ngsAMOVA_genotype_calling_trueMajorMinor:
	# input:
		# bcf=
		# metadata=metadata_file,
	# output:
	# params:
		# outprefix=
	# threads: 4
	# log:
	# run:
		# p = subprocess.Popen(
			# ngsAMOVA
			# + " -i "
			# + input.bcf
			# + " -P "
			# + str(threads)
			# + " -o "
			# + params.outprefix
			# + " -f 'Individual~Region/Population' "
			# + " -doAmova 1 -doEM 1 -doDist 1 --printDistanceMatrix 3"
			# + " --maxEmIter "
			# + str(wildcards.maxEmIter)
			# + " --em-tole 1e-"
			# + str(wildcards.tole)
			# + " -m "
			# + input.metadata,
			# shell=True,
			# stderr=subprocess.PIPE,
			# stdout=subprocess.PIPE,
		# )
		# stdout, stderr = p.communicate()
		# log_file = open(log[0], "w")
		# if p.returncode != 0:
			# if "No shared sites found for pair" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# elif "Total variance is" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# else:
				# print(str(stderr, "utf-8"), file=log_file)
				# raise Exception("ngsAMOVA failed with error: " + str(stderr))
		# else:
			# log_file.write(str(stderr, "utf-8"))
		# log_file.close()
#
#
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_likelihood/ngsAMOVA_maxIter{maxEmIter}_tole{tole}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# params:
		# outprefix="sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_likelihood/ngsAMOVA_maxIter{maxEmIter}_tole{tole}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	# threads: 4
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/amova_results/genotype_likelihood/ngsAMOVA_maxIter{maxEmIter}_tole{tole}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# run:
		# p = subprocess.Popen(
			# ngsAMOVA
			# + " -i "
			# + input.bcf
			# + " -P "
			# + str(threads)
			# + " -o "
			# + params.outprefix
			# + " -f 'Individual~Region/Population' "
			# + " -doAmova 1 -doEM 1 -doDist 1 --printDistanceMatrix 3"
			# + " --maxEmIter "
			# + str(wildcards.maxEmIter)
			# + " --em-tole 1e-"
			# + str(wildcards.tole)
			# + " -m "
			# + input.metadata,
			# shell=True,
			# stderr=subprocess.PIPE,
			# stdout=subprocess.PIPE,
		# )
		# stdout, stderr = p.communicate()
		# log_file = open(log[0], "w")
		# if p.returncode != 0:
			# if "No shared sites found for pair" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# elif "Total variance is" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# else:
				# print(str(stderr, "utf-8"), file=log_file)
				# raise Exception("ngsAMOVA failed with error: " + str(stderr))
		# else:
			# log_file.write(str(stderr, "utf-8"))
		# log_file.close()
#
#
# rule run_ngsAMOVA_genotype_call_2level_postCutoff_doPost_useAncDerFile:
	# input:
		# bcf="sim/{simid}/model_{model_id}/contig_{contig}/data/gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf",
		# tab="sim/{simid}/model_{model_id}/contig_{contig}/data/majorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1.tab",
		# metadata=metadata_file,
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# params:
		# outprefix="sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.log",
	# run:
		# p = subprocess.Popen(
			# ngsAMOVA
			# + " -i "
			# + input.bcf
			# + " --ancderfile "
			# + input.tab
			# + " -o "
			# + params.outprefix
			# + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
			# + input.metadata,
			# shell=True,
			# stderr=subprocess.PIPE,
			# stdout=subprocess.PIPE,
		# )
		# stdout, stderr = p.communicate()
		# log_file = open(log[0], "w")
		# if p.returncode != 0:
			# if "No shared sites found for pair" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# elif "Total variance is" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# else:
				# print(str(stderr, "utf-8"), file=log_file)
				# raise Exception("ngsAMOVA failed with error: " + str(stderr))
		# else:
			# log_file.write(str(stderr, "utf-8"))
		# log_file.close()
#
#
# rule run_ngsAMOVA_genotype_call_2level_postCutoff_doPost_fixedMajorMinor:
	# input:
		# bcf="sim/{simid}/model_{model_id}/contig_{contig}/data/gt/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_postCutoff{postCutoff}_doPost{doPost}.bcf",
		# metadata=metadata_file,
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	# params:
		# outprefix="sim/{simid}/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/amova_results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_fixedMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.log",
	# run:
		# p = subprocess.Popen(
			# ngsAMOVA
			# + " -i "
			# + input.bcf
			# + " -o "
			# + params.outprefix
			# + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
			# + input.metadata,
			# shell=True,
			# stderr=subprocess.PIPE,
			# stdout=subprocess.PIPE,
		# )
		# stdout, stderr = p.communicate()
		# log_file = open(log[0], "w")
		# if p.returncode != 0:
			# if "No shared sites found for pair" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# elif "Total variance is" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# else:
				# print(str(stderr, "utf-8"), file=log_file)
				# raise Exception("ngsAMOVA failed with error: " + str(stderr))
		# else:
			# log_file.write(str(stderr, "utf-8"))
		# log_file.close()
#
#
# rule run_ngsAMOVA_truth:
	# input:
		# vcf="sim/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
		# metadata=metadata_file,
	# output:
		# "sim/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
	# params:
		# outprefix="sim/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}",
	# log:
		# "sim/{simid}/logs/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
	# run:
		# p = subprocess.Popen(
			# ngsAMOVA
			# + " -i "
			# + input.vcf
			# + " -o "
			# + params.outprefix
			# + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
			# + input.metadata,
			# shell=True,
			# stderr=subprocess.PIPE,
			# stdout=subprocess.PIPE,
		# )
		# stdout, stderr = p.communicate()
		# log_file = open(log[0], "w")
		# if p.returncode != 0:
			# if "No shared sites found for pair" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# elif "Total variance is" in str(stderr):
				# shell("echo 0 > " + str(output[0]))
				# print(str(stderr, "utf-8"), file=log_file)
			# else:
				# print(str(stderr, "utf-8"), file=log_file)
				# raise Exception("ngsAMOVA failed with error: " + str(stderr))
		# else:
			# log_file.write(str(stderr, "utf-8"))
		# log_file.close()
#
#
# # ##############################################################################
# # END RULES
# # ##############################################################################
