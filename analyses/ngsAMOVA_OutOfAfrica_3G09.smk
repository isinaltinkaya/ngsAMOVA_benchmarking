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
# ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA/ngsAMOVA"

#version disabling amova get only dist mat
ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/ngsAMOVA_vDisableAmova/ngsAMOVA/ngsAMOVA"

# Average per site depth
# DEPTH=[20,10,5,2,1,0.5,0.2,0.1,0.01]
DEPTH=[10,5,2,1,0.5,0.2,0.1,0.01]

# Number of replicates
n_reps=200

REP=[*range(n_reps)]
# REP=REP[int(config["NN"])]
REP=REP[:20]
# REP=REP[:5]
# REP=REP[:1]


# Set seed
np.random.seed(42)
SEED=np.random.randint(1,2**32-1,n_reps)
# print(SEED)


species=stdpopsim.get_species("HomSap")


exclude_chr_list=['chrX','chrY']
CONTIGID=[c for c in [co.id for co in species.genome.chromosomes] if c not in exclude_chr_list ]

CONTIGID=CONTIGID[config["cc"]]

MODEL=["OutOfAfrica_3G09"]


n_samples_popi=50
ploidy=2
n_pops=3

samples_per_pop=[100]*3

# TOLE=[9,11]
# TOLE=10
TOLE=7


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

# ITERSET=[50,100,150,200,250]
# THRSET=[i for i in range(1,10)]
# THRSET=[1,10]
# THRSET=1


AMOVASET={"noAMOVA":-1,"doAMOVA1":1,"doAMOVA3":3}

rule all:
	input:
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA/windowEst/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}.sfs.csv",
				depth=DEPTH,
				simid=SIMULATION_ID,
				contig=CONTIGID,
				rep=REP,
				model_id=MODEL,
				winsize=[100,50],
				winid=1,
				tole=7,
				doamv=["doAMOVA3"]),
		# expand("simulations/{simid}/model_{model_id}/contig_all/masked_vcfgl_var/{simid}-{model_id}-contig_all-rep{rep}-d{depth}-tole{tole}-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# rep=REP,
				# depth=DEPTH,
				# doamv=["doAMOVA1"]),
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t1-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# rep=REP,
				# depth=DEPTH,
				# doamv=["noAMOVA"],
				# it=ITERSET),
#
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# rep=REP,
				# depth=DEPTH,
				# doamv=["noAMOVA"],
				# it=ITERSET),
#
		# expand("simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_est/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
				# simid=SIMULATION_ID,
				# model_id=MODEL,
				# contig=CONTIGID,
				# tole=TOLE,
				# rep=REP,
				# depth=DEPTH,
				# doamv=["noAMOVA"],
				# it=ITERSET),
#
# rule benchmark_run_ngsAMOVA_sfs_var_minInd2_100K_maxIterSet:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/masked_vcfgl_var/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t1-{doamv}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t1-{doamv}",
		# doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t1-{doamv}.sfs.csv",
	# threads:
		# 1
	# benchmark:
		# "simulations/{simid}/benchmark221010/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t1-{doamv}.sfs.csv"
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 2> {log}
		# """
#
#
# rule benchmark_run_ngsAMOVA_sfs_var_minInd2_100K_maxIterSet_thread10:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/masked_vcfgl_var/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}",
		# doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
	# threads:
		# 10
	# benchmark:
		# "simulations/{simid}/benchmark221010/model_{model_id}/contig_all/subsample_100K/ngsAmova/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv"
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 2> {log}
		# """
#
#
# rule benchmark_run_ngsAMOVA_sfs_var_minInd2_100K_maxIterSet_thread10_getDelta:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/masked_vcfgl_var/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_est/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_est/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}",
		# doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_est/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
	# threads:
		# 10
	# benchmark:
		# "simulations/{simid}/benchmark221016/model_{model_id}/contig_all/subsample_100K/ngsAMOVA_est/{simid}-{model_id}-contig_all_100K-rep{rep}-d{depth}-i{it}-t10-{doamv}.sfs.csv",
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 2> {log}
		# """
#
#
#
# rule run_ngsAMOVA_sfs_var_minInd2_contigAll_thread20:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_all/masked_vcfgl_var/{simid}-{model_id}-contig_all-rep{rep}-d{depth}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_all/masked_vcfgl_var/{simid}-{model_id}-contig_all-rep{rep}-d{depth}-tole{tole}-{doamv}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_all/masked_vcfgl_var/{simid}-{model_id}-contig_all-rep{rep}-d{depth}-tole{tole}-{doamv}",
		# doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
		# metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_all/masked_vcfgl_var/{simid}-{model_id}-contig_all-rep{rep}-d{depth}-tole{tole}-{doamv}.sfs.csv",
	# threads:
		# 20
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter 300 -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 -tole 1e-{wildcards.tole} -m {params.metadata} 2> {log}
		# """
#
#
# rule run_ngsAMOVA_sfs_var_minInd2_contigWindows_thread20:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/windowBcf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB_win{winid}.bcf",
	# output:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA/windowBcf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}.sfs.csv",
	# params:
		# ngsAMOVA=ngsAMOVA,
		# outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA/windowBcf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}",
		# doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
		# metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	# log:
		# "simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsAMOVA/windowBcf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}.sfs.csv",
	# threads:
		# 20
	# shell:
		# """
		# {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter 300 -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 -tole 1e-{wildcards.tole} -m {params.metadata} 2> {log}
		# """
		# # {params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter {wildcards.it} -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 -tole 1e-7 2> {log}
#



rule run_ngsAMOVA_perContigWindow:
	input:
		"simulations/{simid}/model_{model_id}/contig_{contig}/masked_vcfgl_var/windowBcf/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}.bcf"
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA/windowEst/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}.sfs.csv"
	params:
		ngsAMOVA=ngsAMOVA,
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/ngsAMOVA/windowEst/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}",
		doAMOVA=lambda wildcards: AMOVASET[wildcards.doamv],
		metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/metadata_150inds.tsv",
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsAMOVA/windowEst/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-win{winsize}MB-win{winid}-tole{tole}-{doamv}.sfs.csv"
	threads:
		20
	shell:
		"""
		{params.ngsAMOVA} -in {input} -isSim 1 -P {threads} -maxIter 300 -out {params.outprefix} -doAMOVA {params.doAMOVA} -doDist 1  -minInd 2 -tole 1e-{wildcards.tole} -m {params.metadata} 2> {log}
		"""





