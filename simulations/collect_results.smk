# import numpy as np
import pandas as pd
import sys
import os
import subprocess

from itertools import product, combinations
import csv


## Disable warnings (due to msprime version issues)
import warnings

warnings.filterwarnings("ignore")

###############################################################################
# BEGIN CONFIG
###############################################################################
SIMULATION_ID = "sim_demes_v2"

# paths to programs
ANGSD = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_angsd/angsd/angsd"
vcfgl = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"
ngsAMOVA = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA"
ngsLD = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/ngsLD/ngsLD"

# simulation models
MODELS = ["model1", "model2"]

# average per site depth
# DEPTH = [10, 5, 2, 1, 0.5, 0.2, 0.1]
# 230411
# add 20X and 50X to see how genotype calling and gle converges at very high depth
DEPTH = [0.1, 1, 2, 5, 10, 20, 50]
# DEPTH = [0.1,1]
# DEPTH = [2,5]
# DEPTH = [10,20,50]


# contig id x corresponds to contig length x * 1e6
# CONTIGS = [1]
# CONTIGS = [10]
CONTIGS = [1, 10, 100]

# CONTIGS = [50, 100]
# CONTIGS = [1, 2, 10, 50, 100]
# CONTIGS = [1, 10, 100]
# CONTIGS = [2, 50]

# Number of replicates
n_reps = 200
REP = [*range(n_reps)]
REP = REP[:20]

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

# Number of replicates
n_reps = 200
REP = [*range(n_reps)]
REP = REP[:20]

# Genotype Calling
postCutoffDict = {"095": " -postCutoff 0.95 ", "03": " "}

#
#
# END CONFIG
###############################################################################


###############################################################################
# BEGIN RULES
#
#


rule all:
    input:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_fixedMajorMinor_run2.csv",
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_likelihood_inferredMajorMinor.csv",
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_likelihood_fixedMajorMinor.csv",
        "simulations/" + SIMULATION_ID + "/results/collect_results_truth.csv",
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_inferredMajorMinor.csv",
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_fixedMajorMinor.csv",
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/truth/{simid}-{model_id}-{contig}-rep{rep}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            rep=REP,
        ),
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_n_variable_sites_simulated.csv",


#
# ################################
# # with estimated major/minor
# "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
# "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",

# "simulations/{simid}/model_{model_id}/contig_{contig}/truth_2/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",

# # with fixed true major/minor
# "simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
# "simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
# ################################


rule get_bcfStats_true:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.stats",
    shell:
        """
        bcftools stats {input} > {output}
        """


rule get_n_variable_sites_simulated:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.stats",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/truth/{simid}-{model_id}-{contig}-rep{rep}.n_variable_sites",
    shell:
        """
        printf "{wildcards.model_id},{wildcards.contig},{wildcards.rep},$(cat {input}|grep "^SN"|grep "number of SNPs"|cut -d: -f2|tr -d '\t')\n" > {output}
        """


rule collect_n_variable_sites_simulated:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/truth/{simid}-{model_id}-{contig}-rep{rep}.n_variable_sites",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            rep=REP,
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_n_variable_sites_simulated.csv",
    shell:
        """
        printf "model_id,contig,rep,n_variable_sites\n" > {output}
        cat {input} >> {output}
        """


rule prepare_results_genotype_likelihood_fixedMajorMinor:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
    shell:
        """
        ID="genotype_likelihood_maxIter500_tole5"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},$i\n" 
            done >> {output}
        fi
        """


rule prepare_results_genotype_likelihood_inferredMajorMinor:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
    shell:
        """
        ID="genotype_likelihood_maxIter500_tole5_doMajorMinor1"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},$i\n" 
            done >> {output}
        fi
        """


rule prepare_results_genotype_calling_fixedMajorMinor:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
    shell:
        """
        ID="genotype_calling_postCutoff{wildcards.postCutoff}_doPost{wildcards.doPost}"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},$i\n" 
            done >> {output}
        fi
        """


rule prepare_results_genotype_calling_fixedMajorMinor_run2:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
    shell:
        """
        ID="genotype_calling_postCutoff{wildcards.postCutoff}_doPost{wildcards.doPost}"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},$i\n" 
            done >> {output}
        fi
        """


rule prepare_results_genotype_calling_inferredMajorMinor:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
    shell:
        """
        ID="genotype_calling_postCutoff{wildcards.postCutoff}_doPost{wildcards.doPost}_doMajorMinor1"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},{wildcards.depth},$i\n" 
            done >> {output}
        fi
        """


rule prepare_results_truth:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/truth_2/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/results/truth/{simid}-{model_id}-{contig}-rep{rep}.result",
    shell:
        """
        ID="truth"
        if [[ $(cat {input}) == 0 || $(cat {input}) == "" ]];then
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},truth,Phi,Region_in_Total,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},truth,Phi,Population_in_Region,NA\n" >> {output}
            printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},truth,Phi,Population_in_Total,NA\n" >> {output}
        else
            res=($(cat {input} | grep "^Phi" ))
            for i in "${{res[@]}}"; do
                printf "${{ID}},{wildcards.model_id},{wildcards.contig},{wildcards.rep},truth,$i\n" 
            done >> {output}
        fi
        """


rule collect_results_genotype_likelihood_fixedMajorMinor:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_likelihood_fixedMajorMinor.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


rule collect_results_genotype_likelihood_inferredMajorMinor:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_likelihood_inferredMajorMinor.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


rule collect_results_genotype_calling_fixedMajorMinor:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_fixedMajorMinor.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


rule collect_results_genotype_calling_fixedMajorMinor_run2:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_fixedMajorMinor_run2.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


rule collect_results_genotype_calling_inferredMajorMinor:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/genotype_calling/postCutoff{postCutoff}_doPost{doPost}_doMajorMinor1/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
    output:
        "simulations/"
        + SIMULATION_ID
        + "/results/collect_results_genotype_calling_inferredMajorMinor.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


rule collect_results_truth:
    input:
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/results/truth/{simid}-{model_id}-{contig}-rep{rep}.result",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            rep=REP,
        ),
    output:
        "simulations/" + SIMULATION_ID + "/results/collect_results_truth.csv",
    shell:
        """
        printf "AnalysisType,Model,Contig,Rep,Depth,Measure,Level,Value\n" > {output}
        cat {input} >> {output}
        """


# END RULES
###############################################################################
