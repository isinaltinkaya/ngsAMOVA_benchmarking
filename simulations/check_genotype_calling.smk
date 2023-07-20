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
        # expand(
        #     "simulations/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
        #     simid=SIMULATION_ID,
        #     model_id=MODELS,
        #     contig=CONTIGS,
        #     rep=REP,
        # ),
        # expand(
        #     "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.vcf",
        #     simid=SIMULATION_ID,
        #     model_id=MODELS,
        #     contig=CONTIGS,
        #     rep=REP,
        #     depth=DEPTH,
        #     postCutoff=postCutoffDict.keys(),
        #     doPost=["1", "2"],
        # ),
        expand(
            "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/picard_check_concordance/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-ind{ind}.genotype_concordance_summary_metrics",
            simid=SIMULATION_ID,
            model_id=MODELS,
            contig=CONTIGS,
            rep=REP,
            depth=DEPTH,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
            ind=indv_names,
        ),


rule get_truth_vcf:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d50.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    shell:
        """
        bcftools annotate -x "FORMAT/GL,FORMAT/DP" -o {output} {input}
        """


rule called_bcf_to_vcf_for_picard:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.vcf",
    shell:
        """
        bcftools view {input} > {output}
        """


rule picard_check_genotype_calling_concordance:
    input:
        called="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.vcf",
        truth="simulations/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    params:
        prefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/picard_check_concordance/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-ind{ind}",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/picard_check_concordance/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}-ind{ind}.genotype_concordance_summary_metrics",
    shell:
        """
        java -jar /opt/software/picard/2.27.5/picard.jar GenotypeConcordance TRUTH_VCF={input.truth} CALL_VCF={input.called} OUTPUT={params.prefix} TRUTH_SAMPLE={wildcards.ind} CALL_SAMPLE={wildcards.ind}
        """


# rule check_genotype_calling_concordance_collect:
#     input:
#         expand(
#             simid=SIMULATION_ID,
#             contig=CONTIGS,
#             depth=DEPTH,
#             rep=REP,
#             postCutoff=postCutoffDict.keys(),
#             doPost=["1", "2"],
#         ),
#     output:
#         "simulations/sim_demes_v2/model_{model_id}/genotype_calling/genotype_calling_concordance.csv",
#     shell:
#         """
#         echo "rt,rf,rn,nt,nf,nn,n,model,contig,rep,depth" > {output}
#         cat {input} >> {output}
#         """
