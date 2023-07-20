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
        # # expand(
        # #     "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
        # #     simid=SIMULATION_ID,
        # #     model_id=MODELS,
        # #     contig=CONTIGS,
        # #     rep=REP,
        # #     depth=DEPTH,
        # #     postCutoff=postCutoffDict.keys(),
        # #     doPost=["1", "2"],
        # # ),
        # expand(
            # "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
            # simid=SIMULATION_ID,
            # model_id=MODELS,
            # contig=CONTIGS,
            # rep=REP,
            # depth=DEPTH,
        # ),


# doMajorMinor1_2 --> 2 means version 2 that is fixed 230605
rule run_ngsAMOVA_genotype_call_postCutoff_doPost_useAncDerFile_fixed230605:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
        tab="simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1_ancder.tab",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -o "
            + params.outprefix
            + " --ancderfile "
            + input.tab
            + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            elif "Total variance is" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


# doMajorMinor1_2 --> 2 means version 2 that is fixed 230605
rule run_ngsAMOVA_genotype_likelihood_2level_useAncDerFile_fixed230605:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
        tab="simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1_ancder.tab",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    threads: 4
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1_2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -P "
            + str(threads)
            + " -o "
            + params.outprefix
            + " --ancderfile "
            + input.tab
            + " -f 'Individual~Region/Population' -doAmova 1 -doEM 1 -doDist 1 --printDistanceMatrix 3  --maxEmIter 500 --em-tole 1e-5 -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            elif "Total variance is" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


## 230426
rule angsd_get_beagle_for_ngsLD:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/ngsLD/angsd_get_beagle/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.beagle.gz",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/ngsLD/angsd_get_beagle/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsLD/angsd_get_beagle/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.beagle.gz",
    shell:
        """
        (
            {ANGSD} -doglf 2 -doMajorMinor 1 -vcf-gl {input} -out {params.outprefix}
        ) > {log} 2>&1
        """


rule ngsLD_get_LD:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/ngsLD/angsd_get_beagle/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.beagle.gz",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/ngsLD/ngsLD_get_LD/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.ld",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/ngsLD/ngsLD_get_LD/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.ld",
    threads: 5
    shell:
        """
        zcat {input} |cut -f1| sed -e 1d -e 's/_/\t/g' > {output}.ngsld_pos.tsv;
        (
            {ngsLD} --geno {input} --probs --n_ind 40 --n_sites $( expr $(zcat {input} | wc -l) - 1) --n_threads {threads} --pos {output}.ngsld_pos.tsv -out {output}
        ) > {log} 2>&1
        """


###############################################################################


###############################################################################
## 230423
rule angsd_doMajorMinor:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1.mafs.gz",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1.mafs.gz",
    shell:
        """
        ({ANGSD} -vcf-gl {input} -doMajorMinor 1 -doMaf 1 -out {params.outprefix}) 2> {log}
        """


rule angsd_doMajorMinor_get_doMajorMinor_ancDerFile:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1.mafs.gz",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1_ancder.tab",
    shell:
        """
        zcat {input} | cut -f1-4|sed 1d > {output}
        """


rule run_ngsAMOVA_genotype_likelihood_2level_useAncDerFile:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
        tab="simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1_ancder.tab",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    threads: 4
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    # shell:
    # """
    # ({ngsAMOVA} -i {input.bcf} -P {threads} -out {params.outprefix} -f "Individual~Region/Population" -doAMOVA 1 -doEM 1 -doDist 1 -m {input.metadata} --printDistanceMatrix 3  --maxEmIter 500 --em-tole 1e-5 --ancDerFile {input.tab} )  2> {log}
    # """
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -P "
            + str(threads)
            + " -o "
            + params.outprefix
            + " --ancderfile "
            + input.tab
            + " -f 'Individual~Region/Population' -doAmova 1 -doEM 1 -doDist 1 --printDistanceMatrix 3  --maxEmIter 500 --em-tole 1e-5 -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            elif "Total variance is" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


rule run_ngsAMOVA_genotype_call_postCutoff_doPost_useAncDerFile:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
        tab="simulations/{simid}/model_{model_id}/contig_{contig}/angsd_doMajorMinor/doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_doMajorMinor1_ancder.tab",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_doMajorMinor1/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    # shell:
    # """
    # ({ngsAMOVA} -i {input.bcf} -out {params.outprefix} -f "Individual~Region/Population" -doAmova 1 -doDist 2 -m {input.metadata} --printDistanceMatrix 3  --ancDerFile {input.tab} )  2> {log}
    # """
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -o "
            + params.outprefix
            + " --ancderfile "
            + input.tab
            + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            elif "Total variance is" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


###############################################################################
### 230305 rerun all ngsAMOVA
# run ngsAMOVA on results from vcfgl
rule run_ngsAMOVA_genotype_likelihood_2level:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    threads: 4
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    shell:
        """
        ({ngsAMOVA} -i {input.bcf} -P {threads} -out {params.outprefix} -doAMOVA 1 -doEM 1 -doDist 1 -m {input.metadata} --printDistanceMatrix 3  --maxEmIter 500 --em-tole 1e-5 )  2> {log}
        """


# 230306
# collect genotype likelihood results (tole=1e-5)
rule collect_results_run_ngsAMOVA_genotype_likelihood_2level:
    input:
        expand(
            "simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
            simid=SIMULATION_ID,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
    output:
        phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_phi_{model_id}.csv",
    shell:
        """
            bash {collect_results_gle_phi_tole5} {input} > {output.phi}
        """


rule collect_get_mean_shared_nSites_run_ngsAMOVA_genotype_likelihood_2level:
    input:
        expand(
            "simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
            simid=SIMULATION_ID,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
        ),
    output:
        "simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_mean_shared_nSites_{model_id}.csv",
    shell:
        """
        i=0
        for inamv in {input};do
            if [ ${{i}} -eq 0 ];then
                printf "mean_shared_nSites,Model,Contig,Rep,Depth,Tole,AnalysisType\n"
            fi
            i=$((i+1))
            bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
            infile=$(echo ${{inamv}} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
            val=$(zcat ${{infile}} | cut -d, -f11 | datamash mean 1)
            model=$(echo ${{bname}}| cut -d- -f2)
            contig=$(basename ${{bname}} | cut -d- -f3)
            rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
            depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g'|cut -d_ -f1)
            tole=$(dirname ${{infile}} | rev | cut -d/ -f1 | rev | cut -d_ -f3 | sed 's/tole//g')
            type="GLE"
            printf "${{val}},${{model}},${{contig}},${{rep}},${{depth}},${{tole}},${{type}}\n" 
        done > {output}
        """


# rule collect_get_mean_shared_nSites_run_ngsAMOVA_genotype_likelihood_2level:
# 	input:
# 		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
# 			simid=SIMULATION_ID,
# 			contig=CONTIGS,
# 			depth=DEPTH,
# 			rep=REP),
# 	output:
# 		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_mean_shared_nSites_{model_id}.csv"
# 	shell:
# 		"""
# 		i=0
# 		for inamv in {input};do
# 			if [ ${{i}} -eq 0 ];then
# 				printf "mean_shared_nSites,Model,Contig,Rep,Depth,Tole,AnalysisType\n"
# 			fi
# 			i=$((i+1))
# 			bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
# 			infile=$(echo ${{inamv}} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
# 			val=$(zcat ${{infile}} | cut -d, -f11 | datamash mean 1)
# 			model=$(echo ${{bname}}| cut -d- -f2)
# 			contig=$(basename ${{bname}} | cut -d- -f3)
# 			rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
# 			depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g'|cut -d_ -f1)
# 			tole=$(dirname ${{infile}} | rev | cut -d/ -f1 | rev | cut -d_ -f3 | sed 's/tole//g')
# 			type="GLE"
# 			printf "${{val}},${{model}},${{contig}},${{rep}},${{depth}},${{tole}},${{type}}\n"
# 		done > {output}
# 		"""

###############################################################################

#TODO
rerun this without -sites but with majminest 3gl files 1REF 1ALT for estimated majorminor
rerun this with 10gl non-majminest fixed 1REF ref[0]=A 3ALT alts[0]=T 


###################
# use {postCutoff} {doPost}
rule angsd_genotype_calling_postCutoffvar_doPostvar:
    input:
        vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/rmGT/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_rmGT.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    params:
        cutoffarg=lambda wildcards: str(postCutoffDict[str(wildcards.postCutoff)]),
        prefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    shell:
        """
        ({ANGSD} -doMajorMinor 3 -doMaf 1 -doGeno 31 -vcf-gl {input.vcf} -out {params.prefix} -doBcf 1 -doPost {wildcards.doPost} {params.cutoffarg} )2>{log}
        """


# 230526 rerun genotype calling amova with the latest version
# run ngsAMOVA on results from angsd_genotype_calling_postCutoffvar_doPostvar
rule run_ngsAMOVA_genotype_call_postCutoff_doPost_run2:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA_run2/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.log",
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -o "
            + params.outprefix
            + " -f 'Individual~Region/Population' -doAmova 1 -doDist 2 --printDistanceMatrix 3  --isSim 1 -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            elif "Total variance is" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


# 230306 rerun genotype calling amova
# run ngsAMOVA on results from angsd_genotype_calling_postCutoffvar_doPostvar
rule run_ngsAMOVA_genotype_call_postCutoff_doPost:
    input:
        bcf="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
    run:
        p = subprocess.Popen(
            ngsAMOVA
            + " -i "
            + input.bcf
            + " -o "
            + params.outprefix
            + " -doAmova 1 -doDist 2 --printDistanceMatrix 3  -m "
            + input.metadata,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        stdout, stderr = p.communicate()
        log_file = open(log[0], "w")
        if p.returncode != 0:
            if "No shared sites found for pair" in str(stderr):
                shell("echo 0 > " + str(output[0]))
                print(str(stderr, "utf-8"), file=log_file)
            else:
                print(str(stderr, "utf-8"), file=log_file)
                raise Exception("ngsAMOVA failed with error: " + str(stderr))
        else:
            log_file.write(str(stderr, "utf-8"))
        log_file.close()


collect_results_genotype_call_doPost_postCutoff_phi = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/collect_results/collect_results_genotype_call_doPost_postCutoff_phi.sh"


# 230306
rule collect_results_run_ngsAMOVA_genotype_call_postCutoff_doPost:
    input:
        expand(
            "simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
            simid=SIMULATION_ID,
            contig=CONTIGS,
            depth=DEPTH,
            rep=REP,
            postCutoff=postCutoffDict.keys(),
            doPost=["1", "2"],
        ),
    output:
        phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_phi_{model_id}.csv",
    shell:
        """
        bash {collect_results_genotype_call_doPost_postCutoff_phi} {input} > {output.phi}
        """


#
# rule collect_get_mean_shared_nSites_run_ngsAMOVA_genotype_call_postCutoff_doPost:
# input:
# expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
# simid=SIMULATION_ID,
# contig=CONTIGS,
# depth=DEPTH,
# rep=REP,
# postCutoff=postCutoffDict.keys(),
# doPost=["1","2"]),
# output:
# "simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_shared_nSites_{model_id}.csv",
# shell:
# """
# i=0
# for inamv in {input};do
#
# if [ ${{i}} -eq 0 ];then
# printf "mean_shared_nSites,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
# fi
#
# i=$((i+1))
# bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
# infile=$(echo ${{inamv}} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
# model=$(echo ${{bname}}| cut -d- -f2)
# contig=$(basename ${{bname}} | cut -d- -f3)
# rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
# depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g')
# type="GTC"
# doPost=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
# postCutoff=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')
#
# if [[ ${{postCutoff}} == "03" ]];then
# postCutoff="0.3";
# elif [[ ${{postCutoff}} == "095" ]];then
# postCutoff="0.95";
# else
# exit 1;
# fi
#
# if [[ $(cat ${{inamv}}) == "0" ]];then
# val="NA"
# else
# val=$(zcat ${{infile}} | cut -d, -f11 | datamash mean 1)
# fi
# printf "${{val}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"
#
# done > {output}
# """
# rule collect_get_mean_ngsAMOVA_sites_used_run_ngsAMOVA_genotype_call_postCutoff_doPost:
# input:
# expand("simulations/{simid}/logs/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
# simid=SIMULATION_ID,
# contig=CONTIGS,
# depth=DEPTH,
# rep=REP,
# postCutoff=postCutoffDict.keys(),
# doPost=["1","2"]),
# output:
# "simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_ngsAMOVA_sites_used_{model_id}.csv",
# shell:
# """
# i=0
# for inamv in {input};do
#
# if [ ${{i}} -eq 0 ];then
# printf "nSites_processed,nSites_skipped,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
# fi
#
# i=$((i+1))
# bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
# model=$(echo ${{bname}}| cut -d- -f2)
# contig=$(basename ${{bname}} | cut -d- -f3)
# rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
# depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g')
# doPost=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
# postCutoff=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')
# type="GTC"
#
# if [[ ${{postCutoff}} == "03" ]];then
# postCutoff="0.3";
# elif [[ ${{postCutoff}} == "095" ]];then
# postCutoff="0.95";
# else
# exit 1;
# fi
# if [[ $(cat ${{inamv}} | grep -q "No shared sites") ]];then
# processed="NA"
# skipped="NA"
# else
# processed=$(cat ${{inamv}} |grep "Total number of sites processed"|cut -d: -f2|tr -d " ")
# skipped=$(cat ${{inamv}} |grep "Total number of sites skipped for all individual pairs"|cut -d: -f2|tr -d " ")
# fi
# printf "${{processed}},${{skipped}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"
#
# done > {output}
# """
#


################################################################################
# run 2 level AMOVA with ngsAMOVA
# using bcf using the raw simulated genotypes directly from msprime
# 230306
rule run_ngsAMOVA_truth:
    input:
        vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
    shell:
        """
        (
            {ngsAMOVA} -i {input.vcf} -o {params.outprefix} -doAmova 1 -doDist 2 -m {input.metadata} --printDistanceMatrix 3 -pJGCD 3
        )  2> {log}
        """


# rerun with the latest version to check if results are the same
rule run_ngsAMOVA_truth_may23_prepare_gt_vcf:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d50.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    shell:
        """
        bcftools annotate -x "FORMAT/GL,FORMAT/DP" -O b -o {output} {input}
        """


rule run_ngsAMOVA_truth_may23:
    input:
        vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcf_true/{simid}-{model_id}-{contig}-rep{rep}.vcf",
        metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/truth_2/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
    params:
        outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/truth_2/{simid}-{model_id}-{contig}-rep{rep}",
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/truth_2/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
    shell:
        """
        (
            {ngsAMOVA} -i {input.vcf} -o {params.outprefix} -doAmova 1 -doDist 2 -m {input.metadata} --printDistanceMatrix 3  --ancderfile {input.tab} --formula "Individual~Region/Population"
        )  2> {log}
        """


collect_results_truth_phi = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/collect_results/collect_results_truth_phi.sh"


# rule collect_results_run_ngsAMOVA_truth_march23:
# input:
# expand(
# "simulations/{simid}/model_{{model_id}}/contig_{contig}/truth/{simid}-{{model_id}}-{contig}-rep{rep}.amova.csv",
# simid=SIMULATION_ID,
# contig=CONTIGS,
# rep=REP,
# ),
# output:
# phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_truth_phi_{model_id}.csv",
# shell:
# """
# bash {collect_results_truth_phi} {input} > {output.phi}
# """


rule collect_get_mean_shared_nSites_run_ngsAMOVA_truth:
    input:
        expand(
            "simulations/{simid}/model_{{model_id}}/contig_{contig}/truth/{simid}-{{model_id}}-{contig}-rep{rep}.amova.csv",
            simid=SIMULATION_ID,
            contig=CONTIGS,
            rep=REP,
        ),
    output:
        "simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_truth_mean_shared_nSites_{model_id}.csv",
    shell:
        """
        i=0
        for inamv in {input};do
            if [ ${{i}} -eq 0 ];then
                printf "mean_shared_nSites,Model,Contig,Rep,AnalysisType\n"
            fi
            i=$((i+1))
            bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
            infile=$(echo ${{inamv}} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
            val=$(zcat ${{infile}} | cut -d, -f11 | datamash mean 1)
            model=$(echo ${{bname}}| cut -d- -f2)
            contig=$(basename ${{bname}} | cut -d- -f3)
            rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
            type="TRUTH"
            printf "${{val}},${{model}},${{contig}},${{rep}},${{type}}\n"

        done > {output}
        """


# rule get_mean_shared_nSites_run_ngsAMOVA_truth:
# 	input:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
# 	output:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/truth/stats/{simid}-{model_id}-{contig}-rep{rep}.mean_shared_nSites.csv"
# 	shell:
# 		"""
# 		infile=$(echo {input} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
# 		zcat ${{infile}} | cut -d, -f11 | datamash mean 1 > {output}
# 		"""


# rule collect_get_mean_ngsAMOVA_sites_used_run_ngsAMOVA_truth:
# 	input:
# 		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/truth/{simid}-{{model_id}}-{contig}-rep{rep}.amova.csv",
# 			simid=SIMULATION_ID,
# 			contig=CONTIGS,
# 			rep=REP),
# 	output:
# 		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_ngsAMOVA_sites_used_{model_id}.csv",
# 	shell:
# 		"""
# 		i=0
# 		for inamv in {input};do

# 			if [ ${{i}} -eq 0 ];then
# 				printf "nSites_processed,nSites_skipped,Model,Contig,Rep,AnalysisType\n"
# 			fi

# 			i=$((i+1))
# 			bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
# 			model=$(echo ${{bname}}| cut -d- -f2)
# 			contig=$(basename ${{bname}} | cut -d- -f3)
# 		done > {output}
# 		"""

# 				printf "nSites_processed,nSites_skipped,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
# 			fi

# 			i=$((i+1))
# 			bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
# 			model=$(echo ${{bname}}| cut -d- -f2)
# 			contig=$(basename ${{bname}} | cut -d- -f3)
# 			rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
# 			depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g')
# 			type="GTC"
# 			doPost=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
# 			postCutoff=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')

# 			if [[ ${{postCutoff}} == "03" ]];then
# 				postCutoff="0.3";
# 			elif [[ ${{postCutoff}} == "095" ]];then
# 				postCutoff="0.95";
# 			else
# 				exit 1;
# 			fi
# 			if [[ $(cat ${{inamv}} | grep -q "No shared sites") ]];then
# 				printf "NA,NA,${{model}},${{contig}},${{rep}},${{type}},${{doPost}},${{postCutoff}}\n"
# 				continue;
# 			else
# 				processed=$(cat ${{inamv}} |grep "Total number of sites processed"|cut -d: -f2|tr -d " ")
# 				skipped=$(cat ${{inamv}} |grep "Total number of sites skipped for all individual pairs"|cut -d: -f2|tr -d " ")
# 				printf "${{processed}},${{skipped}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"
# 			fi

# 		done > {output}
# 		"""

# # #  cat simulations/sim_demes_v2/logs/model_model2/contig_100/genotype_likelihood/ngsAMOVA_maxIter500_tole5/sim_demes_v2-model2-100-rep0-d1.amova.csv |grep "Total number of sites processed"|cut -d: -f2|tr -d " "
# # ################################################################################


#### BELOW NOTSURE 230411
# vcfgl="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"

# ###################################################
# # define populations
# #
# ploidy=2

# # define population with diploid individuals
# # <pop_id> : <n_individuals_in_pop>
# DEF_POPS={"B1C1":10,
# "B1C2":10,
# "B2C1":10,
# "B2C2":10}


# #note that msprime>1 and stdpopsim_dev version
# #assumes ploidy=2 by default

# haplo_list=[]
# indv_names=[]


# for i, (key, value) in enumerate(DEF_POPS.items()):
# 	haplo_list.append(value*ploidy)
# for ind in range(value):

# 	# Using PLINK-like format: <Family-ID>_<Individual-ID>
# 	# to store <Population-ID>_<Individual-ID>
# 	indv_names.append(f"pop{key}_ind{str(ind+1)}")

# IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(indv_names,2))


# mutation_rate=1.29e-08
# recombination_rate=1.14856e-08


# rule simulation:
# 	output:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
# 	params:
# 		seedfile="simulations/{simid}/model_{model_id}/contig_{contig}/trees/.seed.{simid}-{model_id}-{contig}-rep{rep}.trees",
# 	run:
# 		model=demes.load(str(wildcards.model_id)+".yaml")
# 		demography=msprime.Demography.from_demes(model)
# 		samples=DEF_POPS
# 		seedval=np.random.randint(2**20) + (100*int(wildcards.rep))  + int(wildcards.contig)
# 		with open(params.seedfile,"w") as seedout:
# 			print("Contig :"+str(wildcards.contig)+", rep:"+str(wildcards.rep),file=seedout)
# 			print("Seed value: "+str(seedval), file=seedout)
# 		sequence_length=int(wildcards.contig)*100000
# 		ts = msprime.sim_ancestry([msprime.SampleSet(n, population=p) for p,n in samples.items()],
# 		demography=demography,
# 		sequence_length=sequence_length,
# 		random_seed=seedval,
# 		recombination_rate=recombination_rate)
# 		mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=seedval)
# 		with open(output[0],"w") as tsout:
# 			mts.dump(tsout)


# ## Convert tree sequence to VCF file
# # Using legacy format to avoid multiple instances of sites
# rule tree_to_vcf:
# 	input:
# 		rules.simulation.output,
# 	output:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
# 	run:
# 		ts=tskit.load(input[0])
# 		with open(output[0],"w") as vcfout:
# 			ts.write_vcf(output=vcfout,
# 				contig_id=str(wildcards.contig),
# 				position_transform="legacy",
# 				individual_names=indv_names)
# 		# ploidy=ploidy,


# ## Using only variable sites
# rule vcf_to_vcfgl_var:
# 	input:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf"
# 	output:
# 		"simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
# 	params:
# 		vcfgl=vcfgl,
# 		prefix="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
# 		mode="b",
# 		error_rate=0.002,
# 	log:
# 		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
# 	shell:
# 		"""
# 		RNG_SEED=$(( {wildcards.rep} + 1 ))
# 		({params.vcfgl} -in {input} -out {params.prefix} \
# 			-depth {wildcards.depth} \
# 			-err {params.error_rate}  \
# 			-mode {params.mode} \
# 			-explode 0 \
# 			-seed ${{RNG_SEED}} ) 2> {log}
# 		"""


###############################################################################
###############################################################################
###############################################################################
### Simulation


rule simulation:
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}.trees",
    params:
        seedfile="simulations/{simid}/model_{model_id}/contig_{contig}/trees/{simid}-{model_id}-{contig}-rep{rep}_simulation_params.txt",
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


## Convert tree sequence to VCF file
# Using legacy format to avoid multiple instances of sites
rule tree_to_vcf:
    input:
        rules.simulation.output,
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    run:
        ts = tskit.load(input[0])
        with open(output[0], "w") as vcfout:
            ts.write_vcf(
                output=vcfout,
                contig_id=str(wildcards.contig),
                position_transform="legacy",
                individual_names=indv_names,
            )

            ## Using only variable sites



rule vcf_to_vcfgl_var:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcf/{simid}-{model_id}-{contig}-rep{rep}.vcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    params:
        vcfgl=vcfgl,
        prefix="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}",
        mode="b",
        error_rate=0.002,
    log:
        "simulations/{simid}/logs/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
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


rule get_vcfStats_ofvcfgl:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/vcf_stats/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.stats",
    shell:
        """
        bcftools stats {input} > {output}
        """


###############################################################################
###############################################################################
###############################################################################
### Genotype calling
##


rule remove_simulatedGTs_before_genotypeCalling:
    input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
    output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/rmGT/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_rmGT.bcf",
    shell:
        """
        bcftools annotate -x FORMAT/GT {input} -O b -o {output}
        """



###############################################################################
# Call fasta for each individuals and get NJ trees

rule angsd_doFasta:
	input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/doFasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	params:
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doFasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	shell:
		"""
        {ANGSD} -doFasta 2 -doCounts 1 -vcf-gl {input} -out {params.outprefix}
		"""


rule angsd_inferMajorMinor:
	input:
        "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	output:
        "simulations/{simid}/model_{model_id}/contig_{contig}/inferMajorMinor/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	params:
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/doFasta/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	shell:
		"""
        {ANGSD} -doMajorMinor 1 -vcf-gl {input} -doBcf 1 -out {params.outprefix}
		"""



###############################################################################
# END RULES
###############################################################################
