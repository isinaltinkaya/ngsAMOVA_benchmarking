import numpy as np
import pandas as pd
import sys
import os
import subprocess
import csv




## Disable warnings (due to msprime version issues)
import warnings
warnings.filterwarnings("ignore")

###############################################################################
# BEGIN CONFIG
###############################################################################
SIMULATION_ID = "sim_demes_v2"

# paths to programs
ANGSD = "/maps/projects/lundbeck/scratch/pfs488/Programs/angsd/angsd"
vcfgl = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_vcfgl/vcfgl/vcfgl"
ngsAMOVA = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA"

# simulation models
MODELS = ["model1", "model2"]

# average per site depth
DEPTH = [10, 5, 2, 1, 0.5, 0.2, 0.1]

# contig id x corresponds to contig length x * 1e6
CONTIGS = [1, 2, 10, 50, 100]

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
		# expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_ngsAMOVA_sites_used_{model_id}.csv",model_id=MODELS),
		# expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_truth_mean_shared_nSites_{model_id}.csv", model_id=MODELS),
		# expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_nSites_used_{model_id}.csv",model_id=MODELS),
		# expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_shared_nSites_{model_id}.csv",model_id=MODELS),
		# expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_mean_shared_nSites_{model_id}.csv",model_id=MODELS),
		# #######################################################################################
		# # truth
		# # expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_truth_phi_{model_id}.csv", model_id=MODELS),
		# # expand("simulations/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
		# # 	simid=SIMULATION_ID,
		# # 	model_id=MODELS,
		# # 	contig=CONTIGS,
		# # 	rep=REP),
		# # #######################################################################################
		# # # genotype calling
		# # expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_phi_{model_id}.csv",model_id=MODELS),
		# # expand("simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
		# # 	simid=SIMULATION_ID,
		# # 	model_id=MODELS,
		# # 	contig=CONTIGS,
		# # 	depth=DEPTH,
		# # 	rep=REP,
		# # 	postCutoff=postCutoffDict.keys(),
		# # 	doPost=["1","2"]),
		# # expand("simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		# # 	simid=SIMULATION_ID,
		# # 	model_id=MODELS,
		# # 	contig=CONTIGS,
		# # 	depth=DEPTH,
		# # 	rep=REP,
		# # 	postCutoff=postCutoffDict.keys(),
		# # 	doPost=["1","2"]),
		# # #######################################################################################
		# # # genotype likelihood
	 	expand("simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_phi_{model_id}.csv",model_id=MODELS),
		# # expand("simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
		# # 	simid=SIMULATION_ID,
		# # 	model_id=MODELS,
		# # 	contig=CONTIGS,
		# # 	depth=DEPTH,
		# # 	rep=REP),



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
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	threads: 4
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.amova.csv",
	shell:
		"""
        ({ngsAMOVA} -i {input.bcf} -P {threads} -out {params.outprefix} -doAMOVA 1 -doEM 1 -doDist 1 -m {input.metadata} -printMatrix 3 -sqDist 1 -pJGCD 3 -maxIter 500 -tole 1e-5 )  2> {log}
		"""


collect_results_gle_phi_tole5 = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/collect_results/collect_results_gle_phi_tole5.sh"
# 230306
# collect genotype likelihood results (tole=1e-5)
rule collect_results_run_ngsAMOVA_genotype_likelihood_2level:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP),
	output:
		phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_phi_{model_id}.csv",
	shell:
		"""
			bash {collect_results_gle_phi_tole5} {input} > {output.phi}
		"""

rule collect_get_mean_shared_nSites_run_ngsAMOVA_genotype_likelihood_2level:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_likelihood/ngsAMOVA_maxIter500_tole5/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP),
	output:
		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_likelihood_mean_shared_nSites_{model_id}.csv"
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
		



# rule collect_get_nSites_used_run_ngsAMOVA_genotype_likelihood_2level:

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

# 230306
rule angsd_index_sites_contig:
	input:
		"simulations/{simid}/resources/{contig}.sites"
	output:
		"simulations/{simid}/resources/{contig}.sites.idx"
	params:
		contig_length=lambda wildcards: int(int(wildcards.contig)*1e6)
	shell:
		"""
		{ANGSD} sites index {input}
		"""


###################
# 230306 rerun genotype calling
# use {postCutoff} {doPost}
rule angsd_genotype_calling_postCutoffvar_doPostvar:
	input:
		vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/rmGT/{simid}-{model_id}-{contig}-rep{rep}-d{depth}_rmGT.bcf",
		sites="simulations/{simid}/resources/{contig}.sites",
		sites_idx="simulations/{simid}/resources/{contig}.sites.idx",
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf"
	params:
		cutoffarg=lambda wildcards: str(postCutoffDict[str(wildcards.postCutoff)]),
		prefix="simulations/{simid}/model_{model_id}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	shell:
		"""
		({ANGSD} -doMajorMinor 3 -doMaf 1 -doGeno 31 -vcf-gl {input.vcf} -sites {input.sites} -out {params.prefix} -doBcf 1 -doPost {wildcards.doPost} {params.cutoffarg} )2>{log}
		"""


# 230306 rerun genotype calling amova
# run ngsAMOVA on results from angsd_call_genotypes_postCutoff095_doPost2
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
		p=subprocess.Popen(ngsAMOVA + " -i "+ input.bcf + " -o "+ params.outprefix + " -doAMOVA 2 -doDist 1 -printMatrix 3 -sqDist 1 -pJGCD 3 -m "+ input.metadata, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		stdout, stderr=p.communicate()
		log_file=open(log[0], "w")
		if p.returncode != 0:
			if "No shared sites found for pair" in str(stderr):
				shell("echo 0 > "+ str(output[0]))
				print(str(stderr,"utf-8"), file=log_file)
			else:
				print(str(stderr,"utf-8"), file=log_file)
				raise Exception("ngsAMOVA failed with error: "+ str(stderr))
		else:
			log_file.write(str(stderr,"utf-8"))
		log_file.close()




collect_results_genotype_call_doPost_postCutoff_phi = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/collect_results/collect_results_genotype_call_doPost_postCutoff_phi.sh"
# 230306
rule collect_results_run_ngsAMOVA_genotype_call_postCutoff_doPost:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP,
			postCutoff=postCutoffDict.keys(),
			doPost=["1","2"]),
	output:
		phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_phi_{model_id}.csv",
	shell:
		"""
		bash {collect_results_genotype_call_doPost_postCutoff_phi} {input} > {output.phi}
		"""



rule collect_get_mean_shared_nSites_run_ngsAMOVA_genotype_call_postCutoff_doPost:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP,
			postCutoff=postCutoffDict.keys(),
			doPost=["1","2"]),
	output:
		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_shared_nSites_{model_id}.csv",
	shell:
		"""
		i=0
		for inamv in {input};do

			if [ ${{i}} -eq 0 ];then
				printf "mean_shared_nSites,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
			fi

			i=$((i+1))
			bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
			infile=$(echo ${{inamv}} | sed 's/.amova.csv//g')".joint_geno_count_dist.csv.bgz"
			model=$(echo ${{bname}}| cut -d- -f2)
			contig=$(basename ${{bname}} | cut -d- -f3)
			rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
			depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g')
			type="GTC"
			doPost=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
			postCutoff=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')

			if [[ ${{postCutoff}} == "03" ]];then
				postCutoff="0.3";
			elif [[ ${{postCutoff}} == "095" ]];then
				postCutoff="0.95";
			else
				exit 1;
			fi

			if [[ $(cat ${{inamv}}) == "0" ]];then
				val="NA"
			else
				val=$(zcat ${{infile}} | cut -d, -f11 | datamash mean 1)
			fi
			printf "${{val}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"

		done > {output}
		"""



rule collect_get_nSites_used_run_ngsAMOVA_genotype_call_postCutoff_doPost:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.mafs.gz",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP,
			postCutoff=postCutoffDict.keys(),
			doPost=["1","2"]),
	output:
		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_nSites_used_{model_id}.csv",
	shell:
		"""
		i=0
		for inmaf in {input};do

			val=$(zcat ${{inmaf}} | awk 'END{{print NR-1}}')

			if [ ${{i}} -eq 0 ];then
				printf "mean_nSites_used,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
			fi

			i=$((i+1))
			bname=$(basename ${{inmaf}}|sed 's/.mafs.gz//g')
			model=$(echo ${{bname}}| cut -d- -f2)
			contig=$(basename ${{bname}} | cut -d- -f3)
			rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
			depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g'|cut -d_ -f1)
			type="GTC"
			doPost=$(dirname ${{inmaf}} | rev | cut -d/ -f1 | rev | cut -d_ -f2 | sed 's/doPost//g')
			postCutoff=$(dirname ${{inmaf}} | rev | cut -d/ -f1 | rev | cut -d_ -f1 | sed 's/postCutoff//g')

			if [[ ${{postCutoff}} == "03" ]];then
				postCutoff="0.3";
			elif [[ ${{postCutoff}} == "095" ]];then
				postCutoff="0.95";
			else
				postCutoff="NA";
			fi

			printf "${{val}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"
		done > {output}
		"""


rule collect_get_mean_ngsAMOVA_sites_used_run_ngsAMOVA_genotype_call_postCutoff_doPost:
	input:
		expand("simulations/{simid}/logs/model_{{model_id}}/contig_{contig}/genotype_calling/postCutoff{postCutoff}_doPost{doPost}/ngsAMOVA/{simid}-{{model_id}}-{contig}-rep{rep}-d{depth}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			depth=DEPTH,
			rep=REP,
			postCutoff=postCutoffDict.keys(),
			doPost=["1","2"]),
	output:
		"simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_genotype_calling_mean_ngsAMOVA_sites_used_{model_id}.csv",
	shell:
		"""
		i=0
		for inamv in {input};do

			if [ ${{i}} -eq 0 ];then
				printf "nSites_processed,nSites_skipped,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
			fi

			i=$((i+1))
			bname=$(basename ${{inamv}}|sed 's/.amova.csv//g')
			model=$(echo ${{bname}}| cut -d- -f2)
			contig=$(basename ${{bname}} | cut -d- -f3)
			rep=$(echo ${{bname}} | cut -d- -f4|sed 's/rep//g')
			depth=$(echo ${{bname}}| cut -d- -f5 | sed 's/d//g')
			doPost=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
			postCutoff=$(dirname ${{inamv}} | rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')
			type="GTC"

			if [[ ${{postCutoff}} == "03" ]];then
				postCutoff="0.3";
			elif [[ ${{postCutoff}} == "095" ]];then
				postCutoff="0.95";
			else
				exit 1;
			fi
			if [[ $(cat ${{inamv}} | grep -q "No shared sites") ]];then
				processed="NA"
				skipped="NA"
			else
				processed=$(cat ${{inamv}} |grep "Total number of sites processed"|cut -d: -f2|tr -d " ")
				skipped=$(cat ${{inamv}} |grep "Total number of sites skipped for all individual pairs"|cut -d: -f2|tr -d " ")
			fi
			printf "${{processed}},${{skipped}},${{model}},${{contig}},${{rep}},${{depth}},${{type}},${{doPost}},${{postCutoff}}\n"
		
		done > {output}
		"""

		

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
		outprefix="simulations/{simid}/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}"
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/truth/{simid}-{model_id}-{contig}-rep{rep}.amova.csv",
	shell:
		"""
		(
			{ngsAMOVA} -i {input.vcf} -o {params.outprefix} -doAMOVA 2 -doDist 1 -m {input.metadata} -printMatrix 3 -pJGCD 3 -sqDist 1
		)  2> {log}
		"""




collect_results_truth_phi = "/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/AMOVA_paper_analyses/scripts/collect_results/collect_results_truth_phi.sh"
rule collect_results_run_ngsAMOVA_truth:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/truth/{simid}-{{model_id}}-{contig}-rep{rep}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			rep=REP),
	output:
		phi="simulations/sim_demes_v2/collected_results/collect_march23/sim_demes_v2_truth_phi_{model_id}.csv",
	shell:
		"""
		bash {collect_results_truth_phi} {input} > {output.phi}
		"""


rule collect_get_mean_shared_nSites_run_ngsAMOVA_truth:
	input:
		expand("simulations/{simid}/model_{{model_id}}/contig_{contig}/truth/{simid}-{{model_id}}-{contig}-rep{rep}.amova.csv",
			simid=SIMULATION_ID,
			contig=CONTIGS,
			rep=REP),
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

