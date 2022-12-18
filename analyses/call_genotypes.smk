

## Disable warnings
# There were bunch of warnings due to msprime version issues
import warnings
warnings.filterwarnings('ignore')


###################################################
# CONFIG

SIMULATION_ID="sim_demes_v2"
ANGSD="/maps/projects/lundbeck/scratch/pfs488/Programs/angsd/angsd"

###################################################
# DEPTH 

# average per site depth
DEPTH=[20,10,5,2,1,0.5,0.2,0.1]

###################################################


# Number of replicates
n_reps=200
REP=[*range(n_reps)]
REP=REP[:20]



# END CONFIG
###################################################

CONTIGS=[1,2,10,50,100]

MODELS=["model1","model2"]

rule all:
	input:
		expand("simulations/{simid}/resources/{contig}.sites",
				simid=SIMULATION_ID,
				contig=CONTIGS),
		expand("simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.geno.gz",
				simid=SIMULATION_ID,
				model_id=MODELS,
				contig=CONTIGS,
				depth=DEPTH,
				rep=REP),


# rule generate_reference:
	# output:
		# "simulations/{simid}/resources/{contig}.fa"
	# params:
		# contig_length=lambda wildcards: int(int(wildcards.contig)*1e6)
	# shell:
		# """
		# dd if=/dev/zero bs=1 count={params.contig_length} | tr '\\0' 'A' > {output}
		# sed -i -e '$a\\' {output}
		# sed -i -e '1s/^/>{wildcards.contig}\\n/' {output}
		# """

# rule index_reference:
	# input:
		# "simulations/{simid}/resources/{contig}.fa"
	# output:
		# "simulations/{simid}/resources/{contig}.fa.fai"
	# shell:
		# """
		# samtools faidx {input}
		# """

# rule generate_sites:
	# input:
		# "simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
	# output:
	# shell:
		# """
		# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n'  sim_demes_v2-model2-1-rep0-d1.bcf > sites.txt
		# """

rule generate_sites_contig:
	output:
		"simulations/{simid}/resources/{contig}.sites"
	params:
		contig_length=lambda wildcards: int(int(wildcards.contig)*1e6)
	shell:
		"""
		awk -v X={wildcards.contig} -v XN={params.contig_length} 'BEGIN{{for(c=1;c<XN;c++) printf "%s\\t%d\\tA\\tC\\n",X,c}}' > {output};
		{ANGSD} sites index {output}
		"""



###############################################################################
# ./angsd
# -doMajorMinor 3 # Use pre-specified major and minor
# -sites sites.txt
# -doMaf 1  # Calculate frequencies using fixed major and minor
# -doGeno 31 # 1 + 2 + 4 + 8 + 16 [below]
# -doPost 1 # Estimate the posterior genotype probability based on allele frequency as prior
# -vcf-gl {input.vcf}
###############################################################################
## doGeno
# 1: print out major minor
# 2: print the called genotype as -1,0,1,2 (count of minor)
# 4: print the called genotype as AA, AC, AG, ...
# 8: print all 3 posts (major,major),(major,minor),(minor,minor)
# 16: print the posterior of the called genotype
###############################################################################
rule angsd_call_genotypes:
	input:
		vcf="simulations/{simid}/model_{model_id}/contig_{contig}/vcfgl_var/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.bcf",
		sites="simulations/{simid}/resources/{contig}.sites"
	output:
		"simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.geno.gz"
	log:
		"simulations/{simid}/logs/model_{model_id}/contig_{contig}/call_genotype/{simid}-{model_id}-{contig}-rep{rep}-d{depth}.geno.gz"
	params:
		prefix="simulations/{simid}/model_{model_id}/contig_{contig}/call_genotype/{simid}-{model_id}-{contig}-rep{rep}-d{depth}"
	shell:
		"""
		({ANGSD} -doMajorMinor 3 -doMaf 1 -doGeno 31 -doPost 1 -vcf-gl {input.vcf} -sites {input.sites} -out {params.prefix} -doBcf 1 )2>{log}
		"""

