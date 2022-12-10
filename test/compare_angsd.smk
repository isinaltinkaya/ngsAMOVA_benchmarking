import pandas as pd
import matplotlib.pyplot as plt

SAMPLE = ["popB1C1_ind1", "popB2C1_ind1"]
BCFTOOLS = "/workspaces/test_env/tools/bcftools/bcftools"
ANGSD = "/workspaces/test_env/tools/angsd/angsd"
indat = "sim_demes-model1-1-rep1-d1.bcf"
realSFS = "/workspaces/test_env/tools/angsd/misc/realSFS"
ngsAMOVA="/workspaces/test_env/tools/ngsAMOVA/ngsAMOVA"


rule all:
    input:
        expand("data/{data}.sites.txt.idx", data=[indat]),
        expand("data/{sample}.vcf", sample=SAMPLE),
        expand("results/{sample}.saf.idx", sample=SAMPLE),
        expand("results/{sample1}_{sample2}.sfs",
               sample1=str(SAMPLE[0]),
               sample2=str(SAMPLE[1])),
        expand("results/{sample1}_{sample2}.ngsAMOVA.sfs.csv",
               sample1=str(SAMPLE[0]),
               sample2=str(SAMPLE[1])),
        expand("results/{sample1}_{sample2}.diff.sfs.csv",
                sample1=str(SAMPLE[0]),
                sample2=str(SAMPLE[1])),
        
  


rule bcftools_get_sample:
    input:
        indat,
    output:
        "data/{sample}.vcf",
    shell:
        """
        {BCFTOOLS} view -s {wildcards.sample} {input} > {output}
        """


rule bcftools_get_sites_forANGSD:
    input:
        indat,
    output:
        "data/" + indat + ".sites.txt",
    shell:
        # do not print all ALTs, only the first ALT allele
        """    
        {BCFTOOLS} query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input} | cut -d, -f1 > {output}
        """


rule angsd_index_sites:
    input:
        "data/" + indat + ".sites.txt",
    output:
        "data/" + indat + ".sites.txt.idx",
    shell:
        """
        {ANGSD} sites index {input} -o {output}
        """


rule angsd_doSaf:
    input:
        vcf="data/{sample}.vcf",
        sites="data/" + indat + ".sites.txt",
    output:
        "results/{sample}.saf.idx",
    params:
        prefix="results/{sample}",
    log:
        "logs/{sample}.saf.idx",
    shell:
        """
        (
        {ANGSD} -vcf-gl {input.vcf} -sites {input.sites} -doSaf 5 -doMajorMinor 3 -out {params.prefix}
        ) 2> {log}
        """


rule realSFS_pair:
    input:
        expand("results/{sample}.saf.idx",
                sample=SAMPLE),
    output:
        "results/{sample1}_{sample2}.sfs",
    log:
        "logs/{sample1}_{sample2}.sfs",
    shell:
        """
        (
        {realSFS} -m 0 -cores 1 {input} > {output}
        )2> {log}
        """



rule bcftools_get_pair:
    input:
        indat,
    output:
        "data/"+SAMPLE[0]+"_" + SAMPLE[1]+".vcf",
    params:
        sample1=SAMPLE[0],
        sample2=SAMPLE[1],
    shell:
        """
        {BCFTOOLS} view -s {params.sample1},{params.sample2} {input} > {output}
        """


rule ngsAMOVA:
    input:
        # "data/"+SAMPLE[0]+"_" + SAMPLE[1]+".vcf",
        "data/{sample1}_{sample2}.vcf"
    output:
        # "results/"+SAMPLE[0]+"_" + SAMPLE[1]+".ngsAMOVA.sfs.csv",
        "results/{sample1}_{sample2}.ngsAMOVA.sfs.csv",
    params:
        # prefix="results/"+SAMPLE[0]+"_" + SAMPLE[1]+".ngsAMOVA",
        prefix="results/{sample1}_{sample2}.ngsAMOVA",
    shell:
        """
        {ngsAMOVA} -in {input} -out {params.prefix} -doAMOVA -1 -isSim 1 -minInd 2 -doDist 1
        """




# def read_realSFS_sfs_to_df(sfs_file):
#     with open(sfs_file) as f:
#         return pd.DataFrame([f.read().strip().split(" ")])

# # print(read_realSFS_sfs_to_df("results/"+SAMPLE[0]+"_" + SAMPLE[1]+".sfs"))
# def read_ngsAMOVA_sfs_to_df(sfs_file):
#     with open(sfs_file) as f:
#         #only return cols 3-11
#         return pd.DataFrame([f.read().strip().split(",")[3:12]])
#         # return pd.DataFrame([f.read().strip().split(",")]).)

# # print(read_ngsAMOVA_sfs_to_df("results/"+SAMPLE[0]+"_" + SAMPLE[1]+".ngsAMOVA.sfs.csv"))

# def diff_sfs(realSFS_sfs, ngsAMOVA_sfs):
#     #convert str to float
#     realSFS_sfs = realSFS_sfs.astype(float)
#     ngsAMOVA_sfs = ngsAMOVA_sfs.astype(float)
#     #write to file
#     realSFS_sfs.to_csv("results/"+SAMPLE[0]+"_" + SAMPLE[1]+".diff.sfs.csv", index=False)

# rule get_diff:
#     input:
#         realSFS_sfs="results/{sample1}_{sample2}.sfs",
#         ngsAMOVA_sfs="results/{sample1}_{sample2}.ngsAMOVA.sfs.csv"
#     output: 
#         "results/{sample1}_{sample2}.diff.sfs.csv"
#     run:
#         diff_sfs(read_realSFS_sfs_to_df(input.realSFS_sfs),read_ngsAMOVA_sfs_to_df(input.ngsAMOVA_sfs))


rule plot_diff_R:
    input:
        realSFS_sfs="results/{sample1}_{sample2}.sfs",
        ngsAMOVA_sfs="results/{sample1}_{sample2}.ngsAMOVA.sfs.csv"
        # "results/{sample1}_{sample2}.diff.sfs.csv"
    output:
        "results/{sample1}_{sample2}.diff.sfs.csv"
    shell:
        """
        Rscript plot_test.R {input.realSFS_sfs} {input.ngsAMOVA_sfs} {output}
        """

# use genotypes in vcf file to count SFS
rule angsd_doSaf:
    input:
        vcf="data/{sample}.vcf",
        sites="data/" + indat + ".sites.txt",
    output:
        "results/{sample}.saf.idx",
    params:
        prefix="results/{sample}",
    log:
        "logs/{sample}.saf.idx",
    shell:
        """
        (
            
        ) 2> {log}
        """
