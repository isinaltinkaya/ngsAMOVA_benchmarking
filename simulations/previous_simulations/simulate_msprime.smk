from itertools import product, combinations

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

DEF_POPS= {
		"1" : "1 10",
		"2" : "11 20",
		"3" : "21 30"
		}

rule all:
	input:
		expand("simulations/{sid}_{srlen}MB_rep{rep}/msprime/trees/{sid}_{srlen}MB_rep{rep}.trees",
				rep=REP,
				sid=SIMULATION_ID,
				srlen=REGION_LENGTH),
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
#

rule msprime_simulation:
	output: 
		"simulations/{sid}_{srlen}MB_rep{rep}/msprime/trees/{sid}_{srlen}MB_rep{rep}.trees"
	params:
		rlen=lambda wildcards: int(float(wildcards.srlen)*1e6),
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/msprime/{sid}_{srlen}MB_rep{rep}.ms"
	conda:
		conda_env
	shell:
		"""
		( python3.10 scripts/sim_demography.py {output} 10 10 10 {params.rlen} )2> {log}
		"""

rule msprime2ms:
	output: 
		"simulations/{sid}_{srlen}MB_rep{rep}/msprime/ms/{sid}_{srlen}MB_rep{rep}.ms"
	input: 
		rules.msprime_simulation.output
	conda:
		conda_env
	shell:
		"""
		python3.10 scripts/write_ms.py {input} {output}
		"""




rule ms_to_glf:
	input:
		"simulations/{sid}_{srlen}MB_rep{rep}/msprime/ms/{sid}_{srlen}MB_rep{rep}.ms"
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz"
	params:
		msToGlf="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/msToGlf",
		outbase="simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}",
		seed=42,
		error_rate=0.002,
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz.log"
	shell:
		"""
		({params.msToGlf} -in {input} -out {params.outbase} \
				-singleOut 1 \
				-depth {wildcards.depth} \
				-err {params.error_rate}  \
				-nSites 0 \
				-seed {params.seed}) 2> {log}
		"""

rule split_glf_to_pops:
	input:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz"
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz",
	params:
		pop_size=30,
		pop_range=lambda wildcards: DEF_POPS[wildcards.pop_id],
		sim_glf="simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz",
		outbase="simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}",
		splitgl="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/splitgl",
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz.log",
	shell:
		"""
		({params.splitgl} {params.sim_glf} {params.pop_size} {params.pop_range} > {params.outbase}.glf.gz ) 2> {log}
		"""


rule split_glf_to_inds:
	input:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz"
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz"
	params:
		pop_size=10,
		ind_range="{ind_id} {ind_id}",
		splitgl="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/angsd/misc/splitgl",
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz.log"
	shell:
		"""
		({params.splitgl} {input} {params.pop_size} {params.ind_range} > {output}) 2> {log}
		"""



rule doSaf_perInd:
	input:
		"simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz"
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.saf.idx"
	params:
		outbase="simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}",
		angsd="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/angsd/angsd",
		ref="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/ms/ref.fa",
		fai="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/ms/ref.fa.fai",
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz.log"
	shell:
		"""
		({params.angsd} \
				-glf {input} \
				-doSaf 1 \
				-nInd 1 \
				-ref {params.ref} \
				-fai {params.fai} \
				-isSim 1 \
				-out {params.outbase}) 2> {log}

		"""


rule realSFS_2dsfs_ind2ind:
	input:
		ind1=lambda wildcards: "simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_"+ wildcards.pop_ind_pair.split("-")[0] + ".saf.idx",
		ind2=lambda wildcards: "simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_"+ wildcards.pop_ind_pair.split("-")[1] + ".saf.idx",
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs",
	params:
		realSFS="/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/AMOVA/wd/angsd/misc/realSFS",
	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs.log",
	shell:
		"""
		({params.realSFS} {input.ind1} {input.ind2} > {output} )2>{log}
		"""
