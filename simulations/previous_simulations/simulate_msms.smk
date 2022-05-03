from itertools import product, combinations

POP_IND_ID=list( "_".join(comb) for comb in product(list(map('pop{}'.format,range(1,4))),list(map('ind{}'.format,range(1,11)))))
IND_PAIRS=list("-".join(map(str,comb)) for comb in combinations(POP_IND_ID,2))


DEPTH=[100,20,10,5,2,1,0.5,0.1]

#0.1 to 5
#DEPTH=[100,20,10,0.01]+ [ x / pow(.1, -1) for x in range(1, 50 + 1) ]

simulation_id="sim2"
#simulation_id="sim2"+"_L"+config['region_length']

REGION_LENGTH=range(5,70,5)
REPS=range(1,20,1)

#############################################################################
#### MSMS Simulations
#############################################################################

def calc_theta(L=5130456):
	Ne=7310
	generation_time=25.8
	mu=0.5e-9
	theta=4*Ne*mu*L*generation_time
	return round(theta)

def calc_rho(L):
	Ne=7310
	recombination_rate=1.5e-8
	rho=4*Ne*recombination_rate*L
	return round(rho)


DEF_POPS= {
		"1" : "1 10",
		"2" : "11 20",
		"3" : "21 30"
		}

dadi_params="-I 3 20 20 20 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1"

#############################################################################




rule all:
	input:
		expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/{sid}_{srlen}MB_rep{rep}_d{depth}.glf.gz",
				rep=REPS,
				srlen=REGION_LENGTH,
				sid=simulation_id,
				depth=DEPTH),
		expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/pops/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}.glf.gz",
				rep=REPS,
				srlen=REGION_LENGTH,
				sid=simulation_id,
				depth=DEPTH,
				pop_id=DEF_POPS.keys()),
		expand("simulations/{sid}_{srlen}MB_rep{rep}/glf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.glf.gz",
				rep=REPS,
				srlen=REGION_LENGTH,
				sid=simulation_id,
				depth=DEPTH,
				pop_id=DEF_POPS.keys(),
				ind_id=range(1,11)),
		expand("simulations/{sid}_{srlen}MB_rep{rep}/saf/inds/{sid}_{srlen}MB_rep{rep}_d{depth}_pop{pop_id}_ind{ind_id}.saf.idx",
				rep=REPS,
				srlen=REGION_LENGTH,
				sid=simulation_id,
				depth=DEPTH,
				pop_id=DEF_POPS.keys(),
				ind_id=range(1,11)),
		expand("simulations/{sid}_{srlen}MB_rep{rep}/sfs/ind_pairs/{sid}_{srlen}MB_rep{rep}_d{depth}_indpair_{pop_ind_pair}.sfs",
				rep=REPS,
				srlen=REGION_LENGTH,
				sid=simulation_id,
				depth=DEPTH,
				pop_ind_pair=IND_PAIRS),

rule run_msms:
	output:
		"simulations/{sid}_{srlen}MB_rep{rep}/ms/{sid}_{srlen}MB_rep{rep}.ms"
	params:
		msms_params=dadi_params,
		msms="/willerslev/software/msms/bin/msms",
		nHaplotypes=60,
		nReplicates=1,
		rlen=lambda wildcards: int(float(wildcards.srlen)*1e6),
		theta=lambda wildcards: str(calc_theta(float(wildcards.srlen)*1e6)),
		rho=lambda wildcards: str(calc_rho(float(wildcards.srlen)*1e6)),

	log:
		"simulations/{sid}_{srlen}MB_rep{rep}/logs/ms/{sid}_{srlen}MB_rep{rep}.ms"
	shell:
		"""
		({params.msms} {params.nHaplotypes} {params.nReplicates} -t {params.theta} -r {params.rho} {params.rlen} {params.msms_params} > {output}) 2> {log}
		"""


rule ms_to_glf:
	input:
		"simulations/{sid}_{srlen}MB_rep{rep}/ms/{sid}_{srlen}MB_rep{rep}.ms"
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
