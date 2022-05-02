#!/bin/env python3.10
"""

SOURCE: https://github.com/tszfungc/scripts/blob/39801f1c1e78d72702633ce045b62b5d37d5d723/simulation/sim_demography.py

Function out_of_africa: user manual of msprime.
ooa_gutenkunst2009: original parameters in Gutenkunst 2009 Table 1.
ooa_gravel2011j: params updated according to Gravel 2011. Table 2.

| ID | Pop | 
|----|-----|
| 0  | YRI | 
| 1  | CEU | 
| 2  | CHB | 


"""

	# msprime | msms
	#
import msprime
import math
import sys


def ooa_gutenkunst2009(sample_n_AF, sample_n_EU, sample_n_AS, region_length):

    N_A = 7310
    N_B = 1861
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510

	# convert year to generation
    generation_time = 25
    T_AF = 220e3 / generation_time 
    T_B = 140e3 / generation_time 
    T_EU_AS = 21.2e3 / generation_time 

    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    # Set sample_size of YRI, CEU and CHB you want to output in the current generation
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_n_AF, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=sample_n_EU, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=sample_n_AS, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        msprime.MigrationRateChange(time=T_B, rate=0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dd.print_history()

    # set mutation_rate to you need to genotypes
    return msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        #  length=2.5e8,
        length=region_length,
        recombination_rate=1e-8,
        mutation_rate=1.25e-8
    )

def ooa_gravel2011(sample_n_AF, sample_n_EU, sample_n_AS, region_length):
    N_A = 7310
    N_B = 1861
    N_AF = 14474
    N_EU0 = 1032 
    N_AS0 = 554 
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time 
    T_B = 51e3 / generation_time 
    T_EU_AS = 23e3 / generation_time 
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.0038 
    r_AS = 0.0048 
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5 
    m_AF_EU = 2.5e-5 
    m_AF_AS = 0.78e-5 
    m_EU_AS = 3.11e-5 

    # 0=YRI, 1=CEU and 2=CHB

    # Set sample_size of YRI, CEU and CHB you want to output in the current generation
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_n_AF, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=sample_n_EU, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=sample_n_AS, initial_size=N_AS, growth_rate=r_AS)
    ]

    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]


    demographic_events = [

        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),

        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        msprime.MigrationRateChange(time=T_B, rate=0), 

        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
	
	# print the demographic history
    dd = msprime.DemographyDebugger(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dd.print_history()

    # set mutation_rate to you need to genotypes
    return msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        #  length=2.5e8,
        length=region_length,
        recombination_rate=1e-8,
        mutation_rate=1.25e-8
    )





#  output_file=sys.argv[1]
output_file = open(sys.argv[1],"w")

#  sample_n_AF=int(sys.argv[2])
#  sample_n_EU=int(sys.argv[3])
#  sample_n_AS=int(sys.argv[4])
region_length=int(float(sys.argv[5]))

sample_n_AF=int(sys.argv[2])*2
sample_n_EU=int(sys.argv[3])*2
sample_n_AS=int(sys.argv[4])*2

#  region_length=2.5e8,
#  length=region_length,
ts = ooa_gutenkunst2009(sample_n_AF, sample_n_EU, sample_n_AS,region_length)
ts.dump(output_file)





def ooa_gutenkunst2009(sample_n_AF, sample_n_EU, sample_n_AS, region_length):

	generation_time = 25
	T_OOA = 21.2e3 / generation_time
	T_AMH = 140e3 / generation_time
	T_ANC = 220e3 / generation_time
# We need to work out the starting population sizes based on
# the growth rates provided for these two populations
	r_CEU = 0.004
	r_CHB = 0.0055
	N_CEU = 1000 / math.exp(-r_CEU * T_OOA)
	N_CHB = 510 / math.exp(-r_CHB * T_OOA)
	N_AFR = 12300



	demography = msprime.Demography()
	demography.add_population(
			name="YRI",
			description="Yoruba in Ibadan, Nigeria",
			initial_size=N_AFR,
			)
	demography.add_population(
			name="CEU",
			description=(
				"Utah Residents (CEPH) with Northern and Western European Ancestry"
				),
			initial_size=N_CEU,
			growth_rate=r_CEU,
			)
	demography.add_population(
			name="CHB",
			description="Han Chinese in Beijing, China",
			initial_size=N_CHB,
			growth_rate=r_CHB,
			)
	demography.add_population(
			name="OOA",
			description="Bottleneck out-of-Africa population",
			initial_size=2100,
			)
	demography.add_population(
			name="AMH", description="Anatomically modern humans", initial_size=N_AFR
			)
	demography.add_population(
			name="ANC",
			description="Ancestral equilibrium population",
			initial_size=7300,
			)
#  demography


