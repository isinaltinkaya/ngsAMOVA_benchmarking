

out=simulations/sim_demes_v2/mean_n_variable_sites.csv

printf "model,contig,mean_n_variable_sites\n" > ${out}
for model in "model1" "model2";do
	for contig in 1 10 100;do

		printf "${model},${contig},$(sed 1d simulations/sim_demes_v2/results/collect_results_n_variable_sites_simulated.csv| grep "^${model},${contig},"|cut -d, -f4|datamash mean 1)\n" >> ${out}
	done
done 


