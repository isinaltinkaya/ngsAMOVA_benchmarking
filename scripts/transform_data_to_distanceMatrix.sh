for input in simulations/sim_demes_v2/model_*/contig_*/perIterationDist/sim_demes_v2-model*-*-rep*-d*_tole10.perIterDistances.csv_emIter*.csv;do
	out=$(dirname ${input})/$(basename ${input}).distance_matrix.csv
	cut -d, -f1,3 ${input} |sort -t, -k1|cut -d, -f2|datamash transpose --output-delimiter=, > ${out}
	echo ${out}
done
