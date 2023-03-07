
printf "Measure,Level,Value,Model,Contig,Rep,Depth,AnalysisType\n";
type="GTC"

for infile in "${@}";do
	model=$(basename ${infile%.amova.csv}| cut -d- -f2)
	contig=$(basename ${infile%.amova.csv}| cut -d- -f3)
	rep=$(basename ${infile%.amova.csv}| cut -d- -f4|sed 's/rep//g')
	depth=$(basename ${infile%.amova.csv}| cut -d- -f5|sed 's/d//g')

	if [[ $(cat ${infile}) == "0" ]];then
		printf "Variance_component,Region,NA,${model},${contig},${rep},${depth},${type}\n";
		printf "Variance_component,Population,NA,${model},${contig},${rep},${depth},${type}\n";
		printf "Variance_component,Individual,NA,${model},${contig},${rep},${depth},${type}\n";
	else

		for line in `cat ${infile} | grep "^Variance_component"`;do
			printf "${line},${model},${contig},${rep},${depth},${type}\n";
		done
	fi


done

