

printf "Measure,Level,Value,Model,Contig,Rep,Depth,AnalysisType\n";
type="GTC"

for infile in "${@}";do
	model=$(echo ${infile%.amova.csv}| cut -d- -f2)
	contig=$(echo ${infile%.amova.csv}| cut -d- -f3)
	rep=$(echo ${infile%.amova.csv}| cut -d- -f4|sed 's/rep//g')
	depth=$(echo ${infile%.amova.csv}| cut -d- -f5|sed 's/d//g')
	for line in `cat ${infile} | grep "^Variance_component"`;do
		printf "${line},${model},${contig},${rep},${depth},${type}\n";
	done


done
