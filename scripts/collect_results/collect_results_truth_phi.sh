

printf "Measure,Level,Value,Model,Contig,Rep,AnalysisType\n";
type="Truth"

for infile in "${@}";do
	model=$(echo ${infile%.amova.csv}| cut -d- -f2)
	contig=$(echo ${infile%.amova.csv}| cut -d- -f3)
	rep=$(echo ${infile%.amova.csv}| cut -d- -f4|sed 's/rep//g')
	for line in `cat ${infile} | grep "^Phi"`;do
		printf "${line},${model},${contig},${rep},${type}\n";
	done


done
