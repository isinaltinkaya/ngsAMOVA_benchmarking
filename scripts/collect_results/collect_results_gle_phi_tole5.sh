printf "Measure,Level,Value,Model,Contig,Rep,Depth,Tole,AnalysisType\n";
type="GLE"


for infile in "${@}";do
	model=$(echo ${infile%.amova.csv}| cut -d- -f2)
	contig=$(echo ${infile%.amova.csv}| cut -d- -f3)
	rep=$(echo ${infile%.amova.csv}| cut -d- -f4|sed 's/rep//g')
	depth=$(echo ${infile%.amova.csv}| cut -d- -f5|sed 's/d//g'|cut -d_ -f1)
	tole=$(echo ${infile%.amova.csv} | cut -d- -f5 | cut -d_ -f2 | sed 's/tole//g')
	for line in `cat ${infile} | grep "^Phi"`;do
		printf "${line},${model},${contig},${rep},${depth},${tole},${type}\n";
	done


done
