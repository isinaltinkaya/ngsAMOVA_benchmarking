
printf "Measure,Level,Value,Model,Contig,Rep,Depth,AnalysisType,doPost,postCutoff\n"
type="GTC"

for infile in "${@}";do
	model=$(basename ${infile%.amova.csv}| cut -d- -f2)
	contig=$(basename ${infile%.amova.csv}| cut -d- -f3)
	rep=$(basename ${infile%.amova.csv}| cut -d- -f4 | sed 's/rep//g')
	depth=$(basename ${infile%.amova.csv}| cut -d- -f5 | sed 's/d//g')
	doPost=$(dirname ${infile}| rev | cut -d/ -f2 | rev | cut -d_ -f2 | sed 's/doPost//g')
	postCutoff=$(dirname ${infile}| rev | cut -d/ -f2 | rev | cut -d_ -f1 | sed 's/postCutoff//g')

	if [[ ${postCutoff} == "03" ]];then
		postCutoff="0.3";
	elif [[ ${postCutoff} == "095" ]];then
		postCutoff="0.95";
	else
		exit 1;
	fi


	if [[ $(cat ${infile}) == "0" ]];then
		printf "Phi,Population_in_Region,NoResult,${model},${contig},${rep},${depth},${type},${doPost},${postCutoff}\n";
		printf "Phi,Population_in_Total,NoResult,${model},${contig},${rep},${depth},${type},${doPost},${postCutoff}\n";
		printf "Phi,Region_in_Total,NoResult,${model},${contig},${rep},${depth},${type},${doPost},${postCutoff}\n";
	else

		for line in `cat ${infile} | grep "^Phi"`;do
			printf "${line},${model},${contig},${rep},${depth},${type},${doPost},${postCutoff}\n";
		done
	fi


done
