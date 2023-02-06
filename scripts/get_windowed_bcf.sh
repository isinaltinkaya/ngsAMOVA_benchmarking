contig=$(basename ${1}|cut -d- -f2|cut -d_ -f2)


iter=0
while read -r i j;do

	outfile=${3}_${iter}.bcf
	bcftools view ${2} -r ${contig}:${i}-${j} -O b -o ${outfile} &
	((iter++))

done < "${1}"
