

BED=${1}
FIN=${2}
OUTN=${BED%.bed}

NWIN=$(wc -l ${BED}|cut -d' ' -f1)

#
# for n in $(seq 1 ${NWIN});do
	# WIN=$(sed "${n}q;d" ${BED}|tr '\t' ' ')
	# CHR=$(echo ${WIN}|cut -d ' ' -f1)
	# START=$(echo ${WIN}|cut -d ' ' -f2)
	# END=$(echo ${WIN}|cut -d ' ' -f3)
#
#
	# REG=${CHR}":"${START}"-"${END}
	# WINID="win"${n}
	# echo ${REG} ${OUTN}"_"${WINID}
# done


touch ${FIN}

for n in $(seq 1 ${NWIN});do
	WIN=$(sed "${n}q;d" ${BED})
	WINID="win"${n}
	echo "${WIN}" > ${OUTN}"_"${WINID}

	printf "${WINID}\t" >> ${FIN}
	echo "${WIN}" >> ${FIN}

done



