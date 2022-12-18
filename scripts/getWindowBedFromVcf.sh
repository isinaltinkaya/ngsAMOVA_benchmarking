BCF=${1}
WINSIZE=${2}
BED1=${3}
BED2=${4}

printf "$(bcftools view -H ${BCF} | head -1 | cut -f1-2)\t$(bcftools view -H ${BCF}|tail -1|cut -f2)\n" > ${BED1}
bedtools makewindows -b ${BED1} -w ${WINSIZE} > ${BED2}
