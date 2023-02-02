ngsAMOVA="/maps/projects/lundbeck/scratch/pfs488/AMOVA/runv_ngsAMOVA/ngsAMOVA"

metadata="/maps/projects/lundbeck/scratch/pfs488/AMOVA/simulations/msprime/simulations/sim_demes_v2/metadata_2lvl_with_header.tsv"

while IFS= read -r line; do

    outprefix=${line%.csv.distance_matrix.csv} 
    ${ngsAMOVA} -in_dm ${line} -doEM 0 -doAMOVA 1 -isSim 1 -out ${outprefix} -doDist 1 -printMatrix 0 -sqDist 0 -m ${metadata} --hascolnames 1 &
    echo ${outprefix} >> "${2}"
    

done < "${1}"