#!/usr/bin/env bash

references=( $1 )
samples=( $2 )
sample=${samples[0]}
nuc_channels=( $3 )
mem_channels=( $4 )

#echo ${samples[@]}
#echo ${references[@]}

references_base=()
for ref in ${references[@]}
do
    references_base+=("${ref%_summary.csv}")
done

#echo ${references_base[@]}

#echo ${samples[0]}


matching_ref=()
for ref in ${references_base[@]}
do
    tmp=$( echo "$sample" | grep -z -o "$ref" )
    if [[ $tmp != "" ]]; then
        matching_ref+=("${tmp}_summary.csv")
    fi
done
echo ${matching_ref[@]}


cut -d, -f3,4 --complement "${matching_ref[0]}" > ${sample}


# label 'DNA' channels as '1' for deepcell (nucleus)
for nc in ${nuc_channels[@]}
do
    sed -r -i "s/($nc,)/\11/g" ${sample}
done
# lable specific channels as '2' for deepcell (cytoplasma)
for mc in ${mem_channels[@]}
do
    sed -r -i "s/($mc,)/\12/g" ${sample}
done

