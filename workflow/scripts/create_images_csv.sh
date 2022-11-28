#!/bin/bash



images_csv_header="image,width_px,height_px,num_channels"
echo $images_csv_header > $1
# add images names and dimensions
combined_samples=( "$@" )
unset 'combined_samples[0]'
for sample in ${combined_samples[@]}; do
	basename $sample, >> $1
    #temp_sample=$( echo $sample | sed 's|_seg1234||g' )
    temp_sample=$( echo $sample | sed 's|_mespp_id[0-9][0-9]||g' )
	dims=$( identify -format '%w,%h,%n\n' $temp_sample | head -n 1 )
	truncate -s-1 $1 
	echo  $dims >> $1 
done

#for sample in ${combined_samples[@]}; do
#	basename ${sample}_seg12, >> $1
#	dims=$( identify -format '%w,%h,%n\n' $sample | head -n 1 )
#	truncate -s-1 $1 
#	echo  $dims >> $1 
#done
