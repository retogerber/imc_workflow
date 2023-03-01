#!/bin/bash

while getopts s:p:n:m: flag
do
    case "${flag}" in
        s) summary_panel_in=${OPTARG};;
        p) summary_panel_out=${OPTARG};;
        n) nuc_channels=( ${OPTARG} );;
        m) mem_channels=( ${OPTARG} );;
    esac
done
shift "$(( OPTIND - 1 ))"

# check if file exists
if [ ! -f "$summary_panel_in" ]; then
    echo "File '$summary_panel_in' does not exist" >&2
    exit 1
fi

if [ -z "$summary_panel_out" ]; then
    echo "Variable summary_panel_out (flag -p) is empty!" >&2
    exit 1
fi
if [ -z "$nuc_channels" ]; then
    echo "Variable nuc_channels (flag -n) is empty!" >&2
    exit 1
fi
if [ -z "$mem_channels" ]; then
    echo "Variable mem_channels (flag -m) is empty!" >&2
    exit 1
fi
cut -d, -f3,4 --complement "${summary_panel_in}" > ${summary_panel_out}


# label 'DNA' channels as '1' for deepcell (nucleus)
for nc in ${nuc_channels[@]}
do
    sed -r -i "s/($nc,)/\11/g" ${summary_panel_out}
done
# lable specific channels as '2' for deepcell (cytoplasma)
for mc in ${mem_channels[@]}
do
    sed -r -i "s/($mc,)/\12/g" ${summary_panel_out}
done

