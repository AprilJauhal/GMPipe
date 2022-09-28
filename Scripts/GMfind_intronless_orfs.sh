#!/bin/bash

## Please edit the GMfind_userinput.ctl file to adjust parameters and file locations before running
## To run:
## $ bash path_to/GMPipe/Scripts/GMfind_intronless_orfs.sh path_to/GMfind_userinput.ctl
#-----------------------------------------------------

source $1

echo "CHECKING DIRECTORIES"
timestamp=$(date +%F_%T)
echo $timestamp

if [[ ! -d "${PIPE_PATH}/in/" ]]; then
    mkdir "${PIPE_PATH}/in/"
fi

echo "FINDING ORFs IN GENOME"
timestamp=$(date +%F_%T)
echo $timestamp
Rscript ${SCRIPT_PATH}/ORF_finder.R ${PIPE_PATH} ${ORF_SIZE} ${QUERY_GENOME}
