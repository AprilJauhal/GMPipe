#!/bin/bash

## Please run GMfind_intronless_orfs.bat first to populate Pipe/in with ORF_list.fa and ORF_record.txt for your genome of interest
## Please edit the GMfind_userinput.ctl file to adjust parameters and file locations before running
## Please add the following files to Pipe/in/: master_seq.fa, outgroups.fa, query_seq.fa
##    query_seq.fa: fasta-formatted list of "ingroup" sequences
##    outgroups.fa: fasta-formatted list of known "outgroup" sequences
##    master_seq.fa: diverse subset of sequences from query_seqs.fa including a representatitive from each major subfamily in your ingroup

## To run:
## $ bash path_to/GMPipe/Scripts/GMPipe_finish.sh path_to/GMPipeline_userinput.ctl
## Meant to be run with Snakemake
#-----------------------------------------------------

source $1

timestamp=$(date +%F_%T)
echo "PREPARING FOR TREE TESTING " $timestamp
Rscript ${SCRIPT_PATH}/iqtree_parser.R ${PIPE_PATH}

timestamp=$(date +%F_%T)
echo "EXPORTING RESULTS " $timestamp
if [[ -d ${PIPE_PATH}/out ]]; then
  rm -r ${PIPE_PATH}/out
fi
mkdir ${PIPE_PATH}/out
mkdir ${PIPE_PATH}/out/domain

echo "FORMATING..."
Rscript ${SCRIPT_PATH}/final_formatter.R ${PIPE_PATH}
echo "GENERATING HISTOGRAM..."
Rscript ${SCRIPT_PATH}/histogram.R ${PIPE_PATH}

echo "COPYING FIES TO OUT FOLDER"
cp ${PIPE_PATH}/storage/ML_pass.fa ${PIPE_PATH}/out/PASSING_SEQUENCES.fa
cp ${PIPE_PATH}/storage/ML_fail.fa ${PIPE_PATH}/out/MLtree_fail_sequences.fa
cp ${PIPE_PATH}/storage/ML_undet.fa ${PIPE_PATH}/out/MLtree_undet_sequences.fa
cp ${PIPE_PATH}/storage/undet_by_tree.fa ${PIPE_PATH}/out/undet_by_tree.fa
cp ${PIPE_PATH}/storage/undet_by_constraint.fa ${PIPE_PATH}/out/undet_by_constraint.fa

timestamp=$(date +%F_%T)
echo "PIPELINE FINISHED" $timestamp
