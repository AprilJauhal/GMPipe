#!/bin/bash

## Please run GMfind_intronless_orfs.bat first to populate Pipe/in with ORF_list.fa and ORF_record.txt for your genome of interest
## Please edit the GMfind_userinput.ctl file to adjust parameters and file locations before running
## Please add the following files to Pipe/in/: master_seq.fa, outgroups.fa, query_seq.fa
##    query_seq.fa: fasta-formatted list of "ingroup" sequences
##    outgroups.fa: fasta-formatted list of known "outgroup" sequences
##    master_seq.fa: diverse subset of sequences from query_seqs.fa including a representatitive from each major subfamily in your ingroup

## To run:
## $ bash path_to/GMPipe/Scripts/GMPipe_start.bat path_to/GMPipeline_userinput.ctl

#-----------------------------------------------------

source $1

echo "INITIALIZING GENE MINING PIPELINE"
timestamp=$(date +%F_%T)
echo $timestamp

echo "CHECKING PROGRESS"
timestamp=$(date +%F_%T)
echo $timestamp

if [[ -f ${PIPE_PATH}/storage/hmm/bit_cutoff.txt ]]; then
  echo "CUTOFF FOUND, CONTINUING FROM LAST STOP-POINT"
  source ${PIPE_PATH}/storage/hmm/bit_cutoff.txt
  echo "CUTOFF IS " $CUTOFF
else
  echo "CREATING DIRECTORIES"
  if [[ -d ${PIPE_PATH}/storage/hmm/ ]]; then
    rm -r ${PIPE_PATH}/storage/hmm/pass
    rm -r ${PIPE_PATH}/storage/hmm/
    rm -r ${PIPE_PATH}/storage
  fi
  mkdir ${PIPE_PATH}/storage
  mkdir ${PIPE_PATH}/storage/hmm
  mkdir ${PIPE_PATH}/storage/hmm/pass
  
  echo "ALIGNING SEQUENCES"
  timestamp=$(date +%F_%T)
  echo $timestamp
  ${MUSCLE} -in ${PIPE_PATH}/in/query_seq.fa -out ${PIPE_PATH}/storage/query_seq.afa -quiet
  ${MUSCLE} -in ${PIPE_PATH}/in/master_seq.fa -out ${PIPE_PATH}/storage/master_seq.afa -quiet
  ${MUSCLE} -in ${PIPE_PATH}/in/outgroups.fa -out ${PIPE_PATH}/storage/outgroups.afa -quiet
  ${MUSCLE} -profile -in1 ${PIPE_PATH}/storage/outgroups.afa -in2 ${PIPE_PATH}/storage/master_seq.afa -out ${PIPE_PATH}/storage/reference.afa -quiet
  
  echo "GENERATING HMMER PROFILE"
  timestamp=$(date +%F_%T)
  echo $timestamp
  HMM_THREADS=$(expr $THREADS - 1)
  ${HMMER}/bin/hmmbuild --amino --cpu $HMM_THREADS --informat afa ${PIPE_PATH}/storage/query_seq.hmm ${PIPE_PATH}/storage/query_seq.afa
  
  echo "RUNNING HMMER FOR INGROUPS"
  timestamp=$(date +%F_%T)
  echo $timestamp
  ${HMMER}/bin/hmmsearch --cpu $HMM_THREADS -o in_hmm.out --tblout ${PIPE_PATH}/storage/hmm/ingroup.tblout --domtblout ${PIPE_PATH}/storage/hmm/ingroup.domtblout --noali --notextw --tformat fasta ${PIPE_PATH}/storage/query_seq.hmm ${PIPE_PATH}/in/master_seq.fa
  
  echo "RUNNING HMMER FOR OUTGROUPS"
  timestamp=$(date +%F_%T)
  echo $timestamp
  ${HMMER}/bin/hmmsearch --cpu $HMM_THREADS -o out_hmm.out --tblout ${PIPE_PATH}/storage/hmm/outgroup.tblout --domtblout ${PIPE_PATH}/storage/hmm/outgroup.domtblout --noali --notextw --tformat fasta ${PIPE_PATH}/storage/query_seq.hmm ${PIPE_PATH}/in/outgroups.fa
    
  echo "CALCULATING CUTOFF"
  timestamp=$(date +%F_%T)
  echo $timestamp 
  Rscript ${SCRIPT_PATH}/hmm_cutoff.R ${PIPE_PATH}
  
  echo "RETRIEVING CUTOFF"
  source ${PIPE_PATH}/storage/hmm/bit_cutoff.txt
  echo "CUTOFF IS " $CUTOFF

fi
  
if [[ -f ${PIPE_PATH}/storage/HMMER_PASS.fa ]]; then
  echo "HMMER RESULTS FOUND, CONTINUING FROM LAST STOP POINT"
else 
  echo "CREATING DIRECTORIES"
  if [[ -d ${PIPE_PATH}/storage/opt_ML ]]; then
    rm -r ${PIPE_PATH}/storage/opt_ML
    rm -r ${PIPE_PATH}/storage/ML
  fi
  
  mkdir ${PIPE_PATH}/storage/opt_ML
  mkdir ${PIPE_PATH}/storage/opt_ML/concat
  mkdir ${PIPE_PATH}/storage/ML
  mkdir ${PIPE_PATH}/storage/ML/concat

  
  echo "RUNNING HMMER ON PUTATIVE ORFS"
  timestamp=$(date +%F_%T)
  echo $timestamp
  ${HMMER}/bin/hmmsearch --cpu $HMM_THREADS -o hit_hmm.out --tblout ${PIPE_PATH}/storage/hmm/hit.tblout --domtblout ${PIPE_PATH}/storage/hmm/hit.domtblout --noali --notextw --tformat fasta ${PIPE_PATH}/storage/query_seq.hmm ${PIPE_PATH}/in/ORF_list.fa
  
  echo "PROCESSING HMMER RESULTS"
  timestamp=$(date +%F_%T)
  echo $timestamp
  Rscript ${SCRIPT_PATH}/hmm_scorer.R ${PIPE_PATH}
  
fi

if test -f "${PIPE_PATH}/storage/opt_ML/best_model.txt"; then
	echo "IQTREE OPTIMIZATION COMPLETE, CONTINUING FROM LAST STOP POINT"
else

  timestamp=$(date +%F_%T)
  echo "PREPARING FOR TREE OPTIMIZATION " $timestamp
  Rscript ${SCRIPT_PATH}/iqtree_writer1.R ${PIPE_PATH} ${IQTREE}
  	
  timestamp=$(date +%F_%T)
  echo "RUNNING IQTREE OPTIMIZATION  " $timestamp
  parallel -j${THREADS} < ${PIPE_PATH}/storage/opt_ML/iq_list1.cmd
  	
  timestamp=$(date +%F_%T)
  echo "CONCATENATING FILES  " $timestamp
  	
  cat ${PIPE_PATH}/storage/opt_ML/main_unconstr1.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr2.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr3.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr4.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr5.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr6.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr7.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr8.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr9.treefile ${PIPE_PATH}/storage/opt_ML/main_unconstr10.treefile ${PIPE_PATH}/storage/opt_ML/main_constr1.treefile ${PIPE_PATH}/storage/opt_ML/main_constr2.treefile ${PIPE_PATH}/storage/opt_ML/main_constr3.treefile  ${PIPE_PATH}/storage/opt_ML/main_constr4.treefile ${PIPE_PATH}/storage/opt_ML/main_constr5.treefile ${PIPE_PATH}/storage/opt_ML/main_constr6.treefile ${PIPE_PATH}/storage/opt_ML/main_constr7.treefile ${PIPE_PATH}/storage/opt_ML/main_constr8.treefile ${PIPE_PATH}/storage/opt_ML/main_constr9.treefile ${PIPE_PATH}/storage/opt_ML/main_constr10.treefile > ${PIPE_PATH}/storage/opt_ML/concat/main.treelst
  	
  timestamp=$(date +%F_%T)
  echo "COMPARING RESULTS  " $timestamp
  ${IQTREE} -s ${PIPE_PATH}/storage/reference.afa -m MFP -msub nuclear -redo -z ${PIPE_PATH}/storage/opt_ML/concat/main.treelst -n 0 -zb 10000 -au -pre ${PIPE_PATH}/storage/opt_ML/concat/main_comparison

  timestamp=$(date +%F_%T)
  echo "PREPARING FOR TREE TESTING " $timestamp
  Rscript ${SCRIPT_PATH}/iqtree_writer2.R ${PIPE_PATH}

fi;
