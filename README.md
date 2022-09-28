What is GMPipe?

Before running GMPipe:
1) Using a job scheduler: 
GMPipe includes a thorough (and resource-intensive) phylogenetic clustering step using ML trees and is thus designed to be used on a computer cluster with a job scheduler like OpenLava or Slurm. Any job scheduler will work as long as you have the ability to submit all relevant commands as a single command line input, as this is required for Snakemake to submit jobs to your cluster. If you haven't already, please familiarize yourself with how to submit jobs with your cluster's specific job scheduler before running GMPipe. 

2) Setting up the input files: 
GMPipe requires the following set of protien sequences:
  1) "query_seq.fa:": comprehensive set of sequences from your protein family of interest (eg. olfactory receptors) from either a closely related species or a variety of species
  2) "outgroup.fa" smaller set of representative outgroup sequences representing the diversity of the closest homologs outside of the family of interest.  For ORs, this would be a variety of non-OR GPCR sequences.  
  3) master_seq.fa: diverse subset of sequences from query_seqs.fa including representatives from the major branches of a tree of the query_seq.fa sequences.  Ideally, should contain a similar number of sequences to outgroup.fa. 
  NOTE: the time that GMPipe takes to run is related to the number of sequences in the outgroup.fa and master_seq.fa files--a diverse list is better than a long list here.  Conversely, the length of query_seq.fa doesn't significanntly affect the run time for GMPipe, so it is ok to submit a long list (although quality is still important as query_seq.fa is used for HMM-generation). While GMPipe does not accomoadate user-submitted HMM profiles, you can submit the sequences used to build a reference HMM profile from INTERPRO, PFAM, etc. as your query_seq.fa set.

How to run:
1) Setting up the initialization file and inputs
2) Preparing protein sequence files
3) Running GMPipe_start.sh 
This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

To run:
$ bash <path_to>/GMPipe/Scripts/GMPipe_start.sh <path_to>/GMPipeline_userinput.ctl

4) Running Snakemake

5) Running GMPipe_results.sh
This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

To run:

$ bash <path_to>/GMPipe/Scripts/GMPipe_results.sh <path_to>/GMPipeline_userinput.ctl

6) Manual checks
7) Reccomended additional steps

Notes
Theoretically the three main scripts (how to run steps 3-5) can be run together through a master script, although the logistics of this may depend on your job scheduler. 
