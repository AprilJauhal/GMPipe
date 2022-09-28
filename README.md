What is GMPipe?

Before running GMPipe:
1) Using a job scheduler: 
 -GMPipe includes a thorough (and resource-intensive) phylogenetic clustering step using ML trees and is thus designed to be used on a computer cluster with a job scheduler like OpenLava or Slurm. Any job scheduler will work as long as you have the ability to submit all relevant commands as a single command line input, as this is required for Snakemake to submit jobs to your cluster. If you haven't already, please familiarize yourself with how to submit jobs with your cluster's specific job scheduler before running GMPipe. 

2) Dependencies
Snakemake
MUSCLE (version 3.8.1551)
HMMER= (version 3.1b2)
IQTREE (version 1.6.11)

While it is possible other versions may be compatible, GMPipe has only been tested and validated with these versions.

3) Setting up the input files: 
GMPipe requires the following set of protien sequences:
-"query_seq.fa:": comprehensive set of sequences from your protein family of interest (eg. olfactory receptors) from either a closely related species or a variety of species
-"outgroup.fa" smaller set of representative outgroup sequences representing the diversity of the closest homologs outside of the family of interest.  For ORs, this would be a variety of non-OR GPCR sequences.  
-master_seq.fa: diverse subset of sequences from query_seqs.fa including representatives from the major branches of a tree of the query_seq.fa sequences.  Ideally, should contain a similar number of sequences to outgroup.fa. 
-NOTE: the time that GMPipe takes to run is related to the number of sequences in the outgroup.fa and master_seq.fa files--a diverse list is better than a long list here.  Conversely, the length of query_seq.fa doesn't significanntly affect the run time for GMPipe, so it is ok to submit a long list (although quality is still important as query_seq.fa is used for HMM-generation). While GMPipe does not accomoadate user-submitted HMM profiles, you can submit the sequences used to build a reference HMM profile from INTERPRO, PFAM, etc. as your query_seq.fa set.

How to run:
1) Setting up the initialization file and inputs
2) Preparing protein sequence files
3) Running GMPipe_start.sh 
This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

To run:
$ bash <path_to>/GMPipe/Scripts/GMPipe_start.sh <path_to>/GMPipeline_userinput.ctl

4) Running Snakemake
IMPORTANT: this step must be called from the same directory as the "Snakefile" in the GMPipe directory. The "directory" below refers to the working directory that you selected in GMPipeline_userinput.ctl

To run: 
$ conda activate snakemake
$ snakemake --unlock --directory <directory>
$ snakemake --directory <directory> --jobs <max jobs Snakemake submits at once> --wait-for-files -w 500 --cluster <cluster commands here>
 
 Cluster commands just need to be simple single-threaded jobs, however, all elements required to submit a job to your cluster need to be included in this line. For ease of troubleshooting, please pick a different job name than you use for the snakemake command above.
 Example of cluster commands for OpenLava:
 'bsub -J snakemake -n 1 -o Pipe_snakemake.stdout -e Pipe_snakemake.stderr'"
 
 NOTE: This step often stalls when it is nearly complete (due to Snakemake thinking that all of the jobs have completed due to rounding errors). If the main snakemake job is still running but not submitting more jobs, it may have stalled. If this happens, cancel and resubmit the main Snakemake job (it will resubmitd the required jobs without needing to redo previous steps). 

5) Running GMPipe_results.sh
This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

To run:

$ bash <path_to>/GMPipe/Scripts/GMPipe_results.sh <path_to>/GMPipeline_userinput.ctl

6) Manual checks
7) Reccomended additional steps

Additional notes:
Theoretically the three main scripts (how to run steps 3-5) can be run together through a master script, although the logistics of this may depend on your job scheduler. However, in my experience running the jobs as separate commands makes it easier to figure out if/where a step has stalled.
