# GMPipe: A gene mining pipeline for identifying olfactory recepter repertoires and other multicopy gene families from genomes

## What is GMPipe?

GMPipe is a gene mining pipeline designed to maximize the identification of multicopy gene family homologs in genome assemblies (for intronless genes) or lists of protein sequences. GMPipe has been validated for mining olfactory receptor genes from vertebrate genomes. This algorithm uses custom ingroup and outgroup sequence lists to distinguish between members of a gene family of interest from close off-target relatives. Overall, GMPipe uses HMMER as an initial screening step, paired with phylogenetic clustering to test whether passing sequences are more closly related to "ingroup" than "outgroup" sequences. The result is a pipeline with a low false negative rate, with minimal impact on the false positive rate, especially when paired with additional family-specific filtering steps.  

## Before running GMPipe:  
1) Using a job scheduler:   

    GMPipe includes a thorough (and resource-intensive) phylogenetic clustering step using ML trees and is thus designed to be used on a computer cluster with a job scheduler like OpenLava or Slurm. Any job scheduler will work as long as you have the ability to submit all relevant commands as a single command line input, as this is required for Snakemake to submit jobs to your cluster. If you haven't already, please familiarize yourself with how to submit jobs with your cluster's specific job scheduler before running GMPipe. 

2) Dependencies:  
    * Snakemake (version 5.1.4)
    * MUSCLE (version 3.8.1551)
    * HMMER (version 3.1b2)
    * IQTREE (version 1.6.11) 
    * R (version 4.0.0); packages: 
        * Biostrings
        * stringr
        * dplyr
        * ORFik
        * seqinr
        * ape
        * ggplot2
        * ggpubr

    While it is possible other versions may be compatible, GMPipe has only been tested and validated with these versions.

3) Setting up the input files:  

    GMPipe requires the following set of protien sequences:  
    * "query_seq.fa:": comprehensive set of sequences from your protein family of interest (eg. olfactory receptors) from either a closely related species or a variety of species
    * "outgroup.fa" smaller set of representative outgroup sequences representing the diversity of the closest homologs outside of the family of interest.  For ORs, this would be a variety of non-OR GPCR sequences.  
    * "master_seq.fa": diverse subset of sequences from query_seqs.fa including representatives from the major branches of a tree of the query_seq.fa sequences.  Ideally, should contain a similar number of sequences to outgroup.fa. 
    * Please remove the following special characters from filenames for all of the above files (as well as any spaces or tabs): \` \[ \] \{ \} \(\) \< \> \# \% \& \+ \\ \$ \^ \| \~ \* 
    * NOTE: the time that GMPipe takes to run is related to the number of sequences in the outgroup.fa and master_seq.fa files--a diverse list is better than a long list here.  Conversely, the length of query_seq.fa doesn't significantly affect the run time for GMPipe, so it is ok to submit a long list (although quality is still important as query_seq.fa is used for HMM-generation). While GMPipe does not accomoadate user-submitted HMM profiles, you can submit the sequences used to build a reference HMM profile from INTERPRO, PFAM, etc. as your query_seq.fa set.

## How to run:  
1) Setting up the initialization file and inputs:  

    For each run, you will need to create a directory that has the following:  
    * An "in" folder with the following files (see "Setting up the input files" section above):
        * query_seq.fa
        * outgroup.fa
        * master_seq.fa
    * GMPipeline_userinput.ctl (see below):
    
    Each GMPipe run requires its own "GMPipeline_userinput.ctl" control file.  The file can be located anywhere you want, but to keep things organized I reccomend keeping it in the same "PIPE_PATH." This file includes various settings and paths to dependencies for running GMPipe Scripts. 
    
    An example of the GMPipeline_userinput.ctl file is shown below and is also included in the repository.

    >#Set minimum ORF length (AA)  
    >ORF_SIZE=250 # Minimum ORF size-cutoff threshold; default for Olfactory Receptors (only applies to GMfind_intronless_orfs.sh script)  
    >#Full path to query genome  
    >QUERY_GENOME=[path_to/query_genome.fasta] #(only applies to GMfind_intronless_orfs.sh script)  
    >
    >#Set max threads (for GMPipe_start.sh), recommendation: between 1 and 20 threads  
    >THREADS=20  
    >
    >#Path to output folder (containing this file)  
    >PIPE_PATH=[path_to/run_name] # Custom directory for your GMPipe run, directory needs to exist before running script  
    >#Path to GMPipe/Scripts folder  
    >SCRIPT_PATH=[path_to]/GMPipe/Scripts  
    >
    >#Full paths to dependencies
    >MUSCLE=[path to/muscle-3.8.1551/bin/muscle]  
    >HMMER=[path to/hmmer-3.1b2-linux-intel-x86_64]  
    >IQTREE=[path to/iqtree-1.6.11/bin/iqtree]  
    >
    >#Load R module below, if necessary [Example: module load R/3.6.1]:
    >  

2) Preparing protein sequence files:  

    GMPipe can either screen a specified genome by identifing all of the longest ORF sequences above a particular threshold (will only identify intronless ORFS), or a list of provided protein sequences.   
    
    OPTION 1): Identifying intronless ORFs in a genome

    To run:  
    `$ bash [path_to]/GMPipe/Scripts/GMfind_intronless_orfs.sh [path_to]/GMPipeline_userinput.ctl`
    
    Note that this method will not identify: pseudogenes, fragmented genes, genes with introns
    
    OPTION 2) User-supplied sequence list
    
    Create a fasta (amino acid) file in the "in" folder with the name: ORF_list.fa
    
    NOTE: When identifying these sequences, keep in mind that methods such as tBLASTn that screen based on global pairwise comparisons may miss potential homologs and thus reduce the false negative rate for GMPipe.    

3) Running GMPipe_start.sh:  
    _This step screens sequences based on a HMM profile of the ingroups.fa sequences and optimizes the tree-building steps_  
    
    This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

    To run:  
    `$ bash [path_to]/GMPipe/Scripts/GMPipe_start.sh [path_to]/GMPipeline_userinput.ctl`

4) Running Snakemake:  
    _This step generates and parses trees for each query sequence with ingroup_reps.fa and outgroups.fa sequences_  
    
    IMPORTANT: this step must be called from the same directory as the "Snakefile" in the GMPipe directory. The "directory" below refers to the working directory that you selected in GMPipeline_userinput.ctl

    To run:  
    ```
    $ conda activate snakemake  
    $ snakemake --unlock --directory [directory]  
    $ snakemake --directory [directory] --jobs [max jobs Snakemake submits at once] --wait-for-files -w 500 --cluster [cluster commands here]
    ```
 
    Cluster commands just need to be simple single-threaded jobs, however, all elements required to submit a job to your cluster need to be included in this line. For ease of troubleshooting, please pick a different job name than you use for the snakemake command above.  
    Example of cluster commands for OpenLava:  
    'bsub -J snakemake -n 1 -o Pipe_snakemake.stdout -e Pipe_snakemake.stderr'
    
    The "--jobs" flag for Snakemake allows you to set the number of single threaded batch jobs that Snakemake submits to your cluster at a single time.  The larger the number, the faster it will run. If you have space availiable on your cluster. The more jobs you submit at a time, the faster this step will proceed.
 
    NOTE: This step often stalls when it is nearly complete (due to Snakemake thinking that all of the jobs have completed due to rounding errors). If the main snakemake job is still running but not submitting more jobs, it may have stalled. If this happens, cancel and resubmit the main Snakemake job (it will resubmit the required jobs without needing to redo previous steps). 

5) Running GMPipe_results.sh:  
    _This step sorts, compiles, and analyzes results_  
    
    This script is designed to be submitted to a job scheduler. Be sure to specify the same number of threads as speciied in the GMPipeline_userinput.ctl file. 
While optional for this particular step, GMPipe scripts are designed to run in the same directory as the Snakefile.

    To run:
    `$ bash [path_to]/GMPipe/Scripts/GMPipe_results.sh [path_to]/GMPipeline_userinput.ctl`
    
6) Manual checks:  

    In the out folder, find the names of the sequences flagged as "undetermined" from the following files:
    * undet_by_constraint.fa (sequences flagged based on tree statistics)
    * MLtree_undet_sequences.fa (sequences flagged based on ambiguous clustering)
    For each sequence, go to its corresponding folder in storage/ML and examine the unconstr_##.iqtree files to visually check whether the sequence clusutered with ingroups or outgroups, or perform other checks if desired.
    
    Also examine the "bitscore_hist.png" file.  This image contains two charts that show the distributions of HMM scores. The first shows the distribution of scores for representative ingroups, outgroups, and HMM-passing query sequences that either passed, failed, or were flagged as undetermined for the tree tests.  The tree-passing sequences are hidden for the second grapth to allow for better visualization of the other categories. You should expect to see:
        1) Clear separation of ingroup and outgroup sequences (GMPipe will fail if the distributions are overlapping, but if the distributions are nearly overlapping this may affect the results.)
        2) Clear separation of tree-passing and tree-failing sequences (If you consistently observe sequences with high bit scores that failed the tree test, this could be a sign that the ingroup HMM profile is not specific enough to your family of interest)
        3) Clear separation of outgroup and tree-passing sequences. (The default bit-score cut-off is the value of the highest-scoring outgroup sequence. If you would like to select a more stringent cut-off after viewing the histogram results you can filter the FINAL_RESULTS.txt file based on the bit-score).

7) Reccomended additional steps:  

    GMPipe is designed to cast as wide of a net as possible, while providing specificity through phylogenetic comparisions with representative ingroup and outgroup sequences. This is based on the philosophy that it is easier to add additional filtering steps to narrow a screen, but difficult to make a screeen more broad after it has been run.  
    
    A simple strategy is to align your sequences to a known reference sequence, and check whether there are any ambiguous residues or gaps in the most conserved motifs in the protein.  An advantage of this strategy is that it requires no a priori knowledge of which mutations will or won't inactivate the protein. For Olfactory Receptors this may include the following motifs for example:
    * GN
    * MAYDRYVAIC
    * KAFTCASH
    * PMLNPFIY
    For Olfactory Receptors, including both the GN and PMLNPFIY motifs can also be used to define a "full-length" sequence, as these motifs are found on TM1 and TM7, respectively.  

Additional notes:  
Theoretically the three main scripts (how to run steps 3-5) can be run together through a master script, although the logistics of this may depend on your job scheduler. However, in my experience running the jobs as separate commands makes it easier to figure out if/where a step has stalled. 

Troubleshooting:
    * If the snakemake job is still running but not generating additional jobs, end and restart the job
    * GMPipe is designed to stop prematurely if ingroup and outgroup sequences are too similar. If this happens, check the log file for errors such as "ERROR: INGROUP AND OUTGROUP HMMER BIT SCORES TOO SIMILAR", if this occurs, check and redesign outgroup.fa, ingroup.fa and master_seq.fa sequences. Make sure outgroup sequences are phyllogenetically distinct from ingroup sequences. The master_seq.fa sequences should represent a subset of ingroup sequences.  
    * Ensure that the "snakemake" file is in the same folder that the commands are being called from (which is not necessarily the same directory where the analysis is done). This file can be duplicated and placed into individual directories for each run but should not be edited. 
    * Make sure that you have removed special characters from the sequence names (if you are still having trouble you can try removing stop codons from the sequences themselves as well): \` \[ \] \{ \} \(\) \< \> \# \% \& \+ \\ \$ \^ \| \~ \* 
    * If you are encountering issues running GMPipe, I reccomend trying the tutorial exercise to test if there are any issues with dependencies, etc.
