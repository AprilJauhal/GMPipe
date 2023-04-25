# GMPipe Tutorial

## Overview
This is a basic tutorial for GMPipe. This can be used to preview the program and check if it meets your needs, as well as for troubleshooting. In this exercise, we will extract olfactory receptor sequences from the zebrafish (GRCz10) refseq genome from ncbi, using a set of reedfish sequences as an input (originally from: Policarpo, M., Bemis, K. E., Laurenti, P., Legendre, L., Sandoz, J. C., Rétaux, S., & Casane, D. (2022). Coevolution of the olfactory organ and its receptor repertoire in ray-finned fishes. BMC biology, 20(1), 195. https://doi.org/10.1186/s12915-022-01397-x). For this tutorial, the ORF_list.fa is provided, although it is possible to generate your own set using the GMfind_intronless_orfs.sh scripts and by downloading and unzipping the genome (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.5_GRCz10/GCF_000002035.5_GRCz10_genomic.fna.gz) and then specifying its path in the control file GMPipeline_userinput.ctl. This tutorial will assume that you already have input files prepared. Please refer to the README file for more information on dependencies and file setup.

## Step 1) Setting the control file
Each run should have its own GMPipeline_userinput.ctl file. Please keep the filename the same, but edit the file contents to fit your parameters. Here is an example file below:
```
    #Set minimum ORF length (AA)  
    ORF_SIZE=250 # Minimum ORF size-cutoff threshold; default for Olfactory Receptors (only applies to GMfind_intronless_orfs.sh script)  
    #Full path to query genome  
    QUERY_GENOME=[path_to/query_genome.fasta] #(only applies to GMfind_intronless_orfs.sh script)  
    
    #Set max threads (for GMPipe_start.sh), recommendation: between 1 and 20 threads  
    THREADS=20  
    
    #Path to output folder (containing this file)  
    PIPE_PATH=[path_to/run_name] # Custom directory for your GMPipe run, directory needs to exist before running script  
    #Path to GMPipe/Scripts folder  
    SCRIPT_PATH=[path_to]/GMPipe/Scripts  
    
    #Full paths to dependencies
    MUSCLE=[path to/muscle-3.8.1551/bin/muscle]  
    HMMER=[path to/hmmer-3.1b2-linux-intel-x86_64]  
    IQTREE=[path to/iqtree-1.6.11/bin/iqtree]  
    
    #Load R module below, if necessary [Example: module load R/3.6.1]:
```     
The first two parameters are specific to the GMfind_intronless_orfs.sh script. The ORF_SIZE parameter refers to the minimum ORF size cutoff (AA) for mining ORFs from the genome. The reccomended length for olfactory receptors is 250bp but this can be adjusted if desired. The QUERY_GENOME path is the path to the genome used for ORF extraction.  

Next, we can set the maximum number of threads to be used for GMPipe_start.sh. 20 threads is a good starting point if availiable; the threads make a bigger difference for the snakemake pipeline and are set separately. We also set the paths to the directory where the pipeine will be run (should contain the "in" folder), as well as the path to find the GMPipe scripts (the snakemake is kept separately--in the same directory where the run is executed from). 

Here we can also set the full paths to the dependencies, and if necessary, you can load an R module here. 

## Step 2) GMPipe_Start.sh
This is the initialization script for GMPipe. It generates a HMMER profile for ingroup sequences, generates a bit-score cutoff value based on outgroup sequence scores, and determines the IQTREE substitution model for use in the main pipeline. We are going to assume you are running the tutorial example. The command to run the pipeline is written below and assumes you are executing it from the "GMPipe" directory. However, it is reccomended to execute this command as an HPC job, and some examples are also listed below:

To run: 
```
$ bash Scripts/GMPipe_start.sh tutorial/GMPipeline_userinput.ctl
```
To run with OpenLava (assuming 20 threads set in control file): 
```
$ bsub -J GM_start -n 20 -o tutorial/logs/start.stdout -e tutorial/logs/start.stderr' \
$ "bash Scripts/GMPipe_start.sh tutorial/GMPipeline_userinput.ctl"
```
To run with slurm (assuming 2 days wall-time allotment, but the script should take less than a few hours):
```
$ sbatch --job-name=GMstart-test --cpus-per-task=20 --nodes=1 --mem=1000 --time=02-00:00:00 --output="tutorial/logs/GMstart.%j.out" \
$ --wrap="bash Scripts/GMPipe_start.sh tutorial/GMPipeline_userinput.ctl"
```

After running, there should be a "storage" folder, which should have the following contents:

* HMMER_FAIL.fa   
* HMMER_PASS.fa   
* master_seq.afa  
* outgroups.afa
* reference.afa
* query_seq.hmm
* query_seq.afa
* opt_ML/
* ML/
* hmm/ (if troubleshooting, make sure this contains a "bit_cutoff.txt" file)
    
If any of these sequences are missing, check the log. 

## STEP 3) snakemake
This is the part where the snakemake pipeline is called. In order to work, the snakemake command must be called from the same directory that contains the snakefile. The snakefile can be moved or copied but should not be edited or renamed. The snakemake command creates a snakemake job that calls other jobs from your HPC job scheduler, and some examples are shown below. 

This snakemake pipeline executes as a single-threaded command that calls additional single-threaded commands. The maximum number of cammands is set by the "-w" flag, which is set to 500 in the examples. More commands will allow the job to complete faster, and fewer threads are less intensive on the cluster but take longer. Our tutorial example should run relatively fast. At 200 jobs it should only take a few hours and realistically the number of simultaneous jobs will depend on how busy the clustuer is. Each job is fairly short but thousands are generated by snakemake, which prevents the pipeline from hogging too many computational resources vs. longer individual jobs that use a large number of threads.

Note that when selecting wall time for a slurm command, both the initialization script and snakemake pipelines are restartable if interrupted. Conservative wall times are more convenient but not strictly necessary.

To prepare to run (assuming snakemake is installed as a module):  
```
$ module load conda
$ conda activate snakemake
$ conda activate snakemake  
$ snakemake --unlock --directory [directory]  
```
To run (base command): 
```
$ snakemake --directory [directory] --jobs [max jobs Snakemake submits at once] --wait-for-files -w 500 --cluster [cluster commands here]
```
To run with openlava:
```
$ bsub -J GM_snake -n 1 -o tutorial/logs/GM_snake.stdout -e tutorial/logs/GM_snake.stderr' \
$ "snakemake --directory [directory] --jobs [max jobs Snakemake submits at once] --wait-for-files -w 500 \
$ --cluster 'bsub -J snakemake -n 1 -o Pipe_snakemake.stdout -e Pipe_snakemake.stderr'"
```
To run with slurm (run time will vary based on number of threads, but 10 days is a highly conservative estimate):
```
$ sbatch --job-name=GM_snake --cpus-per-task=1 --nodes=1 --mem=1000 --time=10-00:00:00 --output="tutorial/logs/GMsnake.%j.out" \
$ --wrap="snakemake --directory [directory] --jobs [max jobs Snakemake submits at once] --wait-for-files -w 500 \
$ --cluster "sbatch -J zebra-snake_11 -c 1 -N 1 --mem 1000 --time 1-00:00:00 --output \ 
$ '/powerplant/workspace/hraaxe/GMPipe_2022/test_zebrafish_GRCz11/logs/snake_subjob.%j.out'"
```
Once this part has finished running, there should be an out folder that contains a variety of files, including "PASSING_SEQUENCES.fa" and assuming that you used the example files or extracted ORFs using the GMfind_intronless_orfs.sh and that all the R libraries loaded properly, there should also be file called "bitscore_hist.png" that contains a histogram which can be used to check if the cutoff score used was appropriate. 

If these files do not appear, check if the snakemake job is still running. If it is the only job running it may have stalled. Cancel the job and restart. If this still does not work, check the logs for the main snakemake job for clues. If the pipeline has successfully completed, it should say that it has completed job "0".

STEP 4) Data examination and follow-up

In the output folder there should be the following files: 

* PASSING_SEQUENCES.fa = Sequences that passed both tests with no issues (but should still be aligned and manually inspected in another program)
* MLtree_fail_sequences.fa = Sequences that failed due to clustering with outgroup sequences instead of ingroup sequences  
* MLtree_undet_sequences.fa  = All sequences that were flagged for manual inspection  
* undet_by_constraint.fa = Sequences flagged for inspection specifically based on not consistently grouping with outgroups 
* undet_by_tree.fa = Sequences flagged for inspection specifically based on constrained trees not consistently passing tree tests (pAU)
* outgroup_bitscores.txt = set of all bit_scores from outgroup sequences
* ingroup_bitscores.txt = set of all bit_scores from master_seq sequences
* passing_hit.tbl = hmmer table report of passing hits
* domain/ = directory containing hmmer domain out tables for hits, outgroups, and master_seq sequences (for those interested in hmmer details)

These files may also appear if an ORF record was included in the input folder
    * FINAL_REPORT.txt = Report of all hmmer-passing ORFs. Shows their nucleotide and amino acid sequences, scaffold accession number and position, bit-score and outcome of hmmer and tree tests. 
    * bitscore_hist.png = Histogram file displaying hmmer bit scores for ingroup vs. outgroup sequences overlaid with hmmer_hits sorted based on their outcomes in the tree test.   

If most of the files mentioned above are present (some may be absent if they were empty but check to ensure PASSING_SEQUENCES.fa is present), the next step is to examine some of the output files. 

I highly reccomend examining the unconstrained trees for the undetermined sequences. (Note that these sequences are not included in PASSING_SEQUENCES.fa and if they pass manual inspection you will need to generate a combined list yourself). Essentially, you want to make sure that the hits clustered with ingroups as opposed to outgroups. The constrained trees are easier to automatically parse but may not represent the highest scoring solution, so it is a good idea to manually inspect the unconstrained trees that are flagged. If hits cluster with ingroups instead of outroups for at least 9 out of 10 trees, you can count them as passing sequences. You can identify the paths to the unconstrained trees with the following code:
```
    $ undet_hits=($( grep ">" test_zebrafish/out/MLtree_undet_sequences.fa | sed 's/>//g' ))
    $ for undet_hit in "${undet_hits[@]}"
    $ do
        $ ls test_zebrafish/storage/ML/${undet_hit}/${undet_hit}_unconstr**.iqtree
    $ done
```

Ideally, you should also check the histogram file to see if there was appropriate separation of ingroup and outgroup sequences. The graph on the top has all hmmer-hit scores as well as the outgroup and master_seq (ingroup) scores, and the graph on the bottom is the same but lacks the tree-passing scores for ease of clarity. The most important thing to check for is that there is a large gap between the ingroup and outgroup sequences. In the tutorial example, there is a large gap. However, if there wasn't it would be a sign that the GMPipe input sequences might need adjusting in order to minimize false positives and negatives. If this is the case, make sure that the outgroup sequnces are phyllogenetically distinct from the ingroup sequences. This could happen if the ingroup is technically a subgroup fo teh outgroups and lacks
distinct features or motifs, in which case, you might need to rethink whether separation is possible without setting a high bits-core cut-off. You may also want to consider retroactively raising the bit-score cut-off if there is poor separation between outgroupand tree-passing sequences, depending on the goals of your analysis. The tutorial example is a borderline case. If the goal is obtaining a stringent list of homologs, a stronger bit-score cut-off should be considered. However, if it is possible that your ingroup list of sequences may not represent the full diversity of the family, then there could also be a case to keep these sequences, especially if they pass additional tests following the pipeline. In the case of GMPipe benchmarking, I retained these sequences to get a more accurate estimation of the false positive rate. 


