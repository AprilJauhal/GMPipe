#Set minimum ORF length (AA)
ORF_SIZE=250 # Minimum ORF size-cutoff threshold; default for Olfactory Receptors only applies to GMfind_intronless_orfs.sh script
#Full path to query genome
QUERY_GENOME=<path_to/query_genome.fasta> #only applies to GMfind_intronless_orfs.sh script

#Set max threads (for GMPipe_start.sh), recommendation: between 1 and 20 threads
THREADS=20

#Path to output folder (containing this file)
PIPE_PATH=[path_to/run_name] # Custom directory for your GMPipe run, directory needs to exist before running script
#Path to GMPipe/Scripts folder
SCRIPT_PATH=<path_to>/GMPipe/Scripts

#Full paths to dependencies (examples shown for reccomended versions)
MUSCLE=<path_to/muscle-3.8.1551/bin/muscle>
HMMER=<path_to/hmmer-3.1b2-linux-intel-x86_64>
IQTREE=<path_to/iqtree-1.6.11/bin/iqtree>

#Load R module below, if necessary [Example: module load R/3.6.1]:
