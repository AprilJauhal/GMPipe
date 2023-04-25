# To run: Rscript path_to/iqtree_writer1.R PIPE_PATH IQTREE_path
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose is to generate IQTREE commands for optimization trees built with just master_seq and 
# outgroup sequences to determine optimal settings and ensure trees look correct 

library(stringr)
library(Biostrings)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
IQ_PATH <- commandArgs(trailingOnly=TRUE)[2]
master_afa <- paste0(PIPE_PATH, "/storage/reference.afa")
in_list <- paste0(PIPE_PATH, "/in/master_seq.fa")
out_list <- paste0(PIPE_PATH, "/in/outgroups.fa")
MLpath <- paste0(PIPE_PATH, "/storage/opt_ML/")

# Setting constraints for phylogenetically constrained trees
in_names <- str_flatten(str_remove_all(names(readAAStringSet(in_list)), " "), collapse=",")
out_names <- str_flatten(str_remove_all(names(readAAStringSet(out_list)), " "), collapse=",")
main_constraint <- paste0("((", in_names, "),", out_names, ");")
write(main_constraint, file=paste0(MLpath, "main.constr"))

iq_list1 <- vector()
# writing IQTREE commands tor unconstrained trees [with model_finder]
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr1")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr2")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr3")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr4")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr5")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr6")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr7")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr8")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr9")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_unconstr10")

# writing IQTREE commands tor constrained trees: ((ingroup), outgroup) [with model_finder]
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr1")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr2")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr3")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr4")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr5")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr6")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr7")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr8")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr9")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PIPE_PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PIPE_PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PIPE_PATH, "/storage/opt_ML/main_constr10")

# Writing IQTREE commands
write(iq_list1, file=paste0(MLpath, "iq_list1.cmd"))
