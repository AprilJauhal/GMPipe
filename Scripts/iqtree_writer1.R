# To run: Rscript path_to/iqtree_writer1.R PIPE_PATH IQTREE_path
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own

library(stringr)
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
IQ_PATH <- commandArgs(trailingOnly=TRUE)[2]
master_afa <- paste0(PIPE_PATH, "/storage/reference.afa")
in_list <- paste0(PIPE_PATH, "/in/ingroup_reps.fa")
out_list <- paste0(PIPE_PATH, "/in/outgroups.fa")
MLpath <- paste0(PIPE_PATH, "/storage/opt_ML/")

in_names <- str_remove_all(grep(">", readLines(in_list), value=TRUE), ">")
in_names <- str_flatten(str_remove_all(in_names, " "), collapse=",")
out_names <- str_remove_all(grep(">", readLines(out_list), value=TRUE), ">")
out_names <- str_flatten(str_remove_all(out_names, " "), collapse=",")
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

write(iq_list1, file=paste0(MLpath, "iq_list1.cmd"))

