library(stringr)
PATH <- commandArgs(trailingOnly=TRUE)[1]
IQ_PATH <- commandArgs(trailingOnly=TRUE)[2]
master_afa <- paste0(PATH, "/storage/reference.afa")
in_list <- paste0(PATH, "/in/master_seq.fa")
out_list <- paste0(PATH, "/in/outgroups.fa")
MLpath <- paste0(PATH, "/storage/opt_ML/")

in_names <- str_remove_all(grep(">", readLines(in_list), value=TRUE), ">")
in_names <- str_flatten(str_remove_all(in_names, " "), collapse=",")
out_names <- str_remove_all(grep(">", readLines(out_list), value=TRUE), ">")
out_names <- str_flatten(str_remove_all(out_names, " "), collapse=",")
main_constraint <- paste0("((", in_names, "),", out_names, ");")

write(main_constraint, file=paste0(MLpath, "main.constr"))

iq_list1 <- vector()
#unconstrained
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr1")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr2")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr3")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr4")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr5")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr6")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr7")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr8")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr9")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_unconstr10")

#constrained
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr1")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr2")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr3")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr4")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr5")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr6")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr7")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr8")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr9")
iq_list1[length(iq_list1)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/reference.afa", " -m MFP -msub nuclear -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/opt_ML/main_constr10")

write(iq_list1, file=paste0(MLpath, "iq_list1.cmd"))

