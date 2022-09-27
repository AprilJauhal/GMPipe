library(stringr)
PATH <- commandArgs(trailingOnly=TRUE)[1]
IQ_PATH <- commandArgs(trailingOnly=TRUE)[2]
MUSCLE_PATH <- commandArgs(trailingOnly=TRUE)[3]
master_afa <- paste0(PATH, "/storage/reference.afa")
in_list <- paste0(PATH, "/in/master_seq.fa")
out_list <- paste0(PATH, "/in/outgroups.fa")
query_dir <- dir(paste0(PATH, "/storage/hmm/pass"))

in_names <- str_remove_all(grep(">", readLines(in_list), value=TRUE), ">")
in_names <- str_flatten(str_remove_all(in_names, " "), collapse=",")
out_names <- str_remove_all(grep(">", readLines(out_list), value=TRUE), ">")
out_names <- str_flatten(str_remove_all(out_names, " "), collapse=",")

n=1
best_models <- vector()
while (n<=10){
 model_file <- paste0(PATH, "/storage/opt_ML/main_unconstr", n, ".log")
 best_line <- grep("Best-fit", readLines(model_file), value=TRUE)
 best_models[length(best_models)+1] <- str_split(str_split(best_line, ": ")[[1]][2], " ")[[1]][1]
 n=n+1
}
best_model <- names(sort(table(best_models),decreasing=TRUE))[1]

#Parsing comparison file
comparison_lines <- readLines(paste0(PATH, "/storage/opt_ML/concat/main_comparison.iqtree"))
AU_index <- grep("p-AU", comparison_lines)[1]
index <- AU_index+2
unconstr_scores <- ""
while (index <= AU_index+11) {
  unconstr_line <- str_remove_all(comparison_lines[index], " ")
  unconstr_scores <- paste0(unconstr_scores, str_sub(unconstr_line, start = -1))
  index <- index+1
}
constr_scores <- ""
while (index <= AU_index+21) {
  constr_line <- str_remove_all(comparison_lines[index], " ")
  constr_scores <- paste0(constr_scores, str_sub(constr_line, start = -1))
  index <- index+1
}

unconstr_pass <- str_count(unconstr_scores, "[+]")
constr_pass <- str_count(constr_scores, "[+]")
if(isTRUE(constr_pass<9) | isTRUE(unconstr_pass<9)) warning('ingroups and outgroups not significantly separated in tree, please adjust reference sequences')

iq_list2 <- vector()
iq_list3 <- vector()
muscle_list <- vector()
for (query in query_dir) {
   query_name <- str_remove_all(query, ".fa")
   if (isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr1.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr2.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr3.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr4.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr5.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr6.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr7.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr8.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr9.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_unconstr10.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr1.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr2.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr3.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr4.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr5.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr6.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr7.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr8.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr9.iqtree"))) |
       isFALSE(file.exists(paste0(PATH, "/storage/ML/", query_name, "_constr10.iqtree")))) {
  
  #alignment file
  muscle_list[length(muscle_list)+1]<- paste0(MUSCLE_PATH, " -profile -in1 ", PATH, "/storage/reference.afa -in2 ", PATH, "/storage/hmm/pass/", query_name, ".fa -out ", PATH, "/storage/ML/", query_name, ".afa",  " -quiet")
 
  #unconstrained
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr1")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr2")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr3")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr4")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr5")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr6")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr7")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr8")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr9")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_unconstr10")
  #ingroup constraint
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr1")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr2")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr3")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr4")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr5")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr6")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr7")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr8")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr9")
  iq_list2[length(iq_list2)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -g ", PATH, "/storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre ", PATH, "/storage/ML/", query_name, "_constr10")
   }
  #comparison
  if (isFALSE(file.exists(paste0(PATH, "/storage/ML/concat/", query_name, "_comparison.iqtree")))) {
  iq_list3[length(iq_list3)+1] <- paste0("cat ", paste0(PATH, "/storage/ML/", query_name, "_unconstr1.treefile", " "), 
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr2.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr3.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr4.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr5.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr6.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr7.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr8.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr9.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_unconstr10.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr1.treefile", " "), 
                                         paste0(PATH, "/storage/ML/", query_name, "_constr2.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr3.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr4.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr5.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr6.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr7.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr8.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr9.treefile", " "),
                                         paste0(PATH, "/storage/ML/", query_name, "_constr10.treefile", " "),
                                         "> ", paste0(PATH, "/storage/ML/concat/", query_name, ".treelst")
                                         )
  iq_list3[length(iq_list3)+1] <- paste0(IQ_PATH, " -s ", PATH, "/storage/ML/", query_name, ".afa", " -m ", best_model, " -redo -z ", PATH, "/storage/ML/concat/", query_name, ".treelst -n 0 -zb 10000 -au -pre ", PATH, "/storage/ML/concat/", query_name, "_comparison")
  }}

iq_list2 <- c(muscle_list, iq_list2)

write(iq_list2, file=paste0(PATH, "/storage/ML/iq_list2.cmd"))
write(iq_list3, file=paste0(PATH, "/storage/ML/iq_list3.cmd"))
