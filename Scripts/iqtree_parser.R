# To run: Rscript path_to/iqtree_parser.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose: to parse IQTREE results to determine if hits phylogenetically cluster with ingroup 
# vs. outgroup sequences

library(stringr)
library(ape)
library(Biostrings)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
master_afa <- paste0(PIPE_PATH, "/storage/reference.afa")
test_dir <- dir(paste0(PIPE_PATH, "/storage/ML/"))
test_path <- (paste0(PIPE_PATH, "/storage/ML/"))
ingroups_path <- paste0(PIPE_PATH, "/in/master_seq.fa")
outgroups_path <- paste0(PIPE_PATH, "/in/outgroups.fa")
query_dir <- dir(paste0(PIPE_PATH, "/storage/hmm/pass/"))
query_path <- paste0(PIPE_PATH, "/storage/hmm/pass/")
tree_report_path <- paste0(PIPE_PATH, "/storage/final_report.txt")
pass_path <- paste0(PIPE_PATH, "/storage/ML_pass.fa")
fail_path <- paste0(PIPE_PATH, "/storage/ML_fail.fa")
undet_path <- paste0(PIPE_PATH, "/storage/ML_undet.fa")
undet_notes1 <-  paste0(PIPE_PATH, "/storage/undet_by_tree.fa")
undet_notes2 <-  paste0(PIPE_PATH, "/storage/undet_by_constraint.fa")
ORF_list_path <- (paste0(PIPE_PATH, "/in/ORF_list.fa"))
ORF_record_path <- paste0(PIPE_PATH, "/storage/hmm/hit_report.txt")

# Reading ORF hit record (if applicable)
if (file.exists(ORF_record_path)) {
    ORF_hit_record <- read.table(ORF_record_path)
}

# Reading ORF list
ORF_list <- readAAStringSet(ORF_list_path)

# Extracting gene names from input files
ingroups <- str_squish(names(readAAStringSet(ingroups_path)))
outgroups <- str_squish(names(readAAStringSet(outgroups_path)))

###FUNCTIONS###
iqtree_parser <- function(tr, ingroup, outgroup, query_name) {
    
    #Counting Nodes
    Node_num <- tr$Nnode
    Max_num <- max(tr$edge)
    Tips_num <- Max_num-Node_num
    
    #Parsing branches from tree
    N=Tips_num+1
    if (exists("node_df")) {rm(node_df)}
    while (N <= Max_num ){
        clade_vector <- extract.clade(tr, node=N)$tip.label
        clade_contents <- str_flatten(clade_vector, ";")

        if(query_name %in% clade_vector) {
            node_row <- data.frame(Node=N, Contents=clade_contents, Count=length(clade_vector))
        }

        if (exists("node_df")) {node_df <- rbind(node_df, node_row)}
        if (isFALSE(exists("node_df"))) {node_df <- node_row}

        N=N+1
    }
    
    #Checking contents of smallest branch with query_name
    smallest_node_row <- subset(node_df, Count==min(as.numeric(node_df[["Count"]])))[1,]
    smallest <- str_split(as.character(smallest_node_row$Contents), ";")[[1]]
    incheck <- "no"
    outcheck <- "no"
    for (each in smallest){
      if (each==query_name){}
      if (each %in% ingroup) {incheck <- "yes"}
      if (each %in% outgroup) {outcheck <- "yes"}
    }
    if (incheck=="yes" & outcheck=="no") {result <- "ingroup"}
    if (incheck=="no" & outcheck=="yes") {result <- "outgroup"}
    if (incheck=="yes" & outcheck=="yes") {result <- "undetermined"}
    if (incheck=="no" & outcheck=="no") {result <- "undetermined"}
    
    return(result)
}
###END FUNCTIONS###

ORF_pass <- AAStringSet()
ORF_fail <- AAStringSet()
ORF_undet <- AAStringSet()
undet_ORFs1 <- AAStringSet()
undet_ORFs2 <- AAStringSet()

# Parsing IQTREE sequences for each query in query directory
for (each_query in query_dir){
  query_name <- str_remove(each_query, ".fa")
  each_test_dir <- dir(paste0(test_path, query_name))
  # Verifying presence of all expected IQTREE files
  if (paste0(query_name, "_constr01.treefile") %in% each_test_dir &
      paste0(query_name, "_constr02.treefile") %in% each_test_dir &
      paste0(query_name, "_constr03.treefile") %in% each_test_dir &
      paste0(query_name, "_constr04.treefile") %in% each_test_dir &
      paste0(query_name, "_constr05.treefile") %in% each_test_dir &
      paste0(query_name, "_constr06.treefile") %in% each_test_dir &
      paste0(query_name, "_constr07.treefile") %in% each_test_dir &
      paste0(query_name, "_constr08.treefile") %in% each_test_dir &
      paste0(query_name, "_constr09.treefile") %in% each_test_dir &
      paste0(query_name, "_constr10.treefile") %in% each_test_dir) {
    # Determining if each ORF clusters with ingroups using the iqtree_parser function
    n=1
    in_counter <- 0
    out_counter <- 0
    undet_counter <- 0
    while (n<=10){
      if (n<10) {treename <- paste0(query_name, "_constr0", n, ".treefile")}
      if (n==10) {treename <- paste0(query_name, "_constr", n, ".treefile")}
      tree <- read.tree(paste0(test_path, query_name, "/", treename))
      result <- iqtree_parser(tree, ingroups, outgroups, query_name)
      if (result=="ingroup") {in_counter <- in_counter+1}
      if (result=="outgroup") {out_counter <- out_counter+1}
      if (result=="undetermined") {undet_counter <- undet_counter+1}
      n=n+1
    }
    if (in_counter>=9){final_result <- "pass"
    } else if (out_counter>=9){final_result <- "fail"
    } else {final_result <- "undet"
            undet_ORFs1 <- AAStringSet(c(undet_ORFs1, ORF_list[names(ORF_list)==query_name]))
           }
    
    # Saving results (if there is a hit report)
    if (file.exists(ORF_record_path)) {
        new_result <- subset(ORF_hit_record, ORF_name==query_name)
        new_result$tree <- final_result
        if (exists("tree_report")){tree_report <- rbind(tree_report, new_result)} 
        if (isFALSE(exists("tree_report"))) {tree_report <- new_result}
    }
    
    # Parsing comparison file for tree statistics (pAU-test results)
    comparison_lines <- readLines(paste0(PIPE_PATH, "/storage/ML/concat/", query_name, "/", query_name, "_comparison.iqtree"))
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

    # Flagging ORFs as undetermined based on tree statistics
    unconstr_pass <- str_count(unconstr_scores, "[+]")
    constr_pass <- str_count(constr_scores, "[+]")
    if(isTRUE(constr_pass<9) | isTRUE(unconstr_pass<9)) {
        undet_ORFs2 <- AAStringSet(c(undet_ORFs2, ORF_list[names(ORF_list)==query_name]))
        final_result <- "undet" 
    }
    
    # Adding ORF sequence to passing, failing, and undetermined lists
    ORF <- ORF_list[names(ORF_list)==query_name]
    if (final_result=="pass"){ORF_pass <- AAStringSet(c(ORF_pass, ORF))}
    if (final_result=="fail"){ORF_fail <- AAStringSet(c(ORF_fail, ORF))}
    if (final_result=="undet"){ORF_undet <- AAStringSet(c(ORF_undet, ORF))}
  }
  # Stopping script in the absence of IQTREE files
  else {
      stop('Error: iqtree_parser.R could not find treefiles')}
}

# Adding tree results to hit_record as a "tree_report" (if applicable)
if (file.exists(ORF_record_path)) {
    hit_row <- 1
    while (hit_row<=nrow(ORF_hit_record)){
      if (isFALSE(as.character(ORF_hit_record$Hit[hit_row]) %in% tree_report$Hit)) {
        new_result <- ORF_hit_record[hit_row,]
        new_result$tree <- "NA"
        tree_report <- rbind(tree_report, new_result)
        }
      hit_row <- hit_row+1
    }
    write.table(tree_report, tree_report_path)
}

# Writing output files
writeXStringSet(ORF_pass, file = pass_path)
writeXStringSet(ORF_fail, file = fail_path)
writeXStringSet(ORF_undet, file = undet_path)
writeXStringSet(undet_ORFs1, file=undet_notes1)
writeXStringSet(undet_ORFs2, file=undet_notes2)
