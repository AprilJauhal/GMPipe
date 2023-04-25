# To run: Rscript path_to/histogram.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose is to generate a histogram of bitscores to verify sufficient separation of on-target 
# and off-target sequences

library(ggplot2)
library(ggpubr)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
ingroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/ingroup.tbl")
outgroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/outgroup.tbl")
histo_path <- paste0(PIPE_PATH, "/out/bitscore_hist.png")
final_report_path <- paste0(PIPE_PATH, "/out/FINAL_REPORT.txt")

# Generating histogram if ORF record exists
if (file.exists(final_report_path)) {
    # Reading tables
    final_report <- read.table(final_report_path, header = TRUE)
    ingroup_tbl <- read.table(ingroup_tbl_path, header = TRUE)
    outgroup_tbl <- read.table(outgroup_tbl_path, header = TRUE)

    # Reading passing, failing, undetermined bitscores from final report
    tree_pass_score <- data.frame("bit_score" = subset(final_report, tree=="pass")$bit_score)
    tree_fail_score <- data.frame("bit_score" = subset(final_report, tree=="fail")$bit_score)
    tree_undet_score <- data.frame("bit_score" = subset(final_report, tree=="undet")$bit_score)
    ingroup_score <- data.frame("bit_score" = ingroup_tbl$score)
    outgroup_score <- data.frame("bit_score" = outgroup_tbl$score)

    # Writing bitscores from ingroup and outgroup sequences
    write(as.vector(ingroup_score$bit_score), file=paste0(PIPE_PATH, "/out/ingroup_bitscores.txt"), ncolumns=1)
    write(as.vector(outgroup_score$bit_score), file=paste0(PIPE_PATH, "/out/outgroup_bitscores.txt"), ncolumns=1)

    # Adding NA values as needed to make datasets the same size
    max_len <- max(c(nrow(tree_pass_score),
                     nrow(tree_fail_score),
                     nrow(tree_undet_score),
                     nrow(ingroup_score),
                     nrow(outgroup_score)))
    while (nrow(tree_pass_score)<max_len){
      tree_pass_score <- rbind(tree_pass_score, NA)}
    while (nrow(tree_fail_score)<max_len){
      tree_fail_score <- rbind(tree_fail_score, NA)}
    while (nrow(tree_undet_score)<max_len){
      tree_undet_score <- rbind(tree_undet_score, NA)}
    while (nrow(ingroup_score)<max_len){
      ingroup_score <- rbind(ingroup_score, NA)}
    while (nrow(outgroup_score)<max_len){
      outgroup_score <- rbind(outgroup_score, NA)}

    # Combining data into a dataframe 
    histo_table <- data.frame("tree_pass" = tree_pass_score,
                  "tree_fail" = tree_fail_score, 
                  "tree_undet" = tree_undet_score, 
                  "ingroup" = ingroup_score, 
                  "outgroup" = outgroup_score)
    colnames(histo_table) <- c("tree_pass", "tree_fail", "tree_undet", "ingroup", "outgroup")

    # Generating a histogram of all sequences
    graph1 <- gghistogram(histo_table,
                          x = c("tree_fail", "tree_undet",
                                "ingroup", "outgroup", "tree_pass"),
                          y = "..count..",
                          color = ".x.", fill = ".x.",     
                          merge = TRUE,                    
                          xlab = "bit_score", 
                          palette = "npg"                  
    )
    # Generating a histogram without tree_pass sequences for clarity
    graph2 <- gghistogram(histo_table,
                          x = c("tree_fail", "tree_undet",
                                "ingroup", "outgroup"),
                          y = "..count..",
                          color = ".x.", fill = ".x.",     
                          merge = TRUE,                    
                          xlab = "bit_score", 
                          palette = "npg"                  
    )

    # Combining histograms into a single image
    ggarrange(graph1, graph2, ncol = 1) %>%
      ggexport(filename = histo_path)
} else {print("no histogram generated due to lack of ORF_report in input files")}
