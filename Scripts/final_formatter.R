# To run: Rscript path_to/final_formatter.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose: to retrieve and report hmmer output data from sequences that score above bit-cutoff

library(stringr)
library(dplyr)
library(Biostrings)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
final_report_path <- paste0(PIPE_PATH, "/storage/final_report.txt")
hit_tbl_path <- paste0(PIPE_PATH, "/storage/hmm/hit.tbl")
hit_domtblout_path <- paste0(PIPE_PATH, "/storage/hmm/hit.domtblout")
in_domtblout_path <- paste0(PIPE_PATH, "/storage/hmm/ingroup.domtblout")
out_domtblout_path <- paste0(PIPE_PATH, "/storage/hmm/outgroup.domtblout")
final_report_scores_path <- paste0(PIPE_PATH, "/out/FINAL_REPORT.txt")
hit_tbl_cutoff_path <-paste0(PIPE_PATH, "/out/passing_hit.tbl")
hit_domtbl_cutoff_path <-paste0(PIPE_PATH, "/out/domain/passing_hit.domtbl")
hit_domtbl_path <-paste0(PIPE_PATH, "/out/domain/total_hit.domtbl")
in_domtbl_path <-paste0(PIPE_PATH, "/out/domain/in.domtbl")
pass_path <- paste0(PIPE_PATH, "/storage/HMMER_PASS.fa")

# Reading tables
hit_tbl <- read.table(hit_tbl_path, header = TRUE)
hit_domtblout <- read.table(hit_domtblout_path, header = TRUE)
hit_dom_tblout <- readLines(hit_domtblout_path)[-1]

# Parsing dom_tblout hmmer output table for hits
hit_dom_tblout[1] <- str_remove(hit_dom_tblout[1], "# ")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "#", "num_domains")
hit_dom_tblout <- hit_dom_tblout[-2]
hit_dom_tblout <- hit_dom_tblout[-grep("#", hit_dom_tblout)]
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "target name", "target_name")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "query name", "query_name")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "description of target", "description_of_target")
write(hit_dom_tblout, file = hit_domtbl_path)
hit_domtbl <- read.table(hit_domtbl_path, header = TRUE)

# Parsing dom_tblout hmmer output table for ingroup sequences
in_dom_tblout <- readLines(in_domtblout_path)[-1]
in_dom_tblout[1] <- str_remove(in_dom_tblout[1], "# ")
in_dom_tblout[1] <- str_replace(in_dom_tblout[1], "#", "num_domains")
in_dom_tblout <- in_dom_tblout[-2]
in_dom_tblout <- in_dom_tblout[-grep("#", in_dom_tblout)]
in_dom_tblout[1] <- str_replace(in_dom_tblout[1], "target name", "target_name")
in_dom_tblout[1] <- str_replace(in_dom_tblout[1], "query name", "query_name")
in_dom_tblout[1] <- str_replace(in_dom_tblout[1], "description of target", "description_of_target")
write(in_dom_tblout, file = in_domtbl_path)
in_domtbl <- read.table(in_domtbl_path, header = TRUE)

# Parsing dom_tblout hmmer output table for outgroup sequences
out_domtbl_path <-paste0(PIPE_PATH, "/out/domain/out.domtbl")
out_dom_tblout <- readLines(out_domtblout_path)[-1]
out_dom_tblout[1] <- str_remove(out_dom_tblout[1], "# ")
out_dom_tblout[1] <- str_replace(out_dom_tblout[1], "#", "num_domains")
out_dom_tblout <- out_dom_tblout[-2]
out_dom_tblout <- out_dom_tblout[-grep("#", out_dom_tblout)]
out_dom_tblout[1] <- str_replace(out_dom_tblout[1], "target name", "target_name")
out_dom_tblout[1] <- str_replace(out_dom_tblout[1], "query name", "query_name")
out_dom_tblout[1] <- str_replace(out_dom_tblout[1], "description of target", "description_of_target")
write(out_dom_tblout, file = out_domtbl_path)
out_domtbl <- read.table(out_domtbl_path, header = TRUE)

# Gathering data for hmmer-passing sequences
final_list <- names(readAAStringSet(pass_path))
hit_tbl_cutoff <- subset(hit_tbl, target_name %in% final_list)
hit_domtbl_cutoff <- subset(hit_domtbl, target_name %in% final_list)

# Assembling final report (if ORF record provided)
if (file.exists(final_report_path)) {
    # Verifying and organizing data 
    final_report <- read.table(final_report_path, header = TRUE)
    if (nrow(final_report)!=nrow(hit_tbl_cutoff)){warning("something is wrong with hmmer score tables")}
    hit_tbl_cutoff %>% arrange(target_name)
    final_report %>% arrange(ORF_name)

    # Extracting final bitscores and generating report
    final_report$bit_score <- hit_tbl_cutoff$score
    write.table(final_report, file=final_report_scores_path)
}

# Writing output files
write.table(hit_tbl_cutoff, file=hit_tbl_cutoff_path)
write.table(hit_domtbl_cutoff, file=hit_domtbl_cutoff_path)
