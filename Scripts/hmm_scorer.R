# To run: Rscript path_to/hmm_scorer.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own

library(stringr)

PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1] 
ORF_path <- paste0(PIPE_PATH, "/in/ORF_list.fa")
ORF_record_path <- paste0(PIPE_PATH, "/in/ORF_record.txt")
ORF_hit_record_path <- paste0(PIPE_PATH, "/storage/hmm/hit_report.txt")
hit_list_path <- paste0(PIPE_PATH, "/storage/HMMER_PASS.fa")
miss_list_path <- paste0(PIPE_PATH, "/storage/HMMER_FAIL.fa")
hit_dir_path <- paste0(PIPE_PATH, "/storage/hmm/pass/")
cutoff_path <- paste0(PIPE_PATH, "/storage/hmm/bit_cutoff.txt")
hit_tblout_path <-paste0(PIPE_PATH, "/storage/hmm/hit.tblout")

cutoff <- as.numeric(str_remove_all(readLines(cutoff_path), "[A-Z]+="))
ORF_record <- read.table(ORF_record_path, header = TRUE)

# parsing hmmsearch results for each ORF seqeuences vs. ingroup HMM profile
hit_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/hit.tbl")
hit_tblout <- readLines(hit_tblout_path)[-1]
hit_tblout[1] <- str_remove(hit_tblout[1], "# ")
hit_tblout <- hit_tblout[-2]
hit_tblout <- hit_tblout[-grep("#", hit_tblout)]
hit_tblout[1] <- str_replace(hit_tblout[1], "target name", "target_name")
hit_tblout[1] <- str_replace(hit_tblout[1], "query name", "query_name")
hit_tblout[1] <- str_replace(hit_tblout[1], "description of target", "description_of_target")
write(hit_tblout, ncolumns = 1, file = hit_tbl_path)
hit_tbl <- read.table(hit_tbl_path, header = TRUE)

meets_cutoff <- as.vector(subset(hit_tbl, score >= cutoff)$target_name)
fails_cutoff <- as.vector(subset(hit_tbl, score < cutoff)$target_name)

ORF_hit_record <- subset(ORF_record, ORF_name %in% meets_cutoff)
ORF_hit_names <- as.vector(ORF_hit_record$ORF_name)
ORF_hit_seq <- str_replace(as.vector(ORF_hit_record$Sequence), "[*]", "")

ORF_miss_record <- subset(ORF_record, ORF_name %in% fails_cutoff)
ORF_miss_names <- as.vector(ORF_miss_record$ORF_name)
ORF_miss_seq <- str_replace(as.vector(ORF_miss_record$Sequence), "[*]", "")

# Retrieving sequences for hits
hit_index=1
hit_list <- vector()
while (hit_index<=length(ORF_hit_names)){
  hit_list[length(hit_list)+1] <- paste0(">", ORF_hit_names[hit_index])
  hit_list[length(hit_list)+1] <- ORF_hit_seq[hit_index]
  
  write(c(paste0(">", ORF_hit_names[hit_index]), ORF_hit_seq[hit_index]), 
        file = paste0(hit_dir_path, ORF_hit_names[hit_index], ".fa"), 
        ncolumns = 1)
  
  hit_index=hit_index+1
}

# Retrieving sequences for misses
miss_index=1
miss_list <- vector()
while (miss_index<=length(ORF_miss_names)){
  miss_list[length(miss_list)+1] <- paste0(">", ORF_miss_names[miss_index])
  miss_list[length(miss_list)+1] <- ORF_miss_seq[miss_index]
  miss_index=miss_index+1
}

write(hit_list, file=hit_list_path, ncolumns = 1)
write(miss_list, file=miss_list_path, ncolumns = 1)
write.table(ORF_hit_record, file = ORF_hit_record_path)



