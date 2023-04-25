# To run: Rscript path_to/hmm_scorer.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose is to identify and prepare hmmer-passing sequences for further analysis

library(stringr)
library(Biostrings)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1] 
ORF_path <- paste0(PIPE_PATH, "/in/ORF_list.fa")
ORF_record_path <- paste0(PIPE_PATH, "/in/ORF_record.txt")
ORF_hit_record_path <- paste0(PIPE_PATH, "/storage/hmm/hit_report.txt")
hit_list_path <- paste0(PIPE_PATH, "/storage/HMMER_PASS.fa")
miss_list_path <- paste0(PIPE_PATH, "/storage/HMMER_FAIL.fa")
hit_dir_path <- paste0(PIPE_PATH, "/storage/hmm/pass/")
cutoff_path <- paste0(PIPE_PATH, "/storage/hmm/bit_cutoff.txt")
hit_tblout_path <-paste0(PIPE_PATH, "/storage/hmm/hit.tblout")
hit_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/hit.tbl")

# Reading cutoff value
cutoff <- as.numeric(str_remove_all(readLines(cutoff_path), "[A-Z]+="))

# Reading ORF fasta list
ORF_seqs <- readAAStringSet(ORF_path)
ORF_seqs <- AAStringSet(gsub("*", "", ORF_seqs))

# Parsing hmmsearch results for each ORF seqeuences vs. ingroup HMM profile
hit_tbl_header <- strsplit(str_replace_all(str_replace_all(str_replace_all(str_remove(readLines(hit_tblout_path)[2], "#")
                                                                                  , "description of target", "description_of_target"),
                                                                  "target name", "target_name"),
                                                  "query name", "query_name"), "\\s+")[[1]][-1]
hit_tbl <- read.table(hit_tblout_path)
colnames(hit_tbl) <- hit_tbl_header
write.table(hit_tbl, file = hit_tbl_path)

# Splitting hit-table based on score cutoff into passing and failing sequences
meets_cutoff <- as.vector(subset(hit_tbl, score >= cutoff)$target_name)
fails_cutoff <- as.vector(subset(hit_tbl, score < cutoff)$target_name)

# Extracting data for hmmer-passing and failing sequences
ORF_hit_seq <- ORF_seqs[names(ORF_seqs) %in% meets_cutoff]
ORF_miss_seq <- ORF_seqs[names(ORF_seqs) %in% fails_cutoff]

# Writing individual sequence files for hits
hit_index=1
while (hit_index <= length(ORF_hit_seq)){  
  writeXStringSet(ORF_hit_seq[hit_index], file = paste0(hit_dir_path, names(ORF_hit_seq[hit_index]), ".fa"))
  hit_index=hit_index+1
}

# Writing output files
writeXStringSet(ORF_hit_seq, file=hit_list_path)
writeXStringSet(ORF_miss_seq, file=miss_list_path)

# Updating ORF_record (if applicable)
if (file.exists(ORF_record_path)) {
    ORF_record <- read.table(ORF_record_path, header = TRUE)
    ORF_hit_record <- subset(ORF_record, ORF_name %in% meets_cutoff)
    write.table(ORF_hit_record, file = ORF_hit_record_path)
} 
