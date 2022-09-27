library(stringr)

PATH <- commandArgs(trailingOnly=TRUE)
ORF_path <- paste0(PATH, "/in/ORF_list.fa")
ORF_record_path <- paste0(PATH, "/in/ORF_record.txt")
ORF_hit_record_path <- paste0(PATH, "/storage/hmm/hit_report.txt")
hit_list_path <- paste0(PATH, "/storage/HMMER_PASS.fa")
hit_dir_path <- paste0(PATH, "/storage/hmm/pass/")
cutoff_path <- paste0(PATH, "/storage/hmm/bit_cutoff.txt")
hit_tblout_path <-paste0(PATH, "/storage/hmm/hit.tblout")

cutoff <- as.numeric(str_remove_all(readLines(cutoff_path), "[A-Z]+="))
ORF_record <- read.table(ORF_record_path, header = TRUE)

hit_tbl_path <-paste0(PATH, "/storage/hmm/hit.tbl")
hit_tblout <- readLines(hit_tblout_path)[-1]
hit_tblout[1] <- str_remove(hit_tblout[1], "# ")
hit_tblout <- hit_tblout[-2]
hit_tblout <- hit_tblout[-grep("#", hit_tblout)]
hit_tblout[1] <- str_replace(hit_tblout[1], "target name", "target_name")
hit_tblout[1] <- str_replace(hit_tblout[1], "query name", "query_name")
hit_tblout[1] <- str_replace(hit_tblout[1], "description of target", "description_of_target")
write(hit_tblout, ncolumns = 1, file = hit_tbl_path)
hit_tbl <- read.table(hit_tbl_path, header = TRUE)

meets_cutoff <- as.vector(subset(hit_tbl, score > cutoff)$target_name)

ORF_hit_record <- subset(ORF_record, ORF_name %in% meets_cutoff)
ORF_hit_names <- as.vector(ORF_hit_record$ORF_name)
ORF_hit_seq <- str_replace(as.vector(ORF_hit_record$Sequence), "[*]", "")

index=1
hit_list <- vector()
while (index<=length(ORF_hit_names)){
  hit_list[length(hit_list)+1] <- paste0(">", ORF_hit_names[index])
  hit_list[length(hit_list)+1] <- ORF_hit_seq[index]
  
  write(c(paste0(">", ORF_hit_names[index]), ORF_hit_seq[index]), 
        file = paste0(hit_dir_path, ORF_hit_names[index], ".fa"), 
        ncolumns = 1)
  
  index=index+1
}

write(hit_list, file=hit_list_path, ncolumns = 1)
write.table(ORF_hit_record, file = ORF_hit_record_path)



