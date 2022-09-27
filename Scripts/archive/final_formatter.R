library(stringr)

PATH <- commandArgs(trailingOnly=TRUE)
final_report_path <- paste0(PATH, "/storage/final_report.txt")
hit_tbl_path <- paste0(PATH, "/storage/hmm/hit.tbl")
hit_domtblout_path <- paste0(PATH, "/storage/hmm/hit.domtblout")
in_domtblout_path <- paste0(PATH, "/storage/hmm/ingroup.domtblout")
out_domtblout_path <- paste0(PATH, "/storage/hmm/outgroup.domtblout")

final_report_scores_path <- paste0(PATH, "/out/FINAL_REPORT.txt")
hit_tbl_cutoff_path <-paste0(PATH, "/out/passing_hit.tbl")
hit_domtbl_cutoff_path <-paste0(PATH, "/out/domain/passing_hit.domtbl")
final_report <- read.table(final_report_path, header = TRUE)
hit_tbl <- read.table(hit_tbl_path, header = TRUE)
hit_domtblout <- read.table(hit_domtblout_path, header = TRUE)

hit_domtbl_path <-paste0(PATH, "/out/domain/total_hit.domtbl")
hit_dom_tblout <- readLines(hit_domtblout_path)[-1]
hit_dom_tblout[1] <- str_remove(hit_dom_tblout[1], "# ")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "#", "num_domains")
hit_dom_tblout <- hit_dom_tblout[-2]
hit_dom_tblout <- hit_dom_tblout[-grep("#", hit_dom_tblout)]
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "target name", "target_name")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "query name", "query_name")
hit_dom_tblout[1] <- str_replace(hit_dom_tblout[1], "description of target", "description_of_target")
write(hit_dom_tblout, file = hit_domtbl_path)

hit_domtbl <- read.table(hit_domtbl_path, header = TRUE)

in_domtbl_path <-paste0(PATH, "/out/domain/in.domtbl")
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

out_domtbl_path <-paste0(PATH, "/out/domain/out.domtbl")
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

final_list <- as.vector(final_report$ORF_name)
hit_tbl_list <- which(hit_tbl$target_name %in% final_list)
hit_domtbl_list <- which(hit_domtbl$target_name %in% final_list)
hit_tbl_cutoff <- hit_tbl[hit_tbl_list,]
hit_domtbl_cutoff <- hit_domtbl[hit_domtbl_list,]

if (nrow(final_report)!=nrow(hit_tbl_cutoff)){warning("something is wrong with hmmer score tables")}
hit_tbl_cutoff %>% arrange(target_name)
final_report %>% arrange(ORF_name)

final_report$bit_score <- hit_tbl_cutoff$score

write.table(hit_tbl_cutoff, file=hit_tbl_cutoff_path)
write.table(hit_domtbl_cutoff, file=hit_domtbl_cutoff_path)
write.table(final_report, file=final_report_scores_path)

