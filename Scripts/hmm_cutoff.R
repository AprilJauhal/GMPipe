# To run: Rscript path_to/hmm_cutoff.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own
# Purpose: to determine an appropriate cutoff value for hmmer bitscores

library(stringr)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
ingroup_tblout_path <-paste0(PIPE_PATH, "/storage/hmm/ingroup.tblout")
outgroup_tblout_path <- paste0(PIPE_PATH, "/storage/hmm/outgroup.tblout")
cutoff_path <- paste0(PIPE_PATH, "/storage/hmm/bit_cutoff.txt")
ingroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/ingroup.tbl")
outgroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/outgroup.tbl")

# parsing hmmsearch results for outgroup seqeuences vs. ingroup HMM profile
ingroup_tbl_header <- strsplit(str_replace_all(str_replace_all(str_replace_all(str_remove(readLines(ingroup_tblout_path)[2], "#")
                                                                                  , "description of target", "description_of_target"),
                                                                  "target name", "target_name"),
                                                  "query name", "query_name"), "\\s+")[[1]][-1]
ingroup_tbl <- read.table(ingroup_tblout_path)
colnames(ingroup_tbl) <- ingroup_tbl_header
write.table(ingroup_tbl, file = ingroup_tbl_path)

# parsing hmmsearch results for outgroup seqeuences vs. outgroup HMM profile 
outgroup_tbl_header <- strsplit(str_replace_all(str_replace_all(str_replace_all(str_remove(readLines(outgroup_tblout_path)[2], "#")
                                                                                , "description of target", "description_of_target"),
                                                                "target name", "target_name"),
                                                "query name", "query_name"), "\\s+")[[1]][-1]
outgroup_tbl <- read.table(outgroup_tblout_path)
colnames(outgroup_tbl) <- outgroup_tbl_header
write.table(outgroup_tbl, file = outgroup_tbl_path)

# Cut-off determination
min_ingroup_bit <- min(ingroup_tbl$score)
max_outgroup_bit <- max(outgroup_tbl$score)
print(c("smallest ingroup hmmer bitscore: ", min_ingroup_bit))
print(c("largest outgroup hmmer bitscore: ", max_outgroup_bit))
threshold <- round(max_outgroup_bit)

# Verification that outgroup scores for ingroup and outgroup HMM Profiles do not overlap
if (min_ingroup_bit<=max_outgroup_bit) {stop("ERROR: INGROUP AND OUTGROUP HMMER BIT SCORES TOO SIMILAR")}
if (min_ingroup_bit>max_outgroup_bit) {cutoff <- paste0("CUTOFF=", threshold)}

# Report cutoff value
write(cutoff, file = cutoff_path)
