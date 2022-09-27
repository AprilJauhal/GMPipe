# To run: Rscript path_to/hmm_cutoff.R PIPE_PATH
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMPipeline 
# Not designed to be run on its own

library(stringr)

PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
ingroup_tblout_path <-paste0(PIPE_PATH, "/storage/hmm/ingroup.tblout")
outgroup_tblout_path <- paste0(PIPE_PATH, "/storage/hmm/outgroup.tblout")
cutoff_path <- paste0(PIPE_PATH, "/storage/hmm/bit_cutoff.txt")

# parsing hmmsearch results for outgroup seqeuences vs. ingroup HMM profile
ingroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/ingroup.tbl")
ingroup_tblout <- readLines(ingroup_tblout_path)[-1]
ingroup_tblout[1] <- str_remove(ingroup_tblout[1], "# ")
ingroup_tblout <- ingroup_tblout[-2]
ingroup_tblout <- ingroup_tblout[-grep("#", ingroup_tblout)]
ingroup_tblout[1] <- str_replace(ingroup_tblout[1], "target name", "target_name")
ingroup_tblout[1] <- str_replace(ingroup_tblout[1], "query name", "query_name")
ingroup_tblout[1] <- str_replace(ingroup_tblout[1], "description of target", "description_of_target")
write(ingroup_tblout, ncolumns = 1, file = ingroup_tbl_path)
ingroup_tbl <- read.table(ingroup_tbl_path, header = TRUE)

# parsing hmmsearch results for outgroup seqeuences vs. outgroup HMM profile
outgroup_tbl_path <-paste0(PIPE_PATH, "/storage/hmm/outgroup.tbl")
outgroup_tblout <- readLines(outgroup_tblout_path)[-1]
outgroup_tblout[1] <- str_remove(outgroup_tblout[1], "# ")
outgroup_tblout <- outgroup_tblout[-2]
outgroup_tblout <- outgroup_tblout[-grep("#", outgroup_tblout)]
outgroup_tblout[1] <- str_replace(outgroup_tblout[1], "target name", "target_name")
outgroup_tblout[1] <- str_replace(outgroup_tblout[1], "query name", "query_name")
outgroup_tblout[1] <- str_replace(outgroup_tblout[1], "description of target", "description_of_target")
write(outgroup_tblout, ncolumns = 1, file = outgroup_tbl_path)
outgroup_tbl <- read.table(outgroup_tbl_path, header = TRUE)

# Cut-off determination
min_ingroup_bit <- min(ingroup_tbl$score)
max_outgroup_bit <- max(outgroup_tbl$score)
print(c("smallest ingroup hmmer bitscore: ", min_ingroup_bit))
print(c("largest outgroup hmmer bitscore: ", max_outgroup_bit))
threshold <- round(max_outgroup_bit)

# Verification that outgroup scores for ingroup and outgroup HMM Profiles do not overlap
if (min_ingroup_bit<=max_outgroup_bit) {stop("ERROR: INGROUP AND OUTGROUP HMMER BIT SCORES TOO SIMILAR")}
if (min_ingroup_bit>max_outgroup_bit) {cutoff <- paste0("CUTOFF=", threshold)}

write(cutoff, file = cutoff_path)
