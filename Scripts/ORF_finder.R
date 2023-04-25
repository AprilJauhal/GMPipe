# To run: Rscript path_to/ORF_finder.R PIPE_PATH ORF_Size query_genome_path cores
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMfind_intronless_orfs 
# Not designed to be run on its own
# Identifies all ORFs above the minimum orf size set in the control file and returns an ORF list and report

library(Biostrings)
library(stringr)
library(ORFik)
library(seqinr)

# Setting filepaths
PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
ORF_size <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
query_genome_path <- commandArgs(trailingOnly=TRUE)[3]
ORF_path <- paste0(PIPE_PATH, "/in/ORF_list.fa")
ORF_record_path <- paste0(PIPE_PATH, "/in/ORF_record.txt")

# Read and reverse-compement DNA
fwd_DNA <- readDNAStringSet(query_genome_path)
rev_DNA <- reverseComplement(fwd_DNA)
print("genomic sequences read")
# Prepare main stringsets and dataframes
all_ORFs_fwd <- AAStringSet()
all_ORF_records_fwd <- data.frame("ORF_name"=character(), "Accession"=character(), "Start"=numeric(), 
                              "End"=numeric(), "Length"=numeric(), "AA_Sequence"=character(),
                              "NT_Sequence"=character())
all_ORFs_rev <- all_ORFs_fwd
all_ORF_records_rev <- all_ORF_records_fwd
    
#### FUNCTIONS
ORF_extractor <- function(seqs, strand) {
    # Prepare temporary stringsets and dataframes
    all_ORFs_sub <- AAStringSet()
    all_ORF_records_sub <- data.frame("ORF_name"=character(), "Accession"=character(), "Start"=numeric(),
                                      "End"=numeric(), "Length"=numeric(), "AA_Sequence"=character(),
                                      "NT_Sequence"=character())
    
    each_DNA <- seqs
    names(each_DNA) <- strsplit(names(each_DNA), " ")[[1]][1]
    each_name <- names(each_DNA)
    ORF_iRanges <- findORFs(each_DNA, startCodon = "ATG", minimumLength = ORF_size)$'1'
    n=1

    while (n<=length(ORF_iRanges)){
        #Parsing iRanges data to extract ORFs
        if (strand == "F") {
            ORF_strand <- "_fwdORF_"
            each_ORF_range <- ORF_iRanges[n]
            each_ORF_start <- start(each_ORF_range)
            each_ORF_end <- end(each_ORF_range)
            ORF_dna <- substring(toString(each_DNA), each_ORF_start, each_ORF_end)
        } else if (strand == "R" ) {
            ORF_strand <- "_revORF_"
            each_ORF_range <- ORF_iRanges[n]
            each_ORF_start_R <- start(each_ORF_range) # start for the reversed strand
            each_ORF_end_R <- end(each_ORF_range) # end for the reversed strand
            each_ORF_start <- width(each_DNA)+1-each_ORF_start_R # start for the non-reversed strand
            each_ORF_end <- width(each_DNA)+1-each_ORF_end_R # end for the non-reversed strand
            ORF_dna <- substring(toString(each_DNA), each_ORF_start_R, each_ORF_end_R)
        } else {stop("Invalid Strand Entry")}
        
        # Setting counter for name generation
        counter <- n
        while (nchar(counter)<nchar(length(ORF_iRanges))){
          counter <- paste0(0, counter)
        }
        
        # Preparing ORF for sequence list 
        ORF_name <- paste0(each_name, ORF_strand, counter)
        ORF_sequence <- str_flatten(getTrans(getSequence(ORF_dna)))
        ORF <- AAStringSet(ORF_sequence)                
        names(ORF) <- ORF_name
   
        # Preparing ORF data for record
        new_ORF_record <- data.frame("ORF_name" = ORF_name, "Accession" = each_name, "Start" = each_ORF_start, 
                                     "End"=  each_ORF_end, "Length" = width(each_ORF_range), 
                                     "AA_Sequence" = ORF_sequence, "NT_Sequence" = ORF_dna)
        
        # Grouping variables
        all_ORFs_sub <- (c(all_ORFs_sub, ORF))
        all_ORF_records_sub <- rbind(all_ORF_records_sub, new_ORF_record) 
        
        n=n+1
    }
    # Return ORFs and records to be grouped together after parsing is finished
    list(all_ORFs_sub, all_ORF_records_sub)
}
##################
    
# Find ORFs for each fwd sequence in parallel
fwd_results <- list()
for (i in 1:length(fwd_DNA)) {
  print(paste("finding fwd ORFs for:", names(fwd_DNA[i])))
  fwd_results[[i]] <- ORF_extractor(fwd_DNA[i], "F")
}
rev_results <- list()
for (i in 1:length(rev_DNA)) {
  print(paste("finding rev ORFs for:", names(rev_DNA[i])))
  rev_results[[i]] <- ORF_extractor(rev_DNA[i], "R")
}  
    
# Combine results
all_ORFs_fwd <- do.call(c, lapply(fwd_results, "[[", 1))
all_ORF_records_fwd <- do.call(rbind, lapply(fwd_results, "[[", 2))
all_ORFs_rev <- do.call(c, lapply(rev_results, "[[", 1))
all_ORF_records_rev <- do.call(rbind, lapply(rev_results, "[[", 2))
all_ORFs <- c(all_ORFs_fwd, all_ORFs_rev)
all_ORF_records <- rbind(all_ORF_records_fwd, all_ORF_records_rev)

# Writing ORF list and ORF record
writeXStringSet(all_ORFs, file=ORF_path)
write.table(all_ORF_records, file=ORF_record_path)

