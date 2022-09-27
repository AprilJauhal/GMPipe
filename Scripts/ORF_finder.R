# To run: Rscript path_to/ORF_finder.R PIPE_PATH ORF_Size query_genome_path
# where PIPE_PATH is directory containing "storage", and "in" folders for executing GMfind_intronless_orfs 
# Not designed to be run on its own

 library(Biostrings)
 library(stringr)
 library(ORFik)
 library(seqinr)
 
 PIPE_PATH <- commandArgs(trailingOnly=TRUE)[1]
 ORF_size <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
 query_genome_path <- commandArgs(trailingOnly=TRUE)[3]
 ORF_path <- paste0(PIPE_PATH, "/in/ORF_list.fa")
 ORF_record_path <- paste0(PIPE_PATH, "/in/ORF_record.txt")
 
 all_ORFs <- vector()
 all_ORF_records <- c("ORF_name", "Accession", "Start", "End", "Length", "AA_Sequence", "NT_Sequence")
 fwd_DNA <- readDNAStringSet(query_genome_path)
 rev_DNA <- reverseComplement(fwd_DNA)
 
 #Forward 
 n=1
 while(n<=length(fwd_DNA)){
   each_fwd_DNA <- fwd_DNA[n]
   each_fwd_name <- str_split(names(each_fwd_DNA), " ")[[1]][1]
   print(c("finding ORFs for... ", each_fwd_name, " fwd strand"))
   fwd_ORF_iRanges <- findORFs(each_fwd_DNA, startCodon = "ATG", minimumLength = ORF_size)$'1'
   i=1
   while (i<=length(fwd_ORF_iRanges)){
    each_ORF_range <- fwd_ORF_iRanges[i]
    each_ORF_start <- start(each_ORF_range)
    each_ORF_end <- end(each_ORF_range)
    ORF_dna<- substring(toString(each_fwd_DNA), each_ORF_start, each_ORF_end)
    counter <- i
    while (nchar(counter)<nchar(length(fwd_ORF_iRanges))){
      counter <- paste0(0, counter)
    }
    ORF_name <- paste0(each_fwd_name, "_fwdORF_", counter)
    ORF <- str_flatten(getTrans(getSequence(ORF_dna)))
    all_ORFs[length(all_ORFs)+1] <- paste0(">", ORF_name)
    all_ORFs[length(all_ORFs)+1] <- ORF
    all_ORF_records[length(all_ORF_records)+1] <- ORF_name
    all_ORF_records[length(all_ORF_records)+1] <- each_fwd_name
    all_ORF_records[length(all_ORF_records)+1] <- each_ORF_start
    all_ORF_records[length(all_ORF_records)+1] <- each_ORF_end
    all_ORF_records[length(all_ORF_records)+1] <- width(each_ORF_range)
    all_ORF_records[length(all_ORF_records)+1] <- ORF
    all_ORF_records[length(all_ORF_records)+1] <- ORF_dna
    i=i+1
   }
   n=n+1
 }

 #Reverse
 m=1
 while(m<=length(rev_DNA)){
   each_rev_DNA <- rev_DNA[m]
   each_rev_name <- str_split(names(each_rev_DNA), " ")[[1]][1]
   print(c("finding ORFs for... ", each_rev_name, " rev strand"))
   rev_ORF_iRanges <- findORFs(each_rev_DNA, startCodon = "ATG", minimumLength = ORF_size)$'1'
   j=1
   while (j<=length(rev_ORF_iRanges)){
     each_ORF_range <- rev_ORF_iRanges[j]
     #these are the start/end for the reversed strand
     each_ORF_start <- start(each_ORF_range)
     each_ORF_end <- end(each_ORF_range)
     #these are for the non-reversed strand
     real_start <- width(each_rev_DNA)+1-each_ORF_start
     real_end <- width(each_rev_DNA)+1-each_ORF_end
 
     ORF_dna <- substring(toString(each_rev_DNA), each_ORF_start, each_ORF_end)
     counter <- j
     while (nchar(counter)<nchar(length(rev_DNA))){
       counter <- paste0(0, counter)
     }
     ORF_name <- paste0(each_rev_name, "_revORF_", counter)
     ORF <- str_flatten(getTrans(getSequence(ORF_dna)))
     
     all_ORFs[length(all_ORFs)+1] <- paste0(">", ORF_name)
     all_ORFs[length(all_ORFs)+1] <- ORF
     all_ORF_records[length(all_ORF_records)+1] <- ORF_name
     all_ORF_records[length(all_ORF_records)+1] <- each_rev_name
     all_ORF_records[length(all_ORF_records)+1] <- real_start
     all_ORF_records[length(all_ORF_records)+1] <- real_end
     all_ORF_records[length(all_ORF_records)+1] <- width(each_ORF_range)
     all_ORF_records[length(all_ORF_records)+1] <- ORF
     all_ORF_records[length(all_ORF_records)+1] <- ORF_dna       
     j=j+1
   }
   m=m+1
 }
 
 write(all_ORFs, file=ORF_path, ncolumns=1)
 write(all_ORF_records, file=ORF_record_path, ncolumns = 7)
