 library(Biostrings)
 library(stringr)
 library(ORFik)
 library(seqinr)
 
 PATH <- commandArgs(trailingOnly=TRUE)[1]
 ORF_size <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
 query_genome_path <- commandArgs(trailingOnly=TRUE)[3]
 ORF_path <- paste0(PATH, "/in/ORF_list.fa")
 ORF_record_path <- paste0(PATH, "/in/ORF_record.txt")
 
 all_ORFs <- vector()
 all_ORF_records <- c("ORF_name", "Accession", "Start", "End", "Length", "Sequence")
 fwd_DNA <- readDNAStringSet(query_genome_path)
 rev_DNA <- reverseComplement(fwd_DNA)
 
 #Forward 
 n=1
 while(n<=length(fwd_DNA)){
   each_fwd_DNA <- fwd_DNA[n]
   each_fwd_name <- str_split(names(each_fwd_DNA), " ")[[1]][1]
   print(c("finding ORFs for... ", each_fwd_name, " fwd strand"))
   #made it to here
   fwd_ORF_iRanges <- findORFs(each_fwd_DNA, startCodon = "ATG", minimumLength = ORF_size)$'1'
   i=1
   #didn't make it to here
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
 
     ORF_dna<- substring(toString(each_rev_DNA), each_ORF_start, each_ORF_end)
     counter <- j
     while (nchar(counter)<nchar(length(rev_DNA))){
       counter <- paste0(0, counter)
     }
     ORF_name <- paste0(each_rev_name, "_revORF_", counter)
     ORF <- str_flatten(getTrans(getSequence(ORF_dna)))
     
     ##Verification of adjusted indexes for reverse strand
     #print(DNAString(substring(toString(rev_DNA[m]), each_ORF_start, each_ORF_end)))
     #print(reverseComplement(DNAString(substring(toString(fwd_DNA[m]), real_end, real_start))))
     
     all_ORFs[length(all_ORFs)+1] <- paste0(">", ORF_name)
     all_ORFs[length(all_ORFs)+1] <- ORF
     all_ORF_records[length(all_ORF_records)+1] <- ORF_name
     all_ORF_records[length(all_ORF_records)+1] <- each_rev_name
     all_ORF_records[length(all_ORF_records)+1] <- each_ORF_start
     all_ORF_records[length(all_ORF_records)+1] <- each_ORF_end
     all_ORF_records[length(all_ORF_records)+1] <- width(each_ORF_range)
     all_ORF_records[length(all_ORF_records)+1] <- ORF
     j=j+1
   }
   m=m+1
 }
 
 write(all_ORFs, file=ORF_path, ncolumns=1)
 write(all_ORF_records, file=ORF_record_path, ncolumns = 6)
