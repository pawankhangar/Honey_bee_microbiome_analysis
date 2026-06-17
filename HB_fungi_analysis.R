library(DECIPHER)
packageVersion("DECIPHER")

library(dada2)
packageVersion("dada2") 

library(ShortRead)
packageVersion("ShortRead") 

library(Biostrings)
packageVersion("Biostrings")

library(ggplot2)
packageVersion("ggplot2")

library(stringr) # not strictly required but handy
packageVersion("stringr")

library(readr)
packageVersion("readr")

#set seed for reproducibility
set.seed(0998)

setwd("F:/Input_files/")

# Get path of directory that contains all of the reads

path <- "./ITS_data/" 

#List out the content of the directory

list.files(path)

## Store path of forward and reverse reads

fwd_files <- sort(list.files(path, pattern = "R1", full.names = TRUE)) 
rev_files <- sort(list.files(path, pattern = "R2", full.names = TRUE))

#Splitting the name of the samples of the reads

samples = str_extract(basename(fwd_files), "^[^_]+")

samples = samples[-1]

names(fwd_files) <- samples
names(rev_files) <- samples

#==================Quality Profile====================#

# Forward Read quality
plotQualityProfile(fwd_files[1]) + ggtitle("Forward")

# Reverse Read quality
plotQualityProfile(rev_files[1]) + ggtitle("Reverse")


#==================Reads filtering====================#

filt_dir <- file.path(path, "filtered")
if (!dir.exists(filt_dir)) dir.create(filt_dir)

fwd_F <- file.path(filt_dir, basename(fwd_files))
rev_F <- file.path(filt_dir, basename(rev_files))

names(fwd_F) <- samples
names(rev_F) <- samples

filtered_out <- filterAndTrim(
  fwd = fwd_files, 
  filt = fwd_F,
  rev = rev_files,
  filt.rev = rev_F,
  trimLeft = 10,
  truncLen=c(290,210),
  maxEE = c(2, 2), 
  truncQ = 2, 
  rm.phix = TRUE, 
  compress = TRUE, 
  multithread = FALSE
)  

head(filtered_out)


# Forward Read quality
plotQualityProfile(fwd_F[2])

# Reverse Read quality
plotQualityProfile(rev_F[3])

#=============Learn the Error Rates and Infer Sequences===========#

# Forward read estimates

err_fwd <- learnErrors(fwd_F, multithread = TRUE)

# Reverse read estimates
err_rev <- learnErrors(rev_F, multithread = TRUE)

# Plotting the read estimates
plotErrors(err_fwd, nominalQ = TRUE)


# Plotting the read estimates
plotErrors(err_rev, nominalQ=TRUE)



# Infer and denoising the forward reads
dada_fwd <- dada(fwd_F, err = err_fwd, multithread = TRUE, pool="pseudo")


# Infer and denoising the reverse reads
dada_rev <- dada(rev_F, err=err_rev, multithread=TRUE, pool="pseudo")

#Inspecting the returned object
dada_fwd[[1]]


#================Merge Forward and Reverse Reads==================#
#The forward and reverse reads together to form a contig
mergers <- mergePairs(dada_fwd, fwd_F, 
                      dada_rev, rev_F, verbose=TRUE)
head(mergers[[1]])



#================Constructing a Sequence Table====================#
#Seq table is a matrix with each row representing the samples, columns are the various ASVs, 
#and each cell shows the number of that specific ASV within each sample.

seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))  


#===================Removing Chimeras============================#

#dada2 will align each ASV to the other ASVs, 
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)

sum(seqtab_nochim)/sum(seqtab) #


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_nochim)))

#==============Tracking Reads throughout the process================#	

# small function to get the number of sequences
getN <- function(x) sum(getUniques(x))

track <- cbind(
  filtered_out, 
  sapply(dada_fwd, getN), 
  sapply(dada_rev, getN), 
  sapply(mergers, getN), 
  rowSums(seqtab_nochim)
)

colnames(track) <- c("raw", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chim")
rownames(track) <- samples  
head(track)




#making a summary table of all the trimmed, denoised and merged sequences
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_fwd, getN),
                          dada_r=sapply(dada_rev, getN), merged=sapply(mergers, getN),
                          nonchim=rowSums(seqtab_nochim),
                          final_perc_reads_retained=round(rowSums(seqtab_nochim)/filtered_out[,1]*100, 1))

summary_tab

#Prints out the specified  file in the work dir 
write.table(summary_tab, "./output_files/HB_ITS_read_count_track.tsv", quote=FALSE, sep="\t", col.names=NA)

write.table(seqtab_nochim, "./output_files/seqtab_nochim.tsv", quote=FALSE, sep="\t", col.names=NA)


#======================Assigning taxonomy========================#

#Link to download the reference the dataset https://zenodo.org/record/4587955

unite.ref <- "./ref/ITS_1.fasta"

taxa_fun <- assignTaxonomy(seqtab_nochim, unite.ref, multithread=TRUE, tryRC = TRUE)

write.table(taxa_fun, "./Hb_fun_seq.tsv", quote=FALSE, sep="\t", col.names=NA)

#############################################################################################
#======================Generating the output files========================#


asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")

for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "./output_files/HB_fun_ASVs.fa")

# count table:
asv_tab <- t(seqtab_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "./output_files/HB_fun_ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

##DADA2 pipeline is done here moving onto phyloseq for analysising and plotting of data..

save.image(file = "./output/HB_fun_dada2.RData")
