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

library(stringr) 
packageVersion("stringr")

library(readr)
packageVersion("readr")

#set seed for reproducibility
set.seed(2209)

setwd("F:/Input_files/")

# Get path of directory that contains all of the reads

path <- "./trim_seq" 

#List out the content of the directory

list.files(path)

## Store path of forward and reverse reads

fwd_files <- sort(list.files(path, pattern = "R1", full.names = TRUE)) 
rev_files <- sort(list.files(path, pattern = "R2", full.names = TRUE))

#Splitting the name of the samples of the reads

samples = str_extract(basename(fwd_files), "^[^_]+")


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

fwd_filt <- file.path(filt_dir, basename(fwd_files))
rev_filt <- file.path(filt_dir, basename(rev_files))

names(fwd_filt) <- samples
names(rev_filt) <- samples

filtered_out <- filterAndTrim(
  fwd = fwd_files, 
  filt = fwd_filt,
  rev = rev_files,
  filt.rev = rev_filt,
  truncLen=c(270,220),
  maxEE = c(2, 2), 
  truncQ = 2, 
  rm.phix = TRUE, 
  compress = TRUE, 
  multithread = FALSE
) 


head(filtered_out)


# Forward Read quality
plotQualityProfile(fwd_filt[2])

# Reverse Read quality
plotQualityProfile(rev_filt[2])

#=============Learning the Error Rates and Infer Sequences===========#

# Forward read estimates

err_fwd <- learnErrors(fwd_filt, multithread = TRUE)

# Reverse read estimates
err_rev <- learnErrors(rev_filt, multithread = TRUE)

# Plotting the read estimates
plotErrors(err_fwd, nominalQ = TRUE)


# Plotting the read estimates
plotErrors(err_rev, nominalQ=TRUE)



# Infer and denoising the forward reads
dada_fwd <- dada(fwd_filt, err = err_fwd, multithread = TRUE, pool="pseudo")


# Infer and denoising the reverse reads
dada_rev <- dada(rev_filt, err=err_rev, multithread=TRUE, pool="pseudo")

#Inspecting the returned object
dada_fwd[[2]]


#================Merge Forward and Reverse Reads==================#
#The forward and reverse reads together to form a contig
merging_reads <- mergePairs(dada_fwd, fwd_filt, 
                              dada_rev, rev_filt, verbose=TRUE)
head(merging_reads[[1]])



#================Constructing a Sequence Table====================#
seqtab <- makeSequenceTable(merging_reads)
dim(seqtab) 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))  

saveRDS(seqtab, "./seqtab.rds")
#===================Removing Chimeras============================#

#dada2 will align each ASV to the other ASVs, 
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab_nochim)

sum(seqtab_nochim)/sum(seqtab) 

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab_nochim)))

#==============Tracking Reads throughout the process================#	

# small function to get the number of sequences
getN <- function(x) sum(getUniques(x))

#making a summary table of all the trimmed, denoised and merged sequences
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_fwd, getN),
                          dada_r=sapply(dada_rev, getN), merged=sapply(mergering_reads, getN),
                          nonchim=rowSums(seqtab_nochim),
                          final_perc_reads_retained=round(rowSums(seqtab_nochim)/filtered_out[,1]*100, 1))

summary_tab

#Prints out the specified  file in the work dir 
write.table(summary_tab, "./output_files/Track_read_count.tsv", quote=FALSE, sep="\t", col.names=NA)


#======================Assigning taxonomy========================#

#Link to download the reference the dataset https://zenodo.org/record/4587955

Silva <- "./ref_data/silva_nr99_v138.1_train_set.fa.gz"

Silva_w_sp <- "./ref_data/silva_species_assignment_v138.1.fa.gz"

RDP <- "./ref_data/rdp_19_toGenus_trainset.fa.gz"

GTDB <- "./ref_data/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz" 


seqs <- getSequences(seqtab_nochim)


taxa_silva <- assignTaxonomy(seqtab_nochim, Silva, minBoot = 80, tryRC=TRUE, multithread=TRUE)
taxa_rdp <- assignTaxonomy(seqtab_nochim, RDP, minBoot = 80, tryRC=TRUE, multithread=TRUE)
taxa_gtdb <- assignTaxonomy(seqtab_nochim, GTDB, minBoot = 80, tryRC=TRUE, multithread=TRUE)


Silva_sp <- addSpecies(taxa_sil, Silva_w_sp)
spec_silva <- assignSpecies(seqs, Silva_w_sp, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)


#Prints out the specified  file in the work dir 
write.table(taxa_sil, "./output_files/Hb_bac_silva_seq.tsv", quote=FALSE, sep="\t", col.names=NA)
write.table(taxa_rdp, "./output_files/Hb_bac_rdp_seq.tsv", quote=FALSE, sep="\t", col.names=NA)
write.table(taxa_gtdb, "./output_files/Hb_bac_gtdb_seq.tsv", quote=FALSE, sep="\t", col.names=NA)

#############################################################################################
#======================Generating the output files========================#


asv_seqs <- colnames(seqtab_nochim)
asv_headers <- vector(dim(seqtab_nochim)[2], mode="character")

for (i in 1:dim(seqtab_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "./output_files/HB_bacASVs.fa")

# count table:
asv_tab <- t(seqtab_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "./output_files/HB_bac_counts.tsv", sep="\t", quote=F, col.names=NA)

##DADA2 pipeline is done here moving onto phyloseq for analysising and plotting of data..

save.image(file = "./output_files/HB_bac_dada2.RData")

