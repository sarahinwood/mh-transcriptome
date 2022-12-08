#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)
library(seqinr)

###########
# GLOBALS #
###########

mh_nr_blastx_file <- snakemake@input[["mh_nr_blastx_file"]]
virus_taxids_file <- snakemake@input[["virus_taxids_file"]]
tx_lengths_file <- snakemake@input[["tx_lengths_file"]]
fasta_file <- snakemake@input[["fasta_file"]]

########
# MAIN #
########

##mh results
mh_nr_blastx <- fread(mh_nr_blastx_file)
virus_taxids <- fread(virus_taxids_file, header=FALSE)
tx_lengths <- fread(tx_lengths_file)
tx_lengths <- tx_lengths[,c(1,3)]

setnames(mh_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mh_nr_blastx, transcript_id, evalue, -bit_score)

##extract result with lowest evalue for each peptide (use bit-score if tie)
mh_min_evalues <- mh_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations using taxids
viral_hits <- subset(mh_min_evalues, taxid %in% virus_taxids$V1)
viral_hits_lengths <- merge(viral_hits, tx_lengths, by="transcript_id")
fwrite(viral_hits_lengths, snakemake@output[["best_viral_hits"]])

##PID above 90
pid90 <- subset(viral_hits_lengths, viral_hits_lengths$`%_identical_matches`>90)
# subset fasta for only pid90 transcripts
fasta <- read.fasta(file=fasta_file, as.string=T)
pid90_fasta <- fasta[names(fasta) %in% pid90$transcript_id]
write.fasta(pid90_fasta, names=names(pid90_fasta), file.out=snakemake@output[["pid90_fasta"]])

#write log
sessionInfo()