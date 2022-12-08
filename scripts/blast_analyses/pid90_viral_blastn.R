library(data.table)

blastn <- fread("output/recip_blast/viral_nr_blastx/pid_90_blastn.outfmt6")
transcript_lengths <- fread("output/trinity_abundance/RSEM.isoforms.results")
transcript_lengths <- transcript_lengths[,c(1,3)]

setnames(blastn, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("transcript_id", "nt_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blastn, transcript_id, evalue, -bit_score)

min_evalues <- blastn[,.SD[which.min(evalue)], by=transcript_id]
min_evalues_lengths <- merge(min_evalues, transcript_lengths, by="transcript_id")
fwrite(min_evalues_lengths, "output/recip_blast/viral_nr_blastx/virus_blastn_best.csv")

dwv <- subset(min_evalues_lengths, grepl("Deformed wing", min_evalues$annotation))
moku <- subset(min_evalues_lengths, grepl("Moku", min_evalues$annotation))
iflavirus_best_hits <- full_join(dwv, moku)
fwrite(iflavirus_best_hits, "output/recip_blast/viral_nr_blastx/iflavirus_blastn_best.csv")

##TPM values
trinity_abundance <- fread("output/trinity_abundance/RSEM.isoforms.results")

dwv_virus_abundance <- subset(trinity_abundance, transcript_id %in% dwv$transcript_id)
moku_virus_abundance <- subset(trinity_abundance, transcript_id %in% moku$transcript_id)
