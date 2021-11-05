library(data.table)
library(dplyr)

##blast results
viral_genes_blast_res <- fread("output/MhV_blast/MhV_v_transcriptome/blastn.outfmt6")
setnames(viral_genes_blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("Prodigal_ID", "Trinity_ID", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##prodigal info & trinotate
trinotate <- fread("output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")
recip_viral_genes <- fread("output/recip_blast/nr_blastx/viral_annots.csv")
prodigal_annots <- fread("/Volumes/archive/deardenlab/sarahinwood/mh_projects/mh-viral-genome/output/prodigal_blast/blastp_gff.csv")
id_annot <- prodigal_annots[,c(4,3,13)]
setnames(id_annot, old=c("annotation", "prodigal_nt_id"), new=c("Prodigal_blast_annot", "Prodigal_ID"))

##count number of hits per MhV gene - all results
MhV_hit_counts_all <- viral_genes_blast_res[,length(unique(Trinity_ID)), by=Prodigal_ID]
setnames(MhV_hit_counts_all, old=c("V1"), new=c("all_hits_counts"))

##order blast res so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_genes_blast_res, Prodigal_ID, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- viral_genes_blast_res[,.SD[which.min(evalue)], by=Prodigal_ID]
trinity_counts_all <- viral_genes_blast_res[,length(unique(Prodigal_ID)), by=Trinity_ID]

##merge with recip blast res
trinity_counts_annots <- merge(trinity_counts_all, recip_viral_genes, by.x="Trinity_ID", by.y="transcript_id", all.x=TRUE)
trinity_counts_annots <- trinity_counts_annots[,c(1,2,17,4,12,14,15,16)]
##get rid of BRO hits
non_bro_counts <- subset(trinity_counts_annots, !grepl("bro", annotation, ignore.case=TRUE))
non_bro_counts_hits <- merge(non_bro_counts, (viral_genes_blast_res[,c(1,2)]), all.x=TRUE)
##table of transcripts, MhV gene hits for transcripts with 2+ MhV hits
transcript_v_MhVgene <- subset((non_bro_counts_hits[,c(1,9,2,3)]), V1>1)
transcript_v_MhVgene <- merge(transcript_v_MhVgene, id_annot)
setorder(transcript_v_MhVgene, -V1, Trinity_ID)
fwrite(transcript_v_MhVgene, "output/MhV_blast/MhV_v_transcriptome/transcripts_multiMhV_hits.csv")

##merge with trinotate results
min_evalues_annots <- merge(min_evalues, trinotate, by.x="Trinity_ID", by.y="transcript_id")
min_evalues_all_annots <- merge(min_evalues_annots, id_annot, by="Prodigal_ID")
##write table
fwrite(min_evalues_all_annots, "output/MhV_blast/MhV_v_transcriptome/viral_genes_best_hits.csv")

##compare to recip blast results for LbFV/DaFV
LbFV_DaFV_genes <- fread("data/mh-transcriptome/output/recip_blast/nr_blastx/LbFV_DaFV_annots.csv")
##overlap
recip_prodigal_overlap <- intersect(min_evalues_all_annots$Trinity_ID, LbFV_DaFV_genes$transcript_id)
##merge with prodigal blast res
prodigal_recip_annots <- merge(min_evalues_all_annots, LbFV_DaFV_genes, by.x="Trinity_ID", by.y="transcript_id", all.x=TRUE)
prodigal_recip_annots_simple <- prodigal_recip_annots[,c(1,2,3,42,4,32,16,20,21,24)]
fwrite(prodigal_recip_annots_simple, "output/blast/MhV_v_transcriptome/viral_genes_recip_DaFV_LbFV.csv")

##compare to recip blast results for all viral hits
##overlap
recipv_prodigal_overlap <- intersect(min_evalues_all_annots$Trinity_ID, recip_viral_genes$transcript_id)
##merge with prodigal blast res
prodigal_recipv_annots <- merge(min_evalues_all_annots, recip_viral_genes, by.x="Trinity_ID", by.y="transcript_id", all.x=TRUE)
prodigal_recipv_annots_simple <- prodigal_recipv_annots[,c(1,2,3,42,4,32,16,20,21,24)]
fwrite(prodigal_recipv_annots_simple, "output/blast/MhV_v_transcriptome/viral_genes_recip_all_viral.csv")


