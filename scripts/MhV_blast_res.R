library(data.table)
library(dplyr)

viral_genes_blast_res <- fread("output/MhV_blast/MhV_v_transcriptome/blastn.outfmt6")
prodigal_annots <- fread("data/Mh_prodigal/blastp_gff.csv")
id_annot <- prodigal_annots[,c(3,12)]
setnames(id_annot, old=c("annotation", "prodigal_nt_id"), new=c("Prodigal_blast_annot", "Prodigal_ID"))
trinotate <- fread("output/trinotate/sorted/longest_isoform_annots.csv", na.strings = ".")

setnames(viral_genes_blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("Prodigal_ID", "Trinity_ID", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(viral_genes_blast_res, Prodigal_ID, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- viral_genes_blast_res[,.SD[which.min(evalue)], by=Prodigal_ID]

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
recip_viral_genes <- fread("data/mh-transcriptome/output/recip_blast/nr_blastx/viral_annots.csv")
##overlap
recipv_prodigal_overlap <- intersect(min_evalues_all_annots$Trinity_ID, recip_viral_genes$transcript_id)
##merge with prodigal blast res
prodigal_recipv_annots <- merge(min_evalues_all_annots, recip_viral_genes, by.x="Trinity_ID", by.y="transcript_id", all.x=TRUE)
prodigal_recipv_annots_simple <- prodigal_recipv_annots[,c(1,2,3,42,4,32,16,20,21,24)]
fwrite(prodigal_recipv_annots_simple, "output/blast/MhV_v_transcriptome/viral_genes_recip_all_viral.csv")


