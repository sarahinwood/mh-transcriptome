library("data.table")
library("dplyr")
library("ggplot2")

##ASW results
asw_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")
setnames(asw_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(asw_nr_blastx, transcript_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
asw_min_evalues <- asw_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations
asw_virus <- dplyr::filter(asw_min_evalues, grepl('virus', annotation))
fwrite(asw_virus, "output/recip_blast/nr_blastx/viral_annots.csv")
transcript_id_annot <- asw_virus[,c(1,13)]
fwrite(transcript_id_annot, "output/recip_blast/nr_blastx/viral_transcripts_annots.csv")

##sort of virus annots from trinotate WITHOUT transposon hits
sorted_trinotate <- fread('output/trinotate/viral/viral_annots.csv')
##sum of each viral taxa
blastx.results <- sorted_trinotate[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
genes.per.taxa <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
virus_taxa <- dplyr::filter(genes.per.taxa, grepl('Viruses', V7))
fwrite(virus_taxa, "output/trinotate/viral/sorted_trinotate_viral_taxa.csv")
