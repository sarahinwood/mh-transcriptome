library(data.table)
library(dplyr)

##mh results
mh_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")

setnames(mh_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mh_nr_blastx, transcript_id, evalue, -bit_score)
##ID all LBFV hits (even if better other hit)
Mh_allhits_lbfv <- dplyr::filter(mh_nr_blastx, grepl('Leptopilina boulardi filamentous virus', annotation))
Mh_allhits_lbfv$split_annotation <- tstrsplit(Mh_allhits_lbfv$annotation, "<>", keep=c(1))
Mh_allhits_lbfv_unique <- (unique(Mh_allhits_lbfv$split_annotation))

##am I losing viral hits in this filtering
##TRINITY_DN107665_c0_g1_i1 doesn't have a viral hit in minevalue but blastx is DaFV
##losing to evalues above blast threshold

##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
mh_min_evalues <- mh_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations
mh_virus <- dplyr::filter(mh_min_evalues, grepl('virus', annotation))
mh_virus <- dplyr::filter(mh_virus, !grepl('transposon', annotation))
fwrite(mh_virus, "output/recip_blast/nr_blastx/viral_annots.csv")

LbFV_annots <- dplyr::filter(mh_virus, grepl('Leptopilina boulardi filamentous virus', annotation))
DaFV_annots <- dplyr::filter(mh_virus, grepl('Drosophila-associated filamentous virus', annotation))
LbFV_DaFV <- full_join(LbFV_annots, DaFV_annots)
fwrite(LbFV_DaFV, "output/recip_blast/nr_blastx/LbFV_DaFV_annots.csv")

##merge with viral annots from trinotate
virus_annots <- fread("output/trinotate/viral/genes_viral_annots.csv")
all_virus <- merge(virus_annots, mh_virus, by="transcript_id", all.x=TRUE, all.y=TRUE)
##list of all viral genes
list_viral_transcripts <- all_virus[,c(1)]
##get trinotate info for all viral genes
trinotate.min.eval <- fread('output/trinotate/sorted/longest_isoform_annots.csv', na.strings = ".")
all_virus_trinotate <- merge(list_viral_transcripts, trinotate.min.eval, by="transcript_id", all.x=TRUE)
##merge with blast annots
all_viral_annots <- merge(all_virus_trinotate, mh_virus, by="transcript_id", all.x=TRUE)
fwrite(all_viral_annots, "output/trinotate/all_viral_annots.csv")

transcript_id_annot <- mh_virus[,c(1,13)]
fwrite(transcript_id_annot, "output/recip_blast/nr_blastx/gene_vs_annot.csv")
##search uniprot taxonomy for genus
viral_genera <- fread("output/recip_blast/nr_blastx/genus_gene_vs_annot.csv")
viral_genes.per.genus <- viral_genera[,length(unique(transcript_id)), by=Genus]
fwrite(viral_genes.per.genus, "output/recip_blast/nr_blastx/genes_per_viral_genera.csv")
