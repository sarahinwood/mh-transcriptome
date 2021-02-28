library("data.table")
library("dplyr")

##mh results
mh_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")
setnames(mh_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mh_nr_blastx, transcript_id, evalue, -bit_score)
##ID all LBFV hits (even if better other hit)
Mh_allhits_lbfv <- dplyr::filter(mh_nr_blastx, grepl('Leptopilina boulardi filamentous virus', annotation))
Mh_allhits_lbfv$split_annotation <- tstrsplit(Mh_allhits_lbfv$annotation, "<>", keep=c(1))
Mh_allhits_lbfv_unique <- (unique(Mh_allhits_lbfv$split_annotation))


##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
mh_min_evalues <- mh_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations
mh_virus <- dplyr::filter(mh_min_evalues, grepl('virus', annotation))
mh_virus <- dplyr::filter(mh_virus, !grepl('transposon', annotation))
fwrite(mh_virus, "output/recip_blast/nr_blastx/viral_annots.csv")

transcript_id_annot <- mh_virus[,c(1,13)]
fwrite(transcript_id_annot, "output/recip_blast/nr_blastx/gene_vs_annot.csv")
##search uniprot taxonomy for genus
viral_genera <- fread("output/recip_blast/nr_blastx/genus_gene_vs_annot.csv")
viral_genes.per.genus <- viral_genera[,length(unique(transcript_id)), by=Genus]
fwrite(viral_genes.per.genus, "output/recip_blast/nr_blastx/genes_per_viral_genera.csv")
