library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

##mh results
mh_nr_blastx <- fread("output/recip_blast/nr_blastx/nr_blastx.outfmt3")
virus_taxids <- fread("data/taxids/species_virus_taxids.txt", header=FALSE)

setnames(mh_nr_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("transcript_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(mh_nr_blastx, transcript_id, evalue, -bit_score)

##am I losing viral hits in this filtering
##TRINITY_DN107665_c0_g1_i1 doesn't have a viral hit in minevalue but blastx is DaFV
##losing to evalues above blast threshold

##extract result with lowest evalue for each peptide (use bit-score if tie)
mh_min_evalues <- mh_nr_blastx[,.SD[which.min(evalue)], by=transcript_id]
##filter out virus annotations using taxids
viral_hits <- subset(mh_min_evalues, taxid %in% virus_taxids$V1)
fwrite(viral_hits, "output/recip_blast/nr_blastx/viral_annots.csv")

LbFV_annots <- dplyr::filter(mh_virus, grepl('Leptopilina boulardi filamentous virus', annotation))
DaFV_annots <- dplyr::filter(mh_virus, grepl('Drosophila-associated filamentous virus', annotation))
LbFV_DaFV <- full_join(LbFV_annots, DaFV_annots)
fwrite(LbFV_DaFV, "output/recip_blast/nr_blastx/LbFV_DaFV_annots.csv")

##ID all LBFV hits (even if better other hit) - there are 4 that have a hit to something else
Mh_allhits_lbfv <- dplyr::filter(mh_nr_blastx, grepl('Leptopilina boulardi filamentous virus', annotation))
Mh_allhits_lbfv$split_annotation <- tstrsplit(Mh_allhits_lbfv$annotation, "<>", keep=c(1))
Mh_allhits_lbfv_unique <- (unique(Mh_allhits_lbfv$split_annotation))
best_not_lbfv <- setdiff(Mh_allhits_lbfv$transcript_id, LbFV_annots$transcript_id)
best_not_lbfv_annots <- subset(viral_hits, transcript_id %in% best_not_lbfv)
##2 are to DaFV, two are to bavculo genes (bro virus, IAP parasitoid)

##edited to add viral families and DNA/RNA genome
viral_hits_plot <- fread("output/recip_blast/nr_blastx/viral_annots_plot.csv")
viral_hits_plot$`#gene_id` <- tstrsplit(viral_hits_plot$transcript_id, "_i", keep=c(1))

##merge with viral annots from trinotate
virus_annots_trinotate <- fread("output/trinotate/prokaryotic/viral_genes_plot.csv")
all_virus <- merge(viral_hits_plot, virus_annots_trinotate, by="#gene_id", all.x=TRUE, all.y=TRUE)


##plot recip viral families
##DNA
viral_hits_plot_DNA <- subset(viral_hits_plot, genome=="DNA")
DNA_genes_viral_family <- viral_hits_plot_DNA[,length(unique(`#gene_id`)), by=family]
##plot
ggplot(DNA_genes_viral_family, aes(x=reorder(family, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 21))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of Trinity genes")

##RNA
viral_hits_plot_RNA <- subset(viral_hits_plot, genome=="RNA")
RNA_genes_viral_family <- viral_hits_plot_RNA[,length(unique(`#gene_id`)), by=family]
##plot
ggplot(RNA_genes_viral_family, aes(x=reorder(family, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 21))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of Trinity genes")
