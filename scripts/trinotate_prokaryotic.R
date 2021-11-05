library(data.table)
library(dplyr)
library(ggplot2)

#Get annotation report in right format for venn diagram
annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")
##split blastX results
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
split.first.blastx$e_value <- tstrsplit(split.first.blastx$V5, "E:", keep=c(2))

##split to keep only prokaryotic hits - this is still with a line for every transcript
prokaryotic <- dplyr::filter(split.first.blastx, !grepl('Eukaryota', V7, fixed=TRUE))
prokaryotic$domain <- tstrsplit(prokaryotic$V7, "; ", keep=c(1))
##658 genes, with 837 transcripts have prokaryotic annots
length(unique(prokaryotic$`#gene_id`))

############
# bacteria #
############
#376 genes, 462 transcripts, from 96 different species
bacteria <- dplyr::filter(split.first.blastx, grepl('Bacteria', V7))
length(unique(bacteria$`#gene_id`))
length(unique(bacteria$V7))
##phyla
bacteria$phyla <- tstrsplit(bacteria$V7, "; ", keep=c(2))
fwrite(bacteria, "output/trinotate/prokaryotic/bacterial_genes.csv")

genes_bacterial_phyla <- bacteria[,length(unique(`#gene_id`)), by=phyla]
##plot
ggplot(genes_bacterial_phyla, aes(x=reorder(phyla, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 160))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Bacterial Phyla")+ylab("Number of transcripts")

###########
# viruses #
###########
#264 genes, 340 transcripts, from 41 different species
viruses <- dplyr::filter(split.first.blastx, grepl('Viruses', V7))
length(unique(viruses$`#gene_id`))
length(unique(viruses$V7))
## save csv and manually split to families as tstrsplit getting different classification levels
fwrite(viruses, "output/trinotate/prokaryotic/viral_genes.csv")

viruses_plot <- fread("output/trinotate/prokaryotic/viral_genes_plot.csv")
genes_viral_family <- viruses_plot[,length(unique(`#gene_id`)), by=family]
##plot
ggplot(genes_viral_family, aes(x=reorder(family, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 160))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of transcripts")

###########
# archaea #
###########
#20 genes, 35 transcripts, from 14 different species
archaea <- dplyr::filter(split.first.blastx, grepl('Archaea', V7))
length(unique(archaea$`#gene_id`))
length(unique(archaea$V7))

archaea$phyla <- tstrsplit(archaea$V7, "; ", keep=c(2))
genes_archaeal_phyla <- archaea[,length(unique(`#gene_id`)), by=phyla]
##plot
ggplot(genes_archaeal_phyla, aes(x=reorder(phyla, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Archaeal Phyla")+ylab("Number of transcripts")
