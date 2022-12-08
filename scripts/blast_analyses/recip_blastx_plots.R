library(ggplot2)
library(viridis)
library(data.table)

###############################
## plot recip viral families ##
###############################
##edited to add viral families and DNA/RNA genome
##annotated for plot
viral_hits_plot <- fread("output/recip_blast/viral_nr_blastx/best_viral_hits_plot.csv")

##DNA
viral_hits_plot_DNA <- subset(viral_hits_plot, genome=="DNA")
DNA_genes_viral_family <- viral_hits_plot_DNA[,length(unique(`transcript_id`)), by=family]
##plot
ggplot(DNA_genes_viral_family, aes(x=reorder(family, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 21))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of Trinity genes")

##RNA
viral_hits_plot_RNA <- subset(viral_hits_plot, genome=="RNA")
RNA_genes_viral_family <- viral_hits_plot_RNA[,length(unique(`transcript_id`)), by=family]
##plot
ggplot(RNA_genes_viral_family, aes(x=reorder(family, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  theme_bw()+
  ylim(c(0, 21))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of Trinity genes")

##################
## viral family ##
##################

##counts per family
genes_viral_family <- viral_hits_plot[,length(unique(transcript_id)), by=family]
pid_family <- viral_hits_plot[,mean(`%_identical_matches`), by=family]
family_genes_pid <- merge(genes_viral_family, pid_family, by="family")

##merge with genome type
family_genome_type <- distinct(viral_hits_plot[,c(18,19,20)])
plot_data <- merge(family_genes_pid, family_genome_type, all.x=TRUE, all.y=FALSE)
##round %id to 2sf
plot_data$V1.y <- format(round(plot_data$V1.y, 2), nsmall = 2)
plot_data$pid <- paste(plot_data$V1.y, "%", sep="")
fwrite(plot_data, "output/recip_blast/viral_nr_blastx/viral_annots_plot_data.csv")

##plot
ggplot(plot_data, aes(x=reorder(family, -V1.x), y=V1.x, fill=genome_type))+
  geom_col(alpha=0.8, width=0.9)+
  theme_bw()+
  scale_fill_viridis(discrete=T, name="Genome type")+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic"))+
  xlab("Viral Families")+ylab("Number of Trinity genes")

#################
## genome type ##
#################

##order factors for nicer plot
plot_data$family <- factor(plot_data$family, levels=c("Ascoviridae", "Marseilleviridae", "Mimiviridae", "Podoviridae", "Siphoviridae", "Adintoviridae", "Nudiviridae", "Unclassified DNA virus", "Baculoviridae", "Polydnaviridae", # DNA
                                                      "Narnaviridae", "Nodaviridae", "Luteoviridae", "Picornaviridae", "Solemoviridae", "Iflaviridae", #(+)ssRNA
                                                      "Phenuiviridae", "Rhabdoviridae", #(-)ssRNA
                                                      "Reoviridae", "Totiviridae", #dsRNA
                                                      "Unclassified RNA virus", "Unclassified virus"))
plot_data$genome_type <- factor(plot_data$genome_type, levels=c("dsDNA", "(+)ssRNA", "(-)ssRNA", "dsRNA", "Unclassified RNA", "Unclassified virus"))
##plot - fill for family adds too much info
ggplot(plot_data, aes(x=genome_type, y=V1.x, fill=family))+
  geom_bar(width=0.9, stat="identity")+
  theme_bw()+
  theme(legend.text = element_text(face = c("italic")))+
  scale_fill_viridis(discrete=T, name="Family")+
  xlab("Genome type")+ylab("Number of Trinity genes")+
  scale_x_discrete(labels=c("Unclassified RNA"="Unclassified\nRNA", "Unclassified virus"="Unclassified\nvirus"))

# with %id on bars
ggplot(plot_data, aes(x=genome_type, y=V1.x, fill=family))+
  geom_bar(width=0.9, stat="identity")+
  geom_text(aes(label=pid), position = position_stack(vjust = 0.5), size=2)+
  theme_bw()+
  theme(legend.text = element_text(face = c("italic")))+
  scale_fill_viridis(discrete=T, name="Family")+
  xlab("Genome type")+ylab("Number of Trinity genes")+
  scale_x_discrete(labels=c("Unclassified RNA"="Unclassified\nRNA", "Unclassified virus"="Unclassified\nvirus"))


########################
## certain viral hits ##
########################

##  LbFV
LbFV_annots <- dplyr::filter(viral_hits, grepl('Leptopilina boulardi filamentous virus', annotation))
DaFV_annots <- dplyr::filter(viral_hits, grepl('Drosophila-associated filamentous virus', annotation))
LbFV_DaFV <- full_join(LbFV_annots, DaFV_annots)
fwrite(LbFV_DaFV, "output/recip_blast/nr_blastx/LbFV_DaFV_annots.csv")
##ID all LBFV hits (even if better other hit) - there are 4 that have a hit to something else
Mh_allhits_lbfv <- dplyr::filter(mh_nr_blastx, grepl('Leptopilina boulardi filamentous virus', annotation))
Mh_allhits_lbfv$split_annotation <- tstrsplit(Mh_allhits_lbfv$annotation, "<>", keep=c(1))
Mh_allhits_lbfv_unique <- (unique(Mh_allhits_lbfv$split_annotation))
best_not_lbfv <- setdiff(Mh_allhits_lbfv$transcript_id, LbFV_annots$transcript_id)
best_not_lbfv_annots <- subset(viral_hits, transcript_id %in% best_not_lbfv)
##2 are to DaFV, two are to bavculo genes (bro virus, IAP parasitoid)


## PDV
Mh_PDV <- dplyr::filter(viral_hits_plot, grepl('Polydnaviridae', family))
PDV_gene_hits <- subset(mh_nr_blastx, transcript_id %in% Mh_PDV$transcript_id)
fwrite(PDV_gene_hits, "output/recip_blast/nr_blastx/PDV_all_annots.csv")
notPDV <- subset(PDV_gene_hits, !grepl("ichnovirus|bracovirus", PDV_gene_hits$annotation, ignore.case=T))
##of those hits that are not ichno or braco for these genes, they only have two other virus hits, rest are
###2022 uncharacterised genes, 3177 hypothetical, many from parasitoids and other insects
notPDVunchar <- subset(notPDV, !grepl("uncharacter|hypothetical|GSCOC|unnamed", notPDV$annotation, ignore.case=T))
length(unique(notPDVunchar$transcript_id))
##best hit not polydna/hypo
setorder(notPDVunchar, transcript_id, evalue, -bit_score)
not_PDV_best_hits <- notPDVunchar[,.SD[which.min(evalue)], by=transcript_id]

