library(data.table)
library(ggplot2)

##plot viral genera annots
trinotate_viral_genera <- fread("output/trinotate/viral/genes_per_taxa_viral_genome.csv")
trinotate_viral_genera$annotation_source <- paste("Trinotate")
recip_viral_genera <- fread("output/recip_blast/nr_blastx/genes_per_viral_genera_genome.csv")
recip_viral_genera$annotation_source <- paste("Reciprocal BlastX")

virus_genera <- data.table(full_join(trinotate_viral_genera, recip_viral_genera))
virus_vs_genome <- unique(virus_genera[,c(3,4)])

##plot genera
genera_vs_counts <- virus_genera[,sum((V1)), by=(viral_genera, Genome)]
plot_data <- (merge(genera_vs_counts, virus_vs_genome, all.x=TRUE))
ggplot(plot_data, aes(x=reorder(viral_genera, -V1), y=V1))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col(aes(fill=Genome))+xlab("Viral Genera")+ylab("Number of Annotations")

##plot genome type
genome_vs_counts <- virus_genera[,sum((V1)), by=Genome]
ggplot(genome_vs_counts, aes(x=reorder(Genome, -V1), y=V1))+
  theme_light()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Genome type")+ylab("Number of Annotations")
