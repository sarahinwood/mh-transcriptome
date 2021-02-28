library(data.table)
library(dplyr)
##recip blast viral genes
recip_transcript_id_annot <- fread("output/recip_blast/nr_blastx/gene_vs_annot.csv")
recip_transcript_id_annot$`#gene_id` <- tstrsplit(recip_transcript_id_annot$transcript_id, "_i", keep=c(1))
recip_gene_annot <- recip_transcript_id_annot[,c(3,2)]
recip_gene_annot$annotation <- tstrsplit(recip_gene_annot$annotation, "<>", keep=c(1))
##trinotate viral genes
trinotate_virus <- fread("output/trinotate/viral/genes_viral_annots.csv")
trinotate_virus_id_annot <- trinotate_virus[,c(2,3)]

##all viral genes_v_annots
all_viral <- full_join(recip_gene_annot, trinotate_virus_id_annot)
##expression matrix - filter for viral transcripts
expression_matrix <- fread("output/sample_trinity_abundance/abundance_matrix/all_salmon.gene.TMM.EXPR.matrix")
viral_expression_matrix <- merge(all_viral, expression_matrix, by.x="#gene_id", by.y="V1", all.x=TRUE)
fwrite(viral_expression_matrix, "output/sample_trinity_abundance/viral_expression_matrix.csv")

##transcripts that had hits to viral contigs in genome
transcripts_genome_viral_contigs <- fread("data/blastn_viral_contigs_transcriptome_annots.csv")
transcripts_genome_viral_contigs$contig_id <- tstrsplit(transcripts_genome_viral_contigs$prodigal_nt_id, "_", keep=c(1))
##remove hits to contigs I don't think are viral
transcripts_genome_only_viral_contigs <- data.table(dplyr::filter(transcripts_genome_viral_contigs, !grepl('Scaffold80|Scaffold8624|Scaffold3939|Scaffold28315', contig_id)))
genome_viral_transcripts <- transcripts_genome_only_viral_contigs[,c(4,2,5)]
all_viral_genome_viral <- full_join(all_viral, genome_viral_transcripts)


