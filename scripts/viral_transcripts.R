library(data.table)
library(ggplot2)
library(dplyr)

trinotate.min.eval <- fread('output/trinotate/sorted/longest_isoform_annots.csv', na.strings = ".")
##viral hits
virus_x <- data.table(dplyr::filter(trinotate.min.eval, grepl('Viruses', sprot_Top_BLASTX_hit)))
virus_p <- data.table(dplyr::filter(trinotate.min.eval, grepl('Viruses', sprot_Top_BLASTP_hit)))
virus_annots <- full_join(virus_x, virus_p)
virus_annots <- data.table(dplyr::filter(virus_annots, !grepl('transposon', sprot_Top_BLASTX_hit)))
fwrite(virus_annots, "output/trinotate/viral/genes_viral_annots.csv")



##viral taxa
virus$sprot_Top_BLASTX_hit <- tstrsplit(virus$sprot_Top_BLASTX_hit, "`", keep=1, fixed=TRUE)
split.first.blastx <- virus[,tstrsplit(sprot_Top_BLASTX_hit, "^", fixed=TRUE), by=`#gene_id`]
genes.per.viral.taxa <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
##substitute to keep last string which should be genus
#start of string, 0 or more of any character, 
genes.per.viral.taxa$viral_genera <- data.table(gsub("^.*; ", "", genes.per.viral.taxa$V7))
##have some cases where last string wasn't genus - fix these
genes.per.viral.taxa$viral_genera <- sub("Murine leukemia virus", "Gammaretrovirus", genes.per.viral.taxa$viral_genera)
##would be good to also include whether DNA or RNA virus - search manually
fwrite(genes.per.viral.taxa, "output/trinotate/viral/genes_per_taxa_viral.csv")

##plot viral taxa annots
plot.viral.taxa <- fread("output/trinotate/viral/genes_per_taxa_viral_edited_for_plot.csv")

ggplot(plot.viral.taxa, aes(x=reorder(V7, -V1), y=V1))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  geom_col()+xlab("Viral Genera")+ylab("Number of BlastX Annotations")

##write list of viral transcript IDs to pull out of fasta file
viral_transcripts <- trinotate_virus_annots[,2]
fwrite(list(viral_transcripts), "output/trinotate/viral/viral_transcript_ids.txt")
##^^wont include manual hits e.g. bro genes

##virus transposon transcripts
transposon <- dplyr::filter(trinotate_virus_annots, grepl('transposon', sprot_Top_BLASTX_hit))