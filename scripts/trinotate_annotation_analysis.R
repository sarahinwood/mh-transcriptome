library(data.table)
library(dplyr)
library(VennDiagram)
library(stringr)
library(ggplot2)
library(viridis)
library(gridExtra)

#Get annotation report in right format for venn diagram
annotation.report <- fread('output/trinotate/trinotate/trinotate_annotation_report.txt', na.strings = ".")
pfam <- annotation.report[!is.na(Pfam), unique(`#gene_id`)]
blastx <- annotation.report[!is.na(sprot_Top_BLASTX_hit), unique(`#gene_id`)]
kegg <- annotation.report[!is.na(Kegg), unique(`#gene_id`)]
pfam_go <- annotation.report[!is.na(gene_ontology_Pfam), unique(`#gene_id`)]
bx_go <- annotation.report[!is.na(gene_ontology_BLASTX), unique(`#gene_id`)]
bp_go <- annotation.report[!is.na(gene_ontology_BLASTP), unique(`#gene_id`)]

number.genes <- annotation.report[!is.na(`#gene_id`),length(unique(`#gene_id`))]

##overlap of blast hits for transdecoder predictions
trans_annots <- subset(annotation.report, !is.na(prot_id))
trans_annots <- trans_annots[,c(1,3,7)]
trans_blast <- trans_annots[!is.na(trans_annots$sprot_Top_BLASTX_hit) & !is.na(trans_annots$sprot_Top_BLASTP_hit),]
trans_blast$blastx <- tstrsplit(trans_blast$sprot_Top_BLASTX_hit, "^", keep=c(1), fixed=TRUE)
trans_blast$blastp <- tstrsplit(trans_blast$sprot_Top_BLASTP_hit, "^", keep=c(1), fixed=TRUE)
sum(trans_blast$blastx==trans_blast$blastp)

#Draw Venn Diagram
vd <- venn.diagram(x = list("Pfam"=pfam, "BlastX"=blastx, "Kegg"=kegg), filename=NULL,
                   fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                   main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd)

#Sum of genes with any annotation
long.annotationreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "Pfam", "Kegg"))
any.annotations <- long.annotationreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
any.annotations[,length(unique(`#gene_id`))]
any.annotations[,sum(any_annotations)]

#Sum of genes with blast annotation
long.blastreport <- melt(annotation.report,id.vars = "#gene_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit"))
any.blast <- long.blastreport[,.(any_annotations = any(!is.na(value))),by=`#gene_id`]
any.blast[,length(unique(`#gene_id`))]
any.blast[,sum(any_annotations)]
any_blast_ids <- subset(any.blast, any.blast$any_annotations==TRUE)

#Sum of genes with any Transdecoder predicted protein
transdecoder.report <- melt(annotation.report, id.vars="#gene_id", measure.vars = c("prot_id"))
any.transdecoder <- transdecoder.report[,.(transdecoder_annotation = any(!is.na(value))),by=`#gene_id`]
any.transdecoder[,length(unique(`#gene_id`))]
any.transdecoder[,sum(transdecoder_annotation)]
sum(any.transdecoder$transdecoder_annotation==FALSE)

any_annot_ids <- subset(any.annotations, any.annotations$any_annotations==TRUE)
any_transdecoder_ids <- subset(any.transdecoder, any.transdecoder$transdecoder_annotation==TRUE)

vd3 <- venn.diagram(x = list("Transdecoder"=any_transdecoder_ids$`#gene_id`, "Blast"=any_blast_ids$`#gene_id`), filename=NULL,
                    fill=c("#440154FF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5, height=3, width=6,
                    main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd3)
ggsave(file="output/trinotate/trinotate/transdecoder_blast_venn.svg", plot=vd3, width=6, height=3)

##transdecoder vs blastX
transdecoder <- annotation.report[!is.na(prot_id), unique(`#gene_id`)]
vd4 <- venn.diagram(x = list("Transdecoder"=transdecoder, "BlastX"=blastx), filename=NULL,
                    fill=c("#440154FF", "#FDE725FF"), alpha=0.7, cex = 1, cat.cex=1, lwd=1.5,
                    main=paste("Total Number of Genes = ", number.genes))
grid.newpage()
grid.draw(vd4)

##split blastX results
blastx.results <- annotation.report[!is.na(sprot_Top_BLASTX_hit),.(sprot_Top_BLASTX_hit, `#gene_id`)]
first.blastx.hit <- blastx.results[,tstrsplit(sprot_Top_BLASTX_hit, "`", fixed = TRUE, keep=1), by = `#gene_id`]
split.first.blastx <- first.blastx.hit[,tstrsplit(V1, "^", fixed=TRUE), by=`#gene_id`]
split.first.blastx$e_value <- tstrsplit(split.first.blastx$V5, "E:", keep=c(2))
split.first.blastx$pid <- tstrsplit(split.first.blastx$V4, "%ID", keep=c(1))
gene_pid_eval <- split.first.blastx[!is.na(e_value),.(`#gene_id`, pid, e_value)]
gene_pid_eval$pid <- as.numeric(gene_pid_eval$pid)
gene_pid_eval$e_value <- as.numeric(gene_pid_eval$e_value)

mean(gene_pid_eval$pid)

#pid spread
ggplot(gene_pid_eval, aes(x=pid))+
  geom_histogram(alpha=0.8, fill="#440154FF", colour="#440154FF", binwidth=2, boundary=0)+
  scale_x_continuous(limits=c(0,100))+
  xlab("BlastX sequence identity (%)")+
  ylab("Number of BlastX annotations")+
  theme_bw()

#Annotations per genus
genes.per.genus <- split.first.blastx[,length(unique(`#gene_id`)), by=V7]
setkey(genes.per.genus, V1)
print(genes.per.genus)
fwrite(genes.per.genus, "output/trinotate/genes_per_genus.csv")

#meanwhile in excel sort for genus with most annotations, delete those with low no.
#I'm not interested in, and alter genus name to just genera

#plot annotations per genus - top 25 genera
plot.genes.per.genus <- fread("output/trinotate/genes_per_genus_plot.csv", strip.white=FALSE)
ggplot(plot.genes.per.genus, aes(x=reorder(V7, -V1), y=V1))+
  geom_col(alpha=0.8, fill="#440154FF", colour="#440154FF", width=0.9)+
  ylim(c(0, 7000))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 65, hjust = 1, face = "italic")) +
  xlab("Genus")+ylab("Number of BlastX Annotations")
