#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)

###########
# GLOBALS #
###########

trinotate_file <- snakemake@input[["trinotate_report"]]
longest_isoform_id_list <- snakemake@input[["longest_isoform_ids"]]

trinotate.report <- fread(trinotate_file, na.strings = ".")

##split to keep first blastx hit only
trinotate.sorted <- copy(trinotate.report)
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.report$sprot_Top_BLASTX_hit, "`", fixed=TRUE, keep=c(1))
##split to keep only evalue
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.sorted$blastx_evalue, "^RecName", fixed=TRUE, keep=c(1))
trinotate.sorted$blastx_evalue <- tstrsplit(trinotate.sorted$blastx_evalue, "E:", fixed=TRUE, keep=c(2))
##reorder to get evalue beside blastx annotation
trinotate.sorted <- trinotate.sorted[,c(1,2,3,18,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
##convert evalue column to numeric
trinotate.sorted[,4] <- sapply(trinotate.sorted[,4], as.numeric)
##set order so for each gene the transcript with the lowest evalue hit is at top
setorder(trinotate.sorted, `#gene_id`, `blastx_evalue`, na.last=TRUE)
####extract result with lowest evalue for each gene - what if multiple rows with lowest min?
trinotate.min.eval <- trinotate.sorted[,.SD[which.min(blastx_evalue)], by=`#gene_id`]
##write csv with most sig hit for each gene with annotation
fwrite(trinotate.min.eval, snakemake@output[["best_annot_per_gene"]])

##split to keep annot for longest isoform:
longest_isoform_ids <- fread(longest_isoform_id_list, header=FALSE)
annots_longest_isoform <- merge(trinotate.report, longest_isoform_ids, by.x="transcript_id", by.y="V1", all.x=FALSE, all.y=TRUE)
fwrite(annots_longest_isoform, snakemake@output[["longest_iso_annots"]])

##filter out gene ids for unannotated OR viral annot genes
genes_no_annot <- trinotate.report[is.na(sprot_Top_BLASTX_hit),]
##list of unique gene ids from table of genes with no blastx annot
id_dt_no_annot <- data.table(unique(genes_no_annot$`#gene_id`))
##filter out gene ids for viral annot genes
virus_x <- data.table(dplyr::filter(trinotate.report, grepl('Viruses', sprot_Top_BLASTX_hit)))
virus_p <- data.table(dplyr::filter(trinotate.report, grepl('Viruses', sprot_Top_BLASTP_hit)))
virus_annots_table <- full_join(virus_x, virus_p)
id_dt_viral_annots <- data.table(unique(virus_annots_table$`#gene_id`))
id_dt_viral_or_unann <- full_join(id_dt_no_annot, id_dt_viral_annots)
fwrite(id_dt_viral_or_unann, snakemake@output[["viral_or_unann_transcript_ids"]])
