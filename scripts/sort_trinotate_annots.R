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
signalp_file <- snakemake@input[["signalp"]]
longest_isoform_id_list <- snakemake@input[["longest_isoform_ids"]]

trinotate.report <- fread(trinotate_file, na.strings = ".")
signalp_out <- fread(signalp_file, fill=TRUE, skip=2, col.names=c(
  "name", "Cmax", "pos", "Ymax", "pos", "Smax", "pos", "Smean", "D", "?", "Dmaxcut", "Networks-used"))
longest_isoform_ids <- fread(longest_isoform_id_list, header=FALSE)

##split to keep annot for longest isoform:
annots_longest_isoform <- merge(trinotate.report, longest_isoform_ids, by.x="transcript_id", by.y="V1", all.x=FALSE, all.y=TRUE)

##add in signalp res
signalp_y <- subset(signalp_out, signalp_out$`?`=="Y")
signalp_y$transcript_id <- tstrsplit(signalp_y$name, "::", fixed=TRUE, keep=c(2))
annots_longest_isoform$SignalP <- ifelse(annots_longest_isoform$transcript_id %in% signalp_y$transcript_id, "Y", "N")

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
