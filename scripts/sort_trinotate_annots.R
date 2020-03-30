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

###########
# GLOBALS #
###########

trinotate_file = <- snakemake@input[["trinotate_report"]]

trinotate.report <- fread(trinotate_file)
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

##filter out gene ids for unannotated genes
genes_no_annot <- trinotate.report[is.na(sprot_Top_BLASTX_hit),]
##list of unique gene ids from table of genes with no blastx annot
list_ids_no_annot <- list(unique(genes_no_annot$`#gene_id`))
fwrite(list_ids_no_annot, snakemake@output[["unann_transcript_ids"]])
