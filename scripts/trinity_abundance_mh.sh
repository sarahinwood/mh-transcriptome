#!/usr/bin/env bash

export LD_LIBRARY_PATH=/usr/lib64

 bin/trinity/util/align_and_estimate_abundance.pl \
 	--transcripts output/trinity/Trinity.fasta \
 	--seqType fq \
 	--est_method RSEM \
 	--output_dir output/trinity_abundance \
 	--aln_method bowtie \
 	--prep_reference \
 	--SS_lib_type RF \
 	--thread_count 1 \
 	--trinity_mode \
 	--left output/bbduk_trim/abdo_r1.fq.gz,output/bbduk_trim/head_r1.fq.gz,output/bbduk_trim/sting_r1.fq.gz \
 	--right output/bbduk_trim/abdo_r2.fq.gz,output/bbduk_trim/head_r2.fq.gz,output/bbduk_trim/sting_r1.fq.gz