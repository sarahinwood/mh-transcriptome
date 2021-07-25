#!/usr/bin/env python3
import pathlib2
import peppy

#############
# FUNCTIONS #
#############

##needed to get BUSCO running in new folder
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def get_reads(wildcards):
    input_keys = ['l1r1', 'l2r1', 'l1r2', 'l2r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

###########
# GLOBALS #
###########

bbduk_adapters = '/adapters.fa'

#containers from Tom
bbduk_container = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.3/container-bbmap.bbmap_38.90.sif'
trinotate_container = 'shub://TomHarrop/trinotate_pipeline:v0.0.12'
##program-provided docker containers - latest versions as of 2021-07-09
kraken_container = 'docker://staphb/kraken2:2.1.2-no-db'
trinity_container = 'docker://trinityrnaseq/trinityrnaseq:2.12.0'
busco_container = 'docker://ezlabgva/busco:v5.2.1_cv1'
tidyverse_container = 'docker://rocker/tidyverse:4.1.0'
blast_container= 'docker://ncbi/blast:2.12.0'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'

##old containers
old_busco_container = 'docker://ezlabgva/busco:v4.0.2_cv1'
tom_tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'

#########
# RULES #
#########

rule target:
    input:
        expand('output/busco/{filter}/short_summary.specific.hymenoptera_odb10.{filter}.txt',
                filter=['expression', 'length']),
        #'output/busco/short_summaries/busco_figure.png',
        expand('output/fastqc/{sample}_r{n}_fastqc.zip', sample=all_samples, n=[1,2]),
        'output/fastqc/Mh_venom3_r2_fastqc.zip',
        'output/trinity_stats/stats.txt',
        'output/trinity_stats/xn50.out.txt',
        'output/trinity_stats/bowtie2_alignment_stats.txt',
        expand('output/trinity_stats/isoforms_by_{filter}_bowtie2_alignment_stats.txt',
                filter=['expression', 'length']),
        'output/trinotate/sorted/longest_isoform_annots.csv',
        'output/alt_recip_blast/nr_blastx/nr_blastx.outfmt3',
        expand('output/kraken/kraken_{sample}_out.txt', sample=all_samples),
        'output/kraken/kraken_Mh_venom3_out.txt'


################################################################
##Reciprocal blastx searching for viral annots for unann genes##
################################################################

##need to update using taxid instead
##trial with BLOSUM90 to see whether that increases or decreases hits after having lower e-value threshold
rule recip_nr_blastx:
    input:
        pot_viral_transcripts = 'output/recip_blast/viral_blastx/potential_viral_transcripts.fasta'
    output:
        blastx_res = 'output/alt_recip_blast/nr_blastx/nr_blastx.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/recip_nr_blastx.log'
    shell:
        'blastx '
        '-query {input.pot_viral_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 0.5 '
        '-outfmt "6 std staxids salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_pot_viral_transcripts:
    input:
        length_filtered_transcriptome = 'output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_hit_ids = 'output/recip_blast/viral_blastx/transcripts_viral_hit_ids.txt'
    output:
        pot_viral_transcripts = 'output/recip_blast/viral_blastx/potential_viral_transcripts.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_pot_viral_transcripts.log'
    shell:
        'filterbyname.sh '
        'in={input.length_filtered_transcriptome} '
        'include=t '
        'names={input.transcript_hit_ids} '
        'substring=name '
        'out={output.pot_viral_transcripts} '
        '&> {log}'

rule filter_transcript_ids:
    input:
        blastx_res = 'output/recip_blast/viral_blastx/transcriptome_viral_blastx.outfmt3'
    output:
        transcript_hit_ids = 'output/recip_blast/viral_blastx/transcripts_viral_hit_ids.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/filter_transcript_ids.log'
    script:
        'scripts/recip_viral_blastx_transcript_hit_id_list.R'

rule recip_blastx_viral:
    input:
        unann_transcripts = 'output/trinotate/sorted/unann_transcripts.fasta',
        taxid_list = 'data/taxids/species_virus_taxids.txt'
    output:
        blastx_res = 'output/recip_blast/viral_blastx/transcriptome_viral_blastx.outfmt3'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/recip_blastx_viral.log'
    shell:
        'blastx '
        '-query {input.unann_transcripts} '
        '-db {params.blast_db} '
        '-taxidlist {input.taxid_list} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std staxids salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_unann_transcripts:
    input:
        length_filtered_transcriptome = 'output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        unann_transcript_ids = 'output/trinotate/sorted/ids_genes_no_blastx_or_viral_annot.txt'
    output:
        unann_transcripts = 'output/trinotate/sorted/unann_transcripts.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_unann_transcripts.log'
    shell:
        'filterbyname.sh '
        'in={input.length_filtered_transcriptome} '
        'include=t '
        'names={input.unann_transcript_ids} '
        'substring=name '
        'out={output.unann_transcripts} '
        '&> {log}'

rule sort_trinotate_annots:
    input:
        trinotate_report = 'output/trinotate/trinotate/trinotate_annotation_report.txt',
        longest_isoform_ids = 'output/trinity_filtered_isoforms/isoforms_by_length.txt'
    output:
        viral_or_unann_transcript_ids = 'output/trinotate/sorted/ids_genes_no_blastx_or_viral_annot.txt',
        best_annot_per_gene = 'output/trinotate/sorted/best_annot_per_gene.csv',
        longest_iso_annots = 'output/trinotate/sorted/longest_isoform_annots.csv'
    singularity:
        tidyverse_container
    log:
        'output/logs/sort_trinotate_annots.log'
    script:
        'scripts/sort_trinotate_annots.R'    

#######################################
##Transcriptome assembled & annotated##
#######################################

##annotate transcriptome using separate Trinotate Snakemake workflow
rule trinotate:
    input:
        fasta = 'output/trinity/Trinity.fasta',
        blastdb = 'bin/trinotate_db/uniprot_sprot.pep',
        hmmerdb = 'bin/trinotate_db/Pfam-A.hmm',
        sqldb = 'bin/trinotate_db/Trinotate.sqlite'
    output:
        'output/trinotate/trinotate/trinotate_annotation_report.txt',
        'output/trinotate/trinotate/Trinotate.sqlite'
    params:
        wd = 'output/trinotate'
    threads:
        50
    log:
        'output/logs/trinotate.log'
    singularity:
        trinotate_container
    shell:
        'trinotate_pipeline '
        '--trinity_fasta {input.fasta} '
        '--blast_db {input.blastdb} '
        '--hmmer_db {input.hmmerdb} '
        '--sqlite_db {input.sqldb} '
        '--outdir {params.wd} '
        '--threads {threads} '
        '&> {log}'

##make plot of busco summaries
rule plot_busco:
    input:
        ss = expand('output/busco/{filter_status}/short_summary.specific.hymenoptera_odb10.{filter_status}.txt',
            filter_status=['expression', 'length', 'unfiltered'])
    output:
        busco_plot = 'output/busco/short_summaries/busco_figure.png'
    params:
        ss_dir = 'output/busco/short_summaries'
    singularity:
        busco_container
    threads:
        25
    shell:
        'mkdir {params.ss_dir} & '
        'cp {input.ss} {params.ss_dir}/ & '
        'python3 /busco/scripts/generate_plot.py '
        '-wd {params.ss_dir}'

##busco run on both length-filtered (fasta with longest isoform per gene) and expression-filtered (fasta with most highly expressed isoform per gene) separately
##running on raw Trinity.fasta gives high duplication due to each gene having multiple isoforms
rule filtered_busco:
    input:
        filtered_fasta = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.fasta'
    output:
        'output/busco/{filter}/short_summary.specific.hymenoptera_odb10.{filter}.txt'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco_{filter}.log'))
    params:
        wd = 'output/busco',
        filtered_fasta = lambda wildcards, input: resolve_path(input.filtered_fasta)
    singularity:
        busco_container
    threads:
        25
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.filtered_fasta} '
        '--out {wildcards.filter} '
        '--lineage hymenoptera_odb10 '
        '--cpu {threads} '
        '--augustus_species nasonia '
        '--mode transcriptome '
        '-f '
        '&> {log} '

rule unfiltered_busco:
    input:
        unfiltered_fasta = 'output/trinity/Trinity.fasta'
    output:
        'output/busco/unfiltered/short_summary.specific.hymenoptera_odb10.unfiltered.txt'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco_unfiltered.log'))
    params:
        wd = 'output/busco',
        outdir = 'unfiltered',
        unfiltered_fasta = lambda wildcards, input: resolve_path(input.unfiltered_fasta)
    singularity:
        busco_container
    threads:
        25
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--force '
        '--in {params.unfiltered_fasta} '
        '--out {params.outdir} '
        '--lineage hymenoptera_odb10 '
        '--cpu {threads} '
        '--augustus_species nasonia '
        '--mode transcriptome '
        '-f '
        '&> {log} '

##run on length & expression fil. to see whether this decreases multimapping
##can also use to justify Salmon mapping onto length-filtered
rule fil_bowtie2_alignment_stats:
    input:
        transcriptome = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        alignment_stats = 'output/trinity_stats/isoforms_by_{filter}_bowtie2_alignment_stats.txt'
    params:
        index_basename = 'output/trinity_stats/isoforms_by_{filter}.fasta.index',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    threads:
        50
    singularity:
        trinity_container
    shell:
        'bowtie2-build '
        '{input.transcriptome} '
        '{params.index_basename} || exit 1 ; '
        'bowtie2 '
        '-p 10 '
        '-q '
        '--threads {threads} '
        '-x {params.index_basename} '
        '-1 {params.left} '
        '-2 {params.right} '
        '1> /dev/null 2> {output.alignment_stats}'

##basic read mapping statistics - read representation in assembly
rule bowtie2_alignment_stats:
    input:
        transcriptome = 'output/trinity/Trinity.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        alignment_stats = 'output/trinity_stats/bowtie2_alignment_stats.txt'
    params:
        index_basename = 'output/trinity_stats/Trinity.fasta.index',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    threads:
        50
    singularity:
        trinity_container
    shell:
        'bowtie2-build '
        '{input.transcriptome} '
        '{params.index_basename} || exit 1 ; '
        'bowtie2 '
        '-p 10 '
        '-q '
        '--threads {threads} '
        '-x {params.index_basename} '
        '-1 {params.left} '
        '-2 {params.right} '
        '1> /dev/null 2> {output.alignment_stats}'

##generate fasta files with either longest isoform per gene only OR most highly expressed isoform per gene
rule filter_trinity_isoforms:
    input:
        transcriptome = 'output/trinity/Trinity.fasta',
        isoforms = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.txt'
    output:
        sorted_fasta = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.fasta'
    log:
        'output/logs/filter_isoforms_by_{filter}.log'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input.transcriptome} '
        'include=t '
        'names={input.isoforms} '
        'out={output.sorted_fasta} ' 
        '&> {log}'

##sort isoforms of each transcript based of length or expression level
rule sort_isoforms_r:
    input:
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        expression = 'output/trinity_filtered_isoforms/isoforms_by_expression.txt',
        length = 'output/trinity_filtered_isoforms/isoforms_by_length.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/sort_isoforms_r.log'
    script:
        'scripts/sort_isoforms.R'

##generate exn50 stat - transcriptome equivalent of n50
rule ExN50_stats:
    input:
        abundance = 'output/trinity_abundance/RSEM.isoform.TPM.not_cross_norm',
        transcriptome = 'output/trinity/Trinity.fasta'
    output:
        ExN50_stats = 'output/trinity_stats/xn50.out.txt'
    singularity:
        trinity_container
    log:
        'output/logs/xn50.err.txt'
    shell:
        '/usr/local/bin/util/misc/contig_ExN50_statistic.pl '
        '{input.abundance} '
        '{input.transcriptome} '
        '>{output.ExN50_stats} '
        '2>{log}'

##generate basic assembly stats
rule trinity_stats:
    input:
        transcriptome = 'output/trinity/Trinity.fasta'
    output:
        stats = 'output/trinity_stats/stats.txt'
    singularity:
        trinity_container
    log:
        'output/logs/trinity_stats.log'
    shell:
        '/usr/local/bin/util/TrinityStats.pl '
        '{input.transcriptome} '
        '>{output.stats} '
        '2>{log}'

##generate info for downstream analyses
rule trinity_abundance_to_matrix:
    input:
        gt_map = 'output/trinity/Trinity.fasta.gene_trans_map',
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        'output/trinity_abundance/RSEM.isoform.counts.matrix',
        'output/trinity_abundance/RSEM.isoform.TPM.not_cross_norm'
    params:
        prefix = 'output/trinity_abundance/RSEM'
    singularity:
        trinity_container
    log:
        'output/logs/abundance_estimates_to_matrix.log'
    shell:
        '/usr/local/bin/util/abundance_estimates_to_matrix.pl '
        '--est_method RSEM '
        '--cross_sample_norm none '
        '--out_prefix {params.prefix} '
        '--gene_trans_map {input.gt_map} '
        '{input.abundance} '
        '&> {log}'

##generate info for downstream analyses
rule trinity_abundance:
    input:
        transcripts = 'output/trinity/Trinity.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        'output/trinity_abundance/RSEM.isoforms.results'
    singularity:
        trinity_container
    threads:
        50
    log:
        'output/logs/trinity_abundance.log'
    params:
        outdir = 'output/trinity_abundance',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    shell:
        '/usr/local/bin/util/align_and_estimate_abundance.pl '
        '--transcripts {input.transcripts} '
        '--seqType fq '
        '--est_method RSEM '
        '--aln_method bowtie2 '
        '--output_dir {params.outdir} '
        '--prep_reference '
        '--SS_lib_type RF '
        '--thread_count {threads} '
        '--trinity_mode '
        '--left {params.left} '
        '--right {params.right} '
        '&> {log}'

##assembly
rule Trinity:
    input:
        left = expand('output/bbmerge/{sample}_all_r1.fq.gz', sample=all_samples),
        right = expand('output/bbmerge/{sample}_unmerged_r2.fq.gz', sample=all_samples)
    output:
        'output/trinity/Trinity.fasta',
        'output/trinity/Trinity.fasta.gene_trans_map'
    params:
        outdir = 'output/trinity',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    singularity:
        trinity_container
    threads:
        50
    log:
        'output/logs/trinity.log'
    shell:
        'Trinity '
        '--SS_lib_type RF '
        '--max_memory 300G '
        '--CPU {threads} '
        '--output {params.outdir} '
        '--left {params.left} '
        '--right {params.right} '
        '--seqType fq '
        '&> {log}'

##merge the merged reads and unmerged_r1 reads files for Trinity
rule merge_all_r1_reads:
    input:
        r1 = 'output/bbmerge/{sample}_unmerged_r1.fq.gz',
        merged = 'output/bbmerge/{sample}_merged.fq.gz'
    output:
        joined = 'output/bbmerge/{sample}_all_r1.fq.gz'
    shell:
        'cat {input.r1} {input.merged} > {output.joined}'

##merge overlapping reads (improves assembly quality)
rule bbmerge:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        merged = 'output/bbmerge/{sample}_merged.fq.gz',
        unm1 = 'output/bbmerge/{sample}_unmerged_r1.fq.gz',
        unm2 = 'output/bbmerge/{sample}_unmerged_r2.fq.gz',
        ihist = 'output/bbmerge/{sample}_ihist.txt'
    params:
        adapters = bbduk_adapters
    log:
        'output/logs/bbduk_merge/{sample}.log'
    singularity:
        bbduk_container
    threads:
        20
    shell:
        'bbmerge.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.merged} '
        'outu1={output.unm1} '
        'outu2={output.unm2} '
        'ihist={output.ihist} '
        'verystrict=t '
        'adapters={params.adapters} '
        '&> {log}'

rule kraken_venom3:
    input:
        r1 = 'output/bbduk_trim/Mh_venom3_r1.fq.gz',
        r2 = 'output/bbduk_trim/Mh_venom3_r2.fq.gz',
        db = 'bin/db/kraken_std'
    output:
        out = 'output/kraken/kraken_Mh_venom3_out.txt',
        report = 'output/kraken/reports/kraken_Mh_venom3_report.txt'
    log:
        'output/logs/kraken/kraken_Mh_venom3.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

##kraken db from https://benlangmead.github.io/aws-indexes/k2 as build std was failing
rule kraken:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz',
        db = 'bin/db/kraken_std'
    output:
        out = 'output/kraken/kraken_{sample}_out.txt',
        report = 'output/kraken/reports/kraken_{sample}_report.txt'
    log:
        'output/logs/kraken/kraken_{sample}.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

rule fastqc_venom3:
    input:
        expand('output/bbduk_trim/Mh_venom3_r{n}.fq.gz', n=[1,2])
    output:
        "output/fastqc/Mh_venom3_r{n}_fastqc.zip"
    params:
        outdir = directory('output/fastqc')
    singularity:
        fastqc_container
    shell:
        'fastqc --outdir {params.outdir} {input}'

rule fastqc:
    input:
        expand('output/bbduk_trim/{sample}_r{n}.fq.gz',
            sample=all_samples, n=[1,2])
    output:
        expand('output/fastqc/{sample}_r{n}_fastqc.zip', sample=all_samples, n=[1,2])
    params:
        outdir = directory('output/fastqc')
    singularity:
        fastqc_container
    shell:
        'mkdir -p {params.outdir} ; '
        'fastqc --outdir {params.outdir} {input}'

##trim adapters and low quality seq.
rule bbduk_trim_venom3:
    input:
        r1 = 'output/joined/venom3_r1.fq.gz',
        r2 = 'output/joined/venom3_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim/Mh_venom3_r1.fq.gz',
        r2 = 'output/bbduk_trim/Mh_venom3_r2.fq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/bbduk_trim/venom3.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

##trim adapters and low quality seq.
rule bbduk_trim:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        adapters = bbduk_adapters
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        25
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 ftl=1 ' ##force trim left 1st base (0th) - lib prep kit added T overhang
        '&> {log}'

##need to run for venom3 as sample was excluded from transcriptome
rule cat_reads_venom3:
    input:
        l1r1 = 'data/reads/HM723BCX3/HM723BCX3-6267-15-58-01_S15_L001_R1_001.fastq.gz',
        l2r1 = 'data/reads/HM723BCX3/HM723BCX3-6267-15-58-01_S15_L002_R1_001.fastq.gz',
        l1r2 = 'data/reads/HM723BCX3/HM723BCX3-6267-15-58-01_S15_L001_R2_001.fastq.gz',
        l2r2 = 'data/reads/HM723BCX3/HM723BCX3-6267-15-58-01_S15_L002_R2_001.fastq.gz'
    output: 
        r1 = temp('output/joined/venom3_r1.fq.gz'),
        r2 = temp('output/joined/venom3_r2.fq.gz')
    threads:
        1
    shell:
        'cat {input.l1r1} {input.l2r1} > {output.r1} & '
        'cat {input.l1r2} {input.l2r2} > {output.r2} & '
        'wait'

##cat together any sample files split across multiple lanes
rule cat_reads:
    input:
        unpack(get_reads)
    output: 
        r1 = temp('output/joined/{sample}_r1.fq.gz'),
        r2 = temp('output/joined/{sample}_r2.fq.gz')
    threads:
        1
    shell:
        'cat {input.l1r1} {input.l2r1} > {output.r1} & '
        'cat {input.l1r2} {input.l2r2} > {output.r2} & '
        'wait'