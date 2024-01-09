configfile: "config/snakemake/hg38.config.yml"
workdir: config['workdir']

# include helper functions
include: "rules/common.smk"
include: "rules/make_symlinks.smk"

# preprocessing steps
include: "rules/preprocessing.smk"

# Mapping
include: "rules/minimap2.smk"
include: "rules/winnowmap.smk"
include: "rules/duplomap.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"
include: "rules/pepper_etc.smk"
include: "rules/straglr.smk"
include: "rules/trgt_beta.smk"

# Simulations & Benchmarking
include: "rules/spike_in.smk"
include: "rules/spike_in_svtier1.smk"
include: "rules/spike_in_cmrg.smk"
include: "rules/in_silico_simulations.smk"

# Annotation
include: "rules/repeatmasker.smk"
include: "rules/process_mosaic.smk"

rule all:
    input:
        # tandem repeat genotyping
        expand('output/alignment/hg38/winnowmap/standard/variants/trgt-beta/{specimen}.vcf.gz',  specimen = specimens),
        # alu analysis pipeline rebuild
        expand("analysis/hg38/denovo/files/winnowmap/{setting}/variants/all.filt.vcf.gz", setting = ['standard', 'duplomap']),
        expand("analysis/hg38/denovo/files/winnowmap/{setting}/variants/all.count_repetitive.gff", setting = ['standard', 'duplomap']),
        expand("analysis/hg38/denovo/files/winnowmap/{setting}/variants/all.annotate_features.gff", setting = ['standard', 'duplomap']),
        expand("analysis/hg38/denovo/files/winnowmap/{setting}/repeatmasker/all.filt.alt.fa", setting = ['standard', 'duplomap']),
        expand("analysis/hg38/denovo/files/winnowmap/{setting}/variants/all.annotate_AluY.tsv", setting = ['standard', 'duplomap']),