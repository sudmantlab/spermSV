configfile: "config/snakemake/CHM13.config.yml"
workdir: config['workdir']

# include helper functions
include: "rules/common.smk"

# HiFi QC
include: "rules/HiFiAdapterFilt.smk"

# Mapping
include: "rules/minimap2.smk"
include: "rules/winnowmap.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"
include: "rules/pepper_etc.smk"
include: "rules/straglr.smk"
include: "rules/trgt_beta.smk"