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
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"
include: "rules/pepper_etc.smk"
include: "rules/straglr.smk"
include: "rules/trgt_beta.smk"

# Simulations & Benchmarking
include: "rules/simulations.smk"

# Annotation
include: "rules/repeatmasker.smk"

rule all:
    input:
        expand('output/mapping/hg38/simulations/spike_in/{coverage}/truvari/{svtype}/cumulative/{upper}/summary.json', coverage = ['1X', '2X', '3X', '4X', '5X', '7X', '10X'], svtype = ['INS', 'DEL'], upper = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000]),
        expand(expand('output/mapping/hg38/simulations/spike_in/{{coverage}}/truvari/{{svtype}}/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 5000, 7500, 10000]), svtype = ['INS', 'DEL'], coverage = ['1X', '2X', '3X', '4X', '5X', '7X', '10X']),
        expand('output/mapping/hg38/simulations/control/HG002/truvari/{svtype}/cumulative/{upper}/summary.json', svtype = ['INS', 'DEL'], upper = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000]),
        expand(expand('output/mapping/hg38/simulations/control/HG002/truvari/{{svtype}}/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 5000, 7500, 10000]), svtype = ['INS', 'DEL']),