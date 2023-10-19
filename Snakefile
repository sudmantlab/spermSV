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

# only looking at the cumulative values with different coverage rn because expanding bin wildcards is annoying, and the binned results might not even be that interesting
rule all:
    input:
        expand('output/mapping/hg38/simulations/spike_in/{coverage}/truvari/{svtype}/cumulative/{upper}/summary.json', coverage = ['1X', '2X', '3X', '4X', '5X'], svtype = ['INS', 'DEL'], upper = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000]),
        expand('output/mapping/hg38/simulations/spike_in/1X/truvari/INS/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 500, 7500, 10000]),
        expand('output/mapping/hg38/simulations/spike_in/1X/truvari/DEL/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 500, 7500, 10000]),
        
