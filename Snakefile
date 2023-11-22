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
include: "rules/spike_in.smk"
include: "rules/spike_in_svtier1.smk"
include: "rules/in_silico_simulations.smk"

# Annotation
include: "rules/repeatmasker.smk"

rule all:
    input:
        # for all below: remember to prioritize hg38 mapping > CHM13 mapping until hg38 is complete
        # generate up to date coverage plot
        expand("output/mapping/{refalias}/coverage_stats/plotCoverage/all/coverage_plot.html", refalias = ['hg38']),
        # mosaic calls on actual samples
        expand('output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz', refalias = ['hg38'], specimen = specimens),
        expand('analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa.out.gff', refalias = ['hg38'], analysis = ['alu'], specimen = specimens)
        # # svtier1 benchmarking on hg38 spike-in simulations
        # expand('output/mapping/hg38/simulations/spike_in/{coverage}/truvari/{benchmark_id}/{svtype}/cumulative/{upper}/summary.json', coverage = ['1X', '2X', '3X', '4X', '5X', '7X', '10X'], svtype = ['INS', 'DEL'], upper = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000], benchmark_id = ['svtier1']),
        # expand(expand('output/mapping/hg38/simulations/spike_in/{{coverage}}/truvari/{{benchmark_id}}/{{svtype}}/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 5000, 7500, 10000]), benchmark_id = ['svtier1'], svtype = ['INS', 'DEL'], 
        # coverage = ['1X', '2X', '3X', '4X', '5X', '7X', '10X']),
        # expand('output/mapping/hg38/simulations/control/HG002/truvari/{benchmark_id}/{svtype}/cumulative/{upper}/summary.json', benchmark_id = ['svtier1'], svtype = ['INS', 'DEL'], upper = [100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000]),
        # expand(expand('output/mapping/hg38/simulations/control/HG002/truvari/{{benchmark_id}}/{{svtype}}/bins/{lower}_{upper}/summary.json', zip, lower = [100, 250, 500, 750, 1000, 2500, 5000, 7500], upper = [250, 500, 750, 1000, 2500, 5000, 7500, 10000]), benchmark_id = ['svtier1'], svtype = ['INS', 'DEL'])