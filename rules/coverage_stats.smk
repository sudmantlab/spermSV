rule estimate_coverage:
    input: 
        bam = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output: temp('output/alignment/{refalias}/coverage_stats/{specimen}.default_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input.bam} | genomeCoverageBed -ibam - > {output}'

rule concat_coverage:
    input: 'output/alignment/{refalias}/coverage_stats/{specimen}.default_hist.txt'
    output: 'output/alignment/{refalias}/coverage_stats/{specimen}.default_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """
    
rule estimate_coverage_d:
    input: 
        bam = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output: temp('output/alignment/{refalias}/coverage_stats/{specimen}.d_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input.bam} | genomeCoverageBed -d -ibam - > {output}'

rule concat_coverage_d:
    input: 'output/alignment/{refalias}/coverage_stats/{specimen}.d_hist.txt'
    output: 'output/alignment/{refalias}/coverage_stats/{specimen}.d_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

rule estimate_coverage_bg:
    input: 
        bam = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output: temp('output/alignment/{refalias}/coverage_stats/{specimen}.bg_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input.bam} | genomeCoverageBed -bg -ibam - > {output}'

rule concat_coverage_bg:
    input: 'output/alignment/{refalias}/coverage_stats/{specimen}.bg_hist.txt'
    output: 'output/alignment/{refalias}/coverage_stats/{specimen}.bg_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

rule estimate_coverage_bga:
    input: 
        bam = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output: temp('output/alignment/{refalias}/coverage_stats/{specimen}.bga_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input.bam} | genomeCoverageBed -bga -ibam - > {output}'

rule concat_coverage_bga:
    input: 'output/alignment/{refalias}/coverage_stats/{specimen}.bga_hist.txt'
    output: 'output/alignment/{refalias}/coverage_stats/{specimen}.bga_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

# rule plotcoverage_pdf:
#     input: 
#         bams = expand('output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam', specimen=specimens, allow_missing = True),
#         indices = expand('output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai', specimen=specimens, allow_missing = True)
#     output: 
#         # doesn't specify the rawcount output, due to wildcard conflict in the output naming
#         plot = "output/alignment/{refalias}/{mapper}/standard/coverage_stats/coverage_plot.pdf",
#         rawcounts = "output/alignment/{refalias}/{mapper}/standard/coverage_stats/rawcounts.coverage.tab"
#     conda: "../envs/deeptools.yml"
#     threads: 20
#     params: 
#         format = "pdf",
#         title = "'Coverage - mapping against {refalias}'".format(refalias=config['reference']['alias'])
#     shell: 
#         """
#         plotCoverage -p {threads} \
#         --bamfiles {input.bams} \
#         --plotFile {output.plot} \
#         --plotFileFormat {params.format} \
#         -n 1000000 \
#         --plotTitle {params.title} \
#         --outRawCounts {output.rawcounts} \
#         --ignoreDuplicates \
#         --minMappingQuality 10 
#         """

# rule plotcoverage_all:
#     input: 
#         bams = expand('output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam', specimen=specimens, allow_missing = True),
#         indices = expand('output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai', specimen=specimens, allow_missing = True)
#     output:
#         plot = "output/alignment/{refalias}/{mapper}/standard/coverage_stats/coverage_plot.html",
#         rawcounts = "output/alignment/{refalias}/{mapper}/standard/coverage_stats/rawcounts.coverage.tab"
#     conda: "../envs/deeptools.yml"
#     threads: 10
#     params: 
#         format = "plotly",
#         sample_bp = 1000000
#         title = "'Coverage - mapping against {refalias}'".format(refalias=config['reference']['alias'])
#     shell: 
#         """
#         # forgoes generating raw counts, as generated in pdf rule
#         plotCoverage -p {threads} \
#         --bamfiles {input.bams} \
#         --plotFile {output.plot} \
#         --plotFileFormat {params.format} \
#         -n {params.sample_bp} \
#         --plotTitle {params.title} \
#         --outRawCounts {output.rawcounts} \
#         --ignoreDuplicates \
#         --minMappingQuality 10 
#         """

rule minimap2_coverage_single:
    input:
        bams = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam',
        indices = 'output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output:
        plot = 'output/alignment/{refalias}/{mapper}/standard/coverage_stats/{specimen}/{chr}.coverage.pdf',
        rawcounts = 'output/alignment/{refalias}/{mapper}/standard/coverage_stats/{specimen}/{chr}.rawcounts.tsv'
    conda: 
        "../envs/deeptools.yml"
    threads: 4
    params:
        format = 'pdf',
        sample_bp = 1000000,
        title = "'{specimen} - {chr}'"
    shell:
        """
        plotCoverage -p {threads} \
        --bamfiles {input.bams} \
        --region {wildcards.chr} \
        --plotFile {output.plot} \
        --plotFileFormat {params.format} \
        -n {params.sample_bp} \
        --plotTitle {params.title} \
        --outRawCounts {output.rawcounts} \
        --ignoreDuplicates \
        --minMappingQuality 10 
        """

use rule minimap2_coverage_single as duplomap_coverage_single with:
    input:
        bams = 'output/alignment/{refalias}/minimap2/duplomap/mapped/{specimen}/realigned.bam',
        indices = 'output/alignment/{refalias}/minimap2/duplomap/mapped/{specimen}/realigned.bam.bai'
    output:
        plot = 'output/alignment/{refalias}/minimap2/duplomap/coverage_stats/{specimen}/{chr}.coverage.pdf',
        rawcounts = 'output/alignment/{refalias}/minimap2/duplomap/coverage_stats/{specimen}/{chr}.rawcounts.tsv'