rule estimate_coverage:
    input: 'output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam'
    output: temp('output/mapping/{refalias}/coverage_stats/{specimen}.default_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input} | genomeCoverageBed -ibam - > {output}'

rule concat_coverage:
    input: 'output/mapping/{refalias}/coverage_stats/{specimen}.default_hist.txt'
    output: 'output/mapping/{refalias}/coverage_stats/{specimen}.default_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """
    
rule estimate_coverage_d:
    input: 'output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam'
    output: temp('output/mapping/{refalias}/coverage_stats/{specimen}.d_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input} | genomeCoverageBed -d -ibam - > {output}'

rule concat_coverage_d:
    input: 'output/mapping/{refalias}/coverage_stats/{specimen}.d_hist.txt'
    output: 'output/mapping/{refalias}/coverage_stats/{specimen}.d_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

rule estimate_coverage_bg:
    input: 'output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam'
    output: temp('output/mapping/{refalias}/coverage_stats/{specimen}.bg_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input} | genomeCoverageBed -bg -ibam - > {output}'

rule concat_coverage_bg:
    input: 'output/mapping/{refalias}/coverage_stats/{specimen}.bg_hist.txt'
    output: 'output/mapping/{refalias}/coverage_stats/{specimen}.bg_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

rule estimate_coverage_bga:
    input: 'output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam'
    output: temp('output/mapping/{refalias}/coverage_stats/{specimen}.bga_hist.txt')
    conda: "../envs/deeptools.yml"
    shell: 'samtools view -b {input} | genomeCoverageBed -bga -ibam - > {output}'

rule concat_coverage_bga:
    input: 'output/mapping/{refalias}/coverage_stats/{specimen}.bga_hist.txt'
    output: 'output/mapping/{refalias}/coverage_stats/{specimen}.bga_hist_renamed.txt'
    conda: "../envs/deeptools.yml"
    shell: """
    cat  {input}  | awk '{{print $0, "{wildcards.specimen}"}}' | sed 's| |\t|g' >  {output}
    """

rule plotcoverage_pdf:
    input: expand('output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam', specimen=specimens, allow_missing = True)
    output: 
        # doesn't specify the rawcount output, due to wildcard conflict in the output naming
        plot = "output/mapping/{refalias}/coverage_stats/plotCoverage/all/coverage_plot.pdf",
        rawcounts = "output/mapping/{refalias}/coverage_stats/plotCoverage/all/rawcounts.coverage.tab"
    conda: "../envs/deeptools.yml"
    threads: 20
    params: 
        format = "pdf",
        title = "'Coverage - mapping against T2T CHM13'"
    shell: 
        """
        plotCoverage -p {threads} --bamfiles {input} --plotFile {output.plot} --plotFileFormat {params.format} -n 1000000 --plotTitle {params.title} --outRawCounts {output.rawcounts} --ignoreDuplicates --minMappingQuality 10 
        """

rule plotcoverage_plotly:
    input: expand('output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam', specimen=specimens, allow_missing = True)
    output:
        plot = "output/mapping/{refalias}/coverage_stats/plotCoverage/all/coverage_plot.html",
    conda: "../envs/deeptools.yml"
    threads: 20
    params: 
        format = "plotly",
        title = "'Coverage - mapping against T2T CHM13'"
    shell: 
        """
        # forgoes generating raw counts, as generated in pdf rule
        plotCoverage -p {threads} --bamfiles {input} --plotFile {output.plot} --plotFileFormat {params.format} -n 1000000 --plotTitle {params.title} --ignoreDuplicates --minMappingQuality 10 
        """