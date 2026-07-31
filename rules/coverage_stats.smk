rule hg38_standard_coverage_all:
    input: 
        bams = expand('output/alignment/hg38/minimap2/standard/mapped/{specimen}.sorted.merged.bam', specimen=specimens),
        indices = expand('output/alignment/hg38/minimap2/standard/mapped/{specimen}.sorted.merged.bam.bai', specimen=specimens)
    output:
        plot = "output/alignment/hg38/minimap2/standard/coverage_stats/all.coverage.html",
        rawcounts = "output/alignment/hg38/minimap2/standard/coverage_stats/all.rawcounts.coverage.tab"
    conda: "../envs/deeptools.yml"
    threads: 5
    params: 
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with minimap2 on hg38'",
        labels = specimens
    shell: 
        """
        plotCoverage -p {threads} \
        --bamfiles {input.bams} \
        --plotFile {output.plot} \
        --plotFileFormat {params.format} \
        -n {params.sample_bp} \
        --plotTitle {params.title} \
        --outRawCounts {output.rawcounts} \
        --ignoreDuplicates \
        --minMappingQuality 10 \
        --labels {params.labels}
        """

use rule hg38_standard_coverage_all as hg38_duplomap_coverage_all with:
    input:
        bams = expand('output/alignment/hg38/minimap2/duplomap/mapped/{specimen}/realigned.bam', specimen=specimens),
        indices = expand('output/alignment/hg38/minimap2/duplomap/mapped/{specimen}/realigned.bam.bai', specimen=specimens)
    output:
        plot = "output/alignment/hg38/minimap2/duplomap/coverage_stats/all.coverage.html",
        rawcounts = "output/alignment/hg38/minimap2/duplomap/coverage_stats/all.rawcounts.coverage.tab"
    params: 
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with duplomap on hg38'",
        labels = specimens

use rule hg38_standard_coverage_all as scaffolded_standard_coverage_specimen with:
    input:
        bams = 'output/alignment/{ref}_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam',
        indices = 'output/alignment/{ref}_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output:
        plot = "output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}.coverage.html",
        rawcounts = "output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}.rawcounts.coverage.tab"
    params:
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with minimap2 on {ref}-scaffolded assembly'",
        labels = "{specimen}"

rule hg38_standard_coverage_chr:
    input:
        bam = 'output/alignment/hg38/minimap2/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/hg38/minimap2/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output:
        plot = 'output/alignment/hg38/minimap2/standard/coverage_stats/{specimen}/{chr}.coverage.html',
        rawcounts = 'output/alignment/hg38/minimap2/standard/coverage_stats/{specimen}/{chr}.rawcounts.tsv'
    conda:
        "../envs/deeptools.yml"
    threads: 2
    params:
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with minimap2 on {chr} - {specimen}'",
        region_format = "{chr}"
    shell:
        """
        plotCoverage -p {threads} \
        --bamfiles {input.bam} \
        --region {params.region_format} \
        --plotFile {output.plot} \
        --plotFileFormat {params.format} \
        -n {params.sample_bp} \
        --plotTitle {params.title} \
        --outRawCounts {output.rawcounts} \
        --ignoreDuplicates \
        --minMappingQuality 10 
        """

use rule hg38_standard_coverage_chr as hg38_duplomap_coverage_chr with:
    input:
        bams = 'output/alignment/hg38/minimap2/duplomap/mapped/{specimen}/realigned.bam',
        indices = 'output/alignment/hg38/minimap2/duplomap/mapped/{specimen}/realigned.bam.bai',
    output:
        plot = 'output/alignment/hg38/minimap2/duplomap/coverage_stats/{specimen}/{chr}.coverage.pdf',
        rawcounts = 'output/alignment/hg38/minimap2/duplomap/coverage_stats/{specimen}/{chr}.rawcounts.tsv'
    params:
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with duplomap on {chr} - {specimen}'",
        region_format = "{chr}"

use rule hg38_standard_coverage_chr as scaffolded_standard_coverage_chr with:
    input:
        bam = 'output/alignment/{ref}_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam',
        index = 'output/alignment/{ref}_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam.bai'
    output:
        plot = 'output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}/{chr}.{hap}.coverage.pdf',
        rawcounts = 'output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}/{chr}.{hap}.rawcounts.tsv'
    params:
        format = "plotly",
        sample_bp = 1000000,
        title = "'Coverage with duplomap on {chr} - {specimen}, {hap}'",
        region_format = "{chr}_RagTag_{hap}"