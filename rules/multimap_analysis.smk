rule multimap_analysis:
    input:
        script = "scripts/python/multimap_analysis.py",
        bam = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.bam",
        index = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.bam.bai",
    output:
        bam = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.mapQ_modified.bam",
        readstats = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.mapQ_modified.readstats.json",
        summary = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.mapQ_modified.summary.json",
        modified = "output/alignment/{refalias}/minimap2/standard/mapped/{specimen}.mapQ_modified.reads.txt"
    log:
        "logs/alignment/{refalias}/minimap2/standard/mapped/{specimen}.mapQ_modified.log"
    conda:
        "../envs/sniffles.yml"
    threads: 20
    shell:
        """
        python {input.script} \
        {input.bam} \
        --threads {threads} \
        --output {output.bam} \
        2> >(tee {log} >&2)
        """