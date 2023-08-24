rule samtools_sort:
    input: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/mapping/{refalias}/{mapper}/{setting}/{specimen}/{lane}/{smrtcell}.filt.sorted.bam")
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    input: get_bams_per_sample
    output: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule index_bam:
    input: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    output: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam.bai"
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """