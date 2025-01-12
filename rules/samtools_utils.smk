rule samtools_sort:
    group: "mapping"
    # TODO: Just change all files that have the .filt.sorted suffixes and just Not Do That
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    group: "mapping"
    input: get_temp_bams_per_sample
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule index_bam:
    group: "mapping"
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/mapping.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """
