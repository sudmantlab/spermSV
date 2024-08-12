rule samtools_sort:
    # TODO: Just change all files that have the .filt.sorted suffixes and just Not Do That
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule samtools_haploid_sort:
    input: "output/alignment/self_assembly/{mapper}/standard/mapped/temp/{specimen}.{hap}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/alignment/self_assembly/{mapper}/standard/mapped/temp/{specimen}.{hap}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+',
        hap = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    input: get_temp_bams_per_sample
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule collate_haploid_bams:
    input: get_temp_bams_per_sample
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    output: "output/alignment/self_assembly/{mapper}/standard/mapped/{specimen}.{hap}.sorted.merged.bam"
    threads: 10
    conda: "../envs/mapping.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule index_bam:
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

rule index_haploid_bam:
    input: "output/alignment/self_assembly/{mapper}/standard/mapped/{specimen}.{hap}.sorted.merged.bam"
    output: "output/alignment/self_assembly/{mapper}/standard/mapped/{specimen}.{hap}.sorted.merged.bam.bai"
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+',
        hap = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/mapping.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """
