rule samtools_sort:
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    threads: 20
    conda: "../envs/mapping.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    input: get_bams_per_sample
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    threads: 20
    conda: "../envs/mapping.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

# rule replace_RG:
#     # Creates a more informative read group for the merged BAM files.
#     # Creates a new file and index.
#     # In the future, add a renaming operation to have {specimen}.{smrtcell}.filt.sorted.bam as the output of samtools sort.
#     # This will yield an RG tag that maintains specimen name and smrtcell name in the ID of RG post-merge.
#     input: "output/alignment/{refalias}/{mapper}/{setting}/mapped/{specimen}.sorted.merged.bam"
#     output: "output/alignment/{refalias}/{mapper}/{setting}/mapped/{specimen}.sorted.merged.renamed.bam"
#     threads: 20
#     conda: "../envs/mapping.yml"
#     params:
#         readgroup = config['samtools']['readgroup']
#     shell:
#         """
#         samtools addreplacerg -@ {threads} -r {params.readgroup} {input} -o {output}
#         samtools index -b {output} -@ {threads}
#         """

rule index_bam:
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    threads: 20
    conda: "../envs/mapping.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """