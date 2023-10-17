rule samtools_sort:
    input: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/mapping/{refalias}/{mapper}/{setting}/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    threads: 20
    conda: "../envs/minimap2.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    input: get_bams_per_sample
    output: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    threads: 20
    conda: "../envs/minimap2.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

# rule replace_RG:
#     # Creates a more informative read group for the merged BAM files.
#     # Creates a new file and index.
#     # In the future, add a renaming operation to have {specimen}.{smrtcell}.filt.sorted.bam as the output of samtools sort.
#     # This will yield an RG tag that maintains specimen name and smrtcell name in the ID of RG post-merge.
#     input: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
#     output: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.renamed.bam"
#     threads: 20
#     conda: "../envs/minimap2.yml"
#     params:
#         readgroup = config['samtools']['readgroup']
#     shell:
#         """
#         samtools addreplacerg -@ {threads} -r {params.readgroup} {input} -o {output}
#         samtools index -b {output} -@ {threads}
#         """

rule index_bam:
    input: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    output: "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam.bai"
    threads: 20
    conda: "../envs/minimap2.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """