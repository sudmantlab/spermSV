rule duplomap:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    output:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/psvs.vcf.gz"
    params:
        refgenome = config['reference']['fasta'].strip('.gz'),
        db = config['duplomap']['db'],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/mapping.yml"
    threads: 14
    shell:
        """
        duplomap -@ {threads} \
        -i {input.bam} -d {params.db} \
        -r {params.refgenome} \
        -o {params.outdir} --continue
        """

rule duplomap_index:
    input:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam"
    output:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam.bai"
    threads: 10
    conda: "../envs/mapping.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """