rule duplomap:
    input:
        "output/alignment/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    output:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/psvs.vcf.gz"
    params:
        refgenome = config['reference']['fasta'].strip('.gz'),
        db = config['duplomap']['db'],
        outdir = "output/alignment/{refalias}/{mapper}/{setting}/duplomap"
    conda:
        "../envs/mapping.yml"
    threads: 8
    shell:
        """
        duplomap -@ {threads} \
        -i {input} -d {params.db} \
        -r {params.refgenome} \
        -o {params.outdir}/{wildcards.specimen} --continue
        """

rule duplomap_index:
    input:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam"
    output:
        "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam.bai"
    threads: 8
    conda: "../envs/mapping.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """