rule duplomap:
    input:
        "output/mapping/{refalias}/{mapper}/{setting}/{specimen}.sorted.merged.bam"
    output:
        "output/mapping/{refalias}/{mapper}/{setting}/duplomap/{specimen}/realigned.bam",
        "output/mapping/{refalias}/{mapper}/{setting}/duplomap/{specimen}/psvs.vcf.gz"
    params:
        refgenome = config['reference']['fasta'].strip('.gz'),
        db = config['duplomap']['db'],
        outdir = "output/mapping/{refalias}/{mapper}/{setting}/duplomap"
    conda:
        "../envs/mapping.yml"
    threads: 8
    shell:
        """
        duplomap -i {input} -d {params.db} -r {params.refgenome} -o {params.outdir}/{wildcards.specimen}
        """