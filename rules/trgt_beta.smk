os.environ["RUST_BACKTRACE"] = "full" # set for trgt backtrace

rule trgt_beta:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        "output/alignment/{refalias}/{mapper}/standard/variants/trgt-beta/{specimen}.vcf.gz",
        "output/alignment/{refalias}/{mapper}/standard/variants/trgt-beta/{specimen}.spanning.bam"
    conda: "../envs/trgt.yml"
    threads: 20
    params:
        refgenome = config['reference']['fasta'],
        catalog = config['paths']['trgt']['catalog'],
        binary = config['paths']['trgt']['binary'],
        outdir = "output/alignment/{refalias}/{mapper}/standard/variants/trgt-beta"
    log: 
        "logs/alignment/{refalias}/{mapper}/standard/variants/trgt-beta/{specimen}.log"
    shell:
        """
        mkdir -p {params.outdir}
        {params.binary} --threads {threads} \
        --genome {params.refgenome} \
        --repeats {params.catalog} \
        --reads {input.bam} \
        --output-prefix {params.outdir}/{wildcards.specimen} &> {log}
        """