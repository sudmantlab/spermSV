repo_path = config['paths']['trgt']['repo']

os.environ["RUST_BACKTRACE"] = "full" # set for trgt backtrace
if os.system("echo $PATH | grep /trgt") == 256:
    os.environ["PATH"] += os.pathsep + repo_path

rule trgt_beta:
    input:
        bam = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam.bai"
    output:
        "output/mapping/{refalias}/trgt-beta/{specimen}.vcf.gz",
        "output/mapping/{refalias}/trgt-beta/{specimen}.spanning.bam"
    conda: "../envs/trgt.yml"
    threads: 20
    params:
        binary = "trgt-beta-binary",
        refgenome= config['reference']['fasta'],
        catalog = config['paths']['trgt']['catalog'],
        outdir = "output/mapping/{refalias}/trgt-beta"
    log: 
        "logs/mapping/{refalias}/trgt-beta/{specimen}.log"
    shell:
        """
        mkdir -p {params.outdir}
        {params.binary} --threads {threads} --genome {params.refgenome} --repeats {params.catalog} --reads {input.bam} --output-prefix {params.outdir}/{wildcards.specimen} &> {log}
        """