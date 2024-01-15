rule cmrg_truvari_all:
    # Uses truvari to benchmark CMRG calls.
    # The CMRG benchmark vcf doesn't come with a pre-computed SVLEN field, so all subset analysis needs to be parsed before feeding to truvari.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/all.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/cmrg/all.vcf.gz',
        index = 'output/mapping/hg38/simulations/benchmarks/cmrg/all.vcf.gz.tbi'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/cmrg/all/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} --passonly
        """