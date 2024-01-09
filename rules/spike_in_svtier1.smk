rule svtier1_svtype_set:
    # Given a vcf, creates a vcf for the set of given svtype.
    input:
        vcf='output/mapping/hg38/simulations/{group}/{subgroup}/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz.tbi'
    conda:
        '../envs/truvari.yml'
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVTYPE=="{wildcards.svtype}"' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svtier1_svlen_ins_bins:
    # Takes an vcf containing a set of insertions and creates subset vcfs containing SVs within a given svlen span.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/INS/bins/{lower}_{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/bins/{lower}_{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/bins/{lower}_{upper}.vcf.gz.tbi'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/bins'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN >= {wildcards.lower} & SVLEN <= {wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svtier1_svlen_del_bins:
    # Takes an vcf containing a set of deletions and creates subset vcfs containing SVs within a given svlen span.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/DEL/bins/{lower}_{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/bins/{lower}_{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/bins/{lower}_{upper}.vcf.gz.tbi'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/bins'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN <= -{wildcards.lower} & SVLEN >= -{wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svtier1_svlen_ins_cumulative:
    # Takes an vcf containing a set of insertions and creates subsets by max size, such that sequential size cutoffs can be used to generate cumulative distributions.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/INS/cumulative/{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/cumulative/{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/cumulative/{upper}.vcf.gz.tbi'
    wildcard_constraints:
        upper = '[A-Za-z0-9]+'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/INS/cumulative'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN <= {wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svtier1_svlen_del_cumulative:
    # Takes an vcf containing a set of insertions and creates subsets by max size, such that sequential size cutoffs can be used to generate cumulative distributions.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/DEL/cumulative/{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/cumulative/{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/cumulative/{upper}.vcf.gz.tbi'
    wildcard_constraints:
        upper = '[A-Za-z0-9]+'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/DEL/cumulative'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN >= -{wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svtier1_truvari_all:
    # Uses truvari to benchmark the given svtype subset against the matched svtype subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/all.vcf.gz',
        index = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/all.vcf.gz.tbi',
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/all/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/all"
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} --passonly
        """

use rule svtier1_truvari_all as svtier1_truvari_svlen_bins with:
    # Uses truvari to benchmark the given svtype + svlen binned subset against the matched svtype subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins/{lower}_{upper}.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/bins/{lower}_{upper}.vcf.gz',
        index = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/bins/{lower}_{upper}.vcf.gz.tbi'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/bins/{lower}_{upper}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/bins/{lower}_{upper}"

use rule svtier1_truvari_all as svtier1_truvari_svlen_cumulative with:
    # Uses truvari to benchmark the given svtype + svlen cumulative subset against the matched svtype subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative/{upper}.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/cumulative/{upper}.vcf.gz',
        index = 'output/mapping/hg38/simulations/benchmarks/{benchmark_set}/{svtype}/cumulative/{upper}.vcf.gz.tbi'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/cumulative/{upper}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{benchmark_set}/{svtype}/cumulative/{upper}"