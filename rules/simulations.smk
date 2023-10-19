rule fractional_subsample:
    # Subsamples (with replacement) the parameter-calculated fraction of reads from the config-specified control, in order to achieve the desired spike-in coverage.
    # This fraction is calculated by taking the desired spike-in coverage and dividing it by the config-specified control coverage.
    input:
        "output/mapping/hg38/winnowmap/standard/{control}.sorted.merged.bam"
    output:
        "output/mapping/hg38/simulations/spike_in/{spike_coverage}/{control}_{spike_coverage}.bam"
    conda: "../envs/truvari.yml"
    threads: 20
    params:
        fraction = lambda wildcards, output: np.round(int(wildcards.spike_coverage.strip('X'))/config['simulations']['control']['coverage'], 3),
    shell:
        """
        samtools view -@ {threads} -b -s {params.fraction} {input} > {output}
        """

rule create_spike_in:
    # Merges the specified fraction of "spike-in" sample with the designated "main" file.
    # Right now, the "main" file is hard-coded to sample 894.
    input:
        main = "output/mapping/hg38/winnowmap/standard/894.sorted.merged.bam",
        control = expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/{control}_{spike_coverage}.bam", allow_missing = True, control = config['simulations']['control']['sample'])
    output:
        bam = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam",
        index = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam.bai"
    wildcard_constraints:
        spike_coverage = '[A-Za-z0-9]+'
    conda: "../envs/truvari.yml"
    threads: 20
    shell:
        """
        samtools merge -r -@ {threads} --output-fmt='BAM' {output.bam} {input.main} {input.control}
        samtools index -b -@ {threads} {output.bam}
        """

rule spike_in_sniffles:
    # Mosaic calls on the simulation (894 + HG002 spike-in).
    input:
        bam = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam",
        index = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam.bai"
    output:
        vcf='output/mapping/hg38/simulations/spike_in/{spike_coverage}/all.vcf.gz',
        snf='output/mapping/hg38/simulations/spike_in/{spike_coverage}/all.snf',
        tbi='output/mapping/hg38/simulations/spike_in/{spike_coverage}/all.vcf.gz.tbi'
    wildcard_constraints:
        spike_coverage = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads: 10
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/mapping/hg38/simulation/{spike_coverage}_spike.log"
    benchmark:
        "logs/mapping/hg38/simulation/{spike_coverage}_spike.bench.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {params.refgenome} \
        --tandem-repeats {params.repeats} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """

rule symlink_control:
    # Creates a symlink for the designated control specimen vcf.
    input:
        control_vcf = '/global/scratch/users/stacy-l/spermSV/output/mapping/hg38/sniffles/standard/single_sample/{control}.vcf.gz',
    output:
        control_link = 'output/mapping/hg38/simulations/control/{control}/all.vcf.gz',
    shell:
        """
        ln -s {input.control_vcf} {output.control_link}
        """

rule symlink_benchmark:
    # Creates a symlink for the designated benchmark vcf.
    input:
        svtier1_vcf = config['reference']['benchmarks']['svtier1']['all']
    output:
        svtier1_link = 'output/mapping/hg38/simulations/benchmarks/svtier1/all.vcf.gz'
    shell:
        """
        ln -s {input.svtier1_vcf} {output.svtier1_link}
        """

rule svtype_set:
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

# use rule svtype_set as control_svtype_set with:
#     # Creates a vcf of the set of given svtype from the germline (control) callset.
#     input:
#         vcf = 'output/mapping/hg38/sniffles/standard/single_sample/HG002.vcf.gz'
#     output: 
#         vcf = temp('output/mapping/hg38/simulations/control/HG002/{svtype}/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/control/HG002/{svtype}/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/control/HG002/{svtype}/all.vcf.gz.tbi'

# use rule svtype_set as benchmark_svtype_set with:
#     # Creates a vcf of the set of given svtype from the benchmark callset (below, hard-set to SV Tier 1).
#     input:
#         vcf = config['reference']['benchmarks']['svtier1']
#     output:
#         vcf = temp('output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/all.vcf.gz.tbi'

# use rule spike_in_ins as spike_in_del with:
#     # Creates a vcf of the set of deletions from the simulation callset.
#     output: 
#         vcf = temp('output/mapping/hg38/simulations/spike_in/{spike_coverage}/del/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/del/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/del/all.vcf.gz.tbi'
#     params:
#         svtype = 'DEL'

# use rule spike_in_ins as control_ins with:
#     # Creates a vcf of the set of insertions from the germline (control) callset.
#     input:
#         vcf = 'output/mapping/hg38/sniffles/standard/single_sample/HG002.vcf.gz'
#     output: 
#         vcf = temp('output/mapping/hg38/simulations/control/HG002/ins/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/control/HG002/ins/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/control/HG002/ins/all.vcf.gz.tbi'

# use rule spike_in_ins as control_del with:
#     # Creates a vcf of the set of deletions from the germline (control) callset.
#     input:
#         vcf = 'output/mapping/hg38/sniffles/standard/single_sample/HG002.vcf.gz'
#     output: 
#         vcf = temp('output/mapping/hg38/simulations/control/HG002/del/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/control/HG002/del/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/control/HG002/del/all.vcf.gz.tbi'
#     params:
#         svtype = 'DEL'

# use rule spike_in_ins as benchmark_ins:
#     # Creates a vcf of the set of deletions from the benchmark callset (below, hard-set to SV Tier 1).
#     input:
#         vcf = config['reference']['benchmarks']['svtier1']
#     output:
#         vcf = temp('output/mapping/hg38/simulations/benchmarks/svtier1/ins/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/benchmarks/svtier1/ins/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/benchmarks/svtier1/ins/all.vcf.gz.tbi'
#     params:
#         svtype = 'INS'

# use rule spike_in_ins as benchmark_del:
#     # Creates a vcf of the set of deletions from the benchmark callset (below, hard-set to SV Tier 1).
#     input:
#         vcf = config['reference']['benchmarks']['svtier1']
#     output:
#         vcf = temp('output/mapping/hg38/simulations/benchmarks/svtier1/del/all.vcf'),
#         compressed = 'output/mapping/hg38/simulations/benchmarks/svtier1/del/all.vcf.gz',
#         tbi = 'output/mapping/hg38/simulations/benchmarks/svtier1/del/all.vcf.gz.tbi'
#     params:
#         svtype = 'DEL'

rule svlen_bins:
    # Takes an vcf containing a set of given svtype and creates subset vcfs containing SVs within a given svlen span.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins/{lower}_{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins/{lower}_{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins/{lower}_{upper}.vcf.gz.tbi'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN >= {wildcards.lower} & SVLEN <= {wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule svlen_cumulative:
    # Takes an vcf containing a set of given svtype and creates subsets by max size, such that sequential size cutoffs can be used to generate cumulative distributions.
    input:
        vcf = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative/{upper}.vcf'),
        compressed = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative/{upper}.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative/{upper}.vcf.gz.tbi'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative'
    shell:
        """
        mkdir -p {params.outdir}
        bcftools filter -i 'SVLEN <= {wildcards.upper}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} -c {output.vcf} > {output.compressed}
        tabix {output.compressed}
        """

rule truvari_all:
    # Uses truvari to benchmark the given svtype subset against the matched svtype SV Tier 1 subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/all.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/all.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/all/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/all"
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} --passonly
        """

use rule truvari_all as truvari_svlen_bins with:
    # Uses truvari to benchmark the given svtype + svlen binned subset against the matched svtype SV Tier 1 subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/bins/{lower}_{upper}.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/bins/{lower}_{upper}.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/bins/{lower}_{upper}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/bins/{lower}_{upper}"

use rule truvari_all as truvari_svlen_cumulative with:
    # Uses truvari to benchmark the given svtype + svlen cumulative subset against the matched svtype SV Tier 1 subset.
    input: 
        query = 'output/mapping/hg38/simulations/{group}/{subgroup}/{svtype}/cumulative/{upper}.vcf.gz',
        benchmark = 'output/mapping/hg38/simulations/benchmarks/svtier1/{svtype}/cumulative/{upper}.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/cumulative/{upper}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        outdir = "output/mapping/hg38/simulations/{group}/{subgroup}/truvari/{svtype}/cumulative/{upper}"

# snakemake -pj40 --use-conda --conda-frontend mamba output/mapping/hg38/simulations/spike_in/1X/truvari/INS/cumulative/100/summary.json

# use rule spike_in_truvari_ins as spike_in_truvari_del with:
#     # Uses truvari to benchmark all spike-in deletions against the appropriate filtered SV Tier 1 callset.
#     input: 
#         query='output/mapping/hg38/simulations/spike_in/{spike_coverage}/del/all.vcf.gz'
#     output:
#         expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/simulation/del/all/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     params:
#         category = 'simulation',
#         svtype = 'del',
#         set = 'all',
#         benchmark = config['reference']['benchmarks']['del'],
#         outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"

# rule control_truvari_all:
#     # Uses truvari to benchmark the given svtype subset (control, germline) against the matched svtype SV Tier 1 subset.
#     input: 
#         query = 'output/mapping/hg38/simulations/spike_in/HG002/{svtype}/all.vcf'
#     output:
#         expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/germline/{svtype}/all/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 5
#     params:
#         category = 'germline',
#         set = 'all',
#         benchmark = config['reference']['benchmarks']['{svtype}'],
#         outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"
#     shell:
#         """ 
#         mkdir -p {params.outdir}
#         truvari bench -b {params.benchmark} -c {input.query} -o {params.outdir}/{params.category}/{wildcards.svtype}/{params.set} --passonly
#         """

# use rule spike_in_truvari_ins as control_truvari_del with:
#     # Uses truvari to benchmark all germline deletions against the appropriate filtered SV Tier 1 callset.
#     input: 
#         query='output/mapping/hg38/simulations/spike_in/HG002_germline_cutoff_del.vcf.gz'
#     output:
#         expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/germline/del/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     params:
#         category = 'germline',
#         svtype = 'del',
#         set = 'all',
#         benchmark = config['reference']['benchmarks']['del'],
#         outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"


# pressure testing: dry runs and DAGs
# snakemake -pj40 --use-conda --conda-frontend mamba -np --dag output/mapping/hg38/simulations/control/HG002/truvari/ins/bins/1000_2000/summary.json | dot -Tpng > test_dag.png
# snakemake -pj40 --use-conda --conda-frontend mamba -np --dag output/mapping/hg38/simulations/spike_in/1X/truvari/ins/bins/1000_2000/summary.json | dot -Tpng > test_dag.png
# snakemake -pj40 --use-conda --conda-frontend mamba -np --dag output/mapping/hg38/simulations/spike_in/1X/truvari/ins/cumulative/10000/summary.json | dot -Tpng > test_dag.png
