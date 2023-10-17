rule benchmark_subset_ins:
    # Subsets insertions of config-param size cutoff from the benchmark callset.
    # This ensures that svtype + svlen cutoffs are matched between callsets, and calculations of precision & recall are accurate.
    # This rule is mostly for record-keeping and shouldn't *need* to be run, because these files were already manually generated.
    input:
        config['reference']['benchmarks']['all']
    output:
        vcf = temp(config['reference']['benchmarks']['ins'].strip('.gz')),
        compressed = config['reference']['benchmarks']['ins']
    conda:
        '../envs/truvari.yml'
    params:
        svtype = 'INS',
        svlen = config['simulations']['svlen']
    threads: 5
    shell:
        """
        bcftools filter -i 'SVTYPE=={params.svtype} & SVLEN <= {params.svlen}' {input} -o {output.vcf}
        bgzip -@ {threads} {output.vcf}
        tabix {output.compressed}
        """

rule benchmark_subset_del:
    # Subsets deletions of config-param size cutoff from the benchmark callset.
    # This ensures that svtype + svlen cutoffs are matched between callsets, and calculations of precision & recall are accurate.
    # This rule is mostly for record-keeping and shouldn't *need* to be run, because these files were already manually generated.
    input:
        config['reference']['benchmarks']['all']
    output:
        vcf = temp(config['reference']['benchmarks']['del'].strip('.gz')),
        compressed = config['reference']['benchmarks']['del']
    conda:
        '../envs/truvari.yml'
    params:
        svtype = 'DEL',
        svlen = config['simulations']['svlen']
    threads: 5
    shell:
        """
        bcftools filter -i 'SVTYPE=={params.svtype} & SVLEN <= {params.svlen}' {input} -o {output.vcf}
        bgzip -@ {threads} {output.vcf}
        tabix {output.compressed}
        """

rule fractional_subsample:
    # Subsamples (with replacement) the parameter-calculated fraction of reads from HG002, in order to achieve the desired spike-in coverage.
    # This fraction is calculated by taking the desired spike-in coverage and dividing it by the config-specified HG002 coverage.
    input:
        "output/mapping/hg38/winnowmap/standard/HG002.sorted.merged.bam"
    output:
        "output/mapping/hg38/simulations/spike_in/{specimen}_{spike_coverage}.bam"
    conda: "../envs/truvari.yml"
    threads: 20
    params:
        fraction = lambda wildcards, output: np.round(int(wildcards.spike_coverage.strip('X'))/config['simulations']['refcoverage'], 3)
    shell:
        """
        samtools view -@ {threads} -b -s {params.fraction} {input} > {output}
        """

rule create_spike_in:
    # Merges the specified fraction of "spike-in" sample with the designated "main" file.
    # Right now, the "main" file is hard-coded to sample 894.
    input:
        main_sample = "output/mapping/hg38/winnowmap/standard/894.sorted.merged.bam",
        subsampled = "output/mapping/hg38/simulations/spike_in/HG002_{spike_coverage}.bam"
    output:
        bam = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam",
        index = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam.bai"
    conda: "../envs/truvari.yml"
    threads: 20
    shell:
        """
        samtools merge -r -@ {threads} --output-fmt='BAM' {output.bam} {input.main_sample} {input.subsampled}
        samtools index -b -@ {threads} {output.bam}
        """

rule spike_in_sniffles:
    # Mosaic calls on the simulation (894 + HG002 spike-in).
    input:
        bam = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam",
        index = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.bam.bai"
    output:
        vcf='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.vcf.gz',
        snf='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.snf',
        tbi='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.vcf.gz.tbi'
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

rule spike_in_subset_ins:
    # Subsets insertions of config-param size cutoff from the vcf callset.
    input:
        vcf='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked.vcf.gz'
    output:
        vcf = temp('output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_ins.vcf'),
        compressed = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_ins.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_ins.vcf.gz.tbi'
    conda:
        '../envs/truvari.yml'
    params:
        svtype = 'INS',
        svlen = config['simulations']['svlen']
    threads: 10
    shell:
        """
        bcftools filter -i 'SVTYPE=={params.svtype} & SVLEN <= {params.svlen}' {input.vcf} -o {output.vcf}
        bgzip -@ {threads} {output.vcf}
        tabix {output.compressed}
        """

use rule spike_in_subset_ins as spike_in_subset_del with:
    output: 
        vcf = temp('output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_del.vcf'),
        compressed = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_del.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_del.vcf.gz.tbi'
    params:
        svtype = 'DEL',
        svlen = config['simulations']['svlen']

use rule spike_in_subset_ins as germline_subset_ins with:
    # Takes the HG002 germline callset and applies the same svtype and svlen filters for benchmark comparison.
    input:
        vcf = 'output/mapping/hg38/sniffles/standard/single_sample/HG002.vcf.gz'
    output: 
        vcf = temp('output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_ins.vcf'),
        compressed = 'output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_ins.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_ins.vcf.gz.tbi'

use rule spike_in_subset_ins as germline_subset_del with:
    # Takes the HG002 germline callset and applies the same svtype and svlen filters for benchmark comparison.
    input:
        vcf = 'output/mapping/hg38/sniffles/standard/single_sample/HG002.vcf.gz'
    output: 
        vcf = temp('output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_del.vcf'),
        compressed = 'output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_del.vcf.gz',
        tbi = 'output/mapping/hg38/simulations/spike_in/germline_comp/HG002_germline_cutoff_del.vcf.gz.tbi'
    params:
        svtype = 'DEL',
        svlen = config['simulations']['svlen']

rule spike_in_truvari_ins:
    # Uses truvari to benchmark spike-in insertions against the appropriate filtered SV Tier 1 callset.
    input: 
        query='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_ins.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/simulation/ins/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        category = 'simulation',
        svtype = 'ins',
        benchmark = config['reference']['benchmarks']['ins'],
        outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"
    shell:
        """ 
        mkdir -p {params.outdir}
        truvari bench -b {params.benchmark} -c {input.query} -o {params.outdir}/{params.category}/{params.svtype} --passonly
        """

use rule spike_in_truvari_ins as spike_in_truvari_del with:
    # Uses truvari to benchmark spike-in deletions against the appropriate filtered SV Tier 1 callset.
    input: 
        query='output/mapping/hg38/simulations/spike_in/{spike_coverage}/spiked_cutoff_del.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/simulation/del/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        category = 'simulation',
        svtype = 'del',
        benchmark = config['reference']['benchmarks']['del'],
        outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"

use rule spike_in_truvari_ins as germline_truvari_ins with:
    # Uses truvari to benchmark germline insertions against the appropriate filtered SV Tier 1 callset.
    input: 
        query='output/mapping/hg38/simulations/spike_in/HG002_germline_cutoff_ins.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/germline/ins/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        category = 'germline',
        svtype = 'ins',
        benchmark = config['reference']['benchmarks']['ins'],
        outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"

use rule spike_in_truvari_ins as germline_truvari_del with:
    # Uses truvari to benchmark germline deletions against the appropriate filtered SV Tier 1 callset.
    input: 
        query='output/mapping/hg38/simulations/spike_in/HG002_germline_cutoff_del.vcf.gz'
    output:
        expand("output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari/germline/del/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    params:
        category = 'germline',
        svtype = 'del',
        benchmark = config['reference']['benchmarks']['del'],
        outdir = "output/mapping/hg38/simulations/spike_in/{spike_coverage}/truvari"

# snakemake -pj40 --use-conda --conda-frontend mamba output/mapping/hg38/simulations/spike_in/1X/truvari/simulation/del/summary.json output/mapping/hg38/simulations/spike_in/3X/truvari/simulation/del/summary.json output/mapping/hg38/simulations/spike_in/5X/truvari/simulation/del/summary.json output/mapping/hg38/simulations/spike_in/10X/truvari/simulation/del/summary.json -np