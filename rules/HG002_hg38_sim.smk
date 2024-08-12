
### HG38 SIM ###

use rule minimap2 as hg38_map_all with:
    # This rule maps all phased HG002 reads to the hg38 reference.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/hg38/unsorted/{hap}.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:all_{hap}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule minimap2 as hg38_map_SVs with:
    # This rule maps the phased <=10kB SV-spanning reads to the hg38 reference.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/hg38/unsorted/{hap}_SVs_to_hg38.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:{hap}_SVs_to_hg38\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule merge_SV_all as hg38_hap_1X_merge with:
    # This rule takes all hap1 <=10kB SV-spanning reads mapped to hg38 and merges them with all hap2 reads mapped to hg38.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs_to_hg38.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap2}.bam",
    output:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/1.0_{hap1}_SVs_with_{hap2}.bam"
    conda: "../envs/mapping.yml"
    threads: 10

use rule merge_SV_spike as hg38_hap_fraction_merge with:
    # This rule samples a fraction of the hap1 <=10kB SV-spanning reads mapped to hg38 and merges them with all hap2 reads mapped to hg38.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs_to_hg38.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap2}.bam",
    output:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/{spike}_{hap1}_SVs_with_{hap2}.bam"
    conda: "../envs/mapping.yml"
    threads: 10

use rule sniffles_standard as hg38_spike_in_fraction_germline_calls with:
    # This rule uses the standard Sniffles germline calling mode.
    # The minimum number of reads required to call a variant is dynamically determined from coverage at the region of interest.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads plus all hap2 reads, mapped to hg38.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/{spike}_{hap1}_SVs_with_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/{spike}_{hap1}_SVs_with_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.log"

use rule sniffles_mosaic as hg38_spike_in_fraction_mosaic_calls with:
    # This rule uses the mosaic (low allele frequency, AF) variant calling mode of Sniffles.
    # The minimum number of reads required to call a mosaic variant is set to 1.
    # The minimum AF required to call a mosaic variant is set to 0.
    # The maximum AF for a mosaic variant is 0.2.
    # Germline variants (exceeding 0.2 AF) are not called.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads plus all hap2 reads, mapped to hg38.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/{spike}_{hap1}_SVs_with_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/merged/{spike}_{hap1}_SVs_with_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.log"

rule filter_hap_SVs:
    input:
        script = "scripts/python/filter_hap_SVs.py",
        vcf = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz",
        rnames = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap2}_rnames.txt"
    output:
        vcf = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/hg38/{spike}_{hap1}_SVs_with_{hap2}.{hap2}_filtered.vcf.gz",
        tbi = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/hg38/{spike}_{hap1}_SVs_with_{hap2}.{hap2}_filtered.vcf.gz.tbi",
        json = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/hg38/{spike}_{hap1}_SVs_with_{hap2}.{hap2}_filtered.json"
    threads: 1
    conda: "../envs/sniffles.yml"
    log:
        "alignment/HG002/minimap2/standard/variants/{sniffles_setting}/hg38/{spike}_{hap1}_SVs_with_{hap2}.{hap2}_filtered.log"
    shell:
        """
        python {input.script} {input.rnames} {input.vcf} {output.vcf} {output.json} {log}
        tabix -p vcf {output.vcf} -f
        """

rule truvari_hg38_germline:
    # Query: SVs from a mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hg38, using Sniffles germline
    # Benchmark: The subset of <=10kB hap1 SVs from the HG002 benchmark.
    # When {spike} is 1.0, this determines the set of <=10kB SVs that can be called using only the SV-spanning reads.
    # As {spike} decreases, this determines the proportion of <=10kB SVs called with progressively fewer reads available to support variants.
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{spike}_{hap1}_SVs_with_{hap2}.jl',
        benchmark = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz',
        benchmark_index = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_SVs.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/{spike}_{hap1}_SVs_with_{hap2}/{outfiles}", allow_missing = True,
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

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --includebed {input.regions} \
        --passonly
        """

rule truvari_hg38_mosaic:
    # Query: SVs from a mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hg38, using Sniffles mosaic
    # Benchmark: The subset of <=10kB hap1 SVs from the HG002 benchmark.
    # This is an "uncorrected" evaluation of mosaic mode precision/recall that does not restrict the benchmark to the set of <=10kB SVs that can be called using only the SV-spanning reads.
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.jl',
        benchmark = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz',
        benchmark_index = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_SVs.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/mosaic/{spike}_{hap1}_SVs_with_{hap2}/{outfiles}", allow_missing = True,
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

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --includebed {input.regions} \
        --passonly
        """

use rule truvari_hg38_mosaic as truvari_hg38_mosaic_to_germline_TP with:
    # Query: SVs from a mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hg38, using Sniffles mosaic
    # Benchmark: The true positives from the 1.0 hap1 spike Sniffles germline call, representing the set of <=10kB SVs that can be called using only the SV-spanning reads.
    # This is used to evaluate the "corrected" performance of Sniffles mosaic mode on precision/recall.
    input:
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/hg38/{spike}_{hap1}_SVs_with_{hap2}.jl',
        benchmark = 'output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/1.0_{hap1}_SVs_with_{hap2}/tp-comp.vcf.gz',
        benchmark_index = 'output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/1.0_{hap1}_SVs_with_{hap2}/tp-comp.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_SVs.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/mosaic_to_germline_TP/{spike}_{hap1}_SVs_with_{hap2}/{outfiles}", allow_missing = True,
            outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
