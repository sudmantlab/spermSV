### ARCHIVED OLD RULES ###


use rule minimap2 as self_assembly_SV_mapping with:
    # This rule maps the <=10kB SV-spanning reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/unsorted/{hap1}_SVs_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{hap1}_SVs_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

rule merge_SV_all:
    # This rule takes all hap1 <=10kB SV-spanning reads mapped to hap2 and merges them with all hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/1.0_{hap1}_SVs_to_{hap2}.bam",
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input.bam1} {input.bam2}
        """

rule merge_SV_spike:
    # This rule samples a fraction of the hap1 <=10kB SV-spanning reads mapped to hap2 and merges them with all hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_SVs_to_{hap2}.bam"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '0\.[0-9]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools view -bs {wildcards.spike} {input.bam1} |
        samtools merge -@ {threads} {output} - {input.bam2}
        """

rule merge_all:
    # This rule samples a fraction of *all* hap1 reads mapped to hap2 and merges them with a (1-fraction) amount of hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/subsampled.{spike}_all_{hap1}_to_{hap2}.bam"),
        merged = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_all_to_{hap2}.bam"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
    conda: "../envs/mapping.yml"
    threads: 10
    params:
        bulk = lambda wildcards: str(1 - float(wildcards.spike))
    shell:
        """
        samtools view -@ {threads} -bs {params.bulk} {input.bam2} -o {output.temp}
        samtools view -@ {threads} -bs {wildcards.spike} {input.bam1} | samtools merge -@ {threads} {output.merged} - {output.temp}
        """

use rule sniffles_standard as call_self_assembly_germline_spike with:
    # This rule uses the standard Sniffles germline calling mode.
    # The minimum number of reads required to call a variant is dynamically determined from coverage at the region of interest.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads mapped to hap2.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz.tbi'
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '[01]\.[0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        # simpleRepeat file created using grep filter on RepeatMasker track bed, not directly output from TRF
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.log"

use rule sniffles_mosaic as self_assembly_mosaic_spike with:
    # This rule uses the mosaic (low allele frequency, AF) variant calling mode of Sniffles.
    # The minimum number of reads required to call a mosaic variant is set to 1.
    # The minimum AF required to call a mosaic variant is set to 0.
    # The maximum AF for a mosaic variant is 0.2.
    # Germline variants (exceeding 0.2 AF) are not called.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads mapped to hap2.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz.tbi'
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '[01]\.[0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        # simpleRepeat file created using grep filter on RepeatMasker track bed, not directly output from TRF
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.log"

rule benchmark_liftover_MATERNAL:
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_MATERNAL.dict"
    output:
        vcf = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.vcf.gz",
        tbi = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.vcf.gz.tbi",
        rejected = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.rejected.vcf.gz",
    conda:
        "../envs/gatk.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_MATERNAL.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.MATERNAL.chain"
    shell:
        """
        picard -Xmx3G LiftoverVcf \
        -I {input.vcf} \
        -O {output.vcf} \
        -C {params.chain} \
        --REJECT {output.rejected} \
        --CREATE_INDEX true \
        -R {params.ref} \
        --VERBOSITY ERROR
        """

use rule benchmark_liftover_MATERNAL as benchmark_liftover_PATERNAL with:
    # TODO: Might be able to collapse this into one rule (with MATERNAL) if this doesn't flag the
    # requirement for all input/output to match wildcards.
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_PATERNAL.dict"
    output:
        vcf = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.vcf.gz",
        tbi = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.vcf.gz.tbi",
        rejected = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.rejected.vcf.gz",
    conda:
        "../envs/gatk.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_PATERNAL.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.PATERNAL.chain"

rule benchmark_vcftobed:
    input:
        vcf = "benchmarks/HG002/{filename}.vcf.gz"
    output:
        bed = "benchmarks/HG002/{filename}.bed"
    params:
        vcftobed = "/global/scratch/users/stacy-l/software/ucsc_utilities/vcfToBed"
    shell:
        """
        inputVcf={input.vcf}
        outPrefix=${{inputVcf%".vcf.gz"}}
        {params.vcftobed} {input.vcf} $outPrefix
        """

rule liftover_SV_regions_bed:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs.bed"
    output:
        mapped = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs.liftover.bed",
        unmapped = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs.liftover.unmapped.bed"
    conda:
        "../envs/liftover.yml"
    params:
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap2}.chain"
    shell:
        """
        liftOver {input} {params.chain} {output.mapped} {output.unmapped}
        """

rule truvari_anno_repmask:
    # For some reason it's broken?
    input:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf"
    output:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.repmask.vcf"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda:
        "../envs/truvari.yml"
    params:
        rm_exec = "/global/scratch/users/stacy-l/miniconda3/envs/truvari/bin/RepeatMasker"
    threads: 1
    shell:
        """
        truvari anno repmask -i {input} -o {output} \
        -e {params.rm_exec} \
        -t {threads}
        """

# rule truvari_toplevel_all:
#     input: 
#         # NOTE: We haven't filtered the direct callset on the parameters of <10kB and whatnot so the precision may be very bad. Again.
#         query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz',
#         jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.jl',
#         benchmark = "benchmarks/HG002/sub10kB_{hap1}_to_HG002_{hap2}.vcf.gz",
#         benchmark_index = "benchmarks/HG002/sub10kB_{hap1}_to_HG002_{hap2}.vcf.gz.tbi"
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{hap1}_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     wildcard_constraints:
#         hap1 = '[A-Za-z]+',
#         hap2 = '[A-Za-z]+'
#     conda: "../envs/truvari.yml"
#     threads: 5
#     params:
#         outdir = lambda wildcards, output: os.path.dirname(output[0]),
#     shell:
#         """
#         # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
#         # removing the preemptive outdir frees up the path for truvari to direct outfiles
#         rm -r {params.outdir}

#         # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
#         truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
#         --pctseq 0 \
#         --dup-to-ins \
#         --passonly
#         """

# TODO: Figure out how the hell the region restriction is supposed to work?

rule truvari_toplevel_subset:
    input:
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.jl',
        benchmark = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz',
        benchmark_index = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz.tbi',
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/1.0_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --dup-to-ins \
        --passonly
        """

rule truvari_spike_mosaic:
    # Benchmarks SV calling results from using Sniffles mosaic mode to call variants from:
    # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2
    # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2.
    # This determines the proportion of <=10kB SVs called with progressively fewer reads available to support variants.
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_{set}_to_{hap2}.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_{set}_to_{hap2}.vcf.gz.tbi",
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{spike}_{hap1}_{set}_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
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
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} --pctseq 0 --passonly
        """

rule get_hap_sub10kB_benchmark_SVs:
    # This rule takes in the HG002 structural variant benchmark VCF file, which was created by using dipcall
    # to call SVs on the diploid T2T HG002 assembly against the GRCh38 (hg38) reference genome.
    # It filters the VCF, retaining only contain variants <=10kb that are not INV or BND.
    # It then creates two VCFs, one for each haplotype, containing only the variants that are exclusive to that haplotype.
    # It accomplishes this by filtering for variants that are specifically phased to one of the haplotypes.
    input:
        vcf = 'benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz'
    output:
        maternal = "benchmarks/HG002/sub10kB_MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        maternal_index = "benchmarks/HG002/sub10kB_MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        paternal = "benchmarks/HG002/sub10kB_PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        paternal_index = "benchmarks/HG002/sub10kB_PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p {params.outdir}

        # Get maternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & SVTYPE!="INV" & GT=="1|0"' {input.vcf} -o {output.maternal} -O z9
        tabix {output.maternal}

        # Get paternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & SVTYPE!="INV" & GT=="0|1"' {input.vcf} -o {output.paternal} -O z9
        tabix {output.paternal}
        """

use rule truvari_hg38_germline_all as truvari_hg38_germline_sub10kB with:
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.jl',
        benchmark = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/sub10kB/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),

# use rule truvari_spike_mosaic as truvari_spike_germline with:
#     # Benchmarks SV calling results from using Sniffles germline mode to call variants from:
#     # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2
#     # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2.
#     # This determines the proportion of <=10kB SVs called with progressively fewer reads available to support variants.
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.vcf.gz",
#         jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.jl",
#         benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz",
#         benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz.tbi",
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{spike}_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 5

# use rule truvari_spike_mosaic as truvari_mode_comparison with:
#     # Benchmarks SV calling results, but compares results from two different modes.
#     # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2, called with mosaic mode
#     # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2, called with germline mode
#     # This is experimental and is intended to assess the impact of the different calling modes on SV calling.
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.vcf.gz",
#         jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.jl",
#         benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz",
#         benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz.tbi",
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mode_comparison/{spike}_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 5




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
