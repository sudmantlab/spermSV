rule sniffles_standard:
    input:
        bam = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.vcf.gz',
        snf='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.snf',
        tbi='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.bench.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {params.refgenome} \
        --tandem-repeats {params.repeats} \
        --threads {threads} \
        --mapq {params.mapq} \
        --output-rnames &> {log}
        """

rule sniffles_mosaic:
    # Calls mosaic (somatic) SVs using the --mosaic option.
    # FOR NOW: try config with --minsupport=1; may need to configure minsupport-auto-mult (the coverage minsupport equation)
    # worst case, enable --no-qc flag.
    input:
        bam = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz',
        snf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.snf',
        tbi='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
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
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.bench.log"
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

rule sniffles_mosaic_no_qc:
    # Calls mosaic (somatic) SVs using the --mosaic option, without performing *any* QC filters (--no-qc option).
    input:
        bam = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/winnowmap/standard/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}_no_qc.vcf.gz',
        snf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}_no_qc.snf',
        tbi='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}_no_qc.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}_no_qc.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}_no_qc.bench.log"
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} \
        --reference {params.refgenome} --tandem-repeats {params.repeats} \
        --threads {threads} --mosaic \
        --mapq {params.mapq} \
        --output-rnames \
        --no-qc &> {log}
        """

# rule filter_sniffles:
#     # Filters the output vcf ahead of repeatmasker identification step.
#     input: 
#         vcf='output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.vcf.gz'
#     output: 
#         filtered= temp('output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.unheadered.filtered.vcf'),
#         fa = 'output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.filtered.fa'
#     threads:
#         10
#     run:
#         filter_vcf(input.vcf, output.filtered)

# rule reheader_sniffles:
#     # For some strange reason, the shell command to reheader won't work inside filter_single_vcf. :(
#     input:
#         unfiltered = 'output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.vcf.gz',
#         filtered = 'output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.unheadered.filtered.vcf'
#     output:
#         'output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.filtered.vcf'
#     threads:
#         10
#     shell:
#         """
#         "zcat {input.unfiltered} | head -72 - | cat - {input.filtered} > {output}"
#         """
