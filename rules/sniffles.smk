rule sniffles_standard:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
    benchmark:
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.bench.log"
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
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
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
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"
    benchmark:
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.bench.log"
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
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """

use rule sniffles_mosaic as sniffles_mosaic_duplomap with:
    input:
        bam = "output/alignment/{refalias}/{mapper}/duplomap/{specimen}/realigned.bam",
        index = "output/alignment/{refalias}/{mapper}/duplomap/{specimen}/realigned.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    log:
        "logs/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.log"
    benchmark:
        "logs/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.bench.log"