rule sniffles_standard:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
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
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
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

rule sniffles_standard_hg38_scaffolded:
    """
    The same rule as sniffles_standard, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    and uses the self assembly fasta as a reference.
    """
    input:
        bam = "output/alignment/hg38_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/hg38_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} \
        --mapq {params.mapq} \
        --output-rnames &> {log}
        """

rule sniffles_mosaic_hg38_scaffolded:
    """
    The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    and uses the self assembly fasta as a reference.
    """
    input:
        bam = "output/alignment/hg38_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/hg38_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        # Remove limiter on mapQ due to multimapping
        minsupport = config['sniffles']['minsupport'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/hg38_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"
    shell:
        """
        # Remove limiter on mapQ due to multimapping
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """


use rule sniffles_standard_hg38_scaffolded as sniffles_standard_T2T_scaffolded with:
    input:
        bam = "output/alignment/T2T_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/T2T_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    log:
        "logs/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.log"

use rule sniffles_mosaic_hg38_scaffolded as sniffles_mosaic_T2T_scaffolded with:
    input:
        bam = "output/alignment/T2T_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/T2T_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    log:
        "logs/alignment/T2T_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"

use rule sniffles_mosaic as sniffles_mosaic_duplomap with:
    input:
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        index = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    log:
        "logs/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.log"

rule sniffles_mosaic_scaffolded_qc_all:
    """
    The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    and uses the self assembly fasta as a reference.
    Yields all candidates without filtering.
    """
    input:
        bam = "output/alignment/{scaffolded}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{scaffolded}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/{scaffolded}/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        minsupport = config['sniffles']['minsupport'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.log"
    shell:
        """
        # Yield all candidates, regardless of QC status
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} --mosaic \
        --minsupport=0 \
        --mapq=0 \
        --output-rnames \
        --mosaic-af-min={params.mosaic_af_min} \
        --mosaic-af-max={params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --dev-no-qc &> {log}
        """

use rule sniffles_mosaic_scaffolded_qc_all as sniffles_mosaic_scaffolded_mapq_mod_qc_all with:
    input:
        bam = "output/alignment/{scaffolded}/{mapper}/standard/mapped/{specimen}.sorted.merged.mapQ_modified.bam",
        index = "output/alignment/{scaffolded}/{mapper}/standard/mapped/{specimen}.sorted.merged.mapQ_modified.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/{scaffolded}/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.vcf.gz',
        snf='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.snf',
        tbi='output/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.vcf.gz.tbi'
    log:
        "logs/alignment/{scaffolded}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.log"
