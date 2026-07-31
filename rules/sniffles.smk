### Germline SVs ###

rule sniffles_standard:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+',
        refalias = '[A-Za-z0-9]+',
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


rule sniffles_standard_scaffolded:
    # The same rule as sniffles_standard, except it doesn't use a tandem repeat annotation file
    # and uses the self assembly fasta as a reference.
    input:
        bam = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz',
        snf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.snf',
        tbi='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_standard/{specimen}.log"
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

### Mosaic SVs (low-frequency SVs) ###

rule sniffles_mosaic:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+',
        refalias = '[A-Za-z0-9]+',
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

rule sniffles_mosaic_qc_all:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+',
        refalias = '[A-Za-z0-9]+',
    conda:
        '../envs/sniffles-dev.yml'
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
        "logs/alignment/{refalias}/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.log"
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
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --dev-no-qc &> {log}
        """

use rule sniffles_mosaic_qc_all as sniffles_mosaic_duplomap_qc_all with:
    input:
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        index = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam.bai"
    output:
        vcf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    log:
        "logs/alignment/{refalias}/{mapper}/duplomap/variants/sniffles_mosaic/{specimen}.qc_all.log"

rule sniffles_mosaic_scaffolded:
    # The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file
    # and uses the self assembly fasta as a reference.
    input:
        bam = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz',
        snf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.snf',
        tbi='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    params:
        # Remove limiter on mapQ due to multimapping
        minsupport = config['sniffles']['minsupport'],
        mapq = 0,
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.log"
    shell:
        """
        # Remove limiter on mapQ due to multimapping
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min {params.mosaic_af_min} \
        --mosaic-af-max {params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} &> {log}
        """

rule sniffles_mosaic_scaffolded_qc_all:
    # The same rule as sniffles_mosaic, except it doesn't use a tandem repeat annotation file (only compatible with hg38)
    # and uses the self assembly fasta as a reference.
    # Yields all candidates without filtering.
    input:
        bam = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{ref}_scaffolded/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai",
        fasta = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta"
    output:
        vcf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz',
        snf='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.snf',
        tbi='output/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz.tbi'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        '../envs/sniffles-dev.yml'
    threads:
        10
    params:
        minsupport = 0,
        mapq = 0,
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/{ref}_scaffolded/{mapper}/standard/variants/sniffles_mosaic/{specimen}.qc_all.log"
    shell:
        """
        # Yield all candidates, regardless of QC status
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {input.fasta} \
        --threads {threads} --mosaic \
        --minsupport {params.minsupport} \
        --mapq {params.mapq} \
        --output-rnames \
        --mosaic-af-min={params.mosaic_af_min} \
        --mosaic-af-max={params.mosaic_af_max} \
        --mosaic-qc-strand={params.mosaic_qc_strand} \
        --dev-no-qc &> {log}
        """