os.environ["RUST_BACKTRACE"] = "full" # set for trgt backtrace

rule trgt_repeat_catalog:
    input:
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        index = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    output:
        tempvcf = temp("output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.vcf.gz"),
        tempbam = temp("output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.spanning.bam"),
        vcf = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.sorted.vcf.gz",
        csi = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.sorted.vcf.gz.csi",
        bam = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.sorted.spanning.bam",
        bai = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.sorted.spanning.bam.bai",
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda: "../envs/process_variants.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        catalog = config['paths']['trgt']['repeat_catalog'],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    log: 
        "logs/alignment/{refalias}/{mapper}/standard/variants/trgt/repeat_catalog/{specimen}.log"
    shell:
        """
        mkdir -p {params.outdir}
        trgt --threads {threads} \
        --genome {params.refgenome} \
        --repeats {params.catalog} \
        --reads {input.bam} \
        --output-prefix {params.outdir}/{wildcards.specimen} &> {log}

        # sort + index each file after creation, discard unsorted after complete
        bcftools sort {output.tempvcf} -o {output.vcf}
        bcftools index {output.vcf}
        samtools sort {output.tempbam} -o {output.bam}
        samtools index {output.bam}
        """

use rule trgt_repeat_catalog as trgt_pathogenic with:
    output:
        tempvcf = temp("output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.vcf.gz"),
        tempbam = temp("output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.spanning.bam"),
        vcf = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.sorted.vcf.gz",
        csi = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.sorted.vcf.gz.csi",
        bam = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.sorted.spanning.bam",
        bai = "output/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.sorted.spanning.bam.bai",
    params:
        refgenome = config['reference']['fasta'],
        catalog = config['paths']['trgt']['pathogenic'],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    log: 
        "logs/alignment/{refalias}/{mapper}/standard/variants/trgt/pathogenic/{specimen}.log"