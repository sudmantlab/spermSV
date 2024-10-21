rule extract_inserts:
    input: 
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        fastq = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq"
    output:
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_extracted.tsv",
        merged=temp("output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_merged.txt")
    threads: 10
    resources:
        mem_gb=8
    conda:
        "../envs/somrit.yml"
    params:
        somrit_dir = "code/somrit",
        base_dir = config["workdir"],
    shell:
        """
        cd {params.somrit_dir}
        python somrit.py extract \
        --bam {params.base_dir}/{input.bam} \
        --output-merged {params.base_dir}/{output.merged} \
        --output-tsv {params.base_dir}/{output.tsv} \
        --fastq-file {params.base_dir}/{input.fastq} \
        --threads {threads}
        cd {params.base_dir}
        """

rule bam_to_fastq:
    input:
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
    output:
        fastq =  "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq"
    threads: 5
    conda:
        "../envs/mapping.yml"
    shell:
        "samtools fastq -@ {threads} -c 6 -T '*' {input.bam} -0 {output.fastq}"

rule fastq_index:
    input:
        fastq = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq"
    output:
        fastq_index = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq.fai"
    conda:
        "../envs/mapping.yml"
    shell:
        "samtools faidx {input.fastq}"

rule realign_inserts:
    input:
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_extracted.tsv",
        fastq = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq",
        fastq_index = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq.fai"
    output:
        bam="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned.bam",
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned.tsv"
    threads: 10
    resources:
        mem_gb=20
    conda:
        "../envs/somrit.yml"
    params:
        somrit_dir = "code/somrit",
        base_dir = config["workdir"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        refgenome = config['reference']['fasta'],
        output_prefix = "{specimen}_realigned"
    shell:
        """
        cd {params.somrit_dir}
        python somrit.py realign \
        --bam-list {params.base_dir}/{input.bam} \
        --tsv-list {params.base_dir}/{input.tsv} \
        --fastq-list {params.base_dir}/{input.fastq} \
        --output-dir {params.base_dir}/{params.outdir} \
        --tsv-prefix {params.output_prefix} \
        --bam-prefix {params.output_prefix} \
        --reference-genome {params.refgenome} \
        --threads {threads} \
        --filter-depth \
        --max-insert-size 10000 \
        --max-depth 500 
        cd {params.base_dir}
        """


rule classify_inserts:
    input:
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_extracted.tsv",
        bam="output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        bam_index="output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam.bai",
        realign_tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned.tsv",
        fastq = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.fastq"
    output:
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned_classified.tsv"
    threads: 1
    resources:
        mem_gb=24
    conda:
        "../envs/somrit.yml"
    params:
        somrit_dir = "code/somrit",
        base_dir = config["workdir"],
        repbase="code/somrit/utils/rt_seqs.fa"
    shell:
        """
        cd {params.somrit_dir}
        python somrit.py classify \
        --bam-list {params.base_dir}/{input.bam} \
        --sample-list {wildcards.specimen} \
        --tsv-list {params.base_dir}/{input.tsv} \
        --realign-tsv {params.base_dir}/{input.realign_tsv} \
        --annotation-file {params.repbase} \
        --fastq-list {params.base_dir}/{input.fastq} \
        --output-tsv {params.base_dir}/{output.tsv} 
        cd {params.base_dir}
        """

rule filter_inserts:
    input:
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned_classified.tsv",
        bam="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned.bam",
        bam_index="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned.bam.bai"
    output:
        tsv="output/alignment/{refalias}/{mapper}/somrit/mapped/{specimen}_realigned_classified_filtered.tsv"
    threads: 10
    resources:
        mem_gb=20
    conda:
        "../envs/somrit.yml"
    params:
        somrit_dir = "code/somrit",
        base_dir = config["workdir"],
        refgenome = config['reference']['fasta'],
        repbase="code/somrit/utils/rt_seqs.fa",
        centromeres="code/somrit/utils/hg38_centromeres.tsv",
        telomeres="code/somrit/utils/hg38_telomeres.tsv",
    shell:
        """
        cd {params.somrit_dir}
        python somrit.py filter \
        --threads {threads} \
        --input-tsv {params.base_dir}/{input.tsv} \
        --bam {params.base_dir}/{input.bam} \
        --reference-genome {params.refgenome} \
        --centromeres {params.centromeres} \
        --telomeres {params.telomeres} \
        --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """


rule somrit_index:
    input:
        "output/alignment/{refalias}/{mapper}/somrit/mapped/{prefix}.bam"
    output:
        "output/alignment/{refalias}/{mapper}/somrit/mapped/{prefix}.bam.bai"
    conda:
        "../envs/mapping.yml"
    resources:
        mem_gb=16
    threads: 1
    shell:
        "samtools index -@ {threads} {input}"
