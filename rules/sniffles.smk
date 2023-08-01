import numpy as np
import pandas as pd
import os 

sample_table = pd.read_table("samples.tsv", index_col=False, dtype=str)
all_samples = sample_table['Specimen'].unique()

rule sniffles_single_sample:
    version: subprocess.run(["sniffles --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input:
        bam = "output/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam",
        index = "output/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/sniffles/single_sample/{specimen}.vcf.gz',
        snf='output/sniffles/single_sample/{specimen}.snf',
        tbi='output/sniffles/single_sample/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna"
    log:
        "logs/sniffles/single_sample/{specimen}.log"
    benchmark:
        "logs/sniffles/single_sample/{specimen}.bench.log"
    shell:
        'sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} --reference {params.refgenome} --threads {threads}'

rule sniffles_multi_sample:
    # Note that this will call sniffles multi-input on *all* samples specified in the samples table.
    # If necessary, will build in support later for specifiying a certain subset of samples using a configfile.
    version: subprocess.run(["sniffles --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input: 
        expand('output/sniffles/single_sample/{specimen}.snf', specimen = all_samples)
    output: 
        'output/sniffles/multi_sample/multi_sample.vcf'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna"
    log:
        "logs/sniffles/multi_sample/multi_sample.log"
    benchmark:
        "logs/sniffles/single_sample/multi_sample.bench.log"
    shell:
        'sniffles --input {input} --vcf {output} --reference {params.refgenome} --threads {threads}'