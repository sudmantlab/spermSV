import os
import pandas as pd 

def get_bams_per_sample(wildcards):
    bam_path = "output/minimap2/{setting}/{specimen}/{lane}/{smrtcell}.filt.sorted.bam"
    samples = pd.read_table("samples.tsv", index_col=False, dtype=str)
    samples = samples[samples["Specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [bam_path.format(setting=wildcards.setting, specimen=s[0], lane=s[1], smrtcell = s[2]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

rule minimap2_standard: 
    version: subprocess.run(["minimap2 --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input:
        hifi = "output/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: temp("output/minimap2/standard/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna",
        readgroup = '@RG\\tID:{lane}\\tSM:{specimen}\\tPU:{smrtcell}\\tPL:PACBIO',
        minQ = 10
    conda: "../envs/minimap2.yml"
    threads: 40
    shell: "minimap2 --version && minimap2 -R '{params.readgroup}' -t {threads} -ax map-hifi {params.refgenome} {input.hifi} | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}"

rule minimap2_softclip:
    # Runs minimap2 with -Y setting to soft-clip supplementary alignments.
    # Required by straglr.
    version: subprocess.run(["minimap2 --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input:
        hifi = "output/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: temp("output/minimap2/softclip/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna",
        readgroup = '@RG\\tID:{lane}\\tSM:{specimen}\\tPU:{smrtcell}\\tPL:PACBIO',
        minQ = 10
    conda: "../envs/minimap2.yml"
    threads: 40
    shell: "minimap2 --version && minimap2 -Y -R '{params.readgroup}' -t {threads} -ax map-hifi {params.refgenome} {input.hifi} | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}"

rule samtools_sort:
    version: subprocess.run(["samtools --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input: "output/minimap2/{setting}/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: temp("output/minimap2/{setting}/{specimen}/{lane}/{smrtcell}.filt.sorted.bam")
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    version: subprocess.run(["samtools --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input: get_bams_per_sample
    output: "output/minimap2/{setting}/{specimen}/{specimen}.sorted.merged.bam"
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule index_bam:
    version: subprocess.run(["samtools --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input: "output/minimap2/{setting}/{specimen}/{specimen}.sorted.merged.bam"
    output: "output/minimap2/{setting}/{specimen}/{specimen}.sorted.merged.bam.bai"
    threads: 32
    conda: "../envs/minimap2.yml"
    shell: """
    samtools index -b {input} -@ {threads}
    """