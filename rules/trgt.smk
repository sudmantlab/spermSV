import os
import subprocess

root_path = subprocess.run(["pwd", "-P"], stdout=subprocess.PIPE).stdout.decode('utf-8').strip('\n')
repo_path = "{root_dir}/code/trgt".format(root_dir = root_path)

if os.system("echo $PATH | grep /trgt") == 256:
    os.environ["PATH"] += os.pathsep + repo_path

rule rename_bed:
    # only here for documentation purposes
    # hopefully this never needs to be re-run again...
    input:
        "code/trgt/repeats/repeat_catalog.hg38.bed"
    output:
        "code/trgt/repeats/repeat_catalog.hg38.ncbi.bed"
    params:
        mapping="code/trgt/GRCh38_UCSC2RefSeq.txt" # created by swapping columns in python
    conda: "../envs/trgt.yml"
    threads: 1
    shell:
        """
        cvbio UpdateContigNames \
        -i {input} \
        -o {output} \
        -m {params.mapping} \
        --columns 0 \
        --skip-missing false
        """

rule trgt:
    version: subprocess.run(["code/trgt/trgt-v0.4.0-linux_x86_64 --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input:
        bam = "output/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam",
        index = "output/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam.bai"
    output:
        "output/trgt/{specimen}.vcf.gz",
        "output/trgt/{specimen}.spanning.bam"
    conda: "../envs/trgt.yml"
    threads: 10
    params:
        binary = "trgt-v0.4.0-linux_x86_64",
        refgenome= "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna",
        catalog = "code/trgt/repeats/repeat_catalog.hg38.ncbi.bed",
        outdir = "output/trgt"
    shell:
        """
        {params.binary} --genome {params.refgenome} --repeats {params.catalog} --reads {input.bam} --output-prefix {params.outdir}/{wildcards.specimen}
        """