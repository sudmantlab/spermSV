rule straglr:
    input:
        bam = "output/minimap2/softclip/{specimen}/{specimen}.sorted.merged.bam",
        bai = "output/minimap2/softclip/{specimen}/{specimen}.sorted.merged.bam.bai"
    output:
        "output/straglr/{specimen}.tsv",
        "output/straglr/{specimen}.bed"
    conda:
        "../envs/straglr.yml"
    threads:
        20
    params:
        repo = "code/straglr",
        refgenome = "/global/scratch/users/stacy-l/references/hg38/GCF_000001405.40_GRCh38.p14_genomic.fna",
        out_dir = "output/straglr",
        tmp_dir = "output/straglr/tmp",
    shell:
        """
        python {params.repo}/straglr.py {input.bam} {params.refgenome} {wildcards.specimen} --working_dir {params.out_dir} --tmpdir {params.tmp_dir}
        """