# no support for multi-threading :(
rule svim:
    input:
        bam = "output/hg38_no_alts/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam",
        index = "output/hg38_no_alts/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam.bai"
    output:
        "output/hg38_no_alts/svim/{specimen}/variants.vcf"
        # not explicitly specifying timestamped logfile and debugging files
    params:
        outdir = "output/hg38_no_alts/svim",
        refgenome = "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fa",
    conda: "../envs/svim.yml"
    shell:
        """
        svim alignment {params.outdir}/{wildcards.specimen} {input.bam} {params.refgenome} --verbose --sample {wildcards.specimen}
        """