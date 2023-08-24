rule ref_index:
    input:
        "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fa"
    output:
        "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fai"
    threads: 20
    conda: "../envs/DeepVariant.yml"
    shell:
        """
        samtools faidx {input}
        """

rule DeepVariant:
    version: "1.5.0"
    input:
        bam = "output/hg38_no_alts/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam",
        index = "output/hg38_no_alts/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam.bai"
    params:
        # Make sure to change model to the correct model. 
        # mount_dir should be the highest level working dir, can even be your whole scratch probably
        model="PACBIO",
        refgenome = "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fa",
        refindex = "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fai",
        mount_dir = "/global/scratch/users/stacy-l",
        out_dir = "output/hg38_no_alts/DeepVariant",
        int_dir = "output/hg38_no_alts/DeepVariant/{model}/intermediate",
        cache_dir = "/global/scratch/users/stacy-l/software/singularity",
    output:
        vcf="output/hg38_no_alts/DeepVariant/{model}/{specimen}.vcf.gz",
        gvcf="output/hg38_no_alts/DeepVariant/{model}/{specimen}.g.vcf.gz",
        vcf_index="output/hg38_no_alts/DeepVariant/{model}/{specimen}.vcf.gz.tbi",
        gvcf_index="output/hg38_no_alts/no_altsDeepVariant/{model}/{specimen}.g.vcf.gz.tbi",
        report="output/hg38_no_alts/no_altsDeepVariant/{model}/{specimen}.visual_report.html"
    conda: "../envs/DeepVariant.yml"
    threads: # from FAQ: for best results, set to # of cores on the machine.
        40
    shell: # make sure that there are no spaces after the \ char! that is a newline
        """
        mkdir -p {params.out_dir}
        SINGULARITY_CACHEDIR={params.cache_dir}

        singularity exec --bind {params.mount_dir} \
        docker://google/deepvariant:{version} \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --ref {params.refgenome} \
            --reads {input.bam} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf} \
            --num_shards {threads} \
            --intermediate_results_dir {params.int_dir}/{wildcards.specimen}_int_results/
        """