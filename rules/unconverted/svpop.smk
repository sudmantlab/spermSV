rule svpop_setup:
    # Creates symlinks for svpop run.
    input:
        config = "/global/scratch/users/stacy-l/spermSV/config/svpop/config.json",
        sample_table = "/global/scratch/users/stacy-l/spermSV/config/svpop/samples.tsv",
    output:
        config = "code/svpop/config/config.json",
        sample_table = "code/svpop/config/samples.tsv",
    params:
        config_dir = "code/svpop/config"
    shell:
        """
        mkdir -p {params.config_dir}

        ln -s {input.config} {output.config}
        ln -s {input.sample_table} {output.sample_table}
        """

rule svpop_merge:
    # full documentation on github: https://github.com/EichlerLab/svpop/blob/master/MERGE.md
    # it's... a document
    input:
        config = "code/svpop/config/config.json",
        sample_table = "code/svpop/config/samples.tsv",
    output:
        # temp placeholder
        "output/hg38_no_alts/svpop/placeholder/merge/run.success"
    threads:
        # set max threads avail to resource
        # svpop has thread params per rule execution
        40
    params:
        repo_dir = "code/svpop",
        placeholder_dir = "output/hg38_no_alts/svpop/placeholder",
        merges = ["results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/bed/sv_ins.bed.gz",
        "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/bed/sv_del.bed.gz",
        "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/bed/sv_inv.bed.gz"]
    conda:
        "../envs/svpop.yml"
    shell:
        """
        cd {params.repo_dir}
        snakemake -s Snakefile {params.merges} --cores {threads}
        cd -

        mkdir -p {params.placeholder_dir}
        touch {output}
        """

rule svpop_repeatmasker:
    input:
        config = "code/svpop/config/config.json",
        sample_table = "code/svpop/config/samples.tsv",
    output:
        # temp placeholder
        "output/hg38_no_alts/svpop/placeholder/repeatmasker/run.success"
    threads:
        # set max threads avail to resource
        # svpop has thread params per rule execution
        40
    params:
        repo_dir = "code/svpop",
        placeholder_dir = "output/hg38_no_alts/svpop/placeholder",
        outfiles = ["results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/rmsk/rmsk-table_sv_ins.tsv.gz",
                    "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/rmsk/rmsk-table_sv_del.tsv.gz",
                    "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/rmsk/rmsk-table_sv_inv.tsv.gz"]
    conda:
        "../envs/svpop.yml"
    shell:
        """
        cd {params.repo_dir}
        snakemake -s Snakefile {params.outfiles} --cores {threads}
        cd -

        mkdir -p {params.placeholder_dir}
        touch {output}
        """

rule svpop_trf:
    input:
        config = "code/svpop/config/config.json",
        sample_table = "code/svpop/config/samples.tsv",
    output:
        # temp placeholder
        "output/hg38_no_alts/svpop/placeholder/trf/run.success"
    threads:
        # set max threads avail to resource
        # svpop has thread params per rule execution
        40
    params:
        repo_dir = "code/svpop",
        placeholder_dir = "output/hg38_no_alts/svpop/placeholder",
        outfiles = ["results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/trf/trf-table_sv_ins.tsv.gz",
                    "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/trf/trf-table_sv_del.tsv.gz",
                    "results/variant/sampleset/sniffles_mosaic_t2t/spermSV/all/all/anno/trf/trf-table_sv_inv.tsv.gz"]
    conda:
        "../envs/svpop.yml"
    shell:
        """
        cd {params.repo_dir}
        snakemake -s Snakefile {params.outfiles} --cores {threads}
        cd -

        mkdir -p {params.placeholder_dir}
        touch {output}
        """


# rule svpop_intersect:
#     # full documentation on github: https://github.com/EichlerLab/svpop/blob/master/MERGE.md
#     # it's... a document
#     input:
#         # manually input the desired outfiles... sigh
#         config = "code/svpop/config/config.json",
#         sample_table = "code/svpop/config/samples.tsv"
#     output:
#         # temp placeholder
#         "output/hg38_no_alts/svpop/placeholder/intersect/run.success"
#     threads:
#         # set max threads avail to resource
#         # svpop has thread params per rule execution
#         40
#     params:
#         repo_dir = "code/svpop",
#         placeholder_dir = "output/hg38_no_alts/svpop/placeholder",
#         intersects = ["results/variant/intersect/caller+sniffles_mosaic_hg38+894/caller+svim_hg38+894/szro-match80/all/all/svs_ins/intersect.tsv.gz"] 
#     conda:
#         "../envs/svpop.yml"
#     shell:
#         """
#         cd {params.repo_dir}
#         snakemake -s Snakefile {params.intersects} --cores {threads}
#         cd -

#         mkdir -p {params.placeholder_dir}
#         touch {output}
#         """