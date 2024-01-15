rule pepper_margin_deepvariant_hifi:
    # try to run this through the image build to cache step once to avoid infinite blob DL spam
    # do not alter the output directory structure of outdir/specimen/chr/specimen.region.vcf.gz
    # there are issues with intermediate_files and logs directories not being able to play well w/ each other
    # unless kept in separate directories per chr
    input:
        bam = "output/mapping/{refalias}/minimap2/{setting}/{specimen}.sorted.merged.bam"
    output:
        "output/CHM13/pepper_etc/{specimen}/{region}/{specimen}.{region}.vcf.gz" 
    params:
        mount_dir = config['singularity']['mount_dir'],
        refgenome = config['reference']['fasta'],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        mapQ = config['pepper_etc']['mapQ']
    threads: 4
    conda:
        # borrow existing environment
        "../envs/DeepVariant.yml"
    shell:
        """
        mkdir -p {params.outdir}

        # Run PEPPER-Margin-DeepVariant from docker image
        # heads up, it's a lot of blobs

        singularity exec --bind {params.mount_dir} \
        docker://kishwars/pepper_deepvariant:r0.8 \
        run_pepper_margin_deepvariant call_variant \
        --hifi --phased_output -b {input.bam} -f {params.refgenome} -o {params.outdir} \
        --t {threads} -p {wildcards.specimen}.{wildcards.region} -s {wildcards.specimen} -r {wildcards.region} --pepper_min_mapq {params.mapQ} \
        """