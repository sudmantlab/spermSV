rule exclude_standard:
    # Generates a union exclusion file from .bed files of T2T_CHM13.
    # Note: compared to hg38, there are no gap files for exclusion.
    # Only excludes segdups and simple repeats bedfiles.
    input:
        segdups = config['reference']['annotations']['segdups'],
        simplerepeats = config['reference']['annotations']['segdups']
    conda:
        "../envs/straglr.yml"
    output:
        "output/mapping/{refalias}/straglr/exclude_standard/CHM13_exclude_standard.bed"
    shell:
        """
        # modified (fixed?) from the documentation's suggestion on creating an exclusion bed file
        # only retains first three fields of chr, start, stop 
        # and filters for simple repeats greater than 10000

        awk '$3-$2>=10000' {input.simplerepeats} | cut -f1-3 {input.segdups} - | bedtools sort -i - | bedtools merge -i - -d 1000 > {output}
        """

rule exclude_censat:
    # Generates a union exclusion file from .bed files of T2T_CHM13.
    # Note: compared to hg38, there are no gap files for exclusion.
    # Also excludes centromere + satellite bedfile.
    input:
        segdups = config['reference']['annotations']['segdups'],
        simplerepeats = config['reference']['annotations']['segdups'],
        censat = config['reference']['annotations']['censat']
    conda:
        "../envs/straglr.yml"
    output:
        "output/mapping/{refalias}/straglr/exclude_censat/CHM13_exclude_censat.bed"
    shell:
        """
        # modified (fixed?) from the documentation's suggestion on creating an exclusion bed file
        # only retains first three fields of chr, start, stop 
        # and filters for simple repeats greater than 10000

        awk '$3-$2>=10000' {input.simplerepeats} | cut -f1-3 {input.segdups} {input.censat} - | bedtools sort -i - | bedtools merge -i - -d 1000 > {output}
        """

rule straglr:
    input:
        bam = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam",
        bai = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam.bai",
        exclude = "output/mapping/{refalias}/straglr/{setting}/CHM13_{setting}.bed"
    output:
        "output/mapping/{refalias}/straglr/{setting}/{specimen}.tsv",
        "output/mapping/{refalias}/straglr/{setting}/{specimen}.bed"
    conda:
        "../envs/straglr.yml"
    threads:
        20
    params:
        repo = "code/straglr",
        refgenome = config['reference']['fasta'],
        out_dir = "output/mapping/{refalias}/straglr/{setting}",
        tmp_dir = "output/mapping/{refalias}/straglr/{setting}/tmp"
    shell:
        """
        # example code: genome scan to detect TRs longer than the reference genome by 100bp
        # can make options into params if necessary later
        
        mkdir -p {params.tmp_dir}
        mkdir -p {params.out_dir}

        python {params.repo}/straglr.py {input.bam} {params.refgenome} {params.out_dir}/{wildcards.specimen} \
            --tmpdir {params.tmp_dir} \
            --min_str_len 2 --max_str_len 100 --min_ins_size 100 --genotype_in_size \
            --exclude {input.exclude} --min_support 2 --max_num_clusters 2 --nprocs {threads}
        """