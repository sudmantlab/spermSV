rule filter_mosaic_calls:
    # Reheaders and filters a VCF from sniffles mosaic.
    # 1) Changes the SAMPLE field to a {specimen} named field.
    # 2) Excludes variants with extreme coverage or read support values.
    input:
        vcf = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    conda:
        "../envs/process_variants.yml"
    threads: 1
    params:
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1
    shell:
        """
        echo {wildcards.specimen} | bcftools reheader --threads {threads} -s - {input.vcf} |
        bcftools filter --threads {threads} \
        -i '{params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' \
        - -o {output.filt} --write-index
        """

rule join_mosaic_calls:
    # Joins together filtered calls from all available filtered VCFs.
    input:
        expand('output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
               allow_missing = True, specimen = specimens)
    output:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz",
        index = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz.csi",
    conda:
        "../envs/process_variants.yml"
    threads: 5
    shell:
        """
        bcftools merge --threads {threads} {input} -o {output.vcf} --write-index
        """

rule count_repetitive:
    # Annotates a file with the overlap count and percentage overlap of repetitive or duplicated features 
    # specified in the below bedfiles.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.count_repetitive.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    shell:
        """
        bedtools annotate -names -both -i {input} \
        -files /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_repeatmasker.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_simpleRepeat.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_genomicSuperDups.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_centromeres.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_microsat.bed \
        > {output}
        """

rule annotate_repetitive:
    # Annotates a file with information on repetitive features overlapped in the
    # below bedfiles.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_repetitive.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    shell:
        """
        bedtools intersect -a {input} -wao -names repeatMasker simpleRepeat genomicSuperDups centromeres microsat \
        -b /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_repeatmasker.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_simpleRepeat.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_genomicSuperDups.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_centromeres.bed \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_microsat.bed \
        > {output}
        """

rule svcf_to_gff:
    # Converts an existing Sniffles VCF to a GFF3 format.
    input:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output:
        gff = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    script:
        "../scripts/python/svcf_to_gff.py"

rule annotate_features:
    # Annotates a file with information on curated or putative features overlapped in the below bedfiles.
    # Requires strandedness match 
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_features.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    shell:
        """
        bedtools intersect -a {input} -wao -s -names refseq gencode_v44 \
        -b /global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.ncbiRefSeq.gtf.gz \
        /global/scratch/users/stacy-l/references/hg38_HGSVC/gencode.v44.basic.annotation.gff3.gz \
        > {output}
        """

rule annotate_AluY:
    # Annotates a file with information on RepeatMasker-identified insertions, then subsets to AluY* motif insertions.
    input:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz",
        gff = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa.out.gff"
    output:
        h5 = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/subset.annotate_AluY.h5"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    script:
        "../scripts/python/annotate_AluY.py"

# use rule count_repetitive as count_repetitive_duplomap_minimap2 with:
#     input:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.gff"
#     output:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.count_repetitive.gff"

# use rule annotate_repetitive as annotate_repetitive_duplomap_minimap2 with:
#     input:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.gff"
#     output:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.annotate_repetitive.gff"

# use rule annotate_features as annotate_features_duplomap_minimap2 with:
#     # Annotates a (reduced) GFF3 file with information on curated or putative features overlapped in the
#     # below bedfiles.
#     input:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.gff"
#     output:
#         "analysis/{refalias}/alu/files/duplomap/minimap2/AluY_samples.annotate_features.gff"

# rule annotate_AluY:
#     # Takes the output VCF from sniffles mosaic and annotates with information from
#     # the RepeatMasker GFF output, specifically information on AluY* elements.
#     input:
#         vcf = 'output/alignment/hg38/winnowmap/standard/variants/sniffles_mosaic/{specimen}.vcf.gz'
#         gff = 'analysis/hg38/alu/files/winnowmap_standard/repeatmasker/{specimen}.alt.fa.out.gff'
#     output:
#     params:
#         motif = 'AluY',
#         mosaic_vcf = "",
#         repeatmasker_gff = ""
#     script:
#         "scripts/annotate_alu.py"