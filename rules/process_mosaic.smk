rule preprocess_variants:
    # Takes in the VCF output of sniffles mosaic and preprocesses with the following steps:
    # 1) Prepends the ID field with the sample name (ex. Sniffles2.BND.F38S0 -> 894_Sniffles2.BND.F38S0).
    # 2) Filters out variants with extreme coverage or read support values.
    # 3) Filters out variants on uncertain contigs (_random, chrUn_)
    # 4) Filters out excessive SV lengths (greater than 10000)
    # The ID value passed on to downstream analyses will retain the sample origin.
    input:
        vcf = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    conda:
        "../envs/process_variants.yml"
    params:
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1
    shell:
        """
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY {input.vcf} | \
        bcftools annotate --set-id '{wildcards.specimen}\_%ID' - | \
        bcftools filter - \
        -i 'INFO/SVLEN <= 10000 & {params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' \
        -o {output.filt} --write-index
        """

rule join_filtered_calls:
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
        bcftools concat --threads {threads} -a {input} -o {output.vcf} --write-index
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
        gff = temp("analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff")
    conda:
        "../envs/process_variants.yml"
    threads: 8
    script:
        "../scripts/python/svcf_to_gff.py"

rule annotate_features:
    # Annotates a file with information on curated or putative features overlapped in the below bedfiles.
    # Requires strandedness match.
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

rule write_alt_fasta:
    # Takes a (sniffles) vcf and writes INS alleles to a fasta file for RepeatMasker. 
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output: 
        fa = 'analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 2
    script:
        "../scripts/python/svcf_alt_to_fasta.py"
        
rule analyze_alt_fasta:
    # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
    # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
    input:
        "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa"
    output:
        expand("analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.alt.fa.{suffix}", 
        allow_missing = True, suffix = ['tbl', 'out.gff', 'masked', 'cat'])
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{specimen}.filt.log"
    threads: 10
    conda:
        '../envs/RepeatMasker.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} -dir {params.outdir} {input} &> {log}
        """

rule annotate_repeatmasker:
    # Annotates a file with information on RepeatMasker-identified insertions.
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz",
        out = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa.out"
    output:
        tsv = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_repeatmasker.tsv"
    conda:
        "../envs/process_variants.yml"
    threads: 2
    script:
        "../scripts/python/annotate_repeatmasker.py"