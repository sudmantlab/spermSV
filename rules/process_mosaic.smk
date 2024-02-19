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

rule svcf_to_gff:
    # Converts an existing Sniffles VCF to a GFF3 format.
    input:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output:
        gff = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 4
    script:
        "../scripts/python/svcf_to_gff.py"

rule overlap_repetitive:
    # Creates a file with the overlap count and percentage overlap of repetitive or duplicated features 
    # specified in the below bedfiles.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.overlap_repetitive.gff"
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

rule annotate_genomicSuperDups_bedfile:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_genomicSuperDups.gff"
    conda:
        "../envs/process_variants.yml"
    params:
        name = 'genomicSuperDups',
        database = config['reference']['annotations']['segdups'],
        stranded = '-s' # This is a dummy param that should be default on, only removed (set to "")
        # for bedfiles that do not have strand information.
    threads: 4
    shell:
        """
        bedtools intersect -a {input} {params.stranded} -wao -names {params.name} \
        -b {params.database} \
        > {output}
        """

use rule annotate_genomicSuperDups_bedfile as annotate_repeatMasker_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_repeatMasker.gff"
    params:
        name = 'repeatMasker',
        database = config['reference']['annotations']['repeatmasker'],
        stranded = '-s' # This is a dummy param that should be default on, only removed (set to "")
        # for bedfiles that do not have strand information.

use rule annotate_genomicSuperDups_bedfile as annotate_simpleRepeat_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_simpleRepeat.gff"
    params:
        name = 'simpleRepeat',
        database = config['reference']['annotations']['repeats'],
        stranded = ''

use rule annotate_genomicSuperDups_bedfile as annotate_centromeres_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_centromeres.gff"
    params:
        name = 'centromeres',
        database = config['reference']['annotations']['censat'],
        stranded = '-s' # This is a dummy param that should be default on, only removed (set to "")
        # for bedfiles that do not have strand information.

use rule annotate_genomicSuperDups_bedfile as annotate_microsat_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_microsat.gff"
    params:
        name = 'microsat',
        database = config['reference']['annotations']['microsat'],
        stranded = '-s' # This is a dummy param that should be default on, only removed (set to "")
        # for bedfiles that do not have strand information.

rule annotate_gencode_features:
    # Annotates a file with information on curated or putative features overlapped in Gencode.
    # Requires strandedness match.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_gencode.gff"
    conda:
        "../envs/process_variants.yml"
    threads: 8
    params:
        name = 'gencode_v44',
        database = config['reference']['annotations']['gencode']
    shell:
        """
        bedtools intersect -a {input} -wao -s -names {params.name} \
        -b {params.database} \
        > {output}
        """

use rule annotate_gencode_features as annotate_refseq_features with:
    # Annotates a file with information on curated or putative features overlapped in NCBI Refseq.
    # Requires strandedness match.
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.annotate_refseq.gff"
    params:
        name = 'refseq',
        database = config['reference']['annotations']['refseq']

rule write_alt_fasta:
    # Takes a (sniffles) vcf and writes INS alleles to a fasta file for RepeatMasker. 
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz"
    output: 
        fa = 'analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa'
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    conda:
        "../envs/process_variants.yml"
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
        allow_missing = True, suffix = ['tbl', 'out', 'out.gff', 'masked', 'cat'])
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

rule analyze_repetitive_insertions:
    # Annotates a file with information on RepeatMasker-identified insertions.
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.filt.vcf.gz",
        out = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa.out"
    output:
        tsv = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/all.repetitive_insertions.tsv"
    conda:
        "../envs/process_variants.yml"
    threads: 2
    script:
        "../scripts/python/annotate_repeatmasker.py"

rule annotate_edit_standard:
    # Given a specimen, pulls variants + read names from annotate_repeatmasker output and summarizes 
    # edit distance for read(s) relative to the reference sequence.
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        analysis_utils = "scripts/python/analysis_utils.py",
        bam = "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam",
        tsv = "analysis/{refalias}/denovo/files/{mapper}/standard/variants/all.repetitive_insertions.tsv"
    output:
        tsv = "analysis/{refalias}/denovo/files/{mapper}/standard/variants/{specimen}.annotate_edit_distance.tsv"
    conda:
        "../envs/process_variants.yml"
    threads: 2
    script:
        "../scripts/python/annotate_edit_distance.py"

use rule annotate_edit_standard as annotate_edit_duplomap with:
    input:
        parsing_utils = "scripts/python/parsing_utils.py",
        analysis_utils = "scripts/python/analysis_utils.py",
        bam = "output/alignment/{refalias}/{mapper}/duplomap/mapped/{specimen}/realigned.bam",
        tsv = "analysis/{refalias}/denovo/files/{mapper}/duplomap/variants/all.repetitive_insertions.tsv"
    output:
        tsv = "analysis/{refalias}/denovo/files/{mapper}/duplomap/variants/{specimen}.annotate_edit_distance.tsv"

# rule concat_annotate_edit:
#     # Combine the per-specimen tabular outputs for each setting.
#     input:
#         expand("analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{specimen}.annotate_edit_distance.tsv")
#     output:
#         "analysis/{refalias}/denovo/files/{mapper}/duplomap/variants/all.annotate_edit_distance.tsv"
#     conda:
#         "../envs/process_variants.yml"
#     threads: 2
#     shell:
#         """
#         head -n 1 file1.tsv > combined.tsv
#         for file in *.tsv; do
#             tail -n +2 "$file" >> combined.tsv
#         done
#         """