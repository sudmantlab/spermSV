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
    group: "process_mosaic"
    threads: 1
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

rule preprocess_hg38_scaffolded_variants:
    # Takes in the VCF output of sniffles mosaic and preprocesses with the following steps:
    # 1) Prepends the ID field with the sample name (ex. Sniffles2.BND.F38S0 -> 894_Sniffles2.BND.F38S0).
    # 2) Filters out variants with extreme coverage or read support values.
    # 3) Filters out variants on uncertain contigs (_random, chrUn_)
    # 4) Filters out excessive SV lengths (greater than 10000)
    # The ID value passed on to downstream analyses will retain the sample origin.
    input:
        vcf = 'output/alignment/hg38_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/hg38_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/hg38_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    conda:
        "../envs/process_variants.yml"
    group: "process_mosaic"
    threads: 1
    params:
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1
    shell:
        """
        regions=$(for i in {{1..22}} X Y; do echo -n "chr${{i}}_RagTag_hap1,chr${{i}}_RagTag_hap2,"; done | sed 's/,$//')
        
        bcftools view -r $regions {input.vcf} | \
        bcftools annotate --set-id '{wildcards.specimen}\_%ID' - | \
        bcftools filter - \
        -i 'INFO/SVLEN <= 10000 & {params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' \
        -o {output.filt} --write-index
        """

rule preprocess_T2T_scaffolded_variants:
    # Takes in the VCF output of sniffles mosaic and preprocesses with the following steps:
    # 1) Prepends the ID field with the sample name (ex. Sniffles2.BND.F38S0 -> 894_Sniffles2.BND.F38S0).
    # 2) Filters out variants with extreme coverage or read support values.
    # 3) Filters out variants on uncertain contigs (_random, chrUn_)
    # 4) Filters out excessive SV lengths (greater than 10000)
    # The ID value passed on to downstream analyses will retain the sample origin.
    input:
        vcf = 'output/alignment/T2T_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/T2T_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/T2T_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    conda:
        "../envs/process_variants.yml"
    group: "process_mosaic"
    threads: 1
    params:
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1
    shell:
        """
        regions=$(for i in {{1..22}} X Y; do echo -n "chr${{i}}_RagTag_hap1,chr${{i}}_RagTag_hap2,"; done | sed 's/,$//')
        
        bcftools view -r $regions {input.vcf} | \
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
    group: "process_mosaic"
    threads: 1
    shell:
        """
        bcftools concat --threads {threads} -a {input} -o {output.vcf} --write-index
        """

rule join_filtered_calls_hg38_scaffolded:
    # Joins together filtered calls from all available filtered VCFs.
    input:
        expand('output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
               allow_missing = True, specimen = ['900', '898', 'PD50477f', 'PD50508f', 'PD50519d', 'TSC7237',
               '901', 'TSC6830', '0068', '5619', '895', 'PD50489e', 'PD50521e', '899', '5621',
               'PD50508b', '894', '5980', 'PD50523b', 'PD46180c', 'PD50521b'])
    output:
        vcf = "analysis/hg38_scaffolded/denovo/files/minimap2/standard/variants/all.filt.vcf.gz",
        index = "analysis/hg38_scaffolded/denovo/files/minimap2/standard/variants/all.filt.vcf.gz.csi",
    conda:
        "../envs/process_variants.yml"
    group: "process_mosaic"
    threads: 1
    shell:
        """
        bcftools concat --threads {threads} -a {input} -o {output.vcf} --write-index
        """

rule join_filtered_calls_T2T_scaffolded:
    # Joins together filtered calls from all available filtered VCFs.
    input:
        expand('output/alignment/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
               allow_missing = True, specimen = ['900', '898', 'PD50477f', 'PD50508f', 'PD50519d', 'TSC7237',
               '901', 'TSC6830', '0068', '5619', '895', 'PD50489e', 'PD50521e', '899', '5621',
               'PD50508b', '894', '5980', 'PD50523b', 'HG002', 'PD46180c', 'PD50521b'])
    output:
        vcf = "analysis/T2T_scaffolded/denovo/files/minimap2/standard/variants/all.filt.vcf.gz",
        index = "analysis/T2T_scaffolded/denovo/files/minimap2/standard/variants/all.filt.vcf.gz.csi",
    group: "process_mosaic"
    conda:
        "../envs/process_variants.yml"
    threads: 1
    shell:
        """
        bcftools concat --threads {threads} -a {input} -o {output.vcf} --write-index
        """

rule vcf_to_gff3:
    # Converts an existing Sniffles VCF to a GFF3 format.
    input:
        script = "scripts/python/vcf_to_gff3.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.vcf.gz"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.gff"
    group: "process_mosaic"
    threads: 1
    shell:
        """
        python {input.script} {input.vcf} {output}
        """

rule overlap_repetitive:
    # Creates a file with the overlap count and percentage overlap of repetitive or duplicated features 
    # specified in the below bedfiles.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.overlap_repetitive.tsv"
    group: "process_mosaic"
    conda:
        "../envs/process_variants.yml"
    params:
        repeats = config['reference']['annotations']['repeats'],
        segdups = config['reference']['annotations']['segdups'],
        censat = config['reference']['annotations']['censat'],
        microsat = config['reference']['annotations']['microsat'],
        repeatmasker = config['reference']['annotations']['repeatmasker']
    threads: 1
    shell:
        """
        # Prepare column names for the output CSV file
        echo -e "#CHROM\tPOS\tID\trepeatMasker_count\trepeatMasker_poverlap\tsimpleRepeat_count\tsimpleRepeat_poverlap\tgenomicSuperDup_count\tgenomicSuperDup_poverlap\tcentromeres_count\tcentromeres_poverlap\tmicrosat_count\tmicrosat_poverlap" > {output}
        
        # retain only chrom, pos, ID from input gff
        bedtools annotate -names -both -i {input} \
        -files {params.repeatmasker} {params.repeats} {params.segdups} {params.censat} {params.microsat} \
        | cut -f1,4,9,10-19 - >> {output}

        # Add a new column 'origin' containing the sample name
        awk -F'\t' 'BEGIN {{OFS="\t"}} 
        NR==1 {print $0, "origin"}  # Add header
        NR>1 {
            split($3, a, "_");  # Split the third field by underscore
            print $0, a[1];     # Print all fields and add the first part (sample name)
        }' {output} > {output}.tmp
        
        # Replace the original file with the new one
        mv {output}.tmp {output}
        """

rule annotate_genomicSuperDups_bedfile:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_genomicSuperDups.tsv"
    group: "process_mosaic"
    conda:
        "../envs/process_variants.yml"
    params:
        name = 'genomicSuperDups',
        database = config['reference']['annotations']['segdups']
    threads: 1
    shell:
        """
        echo -e "ID\t#CHROM\tstart\tend\tname\tscore\tstrand\toverlap_bp" > {output}
        
        bedtools intersect \
        -a {input} \
        -s \
        -wao \
        -names {params.name} \
        -b {params.database} \
        | cut -f9-16 - >> {output}

        # Add a new column 'source'
        awk -F'\t' 'BEGIN {{OFS="\t"}} 
        NR==1 {{print "source", $0}}
        NR>1 {{print "{params.name}", $0}}' {output} > {output}.tmp
        
        # Replace the original file with the new one
        mv {output}.tmp {output}
        """

use rule annotate_genomicSuperDups_bedfile as annotate_repeatMasker_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_repeatMasker.tsv"
    group: "process_mosaic"
    params:
        name = 'repeatMasker',
        database = config['reference']['annotations']['repeatmasker']

rule annotate_simpleRepeat_bedfile:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_simpleRepeat.tsv"
    group: "process_mosaic"
    params:
        name = 'simpleRepeat',
        database = config['reference']['annotations']['repeats']
    threads: 1
    shell:
        """
        echo -e "ID\t#CHROM\tstart\tend\tname\tscore\toverlap_bp" > {output}

        bedtools intersect \
        -a {input} \
        -wao \
        -names {params.name} \
        -b {params.database} \
        | cut -f9-15 - >> {output}

        # Add a new column 'source'
        awk -F'\t' 'BEGIN {{OFS="\t"}} 
        NR==1 {{print "source", $0}}
        NR>1 {{print "{params.name}", $0}}' {output} > {output}.tmp
        
        # Replace the original file with the new one
        mv {output}.tmp {output}
        """

rule annotate_centromeres_bedfile:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_centromeres.tsv"
    group: "process_mosaic"
    params:
        name = 'centromeres',
        database = config['reference']['annotations']['censat']
    threads: 1
    shell:
        """
        echo -e "ID\t#CHROM\tstart\tend\tname\toverlap_bp" > {output}

        bedtools intersect \
        -a {input} \
        -wao \
        -names {params.name} \
        -b {params.database} \
        | cut -f9-14 - >> {output}

        # Add a new column 'source'
        awk -F'\t' 'BEGIN {{OFS="\t"}} 
        NR==1 {{print "source", $0}}
        NR>1 {{print "{params.name}", $0}}' {output} > {output}.tmp
        
        # Replace the original file with the new one
        mv {output}.tmp {output}
        """

use rule annotate_centromeres_bedfile as annotate_microsat_bedfile with:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_microsat.tsv"
    group: "process_mosaic"
    params:
        name = 'microsat',
        database = config['reference']['annotations']['microsat'],

rule annotate_gencode_features:
    # Annotates a file with information on curated or putative features overlapped in Gencode.
    # Requires strandedness match.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_gencode.tsv"
    group: "process_mosaic"
    conda:
        "../envs/process_variants.yml"
    threads: 1
    params:
        name = 'gencode_v44',
        database = config['reference']['annotations']['gencode']
    shell:
        """
        echo -e "ID\t#CHROM\tsource\ttype\tstart\tend\tscore\tstrand\tframe\tattribute\toverlap_bp" > {output}

        bedtools intersect \
        -a {input} \
        -wao \
        -s \
        -names {params.name} \
        -b {params.database} \
        | cut -f9-19 \
        | awk -F'\t' '$4 == "gene" && $NF > 0' \
        >> {output}

        # Add a new column 'source'
        awk -F'\t' 'BEGIN {{OFS="\t"}} 
        NR==1 {{print "source", $0}} 
        NR>1 {{print "{params.name}", $0}}' {output} > {output}.tmp

        # Replace the original file with the new one
        mv {output}.tmp {output}
        """

rule get_repetitive_features:
    input:
        genomicSuperDups = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_genomicSuperDups.tsv",
        repeatMasker = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_repeatMasker.tsv",
        simpleRepeat = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_simpleRepeat.tsv",
        centromeres = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_centromeres.tsv",
        microsat = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.annotate_microsat.tsv"
    output:
        unfiltered = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.repetitive_features.unfiltered.tsv",
        filtered = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.repetitive_features.tsv"
    group: "process_mosaic"
    threads: 1
    run:
        import pandas as pd

        genomicSuperDup = pd.read_table(input.genomicSuperDups, index_col=False)
        repeatMasker = pd.read_table(input.repeatMasker, index_col=False)
        simpleRepeat = pd.read_table(input.simpleRepeat, index_col=False)
        centromeres = pd.read_table(input.centromeres, index_col=False)
        microsat = pd.read_table(input.microsat, index_col=False)

        # Concatenate all files
        df = pd.concat([genomicSuperDup, repeatMasker, simpleRepeat, centromeres, microsat])
        df['']
        
        # Save all rows (unfiltered) and annotated rows only
        df.to_csv(output.unfiltered, sep='\t', index=False)
        df.query("start != -1").to_csv(output.filtered, sep='\t', index=False)

rule write_alt_fasta:
    # Takes a (sniffles) vcf and writes INS alleles to a fasta file for RepeatMasker. 
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.vcf.gz"
    output: 
        'analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa'
    group: "process_mosaic"
    threads: 1
    shell:
        """
        zcat {input} | grep -v "^#" | grep "SVTYPE=INS" | awk -F'\t' '$5 != "N" && $5 != "<INS>"' | awk -F'\t' '{{print ">"$3"\\n"$5}}' > {output}
        """
        
rule analyze_alt_fasta:
    # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
    # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
    input:
        "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa"
    output:
        expand("analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.{suffix}", 
        allow_missing = True, suffix = ['tbl', 'out', 'out.gff', 'masked', 'cat'])
    group: "process_mosaic"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{file}.filt.log"
    threads: 5
    conda:
        '../envs/RepeatMasker.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} -dir {params.outdir} {input} &> {log}
        """

rule repeatmasker_out_to_tsv:
    # Converts the RepeatMasker outfile to tsv format.
    input:
        script = "scripts/repeatmasker_out_to_tsv.sh",
        out = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.out"
    output:
        out = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.tsv"
    group: "process_mosaic"
    threads: 1
    shell:
        """
        {input.script} {input.out}
        """

rule annotate_repeatmasker_insertions:
    # Annotates a file with information on RepeatMasker-identified insertions.
    input:
        script = "scripts/python/annotate_vcf_repeatmasker.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.filt.vcf.gz",
        tsv = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.tsv"
    output:
        all_insertions = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.all.repeatmasker_insertions.tsv",
        active = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.active.repeatmasker_insertions.tsv",
        single = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.active.single_type.repeatmasker_insertions.tsv",
        multi = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.active.multi_type.repeatmasker_insertions.tsv",
        low_div = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{file}.active.single_type.low_div.repeatmasker_insertions.tsv"
    group: "process_mosaic"
    wildcard_constraints:
        file = "[A-Za-z0-9]+"
    threads: 1
    shell:
        """
        {input.script} {input.vcf} {input.tsv}
        """

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
    group: "process_mosaic"
    conda:
        "../envs/process_variants.yml"
    threads: 1
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