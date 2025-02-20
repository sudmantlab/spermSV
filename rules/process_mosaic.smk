rule preprocess_variants:
    # Takes in the VCF output of sniffles mosaic and preprocesses with the following steps:
    # 1) Prepends the ID field with the sample name (ex. Sniffles2.BND.F38S0 -> 894_Sniffles2.BND.F38S0).
    # 2) Filters out variants with extreme coverage or read support values.
    # 3) Filters out variants on uncertain contigs (_random, chrUn_)
    # 4) Filters out excessive SV lengths (greater than 10000)
    # The ID value passed on to downstream analyses will retain the sample origin.
    input:
        vcf = 'output/alignment/hg38/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz',
        problematic = '/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_UCSCUnusualRegions.bed',
        centromeres = '/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_centromeres.bed'
    output:
        filt = 'output/alignment/hg38/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/hg38/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    wildcard_constraints:
        refalias = '[A-Za-z0-9]+'
    conda:
        "../envs/process_variants.yml"
    threads: 1
    params:
        region_format = lambda wildcards: "chr${i},",
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1
    shell:
        """
        regions=$(for i in {{1..22}} X Y; do echo -n {params.region_format}; done | sed 's/,$//')
        bcftools view -r $regions {input.vcf} | \
        bcftools annotate --set-id '{wildcards.specimen}\_%ID' - | \
        bcftools filter --soft-filter centromeric --mask-file {input.centromeres} - | \
        bcftools filter --soft-filter problematic --mask-file {input.problematic} - | \
        bcftools filter -i 'INFO/SVLEN <= 10000 & {params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' - \
        -o {output.filt} --write-index
        """

rule preprocess_scaffolded_variants:
    # Preprocesses as above, but relaxes the coverage bounds and removes support read limits due to diploid mapping.
    input:
        vcf = 'output/alignment/{ref}_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/{ref}_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
        index = 'output/alignment/{ref}_scaffolded/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz.csi'
    wildcard_constraints:
        refalias = '[A-Za-z0-9]+'
    conda:
        "../envs/process_variants.yml"
    threads: 1
    params:
        region_format = lambda wildcards: "chr${i}_RagTag_hap1,chr${i}_RagTag_hap2,",
        max_cov = 60,
        min_cov = 20,
        max_freq = 0.1
    shell:
        """
        regions=$(for i in {{1..22}} X Y; do echo -n {params.region_format}; done | sed 's/,$//')
        bcftools view -r $regions {input.vcf} | \
        bcftools annotate --set-id '{wildcards.specimen}\_%ID' - | \
        bcftools filter - \
        -i 'INFO/SVLEN <= 10000 & {params.max_cov} > (DR + DV) & (DR + DV) > {params.min_cov} & (DV <= {params.max_freq} * (DR + DV))' \
        -o {output.filt} --write-index
        """

use rule preprocess_scaffolded_variants as preprocess_pangenome_hg38_variants with:
    input:
        vcf = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.filt.vcf.gz',
        index = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.filt.vcf.gz.csi'
    params:
        region_format = lambda wildcards: "GRCh38#0#chr${i},",
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1

use rule preprocess_scaffolded_variants as preprocess_pangenome_CHM13_variants with:
    input:
        vcf = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.vcf.gz'
    output:
        filt = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.filt.vcf.gz',
        index = 'output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.filt.vcf.gz.csi'
    params:
        region_format = lambda wildcards: "CHM13#0#chr${i},",
        max_cov = 80,
        min_cov = 40,
        max_freq = 0.1

rule join_filtered_calls:
    # Joins together filtered calls from all available filtered VCFs.
    input:
        expand('output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
               allow_missing = True, specimen = specimens)
    output:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-all.filt.vcf.gz",
        index = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-all.filt.vcf.gz.csi",
    wildcard_constraints:
        refalias = '[A-Za-z0-9]+'
    conda:
        "../envs/process_variants.yml"
    threads: 1
    shell:
        """
        bcftools concat --threads {threads} -a {input} -o {output.vcf} --write-index
        """

use rule join_filtered_calls as join_filtered_calls_scaffolded with:
    # Joins together filtered calls from all available filtered VCFs.
    input:
        expand('output/alignment/{ref}_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.filt.vcf.gz',
               allow_missing = True, specimen = [x for x in specimens if x != '900'])
    output:
        vcf = "analysis/{ref}_scaffolded/denovo/files/minimap2/standard/variants/{ref}_scaffolded-all.filt.vcf.gz",
        index = "analysis/{ref}_scaffolded/denovo/files/minimap2/standard/variants/{ref}_scaffolded-all.filt.vcf.gz.csi",

use rule join_filtered_calls as join_filtered_calls_pangenome with:
    input:
        expand('output/alignment/hprc_personalized/variants/hprc-v1.1-mc-chm13.d9/sniffles_mosaic/{ref}/{specimen}.filt.vcf.gz',
               allow_missing = True, ref = ['hg38', 'CHM13'], specimen = [x for x in specimens if x != '900'])
    output:
        vcf = "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.filt.vcf.gz",
        index = "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.filt.vcf.gz.csi",

use rule join_filtered_calls as join_qc_all_filtered_calls with:
    input:
        expand('output/alignment/{refalias}/{mapper}/{setting}/variants/sniffles_mosaic/{specimen}.qc_all.filt.vcf.gz',
               allow_missing = True, specimen = specimens)
    output:
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-all.qc_all.filt.vcf.gz",
        index = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-all.qc_all.filt.vcf.gz.csi",

use rule join_filtered_calls as join_qc_all_filtered_calls_scaffolded with:
    # For some reason, this breaks with the process_variants.yml env but not bcftools.yml??
    input:
        expand('output/alignment/{ref}_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.filt.vcf.gz',
               allow_missing = True, specimen = [x for x in specimens if x != '900'])
    output:
        vcf = "analysis/{ref}_scaffolded/denovo/files/minimap2/standard/variants/{ref}_scaffolded-all.qc_all.filt.vcf.gz",
        index = "analysis/{ref}_scaffolded/denovo/files/minimap2/standard/variants/{ref}_scaffolded-all.qc_all.filt.vcf.gz.csi",

use rule join_filtered_calls as join_qc_all_filtered_calls_pangenome with:
    input:
        expand('output/alignment/hprc_personalized/variants/hprc-v1.1-mc-chm13.d9/sniffles_mosaic/{ref}/{specimen}.qc_all.filt.vcf.gz',
               allow_missing = True, ref = ['hg38', 'CHM13'], specimen = [x for x in specimens if x != '900'])
    output:
        vcf = "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.qc_all.filt.vcf.gz",
        index = "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.qc_all.filt.vcf.gz.csi",

rule vcf_to_gff3:
    # Converts an existing Sniffles VCF to a GFF3 format.
    input:
        script = "scripts/python/vcf_to_gff3.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.vcf.gz"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.gff"
    threads: 1
    shell:
        """
        python {input.script} {input.vcf} {output}
        """

rule overlap_repetitive:
    # Creates a file with the overlap count and percentage overlap of repetitive or duplicated features 
    # specified in the below bedfiles.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.overlap_repetitive.tsv"
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
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_genomicSuperDups.tsv"
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
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_repeatMasker.tsv"
    params:
        name = 'repeatMasker',
        database = config['reference']['annotations']['repeatmasker']

rule annotate_simpleRepeat_bedfile:
    # Annotates a file with information on repetitive features overlapped in the bedfile.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_simpleRepeat.tsv"
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
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_centromeres.tsv"
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
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_microsat.tsv"
    params:
        name = 'microsat',
        database = config['reference']['annotations']['microsat'],

rule annotate_gencode_features:
    # Annotates a file with information on curated or putative features overlapped in Gencode.
    # Requires strandedness match.
    input:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.gff"
    output:
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_gencode.tsv"
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
        genomicSuperDups = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_genomicSuperDups.tsv",
        repeatMasker = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_repeatMasker.tsv",
        simpleRepeat = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_simpleRepeat.tsv",
        centromeres = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_centromeres.tsv",
        microsat = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.annotate_microsat.tsv"
    output:
        unfiltered = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.repetitive_features.unfiltered.tsv",
        filtered = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.repetitive_features.tsv"
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
        "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.vcf.gz"
    output: 
        'analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa'
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
        "analysis/{refalias}/{analysis}/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.out"
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
    threads: 1
    shell:
        """
        {input.script} {input.out}
        """

rule annotate_repeatmasker_insertions:
    # Annotates a file with information on RepeatMasker-identified insertions.
    input:
        script = "scripts/python/annotate_vcf_repeatmasker.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.filt.vcf.gz",
        tsv = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.filt.alt.fa.tsv"
    output:
        all_insertions = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.all.repeatmasker_insertions.tsv",
        active = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.active.repeatmasker_insertions.tsv",
        single = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.active.single_type.repeatmasker_insertions.tsv",
        multi = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.active.multi_type.repeatmasker_insertions.tsv",
        low_div = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.active.single_type.low_div.repeatmasker_insertions.tsv"
    wildcard_constraints:
        file = "[A-Za-z0-9]+"
    threads: 1
    shell:
        """
        {input.script} {input.vcf} {input.tsv}
        """

use rule annotate_repeatmasker_insertions as annotate_repeatmasker_insertions_qc_all with:
    input:
        script = "scripts/python/annotate_vcf_repeatmasker.py",
        vcf = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.filt.vcf.gz",
        tsv = "analysis/{refalias}/denovo/files/{mapper}/{setting}/repeatmasker/{file}.qc_all.filt.alt.fa.tsv"
    output:
        all_insertions = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.all.repeatmasker_insertions.tsv",
        active = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.active.repeatmasker_insertions.tsv",
        single = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.active.single_type.repeatmasker_insertions.tsv",
        multi = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.active.multi_type.repeatmasker_insertions.tsv",
        low_div = "analysis/{refalias}/denovo/files/{mapper}/{setting}/variants/{refalias}-{file}.qc_all.active.single_type.low_div.repeatmasker_insertions.tsv"

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