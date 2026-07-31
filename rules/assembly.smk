rule hifiasm:
    input: get_fastqs_per_sample
    output:
        p_ctg_hap1 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.noseq.gfa",
        ec = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ec.bin",
        ovlp_reverse = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/assembly/hifiasm/{specimen}/{specimen}.asm"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 14
    log: 
        "logs/assembly/hifiasm/{specimen}/{specimen}.asm.log"
    conda: 
        "../envs/HiFiAssembly.yml"
    shell: 
        """
        hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1
        """

rule gfaToFa:
    input: "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.{hap}.p_ctg.gfa"
    output: "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    log: "logs/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa.log"
    threads: 1
    conda: "../envs/HiFiAssembly.yml"
    shell: 
        "gfatools gfa2fa {input} > {output} 2> {log}"

rule ragtag_scaffold:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/assembly_qc.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta_uncompressed'],
        outdir = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}"
    shell:
        """
        mkdir -p {params.outdir}

        ragtag.py scaffold {params.refgenome} \
        {input} \
        -u -w --aligner minimap2 -t {threads} \
        -o {params.outdir}

        # Rename the ragtag output files to include specimen and hap info
        find {params.outdir} -name "ragtag.scaffold.*" \
        -exec sh -c 'for f do dir=$(dirname "$f"); \
        base=$(basename "$f"); suffix=${{base#ragtag.scaffold.}}; \
        mv "$f" "$dir/{wildcards.specimen}.{wildcards.hap}.scaffold.$suffix"; \
        done' sh {{}} +

        # Modify all FASTA headers to include hap information
        sed -i 's/_RagTag$/_RagTag_{wildcards.hap}/' {output.fasta}
        """

rule update_scaffold_agp:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.agp"
    output:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.updated.agp"
    shell:
        """
        awk -F'\t' '{{sub(/_RagTag$/, "_RagTag_{wildcards.hap}", $1);
        print}}' OFS='\t' {input} > {output}
        """


rule cat_scaffolds:
    input:
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/hap1/{specimen}.hap1.scaffold.fasta",
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/hap2/{specimen}.hap2.scaffold.fasta"
    output:
        fa = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta.fai"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 1
    conda:
        "../envs/mapping.yml"
    shell:
        """
        cat {input} > {output.fa}
        samtools faidx {output.fa}
        """

### QC: QUAST (assembly completeness) ###

rule quast_raw:
    # Optional: Run QUAST on the raw assembly fasta without scaffolding.
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.hap1.fa",
        "output/assembly/hifiasm/{specimen}/{specimen}.hap2.fa"
    output:
        "output/assembly/hifiasm/{specimen}/quast/raw/report.html"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/assembly_qc.yml"
    threads: 6
    params:
        outdir = "output/assembly/hifiasm/{specimen}/quast/raw"
    shell:
        """
        quast.py {input} \
        -o {params.outdir} \
        --large --est-ref-size 3100000000 --no-icarus
        """

use rule quast_raw as quast_scaffolded with:
    input:
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/hap1/{specimen}.hap1.scaffold.fasta",
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/hap2/{specimen}.hap2.scaffold.fasta"
    output:
        "output/assembly/hifiasm/{specimen}/quast/{ref}_scaffolded/report.html"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/assembly_qc.yml"
    threads: 6
    params:
        outdir = "output/assembly/hifiasm/{specimen}/quast/{ref}_scaffolded"

### QC: Assembly statistics ###
# TODO: This is dependent on Joana's assembly stats pipeline and other hard-coded paths and should be removed in the final release.

rule run_assembly_stats:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        'output/assembly/assembly_stats/{specimen}.{hap}.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename {wildcards.specimen}.{wildcards.hap} \
        --genomesize {params.genomesize}
        """

rule run_reference_stats:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.fa"
    output:
        'output/assembly/assembly_stats/hg38.ref.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename hg38.ref \
        --genomesize {params.genomesize}
        """

rule run_reference_hap_stats:
    input:
        "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta"
    output:
        'output/assembly/assembly_stats/T2T-HG002.{hap}.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename T2T-HG002.{wildcards.hap} \
        --genomesize {params.genomesize}
        """

rule combine_all_stats:
    input:
        expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = [x for x in specimens if x != '901'], hap = ['hap1', 'hap2'])
        # expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = glob_wildcards("output/assembly/assembly_stats/{specimen}.hap1.tsv").specimen, hap = ['hap1', 'hap2'])
    output:
        'output/assembly/assembly_stats/all.tsv'
    run:
        import pandas as pd
        import glob
        import os

        # Get all matching tsv files at runtime
        input_pattern = "output/assembly/assembly_stats/*.*.tsv"
        all_files = glob.glob(input_pattern)
        
        if not all_files:
            raise ValueError(f"No files found matching the pattern: {input_pattern}")

        # Read and concatenate all files
        df_list = []
        for file in all_files:
            try:
                df = pd.read_csv(file, sep='\t')
                df_list.append(df)
            except pd.errors.EmptyDataError:
                print(f"Warning: {file} is empty and will be skipped.")
        
        if not df_list:
            raise ValueError("No valid data found in any of the input files.")
        
        combined_df = pd.concat(df_list, ignore_index=True)
        
        # Write the combined dataframe to the output file
        combined_df.to_csv(output[0], sep='\t', index=False)

### QC: Flagger (assembly error estimation) ###

rule prepare_flagger_inputs:
    # Flagger requires a BED file that describes the assembly contigs + coordinates.
    input:
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta.fai"
    output:
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.bed})
        
        # Create whole genome bed from fai
        cat {input.fai} | awk '{{print $1"\t0\t"$2}}' > {output.bed}
        
        # Create annotations JSON
        echo '{{"whole_genome": "{output.bed}"}}' > {output.json}
        """

use rule prepare_flagger_inputs as HG002_unscaffolded_prepare with:
    # Use the unscaffolded version of the HG002 assembly as a reference/control.
    input:
        fasta = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa",
        fai = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa.fai"
    output:
        bed = "output/assembly/flagger/HG002/unscaffolded/HG002_whole_genome.bed",
        json = "output/assembly/flagger/HG002/unscaffolded/HG002_annotations_path.json"

rule flagger_bam_to_coverage:
    input:
        bam = "output/alignment/hg38_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam",
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    output:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            bam2cov \
                --bam {input.bam} \
                --output {output.cov} \
                --annotationJson {input.json} \
                --threads {threads} \
                --baselineAnnotation whole_genome
        """

rule run_flagger:
    input:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz",
        alpha_tsv = "config/packages/flagger/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"
    output:
        "output/assembly/flagger/{specimen}/prediction_summary_final.tsv"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """        
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            hmm_flagger \
                --input {input.cov} \
                --outputDir $(dirname {output}) \
                --alphaTsv {input.alpha_tsv} \
                --labelNames Err,Dup,Hap,Col \
                --threads {threads}
        """

rule summarize_flagger:
    input:
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = [x for x in specimens if x != '900'])
    output:
        "output/assembly/flagger/all_prediction_summary_final.tsv"
    threads: 1
    shell:
        """
        echo -e "specimen\tErr\tDup\tHap\tCol\tUnk" > {output}
        find output/assembly/flagger -name "prediction_summary_final.tsv" | while read file; do
            specimen=$(echo "$file" | awk -F'/' '{{print $(NF-1)}}')
            awk -v specimen="$specimen" '$1=="PREDICTION" && $2=="base_level" && $3=="percentage" && $4=="annotation" && $5=="whole_genome" {{print specimen "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}}' "$file" >> {output}
        done
        """

### Optional: RepeatMasker annotation of individual assemblies ###

rule filter_canonical_chromosomes:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta",
        temp_chroms = temp("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/canonical_chroms.txt")
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 1
    conda:
        "../envs/mapping.yml"
    shell:
        """
        samtools faidx {input}
        
        grep "^chr" {input}.fai | grep -v -E "random|chrUn|_h[12]tg" | cut -f1 > {output.temp_chroms}
        
        samtools faidx {input} -r {output.temp_chroms} > {output.fasta}
        """

rule filter_hg38_canonical_chromosomes:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.fa"
    output:
        fasta = "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.fa",
        temp_chroms = temp("/global/scratch/users/stacy-l/references/hg38_HGSVC/canonical_chroms.txt")
    threads: 1
    conda:
        "../envs/mapping.yml"
    shell:
        """
        samtools faidx {input}
        
        awk '$1 ~ /^chr[1-9]$/ || $1 ~ /^chr[1-2][0-9]$/ || $1 ~ /^chr[XYM]$/' {input}.fai | cut -f1 > {output.temp_chroms}
        
        if [ ! -s {output.temp_chroms} ]; then
            echo "No chromosomes matched the pattern" >&2
            exit 1
        fi
        
        samtools faidx {input} -r {output.temp_chroms} > {output.fasta}
        """

rule split_canonical_fasta:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta"
    output:
        expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa", allow_missing = True, chr = chrs) # for all chrs
    conda:
        "../envs/mapping.yml"
    params:
        outdir = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/split_fastas"
    shell:
        """
        mkdir -p {params.outdir}
        awk -v outdir="{params.outdir}" '
        /^>/ {{
            if (file) {{
                close(file)
            }}
            filename = substr($0, 2)
            sub(/_RagTag.*$/, "", filename)
            file = outdir "/" filename ".fa"
            print $0 > file
            next
        }}
        {{ if (file) print > file }}' {input}
        """


rule repeatmasker_per_chr:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/split_fastas/{chr}.fa"
    output:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/per_chr/{chr}.fa.out"
    log:
        "logs/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/{chr}.log"
    conda:
        "../envs/RepeatMasker.yml"
    params:
        engine = config['repeatmasker']['engine'],
        species = config['repeatmasker']['species'],
        outdir = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/per_chr"
    threads: 4
    resources:
        mem_mb = 24000
    shell:
        """
        RepeatMasker -pa {threads} \
            -engine {params.engine} \
            -nocut -gff \
            -species {params.species} \
            -dir {params.outdir} \
            {input} &> {log}
        """

rule repeatmasker_to_bed:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/per_chr/{chr}.fa.out"
    output:
        temp("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/per_chr/{chr}.bed")
    threads: 1
    shell:
        """
        tail -n +4 {input} | \
        awk 'BEGIN{{OFS="\t"}} 
        {{
            print $5, $6-1, $7, $10"#"$11, $1, $9
        }}' > {output}
        """

rule combine_repeatmasker_beds:
    input:
        expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/per_chr/{chr}.bed", allow_missing = True, chr = chrs, hap = ['hap1', 'hap2'])
    output:
        bed = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz",
        tbi = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz.tbi"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/mapping.yml"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | bgzip > {output.bed}
        tabix -p bed {output.bed}
        """

# annotation transfer

rule reverse_hg38_mapping:
    input:
        assembly = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.canonical.fasta",
        ref = "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.fa"
    output:
        "output/assembly/hifiasm/{specimen}/reverse_chain/hg38_scaffolded/{specimen}.{hap}.paf"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    threads: 10
    conda:
        "../envs/mapping.yml"
    shell:
        """
        # Reverse the order of mapping: reference genome is mapped to the assembly
        minimap2 -t {threads} -x asm20 -c --eqx --secondary=no {input.assembly} {input.ref} > {output}
        """

rule reverse_hg38_chain:
    input:
        "output/assembly/hifiasm/{specimen}/reverse_chain/hg38_scaffolded/{specimen}.{hap}.paf"
    output:
        "output/assembly/hifiasm/{specimen}/reverse_chain/hg38_scaffolded/{specimen}.{hap}.chain"
    threads: 1
    shell:
        """
        code/paf2chain/target/release/paf2chain -i {input} > {output}
        """

rule paf2chain:
    input:
        "output/assembly/hifiasm/{specimen}/{scaffolded}/{hap}/{specimen}.{hap}.scaffold.asm.paf"
    output:
        temp("output/assembly/hifiasm/{specimen}/{scaffolded}/{hap}/{specimen}.{hap}.scaffold.chain")
    threads: 1
    shell:
        """
        code/paf2chain/target/release/paf2chain -i {input} > {output}
        """

rule update_chain:
    input:
        script = "scripts/python/update_chain.py",
        chain = "output/assembly/hifiasm/{specimen}/{scaffolded}/{hap}/{specimen}.{hap}.scaffold.chain",
        agp = "output/assembly/hifiasm/{specimen}/{scaffolded}/{hap}/{specimen}.{hap}.scaffold.updated.agp"
    output:
        chain = "output/assembly/hifiasm/{specimen}/{scaffolded}/{hap}/{specimen}.{hap}.scaffold.updated.chain"
    threads: 1
    shell:
        """
        python {input.script} {input.agp} {input.chain} {output.chain}
        """

# rule liftover_hg38_annotations:
#     input:
#         bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/{file}.bed",
#         chain = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.updated.chain"
#     output:
#         mapped = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.bed.gz",
#         mapped_tbi = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.bed.gz.tbi",
#         unmapped = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.unmapped.bed.gz",
#         unmapped_tbi = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.unmapped.bed.gz.tbi"
#     conda:
#         "../envs/liftover.yml"
#     shell:
#         """
#         # Create temp uncompressed bed file paths
#         mapped_bed=$(echo {output.mapped} | sed 's/\.gz$//')
#         unmapped_bed=$(echo {output.unmapped} | sed 's/\.gz$//')

#         # Do not consider beyond first 12 columns of bed spec
#         cols=$(head -n1 {input.bed} | awk '{{print NF}}')
#         if [ "$cols" -gt 12 ]; then
#             liftOver {input.bed} {input.chain} "$mapped_bed" "$unmapped_bed" -bedPlus=12
#         else
#             liftOver {input.bed} {input.chain} "$mapped_bed" "$unmapped_bed"
#         fi

#         echo "Sorting bedfiles"
#         sort -k 1,1 -k2,2n "$mapped_bed" -o "$mapped_bed" 
#         sort -k 1,1 -k2,2n "$unmapped_bed" -o "$unmapped_bed"

#         echo "Gzipping bedfiles"
#         bgzip --force "$mapped_bed"
#         bgzip --force "$unmapped_bed"

#         echo "Indexing bedfiles"
#         tabix {output.mapped}
#         tabix {output.unmapped}
#         """

rule reverse_liftover_alu_coordinates:
    input:
        chain = "output/assembly/hifiasm/{specimen}/reverse_chain/hg38_scaffolded/{specimen}.{hap}.chain",
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/alu_positions.bed"
    output:
        mapped = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/alu_positions.reverse_GRCh38.bed",
        unmapped = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/alu_positions.reverse_GRCh38.unmapped.bed"
    conda:
        "../envs/liftover.yml"
    threads: 1
    shell:
        """
        liftOver {input.positions} {input.chain} {output.mapped} {output.unmapped}
        """

rule summarize_reverse_liftover_alu_coordinates:
    input:
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/alu_positions.reverse_GRCh38.bed", specimen = ['0068', '5619', '5980', '895', 'ASD6058', 'RG6P', 'TSC6830','0EL2', '5621', '894', '899', 'ASD6463', 'SY3L014', 'TSC7237'], hap = ['hap1', 'hap2'])
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output}
        """

rule intersect_reverse_hg38_centromeres:
    input:
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed",
        bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_centromeres.bed"
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.centromeres.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.positions} -b {input.bed} -wb > {output}
        """

rule intersect_reverse_hg38_segdups:
    input:
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed",
        bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_genomicSuperDups.bed"
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.genomicSuperDups.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.positions} -b {input.bed} -wb > {output}
        """

rule intersect_reverse_hg38_microsat:
    input:
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed",
        bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_microsat.bed"
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.microsat.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.positions} -b {input.bed} -wb > {output}
        """

rule intersect_reverse_hg38_gencode:
    input:
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed",
        bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_gencodeV48.bed"
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.gencodeV48.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.positions} -b {input.bed} -wb > {output}
        """

rule intersect_de_novo_repeatmasker:
    input:
        positions = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/alu_positions.bed",
        bed =  "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz",
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/intersected_repeatmasker.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.positions} -b {input.bed} -wb > {output}
        """

rule summarize_de_novo_repeatmasker:
    input:
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/{specimen}/{hap}/intersected_repeatmasker.bed", specimen = ['0068', '5619', '5980', '895', 'ASD6058', 'RG6P', 'TSC6830','0EL2', '5621', '894', '899', 'ASD6463', 'SY3L014', 'TSC7237'], hap = ['hap1', 'hap2'])
    output:
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.de_novo_repeatmasker.bed"
    threads: 1
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n > {output}
        """