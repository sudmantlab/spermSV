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

use rule hifiasm as hifiasm_900 with:
    # a fix for sample 900, which fails from OOM issues from too much data(?)
    input:
        "output/preprocessing/uBAMtoFastq/900/PBmixRevio1489_2_D01_PFTP_30hours_22kbExpressCCSv32hrPE_200pM_HumanSudmant8900_bc2026_CCSExpressIndex/m84139_231209_045828_s4.hifi_reads.bc2026.ccs.fastq.gz",
        "output/preprocessing/uBAMtoFastq/900/PBmixRevio1489_2_C01_PFTP_30hours_22kbExpressCCSv32hrPE_200pM_HumanSudmant8900_bc2026_CCSExpressIndex/m84139_231209_042722_s3.hifi_reads.bc2026.ccs.fastq.gz"
    output: 
        p_ctg_hap1 = "output/assembly/hifiasm/900/900.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/assembly/hifiasm/900/900.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/assembly/hifiasm/900/900.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/assembly/hifiasm/900/900.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/assembly/hifiasm/900/900.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/assembly/hifiasm/900/900.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/assembly/hifiasm/900/900.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/assembly/hifiasm/900/900.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/assembly/hifiasm/900/900.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/assembly/hifiasm/900/900.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/assembly/hifiasm/900/900.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/assembly/hifiasm/900/900.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/assembly/hifiasm/900/900.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/assembly/hifiasm/900/900.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/assembly/hifiasm/900/900.asm.bp.r_utg.noseq.gfa",
        ec = "output/assembly/hifiasm/900/900.asm.ec.bin",
        ovlp_reverse = "output/assembly/hifiasm/900/900.asm.ovlp.reverse.bin",
        ovlp_source = "output/assembly/hifiasm/900/900.asm.ovlp.source.bin"
    params:
        prefix = "output/assembly/hifiasm/900/900.asm"
    log: 
        "logs/assembly/hifiasm/900/900.asm.log"

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

rule quast_raw:
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

rule ragtag_hg38_scaffold:
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

use rule ragtag_hg38_scaffold as ragtag_T2T_scaffold with:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        fasta = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    conda:
        "../envs/assembly_qc.yml"
    threads: 10
    params:
        refgenome = "/global/scratch/users/stacy-l/references/T2T_CHM13/hs1.fa",
        outdir = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{hap}"

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

rule liftover_hg38_annotations:
    input:
        bed = "/global/scratch/users/stacy-l/references/hg38_HGSVC/{file}.bed",
        chain = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.updated.chain"
    output:
        mapped = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.bed.gz",
        mapped_tbi = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.bed.gz.tbi",
        unmapped = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.unmapped.bed.gz",
        unmapped_tbi = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{file}.unmapped.bed.gz.tbi"
    conda:
        "../envs/liftover.yml"
    shell:
        """
        # Create temp uncompressed bed file paths
        mapped_bed=$(echo {output.mapped} | sed 's/\.gz$//')
        unmapped_bed=$(echo {output.unmapped} | sed 's/\.gz$//')

        # Do not consider beyond first 12 columns of bed spec
        cols=$(head -n1 {input.bed} | awk '{{print NF}}')
        if [ "$cols" -gt 12 ]; then
            liftOver {input.bed} {input.chain} "$mapped_bed" "$unmapped_bed" -bedPlus=12
        else
            liftOver {input.bed} {input.chain} "$mapped_bed" "$unmapped_bed"
        fi

        echo "Sorting bedfiles"
        sort -k 1,1 -k2,2n "$mapped_bed" -o "$mapped_bed" 
        sort -k 1,1 -k2,2n "$unmapped_bed" -o "$unmapped_bed"

        echo "Gzipping bedfiles"
        bgzip --force "$mapped_bed"
        bgzip --force "$unmapped_bed"

        echo "Indexing bedfiles"
        tabix {output.mapped}
        tabix {output.unmapped}
        """

use rule liftover_hg38_annotations as liftover_T2T_annotations with:
    input:
        bed = "/global/scratch/users/stacy-l/references/T2T_CHM13/{file}.bed",
        chain = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{hap}/{specimen}.{hap}.scaffold.updated.chain"
    output:
        mapped = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{hap}/{file}.bed",
        unmapped = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{hap}/{file}.unmapped.bed"

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

rule svbyeye_alignment:
    input:
        "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/{specimen}.{hap}.scaffold.fasta"
    output:
        "output/assembly/hifiasm/{specimen}/svbyeye/hg38_scaffolded/{specimen}.{hap}.paf"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+"
    params:
        refgenome = config['reference']['fasta']
    threads: 10
    conda:
        "../envs/mapping.yml"
    shell:
        """
        minimap2 -t {threads} -x asm20 -c --eqx --secondary=no {params.refgenome} {input} > {output}
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
        
        # more precise matching
        awk '$1 ~ /^chr[1-9]$/ || $1 ~ /^chr[1-2][0-9]$/ || $1 ~ /^chr[XYM]$/' {input}.fai | cut -f1 > {output.temp_chroms}
        
        # Check if we got any matches
        if [ ! -s {output.temp_chroms} ]; then
            echo "No chromosomes matched the pattern" >&2
            exit 1
        fi
        
        samtools faidx {input} -r {output.temp_chroms} > {output.fasta}
        """

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

rule create_ref_dict:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.fa"
    output:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.dict"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/gatk"
    shell:
        """
        gatk CreateSequenceDictionary -R {input}
        """

rule liftover_hg38_scaffolded:
    input:
        vcf = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.vcf.gz",
        chain = "output/assembly/hifiasm/{specimen}/reverse_chain/hg38_scaffolded/{specimen}.{hap}.chain",
        ref = "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.fa",
        ref_dict = "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.canonical.dict"
    output:
        vcf = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.{hap}.hg38.vcf.gz",
        reject = "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.{hap}.rejected.vcf.gz"
    log:
        "logs/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/liftover/{specimen}.{hap}.liftover.log"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/gatk"
    params:
        java_opts = "-Xmx24g"
    shell:
        """
        gatk --java-options "{params.java_opts}" LiftoverVcf \
            -I {input.vcf} \
            -O {output.vcf} \
            -C {input.chain} \
            -R {input.ref} \
            --REJECT {output.reject} \
            --RECOVER_SWAPPED_REF_ALT true \
            --WARN_ON_MISSING_CONTIG true \
            --MAX_RECORDS_IN_RAM 100000 &> {log}
        """