localrules: repeatmasker_scaffolded_to_bed

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

rule repeatmasker_scaffolded:
    input:
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta",
    output:
        out = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/repeatmasker/{specimen}.diploid.fasta.out"
    log:
        "logs/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta.out.log"
    conda:
        "../envs/RepeatMasker.yml"
    params:
        engine = config['repeatmasker']['engine'],
        species = config['repeatmasker']['species'],
        outdir = "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/repeatmasker"
    threads: 36
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff -species {params.species} -dir {params.outdir} {input} &> {log}
        """

rule repeatmasker_scaffolded_to_bed:
    input:
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta.out"
    output:
        "output/assembly/hifiasm/{specimen}/{ref}_scaffolded/{specimen}.diploid.fasta.bed"
    threads: 1
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}NR>3 {{if($9=="C"){{strand="-"}}else{{strand="+"}}}}{{print $5,$6-1,$7,$10,".",strand}}' {input} > {output}
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