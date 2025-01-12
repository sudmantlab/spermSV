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


rule minigraph_cactus_pangenome:
    # Need to run snakemake with --use-singularity --singularity-args "--fakeroot --writable-tmpfs"
    # Figure this out later
    input:
        config = "config/packages/minigraph-cactus/samples.txt"
    output:
        "output/assembly/minigraph-cactus/svg.vcf.gz",
        expand("output/assembly/minigraph-cactus/chrom-alignments/{chr}.vg", chr = chrs)
    container:
        "code/cactus/cactus_v2.8.1.sif"
    threads: 40
    shell:
        """
        cactus-pangenome ./js \
        {input.config} \
        --outDir output/assembly/minigraph-cactus \
        --outName svg \
        --reference CHM13 hg38 \
        --vcfReference CHM13 hg38 \
        --refContigs $(for i in $(seq 22); do printf "chr$i "; done ; echo "chrX chrY chrM") \
        --vcf --giraffe --gfa --gbz --vg
        """