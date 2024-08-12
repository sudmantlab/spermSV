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
    threads: 40
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
    log: "logs/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: 
        "gfatools gfa2fa {input} > {output} 2> {log}"

rule quast_raw:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.hap1.fa",
        "output/assembly/hifiasm/{specimen}/{specimen}.hap2.fa"
    output:
        "output/assembly/hifiasm/{specimen}/quast/raw/report.html"
    conda:
        "../envs/assembly_qc.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        outdir = "output/assembly/hifiasm/{specimen}/quast/raw"
    shell:
        """
        quast.py {input} -r {params.refgenome} \
        -o {params.outdir} \
        --large --no-icarus --fragmented
        """

rule ragtag_scaffold:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        expand("output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.{ext}", allow_missing = True, ext = ["agp", "stats", "asm.paf", "asm.paf.log", "confidence.txt", "err"]),
        temp("output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.temp.fasta")
    conda:
        "../envs/assembly_qc.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        outdir = "output/assembly/hifiasm/{specimen}/scaffolded/{hap}"
    shell:
        """
        mkdir -p {params.outdir}
        ragtag.py scaffold {params.refgenome} {input} \
        -u --mm2-params='-t {threads} -x asm5 --eqx --cs' \
        -o {params.outdir}

        for file in {params.outdir}/ragtag.scaffold.*; do
            new_name="${{file/ragtag.scaffold/{wildcards.specimen}.{wildcards.hap}}}"
            mv "$file" "${{new_name/scaffolded\/{wildcards.hap}\//}}"
        done
        """

rule rename_scaffolded_ctgs:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.temp.fasta"
    output:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.fasta"
    shell:
        """
        sed 's/>\(.*\)_RagTag/>\\1_RagTag_{wildcards.hap}/' {input} > {output}
        """

use rule quast_raw as quast_scaffolded with:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.hap1.scaffold.fasta",
        "output/assembly/hifiasm/{specimen}/{specimen}.hap2.scaffold.fasta"
    output:
        "output/assembly/hifiasm/{specimen}/quast/scaffolded/report.html"
    conda:
        "../envs/assembly_qc.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        outdir = "output/assembly/hifiasm/{specimen}/quast/scaffolded"

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