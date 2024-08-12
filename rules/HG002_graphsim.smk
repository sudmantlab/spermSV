rule get_diploid_benchmark_callset:
    # This rule takes in the HG002 structural variant benchmark VCF file, which was created by using dipcall
    # to call SVs on the diploid T2T HG002 assembly against the GRCh38 (hg38) reference genome.
    # It filters the VCF, retaining only contain variants <=10kb that are not BND.
    # Only phased (1|0) or (0|1) variants are retained.
    input:
        vcf = 'benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz'
    output:
        callset = "benchmarks/HG002/sub10kB_diploid_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        index = "benchmarks/HG002/sub10kB_diploid_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p {params.outdir}

        # Get maternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & (GT="1|0" | GT="0|1")' {input.vcf} -o {output.callset} -O z9
        tabix {output.callset}
        """

rule haplograph:
    # Need to run snakemake with --use-singularity --singularity-args "--fakeroot --writable-tmpfs"
    # Figure this out later
    input:
        config = "config/packages/minigraph-cactus/HG002_haplograph.txt"
    output:
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.vcf.gz",
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.gbz",
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.gfa.gz",
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.xg",
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.paf",
        expand("output/assembly/minigraph-cactus/HG002_sim/chrom-alignments/{chr}.vg", chr = chrs)
    container:
        "code/cactus/cactus_v2.8.1.sif"
    threads: 40
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        cactus-pangenome ./js \
        {input.config} \
        --outDir {params.outdir} \
        --outName haplograph \
        --reference hg38 CHM13 \
        --vcfReference hg38 CHM13 \
        --refContigs $(for i in $(seq 22); do printf "chr$i "; done ; echo "chrX chrY chrM") \
        --vcf --giraffe --gfa --gbz --chrom-vg --xg
        """

rule vg_construct:
    # Need to run snakemake with --use-singularity --singularity-args "--fakeroot --writable-tmpfs"
    # Figure this out later
    input:
        callset = "benchmarks/HG002/sub10kB_diploid_GRCh38_HG2-T2TQ100-V1.0.vcf.gz"
    output:
        "output/assembly/vg/HG002_sim/{chromosome}.vg"
    container:
        "code/vg/vg_latest.sif"
    params:
        refgenome = config['reference']['fasta_uncompressed']
    log:
        "logs/assembly/vg/HG002_sim/{chromosome}.log"
    threads: 1
    shell:
        """
        echo "Constructing graph for {wildcards.chromosome}..."
        time vg construct -C -S -R {wildcards.chromosome} -r {params.refgenome} -v {input.callset} -t {threads} -m 32 -p -a > {output} 2>> {log}
        """

rule vg_combine:
    input:
        expand("output/assembly/vg/HG002_sim/{chromosome}.vg", allow_missing = True, chromosome = [f"chr{i}" for i in list(range(1, 23)) + ['X', 'Y']])
    output:
        "output/assembly/vg/HG002_sim/all.vg"
    container:
        "code/vg/vg_latest.sif"
    threads: 1
    shell:
        """
        vg combine {input} > {output}
        """

use rule vg_combine as combine_haplograph with:
    input:
        expand("output/assembly/minigraph-cactus/HG002_sim/chrom-alignments/{chr}.vg", chr = chrs)
    output:
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.vg"
    container:
        "code/vg/vg_latest.sif"
    threads: 1

rule vg_index:
    input:
        vg = "output/assembly/vg/HG002_sim/{name}.vg",
        callset = "benchmarks/HG002/sub10kB_diploid_GRCh38_HG2-T2TQ100-V1.0.vcf.gz"
    output:
        gbwt = "output/assembly/vg/HG002_sim/{name}.gbwt",
        xg = "output/assembly/vg/HG002_sim/{name}.xg"
    threads: 10
    container:
        "code/vg/vg_latest.sif"
    shell:
        """
        vg index -t {threads} -x {output.xg} -G {output.gbwt} -v {input.callset} {input.vg}
        """

rule haplograph_gbwt:
    input:
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.gbz"
    output:
        "output/assembly/minigraph-cactus/HG002_sim/haplograph.gbwt"
    threads: 1
    container:
        "code/vg/vg_latest.sif"
    shell:
        """
        # Extract the GBWT index from a GBZ file.
        vg gbwt  -o {output} -Z {input}
        """

rule hprc_gbwt:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hprc-v1.1-mc-grch38.gbz"
    output:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hprc-v1.1-mc-grch38.gbwt"
    threads: 20
    shell:
        """
        vg gbwt --num-threads {threads} -o {output} -Z {input}
        """

rule sample_fastq:
    input:
        "output/preprocessing/HiFiAdapterFilt/HG002/placeholder_for_HPRC_revio_data/m84039_230117_233243_s1.hifi_reads.default.ccs.filt.fastq.gz",
    output:
        "output/assembly/vg/HG002_sim/temp/sampled.fastq"
    conda:
        "../envs/mapping.yml"
    threads: 1
    params:
        reads = 10000
    shell:
        """
        seqtk sample {input} {params.reads} > {output}
        """
rule vg_sim:
    input:
        sampled = "output/assembly/vg/HG002_sim/temp/sampled.fastq",
        gbwt = "output/assembly/vg/HG002_sim/{name}.gbwt",
        xg = "output/assembly/vg/HG002_sim/{name}.xg"
    output:
        # "output/assembly/vg/HG002_sim/sim/{name}_simulated.gam"
        "output/assembly/vg/HG002_sim/sim/{name}_simulated.fastq"
    container:
        "code/vg/vg_latest.sif"
    params:
        reads = 10, # number of reads
        seed = 42
    threads: 10
    shell:
        """
        # can use an additional parameter to sample from a specific path (haplotype): ie -P HG002#0#chr21#0
        vg sim -t {threads} --random-seed {params.seed} -x {input.xg} -g {input.gbwt} -m HG002 -n {params.reads} -F {input.sampled} > {output}
        """