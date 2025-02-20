localrules: convert_hg38_repeat_bed_hprc, convert_CHM13_repeat_bed_hprc

rule hprc_personalized_gbwt:
    input:
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz"
    output:
        ri = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.ri",
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/gbwt/{prefix}.gbwt.log"
    params:
        container = config["singularity"]["minigraph-cactus"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["minigraph-cactus"]["bind_paths"]]),
        singularity_args = config["singularity"]["minigraph-cactus"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home $(dirname {input.gbz}) \
        {params.container} \
        vg gbwt --num-threads {threads} \
        -p -r {output.ri} -Z {input.gbz}) 2> {log}
        """

rule hprc_personalized_haplotype_sampling:
    input:
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz",
        ri = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.ri",
    output:
        haplotypes = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.hapl"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/haplotypes/{prefix}.haplotypes.log"
    params:
        container = config["singularity"]["minigraph-cactus"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["minigraph-cactus"]["bind_paths"]]),
        singularity_args = config["singularity"]["minigraph-cactus"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home $(dirname {input.gbz}) \
        {params.container} \
        vg haplotypes --threads {threads} \
        -v 2 -H {output.haplotypes} \
        -r {input.ri} {input.gbz}) 2> {log}
        """

def get_fastas_per_sample(wildcards, sample_table = config['sample_table']) -> List:
    """
    A helper function for creating a list of fastq.gz paths for each sample ahead of personalized pangenome construction.
    TODO: The base path should maybe be set in the configfile.
    """
    basepath = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fasta.gz"
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [basepath.format(specimen=s[0], group = s[1], lane=s[2], smrtcell = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

rule convert_fastq_to_fasta:
    input:
        fastq = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    output:
        fasta = temp("output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fasta.gz")
    threads: 16
    shell:
        """
        pigz -dc {input} | parallel --pipe --keep-order -j {threads} --block 250M 'sed -n "1~4s/^@/>/p;2~4p"' | pigz -p {threads} > {output}
        """

rule hprc_personalized_kmer_counting:
    input:
        get_fastas_per_sample
    output:
        fa_list = temp("output/alignment/hprc_personalized/kmers/{specimen}.fa.list"),
        kff = "output/alignment/hprc_personalized/kmers/{specimen}.kff"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/kmers/{specimen}.kmer.log"
    params:
        kmc = config["paths"]["kmc"]["path"],
        outpath = lambda wildcards, output: os.path.splitext(output.kff)[0],
        tmpdir = "output/alignment/hprc_personalized/kmers/{specimen}_tmp",
        mem = 152 # note that this is set via config as well but needs to be used as an arg here
    threads: 16
    shell:
        """
        mkdir -p {params.tmpdir}
        echo "{input}" | tr ' ' '\\n' > {output.fa_list}
        {params.kmc} -k29 -m{params.mem} -sm -t{threads} -v -okff -fm @{output.fa_list} {params.outpath} {params.tmpdir} &> {log}
        rm -r {params.tmpdir}
        """

rule hprc_personalized_graph_construction:
    input:
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz",
        haplotypes = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.hapl",
        kmer = "output/alignment/hprc_personalized/kmers/{specimen}.kff"
    output:
        gbz = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.gbz"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/graphs/{prefix}/{specimen}.graph.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["minigraph-cactus"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["minigraph-cactus"]["bind_paths"]]),
        singularity_args = config["singularity"]["minigraph-cactus"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg haplotypes \
        --threads {threads} \
        -v 2 --include-reference --diploid-sampling \
        -i {input.haplotypes} -k {input.kmer} \
        -g {output.gbz} {input.gbz}) 2> {log}
        """

rule hprc_personalized_graph_index:
    input:
        # download from HPRC
        gbz = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.gbz"
    output:
        "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.dist",
        "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.longread.zipcodes", 
        "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.longread.withzip.min"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/graphs/{prefix}/{specimen}.index.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"],
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home}/$(dirname {input.gbz}) \
        {params.container} \
        vg autoindex --workflow lr-giraffe \
        --threads {threads} \
        --prefix {wildcards.specimen} \
        --gbz $(basename {input.gbz})) 2> {log}
        """

rule hprc_personalized_haplotype_mapping:
    input:
        gbz = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.gbz",
        fastq = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz",
        indices = ["output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.dist",
        "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.longread.zipcodes", 
        "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.longread.withzip.min"]
    output:
        gam = temp("output/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.gam"),
        report = "logs/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.report.tsv"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg giraffe -b hifi \
        -Z {input.gbz} \
        -f {input.fastq} \
        --threads {threads} \
        --report-name {output.report} \
        -p > {output.gam}) 2> {log}
        """

def get_temp_gams_per_sample(wildcards, sample_table = config['sample_table']) -> List[str]:
    """
    A helper function for creating a list of expected per-smrtcell gam outputs per sample to be merged into a single gam.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    basepath = "output/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.gam"
    input_samples = [basepath.format(specimen=s[0], prefix = wildcards.prefix, lane=s[2], smrtcell=s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    return input_samples

rule hprc_personalized_merge_gams:
    input:
        get_temp_gams_per_sample
    output:
        "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.mapped.gam"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    threads: 1
    shell:
        "cat {input} > {output}"

rule hprc_personalized_surject:
    input:
        gam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.mapped.gam",
        gbz = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.gbz"
    output:
        bam = temp("output/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}.surjected.bam")
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    log:
        "logs/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 18
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg surject -b \
        --threads {threads} \
        -x {input.gbz} \
        {input.gam} \
        --prune-low-cplx > {output.bam}) 2> {log}
        """

rule hprc_surjected_bam_sort:
    input:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/temp/{specimen}.surjected.bam"
    output:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam",
        bai = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam.bai"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    threads: 8
    conda:
        "../envs/mapping.yml"
    shell:
        """
        samtools sort -@ {threads} {input.bam} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule hprc_surject_ref_fasta:
    input:
        gbz = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.gbz"
    output:
        fasta = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.{ref}.fasta"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    threads: 8
    log:
        "logs/alignment/hprc_personalized/graphs/{prefix}/{specimen}.{ref}.fasta.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg paths --extract-fasta -x {input.gbz} --paths-by {wildcards.ref} --threads {threads}> {output.fasta}) > {output.fasta} 2> {log}
        """

rule convert_hg38_repeat_bed_hprc:
    # Quick conversion of GRCh38 chr names (chr{n}) to graph-compatible names (GRCh38#0#chr{n})
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/GRCh38_simpleRepeat.bed"
    output:
        "/global/scratch/users/stacy-l/references/HPRC/GRCh38_simpleRepeat.bed"
    shell:
        """
        sed 's/chr/GRCh38#0#chr/' {input} > {output}
        """

rule convert_CHM13_repeat_bed_hprc:
    # Quick conversion of CHM13 chr names (chr{n}) to graph-compatible names (CHM13#0#chr{n})
    input:
        "/global/scratch/users/stacy-l/references/T2T_CHM13/chm13v2.0_simpleRepeat.bed"
    output:
        "/global/scratch/users/stacy-l/references/HPRC/CHM13_simpleRepeat.bed"
    shell:
        """
        sed 's/chr/CHM13#0#chr/' {input} > {output}
        """

rule hprc_sniffles_mosaic_hg38:
    input:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam",
        index = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam.bai",
        fasta = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.GRCh38.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HPRC/GRCh38_simpleRepeat.bed"
    output:
        vcf = "output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.vcf.gz"
    log:
        "logs/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.log"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    threads: 10
    conda:
        "../envs/sniffles-dev.yml"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --reference {input.fasta} \
        --tandem-repeats {input.repeats} \
        --minsupport 0 \
        --mapq 20 \
        --output-rnames \
        --mosaic \
        --mosaic-af-min 0 \
        --mosaic-af-max 0.2 \
        --mosaic-qc-strand False \
        --threads {threads} \
        --dev-no-qc &> {log}
        """

use rule hprc_sniffles_mosaic_hg38 as hprc_sniffles_mosaic_CHM13 with:
    input:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam",
        index = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam.bai",
        fasta = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.CHM13.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HPRC/CHM13_simpleRepeat.bed"
    output:
        vcf = "output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.vcf.gz"
    log:
        "logs/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.log"

rule hprc_sniffles_mosaic_hg38_qc_all:
    input:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam",
        index = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam.bai",
        fasta = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.GRCh38.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HPRC/GRCh38_simpleRepeat.bed"
    output:
        vcf = "output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.qc_all.vcf.gz"
    log:
        "logs/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/hg38/{specimen}.qc_all.log"
    wildcard_constraints:
        specimen = "[A-Za-z0-9]+",
        prefix = "[A-Za-z0-9\-\.]+",
    threads: 10
    conda:
        "../envs/sniffles-dev.yml"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --reference {input.fasta} \
        --tandem-repeats {input.repeats} \
        --minsupport 0 \
        --mapq 20 \
        --output-rnames \
        --mosaic \
        --mosaic-af-min 0 \
        --mosaic-af-max 0.2 \
        --mosaic-qc-strand False \
        --threads {threads} \
        --dev-no-qc &> {log}
        """

use rule hprc_sniffles_mosaic_hg38_qc_all as hprc_sniffles_mosaic_CHM13_qc_all with:
    input:
        bam = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam",
        index = "output/alignment/hprc_personalized/mapped/{prefix}/{specimen}.surjected.bam.bai",
        fasta = "output/alignment/hprc_personalized/graphs/{prefix}/{specimen}.CHM13.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HPRC/CHM13_simpleRepeat.bed"
    output:
        vcf = "output/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.qc_all.vcf.gz"
    log:
        "logs/alignment/hprc_personalized/variants/{prefix}/sniffles_mosaic/CHM13/{specimen}.qc_all.log"

# rule chrom_graphs:
#     # Need to run snakemake with --use-singularity --singularity-args "--fakeroot --writable-tmpfs"
#     # Figure this out later
#     input:
#         config = "config/packages/minigraph-cactus/samples.txt"
#     output:
#         expand("output/assembly/minigraph-cactus/svg.{exts}", exts = ["vcf.gz", "gbz", "gfa.gz", "xg", "paf"]),
#         expand("output/assembly/minigraph-cactus/chrom-alignments/{chr}.vg", chr = chrs)
#     log:
#         "logs/assembly/minigraph-cactus/chrom-alignments/minigraph-cactus.log"
#     params:
#         home = config["workdir"],
#         container = config["singularity"]["minigraph-cactus"]["container"],
#         bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["minigraph-cactus"]["bind_paths"]]),
#         singularity_args = config["singularity"]["minigraph-cactus"]["args"]
#     threads: 56
#     shell:
#         """
#         (singularity exec {params.singularity_args} \
#         {params.bind_paths} \
#         --home $(dirname {input.gbz}) \
#         {params.container} \
#         cactus-pangenome ./js \
#         {input.config} \
#         --outDir output/assembly/minigraph-cactus \
#         --outName svg \
#         --reference CHM13 hg38 \
#         --vcfReference CHM13 hg38 \
#         --refContigs $(for i in $(seq 22); do printf "chr$i "; done ; echo "chrX chrY chrM") \
#         --chrom-vg --maxCores 56 --indexCores 32 --mapCores 8 --alignCores 16 \
#         --vcf --giraffe --gfa --gbz --chrom-vg --xg) 2> {log}
#         """

# rule combine_chrom_graphs:
#     input:
#         expand("output/assembly/minigraph-cactus/chrom-alignments/{chr}.vg", chr = chrs)
#     output:
#         "output/assembly/minigraph-cactus/svg.vg"
#     log:
#         "logs/assembly/minigraph-cactus/chrom-alignments/combine_chrom_graphs.log"
#     params:
#         home = config["workdir"],
#         container = config["singularity"]["vg"]["container"],
#         bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
#         singularity_args = config["singularity"]["vg"]["args"]
#     shell:
#         """
#         (singularity exec {params.singularity_args} \
#         {params.bind_paths} \
#         --home $(dirname {input.gbz}) \
#         {params.container} \
#         vg combine {input} > {output}) 2> {log}
#         """