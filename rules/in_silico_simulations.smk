rule subset_training:
    input:
        "output/preprocessing/uBAMtoFastq/HG002/placeholder_for_HPRC_revio_data/m84039_230117_233243_s1.hifi_reads.default.ccs.fastq.gz"
    output:
        "output/in_silico/badread/training/HG002.fastq.gz"
    threads: 10
    conda:
        "../envs/mapping.yml"
    params:
        seed = 42,
        n = 12000000
    shell:
        """
        seqtk sample -s{params.seed} {input} {params.n} > {output}
        """
        

rule training2paf:
    input:
        "output/in_silico/badread/training/HG002.fastq.gz"
    output:
        "output/in_silico/badread/training/HG002.paf.gz"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        refgenome = config["reference"]["fasta_uncompressed"]
    threads: 14
    shell:
        """
        minimap2 {params.refgenome} {input} -t {threads} -x map-hifi -y -L --eqx -c --cs --MD | gzip > {output}
        """

rule train_error_model:
    input:
        fastq = "output/in_silico/badread/training/HG002.fastq.gz",
        paf = "output/in_silico/badread/training/HG002.paf.gz"
    output:
        "output/in_silico/badread/models/HG002.error.gz"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    threads: 1
    params:
        refgenome = config["reference"]["fasta_uncompressed"]
    log:
        "output/in_silico/badread/logs/HG002.error.log"
    shell:
        """
        badread error_model --reference {params.refgenome} --reads {input.fastq} --alignment {input.paf} --debug | gzip -c > {output} 2>{log}
        """

rule train_qscore_model:
    input:
        fastq = "output/in_silico/badread/training/HG002.fastq.gz",
        paf = "output/in_silico/badread/training/HG002.paf.gz"
    output:
        "output/in_silico/badread/models/HG002.qscore.gz"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    threads: 1
    params:
        refgenome = config["reference"]["fasta_uncompressed"]
    log:
        "output/in_silico/badread/logs/HG002.qscore.log"
    shell:
        """
        badread qscore_model --reference {params.refgenome} --reads {input.fastq} --alignment {input.paf} --debug | gzip -c > {output} 2>{log}
        """

rule generate_dipcall_bed:
    input:
        vcf = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        tbi = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        h1 = "output/in_silico/VISOR/dipcall/dipcall.h1.bed",
        h2 = "output/in_silico/VISOR/dipcall/dipcall.h2.bed",
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        refgenome = config["reference"]["fasta_uncompressed"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    threads: 1
    shell:
        """
        # CMRG: CONTRAC, DUP, SIMPLEDEL, SIMPLEINS, SUBSDEL
        # dipcall: DEL and INS only

        bcftools filter -i "INFO/SVLEN >= 50 & INFO/SVLEN <= 10000" {input.vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN[\t%SAMPLE=%GT]\n' - | grep -w "DEL" | grep "1|0" | grep -v "^#"  | awk 'OFS=FS="\t"''{{print $1, $2, $2+4, "deletion", "None", "0"}}' >> {output.h1}
        bcftools filter -i "INFO/SVLEN >= 50 & INFO/SVLEN <= 10000" {input.vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN[\t%SAMPLE=%GT]\n' - | grep -w "DEL" | grep "0|1" | grep -v "^#"  | awk 'OFS=FS="\t"''{{print $1, $2, $2+4, "deletion", "None", "0"}}' >> {output.h2}
        bcftools filter -i "INFO/SVLEN >= 50 & INFO/SVLEN <= 10000" {input.vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN[\t%SAMPLE=%GT]\t%ALT\n' - | grep -w "INS" | grep "1|0" | grep -v "^#"  | awk 'OFS=FS="\t"''{{print $1, $2, $2+4, "insertion", $6, "0"}}' >> {output.h1}
        bcftools filter -i "INFO/SVLEN >= 50 & INFO/SVLEN <= 10000" {input.vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN[\t%SAMPLE=%GT]\t%ALT\n' - | grep -w "INS" | grep "0|1" | grep -v "^#"  | awk 'OFS=FS="\t"''{{print $1, $2, $2+4, "insertion", $6, "0"}}' >> {output.h2}
        """

rule generate_CMRG_bed:
    input:
        vcf = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        tbi = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi"
    output:
        h1 = "output/in_silico/VISOR/CMRG/CMRG.h1.bed",
        h2 = "output/in_silico/VISOR/CMRG/CMRG.h2.bed",
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        refgenome = config["reference"]["fasta_uncompressed"],
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    threads: 1
    shell:
        """
        # CMRG: CONTRAC, DUP, SIMPLEDEL, SIMPLEINS, SUBSDEL
        # dipcall: DEL and INS only

        bcftools query -f '%CHROM\t%POS\t%INFO/REPTYPE\t%INFO/BREAKSIMLENGTH[\t%SAMPLE=%GT]\n' {input.vcf} | grep "DEL" | grep "1|0" | grep -v "^#"  | awk 'OFS=FS="\t"''{{if ($3 > 50) print $1, $2, $2+4, "deletion", "None", "0"}}' >> {output.h1}
        bcftools query -f '%CHROM\t%POS\t%INFO/REPTYPE\t%INFO/BREAKSIMLENGTH[\t%SAMPLE=%GT]\n' {input.vcf} | grep "DEL" | grep "0|1" | grep -v "^#"  | awk 'OFS=FS="\t"''{{if ($3 > 50) print $1, $2, $2+4, "deletion", "None", "0"}}' >> {output.h2}
        bcftools query -f '%CHROM\t%POS\t%INFO/REPTYPE\t%INFO/BREAKSIMLENGTH[\t%SAMPLE=%GT]\t%ALT\n' {input.vcf} | grep "INS" | grep "1|0" | grep -v "^#"  | awk 'OFS=FS="\t"''{{if ($3 > 50) print $1, $2, $2+4, "insertion", $6, "0"}}' >> {output.h1}
        bcftools query -f '%CHROM\t%POS\t%INFO/REPTYPE\t%INFO/BREAKSIMLENGTH[\t%SAMPLE=%GT]\t%ALT\n' {input.vcf} | grep "INS" | grep "0|1" | grep -v "^#"  | awk 'OFS=FS="\t"''{{if ($3 > 50) print $1, $2, $2+4, "insertion", $6, "0"}}' >> {output.h2}
        """

rule generate_alu_bed:
    input:
        script = "scripts/python/create_hack_bed.py",
        alu_fasta = "output/in_silico/repeats/dfam_AluY_homininae.fasta",
        L1_fasta = "output/in_silico/repeats/dfam_L1_homininae.fasta"
    output:
        "output/in_silico/VISOR/alu/hack_insertions.bed"
    conda:
        "../envs/biopython.yml"
    params:
        maxdims = "output/in_silico/VISOR/hg38.maxdims.tsv",
        alu_count = 1000,
        L1_count = 50,
    threads: 1
    shell:
        """
        python {input.script} \
        --chrom_info {params.maxdim} \
        --alu_count {params.alu_count} \
        --L1_count {params.L1_count} \
        --alu_fasta {input.alu_fasta} \
        --L1_fasta {input.L1_fasta} \
        --output {output} --verbose
        """

rule HACk_benchmark:
    input:
        h1 = "output/in_silico/VISOR/{benchmark}/{benchmark}.h1.bed",
        h2 = "output/in_silico/VISOR/{benchmark}/{benchmark}.h2.bed"
    output:
        "output/in_silico/VISOR/{benchmark}/haps/h1.fa",
        "output/in_silico/VISOR/{benchmark}/haps/h1.fa.fai",
        "output/in_silico/VISOR/{benchmark}/haps/h2.fa",
        "output/in_silico/VISOR/{benchmark}/haps/h2.fa.fai"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        refgenome = config["reference"]["fasta_uncompressed"],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 7
    shell:
        """
        VISOR HACk -g {params.refgenome} -b {input.h1} {input.h2} -o {params.outdir}
        """

rule HACk_alu:
    input:
        "output/in_silico/VISOR/alu/hack_insertions.bed"
    output:
        "output/in_silico/VISOR/alu/haps/h1.fa",
        "output/in_silico/VISOR/alu/haps/h1.fa.fai"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        refgenome = config["reference"]["fasta_uncompressed"],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 7
    shell:
        """
        rm -r {params.outdir} # cleans existing folder to avoid triggering error
        VISOR HACk -g {params.refgenome} -b {input} -o {params.outdir}
        """

rule split_fa:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.fa"
    output:
        expand("output/in_silico/references/hg38/chroms/{chr}/{chr}.fa", chr = chrs)
    params:
        faSplit = "/global/scratch/users/stacy-l/software/ucsc_utilities/faSplit",
        outdir = "output/in_silico/references/hg38/chroms"
    shell:
        """
        mkdir -p {params.outdir}
        {params.faSplit} byname {input} {params.outdir}/

        # Create directories and move files for each chromosome
        for chr in chr{{1..22}} chrX chrY; do
            mkdir -p {params.outdir}/$chr
            mv {params.outdir}/$chr.fa {params.outdir}/$chr/.
        done
        """

rule index_fa:
    input:
        "output/in_silico/{file}.fa"
    output:
        "output/in_silico/{file}.fa.fai"
    conda:
        "../envs/mapping.yml"
    threads: 1
    shell:
        """
        samtools faidx {input}
        """

use rule split_fa as split_synthetic with:
    input:
        "output/in_silico/VISOR/alu/haps/h1.fa"
    output:
        expand("output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa", chr = chrs)
    params:
        faSplit = "/global/scratch/users/stacy-l/software/ucsc_utilities/faSplit",
        outdir = "output/in_silico/VISOR/alu/haps/chroms"

rule generate_benchmark_dims:
    input:
        expand("output/in_silico/VISOR/{benchmark}/haps/{hap}.fa.fai", allow_missing = True, hap = ['h1', 'h2'])
    output:
        haplochroms = "output/in_silico/VISOR/{benchmark}/haplochroms.dim.tsv",
        maxdims = "output/in_silico/VISOR/{benchmark}/maxdims.tsv"
    threads: 1
    shell:
        """
        cat {input} | cut -f1,2 - > {output.haplochroms}
        cat {output.haplochroms} | sort  | awk '$2 > maxvals[$1] {{lines[$1]=$0; maxvals[$1]=$2}} END {{ for (tag in lines) print lines[tag] }}' > {output.maxdims}
        """

use rule generate_benchmark_dims as generate_alu_dims with:
    input:
        "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa.fai"
    output:
        haplochroms = "output/in_silico/VISOR/alu/haps/chroms/{chr}.haplochrom.dim.tsv",
        maxdims = "output/in_silico/VISOR/alu/haps/chroms/{chr}.maxdim.tsv"

rule generate_alu_laser_bed:
    input:
        "output/in_silico/VISOR/alu/haps/chroms/{chr}.maxdim.tsv"
    output:
        "output/in_silico/VISOR/alu/laser/af_{af}/{chr}.bed"
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    params:
        freq = lambda wildcards: float(wildcards.af)
    threads: 1
    shell:
        """
        # assumes uniform spike-in across all chromosomes
        awk 'OFS=FS="\t"''{{print $1, "1", $2, "100.0", "{params.freq}"}}' {input} > {output}
        """

use rule generate_alu_laser_bed as generate_alu_hifi_laser_bed with:
    input:
        "output/in_silico/VISOR/alu/haps/chroms/{chr}.maxdim.tsv"
    output:
        "output/in_silico/VISOR/alu/laser/hifi/af_{af}/{chr}.bed"

use rule generate_alu_laser_bed as generate_alu_newmodel_hifi_laser_bed with:
    input:
        "output/in_silico/VISOR/alu/haps/chroms/{chr}.maxdim.tsv"
    output:
        "output/in_silico/VISOR/alu/laser/newmodel_hifi/af_{af}/{chr}.bed"

rule simulate_alu_hifi_reads:
    input:
        chrom = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa",
        chrom_index = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa.fai",
        ref_chrom = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa",
        ref_index = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa.fai",
        bed = "output/in_silico/VISOR/alu/laser/hifi/af_100/{chr}.bed",
        error_model = "code/Badread/badread/error_models/pacbio2021.gz",
        qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
    output:
        "output/in_silico/VISOR/alu/laser/hifi/af_100/{chr}/sim.srt.bam",
        "output/in_silico/VISOR/alu/laser/hifi/af_100/{chr}/sim.srt.bam.bai"
    params:
        chromdir = lambda wildcards, input: os.path.dirname(input.chrom),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        coverage = 50,
        identity_min = 99,
        identity_max = 100,
        identity_stdev = 0.5,
        length_mean = 16000,
        length_stdev = 2000,
        junk_reads = 0,
        random_reads = 0,
        glitches_rate = 0,
        glitches_size = 0,
        glitches_skip = 0,
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    threads: 4
    log:
        "output/in_silico/VISOR/alu/laser/hifi/af_100/logs/{chr}.log"
    shell:
        """
        rm -r {params.outdir} # prevents existing directory error

        VISOR LASeR \
        -g {input.chrom} \
        -s {params.chromdir} \
        -b {input.bed} \
        -o {params.outdir} \
        --threads {threads} \
        --coverage {params.coverage} \
        --identity_min {params.identity_min} \
        --identity_max {params.identity_max} \
        --identity_stdev {params.identity_stdev} \
        --length_mean {params.length_mean} \
        --length_stdev {params.length_stdev} \
        --junk_reads {params.junk_reads} \
        --random_reads {params.random_reads} \
        --glitches_rate {params.glitches_rate} \
        --glitches_size {params.glitches_size} \
        --glitches_skip {params.glitches_skip} \
        --error_model {input.error_model} \
        --qscore_model {input.qscore_model} \
        --read_type pacbio \
        --tag --fastq --compress > {log}
        """

use rule simulate_alu_hifi_reads as simulate_ref_hifi_reads with:
    input:
        chrom = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa",
        chrom_index = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa.fai",
        ref_chrom = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa",
        ref_index = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa.fai",
        bed = "output/in_silico/VISOR/alu/laser/hifi/af_0/{chr}.bed",
        error_model = "code/Badread/badread/error_models/pacbio2021.gz",
        qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
    output:
        "output/in_silico/VISOR/alu/laser/hifi/af_0/{chr}/sim.srt.bam",
        "output/in_silico/VISOR/alu/laser/hifi/af_0/{chr}/sim.srt.bam.bai"
    log:
        "output/in_silico/VISOR/alu/laser/hifi/af_0/logs/{chr}.log"
    params:
        chromdir = lambda wildcards, input: os.path.dirname(input.chrom),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        coverage = 50,
        identity_min = 99,
        identity_max = 100,
        identity_stdev = 0.5,
        length_mean = 16000,
        length_stdev = 2000,
        junk_reads = 0,
        random_reads = 0,
        glitches_rate = 0,
        glitches_size = 0,
        glitches_skip = 0,

use rule simulate_alu_hifi_reads as simulate_mixed_hifi_reads with:
    input:
        chrom = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa",
        chrom_index = "output/in_silico/VISOR/alu/haps/chroms/{chr}/{chr}.fa.fai",
        ref_chrom = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa",
        ref_index = "output/in_silico/references/hg38/chroms/{chr}/{chr}.fa.fai",
        bed = "output/in_silico/VISOR/alu/laser/hifi/af_8/{chr}.bed",
        error_model = "code/Badread/badread/error_models/pacbio2021.gz",
        qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
    output:
        "output/in_silico/VISOR/alu/laser/hifi/af_8/{chr}/sim.srt.bam",
        "output/in_silico/VISOR/alu/laser/hifi/af_8/{chr}/sim.srt.bam.bai"
    log:
        "output/in_silico/VISOR/alu/laser/hifi/af_8/logs/{chr}.log"
    params:
        chromdir = lambda wildcards, input: os.path.dirname(input.chrom),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        coverage = 50,
        identity_min = 99,
        identity_max = 100,
        identity_stdev = 0.5,
        length_mean = 16000,
        length_stdev = 2000,
        junk_reads = 0,
        random_reads = 0,
        glitches_rate = 0,
        glitches_size = 0,
        glitches_skip = 0,

rule combine_simulated_bams:
    input:
        expand("output/in_silico/VISOR/alu/laser/hifi/{af}/{chr}/sim.srt.bam", allow_missing = True, chr = chrs)
    output:
        "output/in_silico/VISOR/alu/laser/hifi/{af}/all.bam"
    threads: 28
    conda:
        "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
    shell:
        """
        samtools merge -r -@ {threads} --output-fmt='BAM' --write-index {output} {input}
        """

use rule sniffles_mosaic as call_simulated_SVs with:
    input:
        bam = "output/in_silico/VISOR/alu/laser/hifi/{af}/all.bam"
    output:
        vcf='output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.vcf.gz',
        snf='output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.snf',
        tbi='output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        10
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.log"


rule benchmark_simulated_SVs:
    input:
        query = 'output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.vcf.gz',
        query_index = 'output/in_silico/VISOR/alu/laser/hifi/{af}/variants/sniffles_mosaic/callset.vcf.gz.tbi',
        benchmark = '',
        benchmark_index = ''

# rule hg38_full_CMRG_benchmark:
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz",
#         query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz.tbi",
#         jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.jl",
#         benchmark = "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
#         benchmark_index = "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
#         includebed = "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.bed"
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/hg38_full_CMRG_benchmark/{file}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 1
#     params:
#         refgenome = config['reference']['fasta'],
#         outdir = lambda wildcards, output: os.path.dirname(output[0])
#     shell:
#         """
#         # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification?)
#         truvari bench \
#         -f {params.refgenome} \
#         -b {input.benchmark} \
#         -c {input.query} \
#         -o {params.outdir}/bench \
#         -r 1000 \
#         --includebed {input.includebed} \
#         --dup-to-ins \
#         --passonly

#         mv {params.outdir}/bench/* {params.outdir}/.
#         rm -r {params.outdir}/bench
#         """

# rule run_benchmark_laser:
#     # Run in parallel and merge, otherwise performance is prohibitively slow.
#     input:
#         haps = expand("output/in_silico/VISOR/{benchmark}/haps/{haplotype}.fa", allow_missing = True, haplotype = ["h1", "h2"]),
#         bed = "output/in_silico/VISOR/{benchmark}/laser.af_{af}.bed",
#         error_model = "code/Badread/badread/error_models/pacbio2021.gz",
#         qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
#     output:
#         expand("output/in_silico/VISOR/{benchmark}/laser.af_{af}/sim.srt.{ext}", allow_missing = True, ext = ["bam", "bai"])
#     params:
#         hapsdir = lambda wildcards, input: os.path.dirname(input.haps[0]),
#         refgenome = config["reference"]["fasta_uncompressed"],
#         outdir = lambda wildcards, output: os.path.dirname(output[0]),
#         coverage = 50,
#         length_mean = 16000,
#         length_stdev = 2000,
#     conda:
#         "/global/scratch/users/stacy-l/miniconda3/envs/VISOR"
#     threads: 28
#     shell:
#         """
#         VISOR LASeR \
#         -g {params.refgenome} \
#         -s {params.hapsdir} \
#         -b {input.bed} \
#         -o {params.outdir} \
#         --threads {threads} \
#         --coverage {params.coverage} \
#         --length_mean {params.length_mean} \
#         --length_stdev {params.length_stdev} \
#         --error_model {input.error_model} \
#         --qscore_model {input.qscore_model} \
#         --read_type pacbio \
#         --tag --fastq --compress
#         """
        