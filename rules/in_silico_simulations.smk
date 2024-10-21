rule split_index_ref:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.fa"
    output:
        directory("output/in_silico/references/hg38")
    conda:
        "../envs/mapping.yml"
    params:
        faSplit = "/global/scratch/users/stacy-l/software/ucsc_utilities/faSplit"
    shell"
        """
        {params.faSplit} byname {input} {output}/
        """


rule subsample_training:
    input:
        "output/alignment/hg38/minimap2/standard/mapped/HG002.sorted.merged.bam"
    output:
        "output/in_silico/badread/subsampled_HG002.bam"
    conda:
        "../envs/mapping.yml"
    threads: 5
    params:
        seed = 42
    shell:
        """
        samtools view -b -@ {threads} \
        --subsample 0.05 \
        --subsample-seed {params.seed} \
        -o {output} \
        {input}
        """

rule training2fastq:
    input:
        "output/in_silico/badread/subsampled_HG002.bam"
    output:
        "output/in_silico/badread/subsampled_HG002.fastq.gz"
    threads: 5
    conda:
        "../envs/mapping.yml"
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T '*' {input} -0 {output}
        """

rule training2paf:
    input:
        "output/in_silico/badread/subsampled_HG002.fastq.gz"
    output:
        "output/in_silico/badread/subsampled_HG002.paf.gz"
    conda:
        "VISOR"
    params:
        refgenome = config["reference"]["fasta"].strip('.gz')
    threads: 20
    shell:
        """
        minimap2 {params.refgenome} {input} -t {threads} -x map-hifi -y -L --eqx -c --cs --MD | gzip > {output}
        """

rule train_error_model:
    input:
        fastq = "output/in_silico/badread/subsampled_HG002.fastq.gz",
        paf = "output/in_silico/badread/subsampled_HG002.paf.gz"
    output:
        "output/in_silico/badread/subsampled_HG002.error_model"
    conda:
        "VISOR"
    threads: 1
    params:
        refgenome = config["reference"]["fasta"].strip('.gz')
    log:
        "output/in_silico/badread/subsampled_HG002.error_model.log"
    shell:
        """
        badread error_model --reference {params.refgenome} --reads {input.fastq} --alignment {input.paf} --debug > {output} 2>{log}
        """

rule train_qscore_model:
    input:
        fastq = "output/in_silico/badread/subsampled_HG002.fastq.gz",
        paf = "output/in_silico/badread/subsampled_HG002.paf.gz"
    output:
        "output/in_silico/badread/subsampled_HG002.qscore_model"
    conda:
        "VISOR"
    threads: 1
    params:
        refgenome = config["reference"]["fasta"].strip('.gz')
    log:
        "output/in_silico/badread/subsampled_HG002.qscore_model.log"
    shell:
        """
        badread qscore_model --reference {params.refgenome} --reads {input.fastq} --alignment {input.paf} --debug > {output} 2>{log}
        """

rule generate_dipcall_bed:
    input:
        # vcf = "ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz",
        # tbi = "ALL.wgs.integrated_sv_map_v2_GRCh38.20130502.svs.genotypes.vcf.gz.tbi",
        vcf = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        tbi = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        h1 = "output/in_silico/VISOR/dipcall/dipcall.h1.bed",
        h2 = "output/in_silico/VISOR/dipcall/dipcall.h2.bed",
    conda:
        "VISOR"
    params:
        refgenome = config["reference"]["fasta"].strip('.gz'),
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
        "VISOR"
    params:
        refgenome = config["reference"]["fasta"].strip('.gz'),
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
        alu_count = 10,
        L1_count = 1,
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
        "VISOR"
    params:
        refgenome = config["reference"]["fasta"].strip('.gz'),
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
        "VISOR"
    params:
        refgenome = config["reference"]["fasta"].strip('.gz'),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 7
    shell:
        """
        rm -r {params.outdir} # cleans existing folder to avoid triggering error
        VISOR HACk -g {params.refgenome} -b {input} -o {params.outdir}
        """

rule generate_benchmark_dims:
    input:
        expand("output/in_silico/VISOR/{benchmark}/haps/{hap}.fa.fai", allow_missing = True, hap = ['h1', 'h2'])
    output:
        haplochroms = "output/in_silico/VISOR/{benchmark}/haplochroms.dim.tsv",
        maxdims = "output/in_silico/VISOR/{benchmark}/maxdims.tsv"
    conda:
        "VISOR"
    params:
        hapsdir = lambda wildcards, input: os.path.dirname(input[0])
    threads: 1
    shell:
        """
        cat {input} | cut -f1,2 - > {output.haplochroms}
        cat {output.haplochroms} | sort  | awk '$2 > maxvals[$1] {{lines[$1]=$0; maxvals[$1]=$2}} END {{ for (tag in lines) print lines[tag] }}' > {output.maxdims}
        """

rule generate_alu_dims:
    input:
        "output/in_silico/VISOR/alu/haps/h1.fa.fai"
    output:
        haplochroms = "output/in_silico/VISOR/{benchmark}/haplochroms.dim.tsv",
        maxdims = "output/in_silico/VISOR/{benchmark}/maxdims.tsv"
    conda:
        "VISOR"
    params:
        hapsdir = lambda wildcards, input: os.path.dirname(input[0])
    threads: 1
    shell:
        """
        cat {input} | cut -f1,2 - > {output.haplochroms}
        cat {output.haplochroms} | sort  | awk '$2 > maxvals[$1] {{lines[$1]=$0; maxvals[$1]=$2}} END {{ for (tag in lines) print lines[tag] }}' > {output.maxdims}
        """

rule generate_laser_bed:
    input:
        "output/in_silico/VISOR/{benchmark}/maxdims.tsv"
    output:
        "output/in_silico/VISOR/{benchmark}/laser.af_{af}.bed"
    conda:
        "VISOR"
    params:
        normal = lambda wildcards: 100 - float(wildcards.af)
    threads: 1
    shell:
        """
        # assumes uniform spike-in across all chromosomes
        awk 'OFS=FS="\t"''{{print $1, "1", $2, "100.0", "{params.normal}"}}' {input} > {output}
        """

rule run_alu_laser:
    input:
        haps = "output/in_silico/VISOR/alu/haps/h1.fa",
        bed = "output/in_silico/VISOR/alu/laser.af_6.67.bed",
        error_model = "code/Badread/badread/error_models/pacbio2021.gz",
        qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
    output:
        expand("output/in_silico/VISOR/alu/laser.af_6.67/sim.srt.{ext}", allow_missing = True, ext = ["bam", "bai"])
    params:
        hapsdir = lambda wildcards, input: os.path.dirname(input.haps),
        refgenome = config["reference"]["fasta"].strip('.gz'),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        coverage = 2,
        length_mean = 16000,
        length_stdev = 2000,
    conda:
        "VISOR"
    threads: 56
    shell:
        """
        rm -r {params.outdir} # prevents flagged error
        VISOR LASeR \
        -g {params.refgenome} \
        -s {params.hapsdir} \
        -b {input.bed} \
        -o {params.outdir} \
        --threads {threads} \
        --coverage {params.coverage} \
        --length_mean {params.length_mean} \
        --length_stdev {params.length_stdev} \
        --error_model {input.error_model} \
        --qscore_model {input.qscore_model} \
        --read_type pacbio \
        --tag --fastq --compress
        """

rule run_benchmark_laser:
    # Run in parallel and merge, otherwise performance is prohibitively slow.
    input:
        haps = expand("output/in_silico/VISOR/{benchmark}/haps/{haplotype}.fa", allow_missing = True, haplotype = ["h1", "h2"]),
        bed = "output/in_silico/VISOR/{benchmark}/laser.af_{af}.bed",
        error_model = "code/Badread/badread/error_models/pacbio2021.gz",
        qscore_model = "code/Badread/badread/qscore_models/pacbio2021.gz"
    output:
        expand("output/in_silico/VISOR/{benchmark}/laser.af_{af}/sim.srt.{ext}", allow_missing = True, ext = ["bam", "bai"])
    params:
        hapsdir = lambda wildcards, input: os.path.dirname(input.haps[0]),
        refgenome = config["reference"]["fasta"].strip('.gz'),
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        coverage = 50,
        length_mean = 16000,
        length_stdev = 2000,
    conda:
        "VISOR"
    threads: 28
    shell:
        """
        VISOR LASeR \
        -g {params.refgenome} \
        -s {params.hapsdir} \
        -b {input.bed} \
        -o {params.outdir} \
        --threads {threads} \
        --coverage {params.coverage} \
        --length_mean {params.length_mean} \
        --length_stdev {params.length_stdev} \
        --error_model {input.error_model} \
        --qscore_model {input.qscore_model} \
        --read_type pacbio \
        --tag --fastq --compress
        """
        