### shared rules ###
rule split_HG002_reference:
    # This rule separates the diploid T2T HG002 assembly into two haplotype assemblies.
    # The two relevant haplotypes in all downstream operations are MATERNAL and PATERNAL.
    # Each haplotype assembly contains chr1-chr22 autosomes, suffixed with _{hap}.
    # The maternal assembly contains chrX_MATERNAL, and the paternal assembly contains chrY_PATERNAL.
    input:
        "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz"
    output:
        fasta = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta.gz",
        index = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta.gz.fai"
    conda: "../envs/mapping.yml"
    shell:
        """
        zcat {input} | grep '>.*{wildcards.hap}' | sed 's/>//' | seqtk subseq {input} - | bgzip > {output.fasta}
        samtools faidx {output.fasta}
        """

use rule minimap2 as diploid_self_mapping with:
    # Maps the HG002 HiFi reads to the HG002 T2T assembly.
    # The output BAM goes through sorting and merging operations (not shown in this file), ultimately
    # producing the file (output/alignment/HG002/minimap2/standard/mapped/self/diploid/HG002.sorted.merged.bam)
    # used as input for the extract_hap_reads rule.
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: 
        temp("output/alignment/HG002/minimap2/standard/mapped/unsorted/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

rule split_self_mapped_hap:
    # This rule takes in the mapped HiFi reads and splits them into two separate BAM files,
    # one for each haplotype.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/HG002.sorted.merged.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.bam"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools view -H {input} | grep '^@SQ' |awk '$2 ~ /{wildcards.hap}/' | cut -f 2 | sed 's/SN://' | xargs samtools view -@ {threads} --with-header -b {input} -o {output}
        """

rule extract_hap_reads:
    # This rule converts the hap-aligned reads in the BAM and converts them to a FASTQ format.
    # This step is used to phase the reads.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.bam"
    output:
        fastq = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.fastq.gz"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T '*' {input.bam} -0 {output.fastq}
        """

rule create_hap_read_name_lists:
    # This rule creates a list of read names mapped to each haplotype.
    input:
        maternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/MATERNAL.bam",
        paternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/PATERNAL.bam"
    output:
        maternal_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/MATERNAL_rnames.txt",
        paternal_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/PATERNAL_rnames.txt"
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools view -@ {threads} {input.maternal_bam} | cut -f1 | sort | uniq > {output.maternal_reads}
        samtools view -@ {threads} {input.paternal_bam} | cut -f1 | sort | uniq > {output.paternal_reads}
        """

rule get_hap_benchmark_SVs:
    # This rule takes in the HG002 structural variant benchmark VCF file, which was created by using dipcall
    # to call SVs on the diploid T2T HG002 assembly against the GRCh38 (hg38) reference genome.
    # It filters the VCF, retaining only variants that are not INV or BND.
    # It creates two VCFs, one for each haplotype, containing only the variants that are exclusive to that haplotype.
    # It accomplishes this by filtering for variants that are specifically phased to one of the haplotypes.
    input:
        vcf = 'benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz'
    output:
        maternal = "benchmarks/HG002/MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        maternal_index = "benchmarks/HG002/MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        paternal = "benchmarks/HG002/PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        paternal_index = "benchmarks/HG002/PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p {params.outdir}

        # Get maternal callset
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & GT=="1|0"' {input.vcf} -o {output.maternal} -O z9
        tabix {output.maternal}

        # Get paternal callset
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & GT=="0|1"' {input.vcf} -o {output.paternal} -O z9
        tabix {output.paternal}
        """

rule get_hap_sub10kB_benchmark_SVs:
    # This rule takes in the HG002 structural variant benchmark VCF file, which was created by using dipcall
    # to call SVs on the diploid T2T HG002 assembly against the GRCh38 (hg38) reference genome.
    # It filters the VCF, retaining only contain variants <=10kb that are not INV or BND.
    # It then creates two VCFs, one for each haplotype, containing only the variants that are exclusive to that haplotype.
    # It accomplishes this by filtering for variants that are specifically phased to one of the haplotypes.
    input:
        vcf = 'benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz'
    output:
        maternal = "benchmarks/HG002/sub10kB_MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        maternal_index = "benchmarks/HG002/sub10kB_MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        paternal = "benchmarks/HG002/sub10kB_PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        paternal_index = "benchmarks/HG002/sub10kB_PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p {params.outdir}

        # Get maternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & SVTYPE!="INV" & GT=="1|0"' {input.vcf} -o {output.maternal} -O z9
        tabix {output.maternal}

        # Get paternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & SVTYPE!="INV" & GT=="0|1"' {input.vcf} -o {output.paternal} -O z9
        tabix {output.paternal}
        """

use rule sniffles_standard as hap_hg38_germline_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap} phased reads mapped to hg38.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats'],
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/{hap}.log"

rule truvari_hg38_germline_all:
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.jl',
        benchmark = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/{hap}/all/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --dup-to-ins \
        --passonly
        """

use rule truvari_hg38_germline_all as truvari_hg38_germline_sub10kB with:
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.jl',
        benchmark = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/{hap}/sub10kB/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),


rule extract_spanning_sv_reads:
    # This rule takes in the VCF file for a single haplotype, the BAM file of HiFi reads mapped to the haplotype,
    # and the list of read names mapped to the haplotype.
    # It then extracts the reads that completely span the SVs and writes them to a FASTQ file.
    # It also produces a bedfile of regions corresponding to the input vcf coordinates for use in downstream
    # benchmarking.
    # For each read:
    # Calculates the read's span based on its CIGAR string.
    # Checks if the read fully spans any SV in the same chromosome.
    # If a read spans an SV completely, it's written to the output FASTQ file.
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        bam = "output/alignment/hg38/minimap2/standard/mapped/HG002.sorted.merged.bam",
        haplotype_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}_rnames.txt",
        script = "scripts/extract_spanning_sv_reads.sh"
    output:
        fastq = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.fastq.gz",
        bed = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.bed"
    log:
        "logs/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.log"
    conda: "../envs/truvari.yml"
    threads: 5
    shell:
        """
        bash {input.script} {input.vcf} {input.bam} {input.haplotype_reads} {output.fastq} {threads} {wildcards.hap} > {log} 2>&1
        """

use rule samtools_sort as generic_sort with:
    # This rule sorts any BAM file in the designated directory.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/{subdirs}/unsorted/{filename}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/{subdirs}/{filename}.bam"
    wildcard_constraints:
        subdirs = "(?!.*unsorted)(?!.*merged).*"
    conda: "../envs/mapping.yml"
    threads: 5

use rule index_bam as generic_index with:
    # This rule creates an index for any BAM file in the designated directory.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam.bai"
    conda: "../envs/mapping.yml"
    threads: 5

rule generic_samtools_coverage:
    # This rule produces a text/ASCII-based summary of coverage on each chr for any BAM file in the designated
    # directory.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam"
    output:
        report = "output/alignment/HG002/minimap2/standard/mapped/{filename}.coverage.txt"
    conda: "../envs/mapping.yml"
    threads: 1
    shell:
        """
        samtools coverage -o {output.report} {input.bam}
        """

rule vcf2df:
    # This rule converts vcf files to joblib files, which can be read into pandas DataFrames.
    # TODO: Note that this parsing allows skip levels across subdirectories.
    input:
        'output/alignment/HG002/minimap2/standard/variants/{filename}.vcf.gz'
    output:
        'output/alignment/HG002/minimap2/standard/variants/{filename}.jl'
    conda:
        '../envs/truvari.yml'
    threads: 1
    shell:
        """
        truvari vcf2df --info --format {input} {output}
        """

rule diploid_assembly_coverage:
    input:
        bams = expand("output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.bam", hap = ['MATERNAL', 'PATERNAL', 'HG002.sorted.merged']),
        indices = expand("output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.bam.bai", hap = ['MATERNAL', 'PATERNAL', 'HG002.sorted.merged'])
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/coverage.html"
    conda: "../envs/deeptools.yml"
    threads: 20
    params:
        format = "plotly",
        title = "'Phased read alignment coverage to diploid assembly'"
    shell:
        """
        plotCoverage -p {threads} --bamfiles {input.bams} --plotFile {output} --plotFileFormat {params.format} -n 1000000 --plotTitle {params.title} --ignoreDuplicates --minMappingQuality 10 
        """

rule self_assembly_coverage:
    input:
        hap_bams = expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam", allow_missing = True, hap1 = ['MATERNAL', 'PATERNAL']),
        sv_bams = expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/1.0_{hap1}_SVs_to_{hap2}.bam", allow_missing = True, hap1 = ['MATERNAL', 'PATERNAL']),
        hap_indices = expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam.bai", allow_missing = True, hap1 = ['MATERNAL', 'PATERNAL']),
        sv_indices = expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/1.0_{hap1}_SVs_to_{hap2}.bam.bai", allow_missing = True, hap1 = ['MATERNAL', 'PATERNAL']),
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/coverage.html"
    conda: "../envs/deeptools.yml"
    threads: 20
    params:
        format = "plotly",
        title = "'Phased read alignment coverage to self-assembly'"
    shell:
        """
        plotCoverage -p {threads} --bamfiles {input.hap_bams} {input.sv_bams} --plotFile {output} --plotFileFormat {params.format} -n 1000000 --plotTitle {params.title} --ignoreDuplicates --minMappingQuality 10 
        """

use rule self_assembly_coverage as hg38_coverage with:
    input:
        hap_bams = expand("output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam", hap = ['MATERNAL', 'PATERNAL']),
        sv_bams = expand("output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.bam", hap = ['MATERNAL', 'PATERNAL']),
        hap_indices = expand("output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam.bai", hap = ['MATERNAL', 'PATERNAL']),
        sv_indices = expand("output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}_SVs.bam.bai", hap = ['MATERNAL', 'PATERNAL'])
    output:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/coverage.html"
    conda: "../envs/deeptools.yml"
    threads: 20
    params:
        format = "plotly",
        title = "'Phased read alignment coverage to hg38'"

### T2T HG002 SIM ###
use rule minimap2 as self_assembly_mapping with:
    # This rule maps HiFi reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    # For example, this rule can map MATERNAL phased reads to either the MATERNAL assembly or the PATERNAL assembly.
    # We expect optimal mapping for matched haplotype mapping and worse mapping for cross-haplotype mapping.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap1}.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/unsorted/{hap1}_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{hap1}_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule minimap2 as self_assembly_SV_mapping with:
    # This rule maps the <=10kB SV-spanning reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/unsorted/{hap1}_SVs_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{hap1}_SVs_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

rule merge_SV_all:
    # This rule takes all hap1 <=10kB SV-spanning reads mapped to hap2 and merges them with all hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/1.0_{hap1}_SVs_to_{hap2}.bam",
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input.bam1} {input.bam2}
        """

rule merge_SV_spike:
    # This rule samples a fraction of the hap1 <=10kB SV-spanning reads mapped to hap2 and merges them with all hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_SVs_to_{hap2}.bam"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '0\.[0-9]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools view -bs {wildcards.spike} {input.bam1} |
        samtools merge -@ {threads} {output} - {input.bam2}
        """

rule merge_all:
    # This rule samples a fraction of *all* hap1 reads mapped to hap2 and merges them with a (1-fraction) amount of hap2 reads mapped to hap2.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/subsampled.{spike}_all_{hap1}_to_{hap2}.bam"),
        merged = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_all_to_{hap2}.bam"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
    conda: "../envs/mapping.yml"
    threads: 10
    params:
        bulk = lambda wildcards: str(1 - float(wildcards.spike))
    shell:
        """
        samtools view -@ {threads} -bs {params.bulk} {input.bam2} -o {output.temp}
        samtools view -@ {threads} -bs {wildcards.spike} {input.bam1} | samtools merge -@ {threads} {output.merged} - {output.temp}
        """

use rule sniffles_standard as call_self_assembly_germline with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all hap1 reads mapped to hap2.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz.tbi'
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        # simpleRepeat file created using grep filter on RepeatMasker track bed, not directly output from TRF
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.log"

use rule sniffles_standard as call_self_assembly_germline_spike with:
    # This rule uses the standard Sniffles germline calling mode.
    # The minimum number of reads required to call a variant is dynamically determined from coverage at the region of interest.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads mapped to hap2.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz.tbi'
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '[01]\.[0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        # simpleRepeat file created using grep filter on RepeatMasker track bed, not directly output from TRF
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.log"

use rule sniffles_mosaic as self_assembly_mosaic_spike with:
    # This rule uses the mosaic (low allele frequency, AF) variant calling mode of Sniffles.
    # The minimum number of reads required to call a mosaic variant is set to 1.
    # The minimum AF required to call a mosaic variant is set to 0.
    # The maximum AF for a mosaic variant is 0.2.
    # Germline variants (exceeding 0.2 AF) are not called.
    # It calls SVs from a BAM containing some {spike} amount of the hap1 SV-spanning reads mapped to hap2.
    # {spike} ranges from 0.01 (1% of SV-spanning reads) to 1.0 (all extracted SV-spanning reads).
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged/{spike}_{hap1}_{set}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz.tbi'
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+',
        spike = '[01]\.[0-9]+'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        # simpleRepeat file created using grep filter on RepeatMasker track bed, not directly output from TRF
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.log"

rule benchmark_liftover_MATERNAL:
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_MATERNAL.dict"
    output:
        vcf = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.vcf.gz",
        tbi = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.vcf.gz.tbi",
        rejected = "benchmarks/HG002/sub10kB_{hap}_to_HG002_MATERNAL.rejected.vcf.gz",
    conda:
        "../envs/gatk.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_MATERNAL.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.MATERNAL.chain"
    shell:
        """
        picard -Xmx3G LiftoverVcf \
        -I {input.vcf} \
        -O {output.vcf} \
        -C {params.chain} \
        --REJECT {output.rejected} \
        --CREATE_INDEX true \
        -R {params.ref} \
        --VERBOSITY ERROR
        """

use rule benchmark_liftover_MATERNAL as benchmark_liftover_PATERNAL with:
    # TODO: Might be able to collapse this into one rule (with MATERNAL) if this doesn't flag the
    # requirement for all input/output to match wildcards.
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_PATERNAL.dict"
    output:
        vcf = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.vcf.gz",
        tbi = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.vcf.gz.tbi",
        rejected = "benchmarks/HG002/sub10kB_{hap}_to_HG002_PATERNAL.rejected.vcf.gz",
    conda:
        "../envs/gatk.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_PATERNAL.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.PATERNAL.chain"

rule benchmark_vcftobed:
    input:
        vcf = "benchmarks/HG002/{filename}.vcf.gz"
    output:
        bed = "benchmarks/HG002/{filename}.bed"
    params:
        vcftobed = "/global/scratch/users/stacy-l/software/ucsc_utilities/vcfToBed"
    shell:
        """
        inputVcf={input.vcf}
        outPrefix=${{inputVcf%".vcf.gz"}}
        {params.vcftobed} {input.vcf} $outPrefix
        """

rule liftover_SV_regions_bed:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap1}_SVs.bed"
    output:
        mapped = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs.liftover.bed",
        unmapped = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_SVs.liftover.unmapped.bed"
    conda:
        "../envs/liftover.yml"
    params:
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap2}.chain"
    shell:
        """
        liftOver {input} {params.chain} {output.mapped} {output.unmapped}
        """

rule truvari_anno_trf:
# liftOver -bedPlus=3 /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/anno.trf.bed \
# /global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.MATERNAL.chain \
# /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/anno.trf.MATERNAL.bed \
# /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/anno.trf.MATERNAL.unmapped.bed
# anno.trf.MATERNAL.bed | bgzip - > anno.trf.MATERNAL.bed.gz 
# then do a whole python notebook anno explode...
# liftOver -bedPlus=3 /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/annos.MATERNAL.unlifted.bed.gz \
# /global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.MATERNAL.chain \
# /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/annos.MATERNAL.lifted.bed \
# /global/scratch/users/stacy-l/spermSV/benchmarks/HG002/annos.MATERNAL.lifted.unmapped.bed
    input:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/{file}.vcf.gz"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf"),
        anno = "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf.gz"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda:
        "../envs/truvari.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        trf_ref = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/annos.{hap2}.merged.bed.gz",
        trf_exec = "/global/scratch/users/stacy-l/miniconda3/envs/truvari/bin/trf",
    threads: 1
    shell:
        """
        truvari anno trf -i {input} -o {output.temp} \
        -e {params.trf_exec} \
        -r {params.trf_ref} \
        -f {params.ref} \
        -t {threads}

        bgzip {output.temp} -c > {output.anno}
        """

use rule truvari_anno_trf as sniffles_ano_trf with:
    input:
        "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/{file}.vcf.gz"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/annotated/{file}.trf.vcf"),
        anno = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/annotated/{file}.trf.vcf.gz"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda:
        "../envs/truvari.yml"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        trf_ref = "/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/annos.{hap2}.merged.bed.gz",
        trf_exec = "/global/scratch/users/stacy-l/miniconda3/envs/truvari/bin/trf",
    threads: 1

rule truvari_anno_repmask:
    # For some reason it's broken?
    input:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf"
    output:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.repmask.vcf"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda:
        "../envs/truvari.yml"
    params:
        rm_exec = "/global/scratch/users/stacy-l/miniconda3/envs/truvari/bin/RepeatMasker"
    threads: 1
    shell:
        """
        truvari anno repmask -i {input} -o {output} \
        -e {params.rm_exec} \
        -t {threads}
        """

rule truvari_toplevel_all:
    input: 
        # NOTE: We haven't filtered the direct callset on the parameters of <10kB and whatnot so the precision may be very bad. Again.
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.jl',
        benchmark = "benchmarks/HG002/sub10kB_{hap1}_to_HG002_{hap2}.vcf.gz",
        benchmark_index = "benchmarks/HG002/sub10kB_{hap1}_to_HG002_{hap2}.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{hap1}_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --dup-to-ins \
        --passonly
        """

# TODO: Figure out how the hell the region restriction is supposed to work?

rule truvari_toplevel_subset:
    input:
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.jl',
        benchmark = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz',
        benchmark_index = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz.tbi',
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/1.0_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} \
        --pctseq 0 \
        --dup-to-ins \
        --passonly
        """

rule truvari_spike_mosaic:
    # Benchmarks SV calling results from using Sniffles mosaic mode to call variants from:
    # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2
    # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2.
    # This determines the proportion of <=10kB SVs called with progressively fewer reads available to support variants.
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.vcf.gz",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_{set}_to_{hap2}.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_{set}_to_{hap2}.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_{set}_to_{hap2}.vcf.gz.tbi",
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{spike}_{hap1}_{set}_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        # a temp fix for snakemake's behavior of preemptively creating output directories, which truvari does not like at all
        # removing the preemptive outdir frees up the path for truvari to direct outfiles
        rm -r {params.outdir}

        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification)
        truvari bench -b {input.benchmark} -c {input.query} -o {params.outdir} --pctseq 0 --passonly
        """

# use rule truvari_spike_mosaic as truvari_spike_germline with:
#     # Benchmarks SV calling results from using Sniffles germline mode to call variants from:
#     # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2
#     # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2.
#     # This determines the proportion of <=10kB SVs called with progressively fewer reads available to support variants.
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.vcf.gz",
#         jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.jl",
#         benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz",
#         benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz.tbi",
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{spike}_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 5

# use rule truvari_spike_mosaic as truvari_mode_comparison with:
#     # Benchmarks SV calling results, but compares results from two different modes.
#     # Query: A mixture of a fractional amount of <=10kB hap1 SV reads with all hap2 reads, mapped to hap2, called with mosaic mode
#     # Benchmark: All <=10kB hap1 SV reads with all hap2 reads, mapped to hap2, called with germline mode
#     # This is experimental and is intended to assess the impact of the different calling modes on SV calling.
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.vcf.gz",
#         jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{spike}_{hap1}_SVs_to_{hap2}.jl",
#         benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz",
#         benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/1.0_{hap1}_SVs_to_{hap2}.vcf.gz.tbi",
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mode_comparison/{spike}_{hap1}_SVs_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     conda: "../envs/truvari.yml"
#     threads: 5

