use rule minimap2 as minimap2_HG002_benchmark with:
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: 
        temp("output/alignment/HG002/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

rule split_HG002_reference:
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

rule split_mapped_hap:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/HG002.sorted.merged.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/HG002_{hap}.bam"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools view -H {input} | grep '^@SQ' |awk '$2 ~ /{wildcards.hap}/' | cut -f 2 | sed 's/SN://' | xargs samtools view -@ {threads} --with-header -b {input} -o {output}
        """

rule extract_hap_reads:
    """
    Produces fastq files for haplotype-mapped reads per haplotype.
    """
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/HG002_{hap}.bam"
    output:
        fastq = "output/alignment/HG002/minimap2/standard/mapped/fastq/HG002_{hap}.fastq.gz"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T '*' {input.bam} -0 {output.fastq}
        """

rule create_hap_read_name_lists:
    input:
        maternal_bam = "output/alignment/HG002/minimap2/standard/mapped/HG002_MATERNAL.bam",
        paternal_bam = "output/alignment/HG002/minimap2/standard/mapped/HG002_PATERNAL.bam"
    output:
        maternal_reads = "output/alignment/HG002/minimap2/standard/mapped/fastq/MATERNAL_rnames.txt",
        paternal_reads = "output/alignment/HG002/minimap2/standard/mapped/fastq/PATERNAL_rnames.txt"
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools view -@ {threads} {input.maternal_bam} | cut -f1 | sort | uniq > {output.maternal_reads}
        samtools view -@ {threads} {input.paternal_bam} | cut -f1 | sort | uniq > {output.paternal_reads}
        """

rule get_hap_callsets:
    # Subsets the HG002 SV benchmar vcf to only contain variants <10kb that are not BND on each hap.
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
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & GT=="1|0"' {input.vcf} -o {output.maternal} -O z9
        tabix {output.maternal}

        # Get paternal callset
        bcftools filter -i 'SVLEN <= 10000 & SVTYPE!="BND" & GT=="0|1"' {input.vcf} -o {output.paternal} -O z9
        tabix {output.paternal}
        """

rule extract_spanning_sv_reads:
    input:
        vcf = "benchmarks/HG002/sub10kB_{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        bam = "output/alignment/hg38/minimap2/standard/mapped/HG002.sorted.merged.bam",
        haplotype_reads = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_rnames.txt",
        script = "scripts/extract_spanning_sv_reads.sh"
    output:
        fastq = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.fastq.gz",
        bed = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.hg38.sv_regions.bed"
    log:
        "logs/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.log"
    conda: "../envs/truvari.yml"
    threads: 5
    shell:
        """
        bash {input.script} {input.vcf} {input.bam} {input.haplotype_reads} {output.fastq} {threads} {wildcards.hap} > {log} 2>&1
        """

use rule minimap2 as cross_map_hap_SVs with:
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/temp/{hap1}_10kB_SVs_to_{hap2}.unsorted.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{hap1}_10kB_SVs_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule minimap2 as hg38_map_hap_SVs with:
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/temp/{hap}_10kB_SVs_to_hg38.unsorted.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:{hap}_10kB_SVs_to_hg38\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule minimap2 as cross_map_hap_all with:
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/fastq/HG002_{hap1}.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/temp/all_{hap1}_to_{hap2}.unsorted.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:all_{hap1}_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule minimap2 as hg38_map_hap_all with:
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/fastq/HG002_{hap}.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/temp/all_{hap}_to_hg38.unsorted.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:all_{hap}_to_hg38\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule samtools_sort as generic_sort with:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/temp/{filename}.unsorted.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam"
    conda: "../envs/mapping.yml"
    threads: 5

use rule index_bam as generic_index with:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam.bai"
    conda: "../envs/mapping.yml"
    threads: 5

rule samtools_coverage:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/{filename}.bam"
    output:
        report = "output/alignment/HG002/minimap2/standard/mapped/reports/samtools/{filename}.coverage.txt"
    conda: "../envs/mapping.yml"
    threads: 1
    shell:
        """
        samtools coverage -o {output.report} {input.bam}
        """

rule hap_spike_in_1X_merge:
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/{hap1}_10kB_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/all_{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.bam",
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input.bam1} {input.bam2}
        """

rule hap_spike_in_fraction_merge:
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/{hap1}_10kB_SVs_to_{hap2}.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/all_{hap2}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.bam"
    wildcard_constraints:
        spike = '0\.[0-9]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools view -bs {wildcards.spike} {input.bam1} |
        samtools merge -@ {threads} {output} - {input.bam2}
        """

use rule hap_spike_in_1X_merge as hg38_hap_1X_merge with:
    # Serves as the "control", where the max # of hap2 benchmarked SVs should be callable, 
    # and the max # of the <=10kB hap1 benchmarked SVs should be callable.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/{hap1}_10kB_SVs_to_hg38.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/all_{hap2}_to_hg38.bam",
    output:
        "output/alignment/HG002/minimap2/standard/mapped/merged_1.0_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam"
    conda: "../envs/mapping.yml"
    threads: 10

use rule hap_spike_in_fraction_merge as hg38_hap_fraction_merge with:
    # The max # of hap2 benchmarked SVs should be callable, 
    # and a subset of the max # of <=10kB hap1 benchmarked SVs should be callable.
    input:
        bam1 = "output/alignment/HG002/minimap2/standard/mapped/{hap1}_10kB_SVs_to_hg38.bam",
        bam2 = "output/alignment/HG002/minimap2/standard/mapped/all_{hap2}_to_hg38.bam",
    output:
        "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam"
    conda: "../envs/mapping.yml"
    threads: 10

# On hold until I can decide whether or not it's useful to calculate coverage on HG002 T2T coordinates.
# rule gatk_dict:
#     input: "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta"
#     output: 
#         dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.dict"
#     conda: "../envs/gatk.yml"
#     threads: 1
#     shell:
#         """
#         picard CreateSequenceDictionary -R {output} -O {output.dict}
#         """

# rule gatk_SV_coverage:
#     input:
#         bam = "output/alignment/HG002/minimap2/standard/mapped/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.bam",
#         bed = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.sv_regions.bed"
#     output:
#         # amend after the outputs are more clear
#         "output/alignment/HG002/minimap2/standard/mapped/reports/gatk/merged_1.0_{hap1}_10kB_SVs_to_{hap2}"
#     conda: "../envs/gatk.yml"
#     threads: 2
#     resources:
#         mem_mb = 60000
#     params:
#         mem_mb = 58000,
#         refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta"
#     shell:
#         """
#         gatk \
#         --java-options "-Xmx{params.mem_mb}m -XX:ParallelGCThreads={threads}" \
#         DepthOfCoverage \
#         -R {params.refgenome} \
#         -O {output} \
#         -I {input.bam} \
#         -pt readgroup \
#         -L {input.bed}
#         """

use rule sniffles_standard as hap_specific_germline_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/all_{hap1}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/all_{hap1}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.log"

use rule sniffles_standard as hap_hg38_germline_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/all_{hap}_to_hg38.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/all_{hap}_to_hg38.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap}_to_hg38.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap}_to_hg38.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap}_to_hg38.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap}_to_hg38.log"

use rule sniffles_standard as cross_hap_spike_in_fraction_germline_calls with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.log"

use rule sniffles_standard as hg38_spike_in_fraction_germline_calls with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.log"

use rule sniffles_mosaic as cross_hap_spike_in_fraction_mosaic_calls with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.log"

use rule sniffles_mosaic as hg38_spike_in_fraction_mosaic_calls with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.log"

rule vcf2df:
    input:
        'output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/{filename}.vcf.gz'
    output:
        'output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/{filename}.jl'
    conda:
        '../envs/truvari.yml'
    threads: 1
    shell:
        """
        truvari vcf2df --info --format {input} {output}
        """

rule liftover_sv_regions:
    # FLAG FOR DEPRECATION?
    # Uses UCSC liftOver to convert coordinates of the benchmark <=10kB SVs from hg38 to HG002 T2T.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.hg38.sv_regions.bed"
    output:
        bed = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.HG002_T2T.sv_regions.bed",
        unmapped = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.HG002_T2T.sv_regions.unmapped.bed"
    conda:
        '../envs/liftover.yml'
    params:
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap}.chain"
    threads: 1
    shell:
        """
        liftOver {input} {params.chain} {output.bed} {output.unmapped}
        """

rule swap_hap_sv_regions:
    # FLAG FOR DEPRECATION
    # Assuming that the coordinates of the chrs are not significantly different between MATERNAL and PATERNAL,
    # This allows for assessment of hap1 SVs called on hap2 coordinates.
    # chrX and chrY contigs are dropped because they cannot be swapped.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.HG002_T2T.sv_regions.bed"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap}_10kB_SVs_spanning.HG002_T2T.sv_regions.{hap2}_coords.bed"
    threads: 1
    shell:
        """
        if [ "{wildcards.hap}" = "MATERNAL" ]; then
            sed 's/MATERNAL/PATERNAL/g' {input} | grep -v -e 'chrX' > {output}
        else
            sed 's/PATERNAL/MATERNAL/g' {input} | grep -v -e 'chrY' > {output}
        fi
        """

rule truvari_cross_hap_subset_vs_ctrl:
    # Benchmarks the recall/precision of SVs called with germline mode from:
    # query: 10kB hap1 SV reads spiked into hap2
    # benchmark: all hap1 reads spiked into hap2.
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.jl',
        benchmark = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.vcf.gz',
        benchmark_index = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/all_{hap1}_to_{hap2}.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.HG002_T2T.sv_regions.{hap2}_coords.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/T2T/subset_comparison/merged_1.0_{hap1}_10kB_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
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
        --includebed {input.regions} \
        --passonly
        """

rule truvari_hg38_germline:
    # TODO: Figure out if any of this really makes sense oh my god
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.jl',
        benchmark = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz',
        benchmark_index = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.hg38.sv_regions.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
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
        --includebed {input.regions} \
        --passonly
        """

rule truvari_hg38_mosaic:
    # TODO: Figure out if any of this really makes sense oh my god
    input: 
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.jl',
        benchmark = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz',
        benchmark_index = '/global/scratch/users/stacy-l/spermSV/benchmarks/HG002/sub10kB_{hap1}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.hg38.sv_regions.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
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
        --includebed {input.regions} \
        --passonly
        """

use rule truvari_hg38_mosaic as truvari_hg38_mosaic_to_germline_TP with:
    # Uses the true positives from the 1X hap1 spike germline call to evaluate the recall/precision of the mosaic germline call.
    input:
        query = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.vcf.gz',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}_to_hg38.jl',
        benchmark = 'output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/merged_1.0_{hap1}_10kB_SVs_with_{hap2}/tp-comp.vcf.gz',
        benchmark_index = 'output/alignment/HG002/minimap2/standard/variants/truvari/hg38/germline/merged_1.0_{hap1}_10kB_SVs_with_{hap2}/tp-comp.vcf.gz.tbi',
        regions = "output/alignment/HG002/minimap2/standard/mapped/fastq/{hap1}_10kB_SVs_spanning.hg38.sv_regions.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/mosaic_to_germline_TP/merged_{spike}_{hap1}_10kB_SVs_with_{hap2}/{outfiles}", allow_missing = True,
            outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),

rule truvari_spike_mosaic_matched:
    # Benchmarks the recall/precision of SVs called with mosaic mode from:
    # query: A fractional amount of 10kB hap1 SV reads spiked into hap2
    # benchmark: All 10kB hap1 SV reads spiked into hap2.
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz.tbi",
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/T2T/mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
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

use rule truvari_spike_mosaic_matched as truvari_spike_germline_matched with:
    # Benchmarks the recall/precision of SVs called with germline mode from:
    # query: A fractional amount of 10kB hap1 SV reads spiked into hap2
    # benchmark: All 10kB hap1 SV reads spiked into hap2.
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz.tbi",
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/T2T/germline/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5

use rule truvari_spike_mosaic_matched as truvari_mode_comparison with:
    # Benchmarks the recall/precision of 10kB hap1 SV reads spiked into hap2.
    # query: All 10kB hap1 SV reads spiked into hap2, mosaic mode.
    # benchmark: All 10kB hap2 SV reads spiked into hap2, germline mode.
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/merged_1.0_{hap1}_10kB_SVs_to_{hap2}.vcf.gz.tbi",
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/T2T/mode_comparison/merged_{spike}_{hap1}_10kB_SVs_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 5
