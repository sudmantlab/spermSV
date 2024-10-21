# Section 1: Generic rules

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

# Section 2: Read phasing

rule split_T2T_HG002:
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

use rule minimap2 as T2T_self_mapping with:
    # Maps the HG002 HiFi reads to the HG002 T2T assembly.
    # The output BAM goes through sorting and merging operations (not shown in this file), ultimately
    # producing the file (output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam)
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

rule split_T2T_self_mapped_by_hap:
    # This rule takes in the mapped HiFi reads and splits them into two separate BAM files,
    # one for each haplotype.
    input:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/all_{hap}.bam"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools view -H {input} | grep '^@SQ' |awk '$2 ~ /{wildcards.hap}/' | cut -f 2 | sed 's/SN://' | xargs samtools view -@ {threads} --with-header -b {input} -o {output}
        """

rule create_hap_read_name_lists:
    # This rule creates a list of read names mapped to each haplotype.
    input:
        maternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/all_MATERNAL.bam",
        paternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/all_PATERNAL.bam"
    output:
        maternal_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/MATERNAL_rnames.txt",
        paternal_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/PATERNAL_rnames.txt",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        # Extract read names
        samtools view -@ {threads} {input.maternal_bam} | cut -f1 | sort | uniq > {output.maternal_reads}.tmp
        samtools view -@ {threads} {input.paternal_bam} | cut -f1 | sort | uniq > {output.paternal_reads}.tmp
        
        # Find shared read names
        comm -12 {output.maternal_reads}.tmp {output.paternal_reads}.tmp > {output.shared_reads}
        
        # Remove shared reads from maternal and paternal lists
        comm -23 {output.maternal_reads}.tmp {output.shared_reads} > {output.maternal_reads}
        comm -23 {output.paternal_reads}.tmp {output.shared_reads} > {output.paternal_reads}
        
        # Clean up temporary files
        rm {output.maternal_reads}.tmp {output.paternal_reads}.tmp
        """

rule filter_dipcall_benchmark_SVs:
    # This rule takes in the HG002 structural variant benchmark VCF file, which was created by using dipcall
    # to call SVs on the diploid T2T HG002 assembly against the GRCh38 (hg38) reference genome.
    # It filters the VCF, retaining only variants that are:
    # 1) not INV or BND
    # 2) exclusively phased to one of the haplotypes
    # 3) greater/equal to 50 bp (Sniffles default minimum is 50 bp)
    # It creates two VCFs, one for each haplotype, containing only the variants that are exclusive to that haplotype.
    # It accomplishes this by filtering for variants that are specifically phased to one of the haplotypes.
    input:
        vcf = 'benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz'
    output:
        unphased = "benchmarks/HG002/unphased_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        unphased_index = "benchmarks/HG002/unphased_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        maternal = "benchmarks/HG002/MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        maternal_index = "benchmarks/HG002/MATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        paternal = "benchmarks/HG002/PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        paternal_index = "benchmarks/HG002/PATERNAL_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        homozygous = "benchmarks/HG002/homozygous_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        homozygous_index = "benchmarks/HG002/homozygous_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p {params.outdir}
        # Get "unphased" callset (no filter on phase)
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & INFO/SVLEN >= 50' {input.vcf} -o {output.unphased} -O z9
        tabix {output.unphased}

        # Get maternal callset
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & INFO/SVLEN >= 50 & GT=="0|1"' {input.vcf} -o {output.maternal} -O z9
        tabix {output.maternal}

        # Get paternal callset
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & INFO/SVLEN >= 50 & GT=="1|0"' {input.vcf} -o {output.paternal} -O z9
        tabix {output.paternal}

        # Get homozygous callset
        bcftools filter -i 'SVTYPE!="BND" & SVTYPE!="INV" & INFO/SVLEN >= 50 & GT=="1|1"' {input.vcf} -o {output.homozygous} -O z9
        tabix {output.homozygous}
        """

use rule filter_dipcall_benchmark_SVs as filter_CMRG_benchmark_SVs with:
    # This rule uses the Wagner 2022 SV benchmark VCF file, splitting it by haplotype and filtering as above.
    input:
        vcf = "benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0.vcf.gz"
    output:
        unphased = "benchmarks/HG002/unphased_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        unphased_index = "benchmarks/HG002/unphased_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
        maternal = "benchmarks/HG002/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        maternal_index = "benchmarks/HG002/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
        paternal = "benchmarks/HG002/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        paternal_index = "benchmarks/HG002/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
        homozygous = "benchmarks/HG002/homozygous_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        homozygous_index = "benchmarks/HG002/homozygous_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi"
    conda: "../envs/truvari.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])

rule disambiguate_T2T_self_mapped:
    # This rule takes in the BAM of {hap} phased reads mapped to self and removes reads that are mapped to both haplotypes.'
    # We retain the bamfile that contains ambiguous reads for now, because we want to understand the impact of removing ambiguous reads.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/all_{hap}.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.bam",
    conda:
        "../envs/truvari.yml"
    threads: 5
    shell:
        """
        samtools view -h -b -@{threads} -N ^{input.shared_reads} {input.bam} -o {output}
        """

rule unambiguous_T2T_self_mapped_to_fastq:
    # This rule converts the unambiguous hap-aligned reads in the BAM and converts them to a FASTQ format.
    # This step is used to phase the reads.'
    # TODO: We may not need this, as we consider the hg38 mapping and deambiguation separately.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.bam"
    output:
        fastq = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.fastq.gz"
    wildcard_constraints:
        hap = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T '*' {input.bam} -0 {output.fastq}
        """

# Section 3: Mapping and variant calling 

rule unphased_T2T_self_mapped_to_fastq:
    input:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.fastq.gz"
    conda: "../envs/mapping.yml"
    threads: 40
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T '*' {input} -0 {output}
        """

use rule minimap2 as map_unphased_to_hg38 with:
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/hg38/unsorted/unphased.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:unphased\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 40

use rule minimap2 as remap_ambiguous_to_hg38 with:
    # This rule maps all phased (but not unambiguously phased) HG002 reads to the hg38 reference.
    # We don't input the unambiguous FASTQs here because we want to understand the impact of removing ambiguous reads.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/all_{hap}.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/hg38/unsorted/all_{hap}.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = "@RG\\tID:HG002\\tDS:all_{hap}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10

use rule disambiguate_T2T_self_mapped as disambiguate_hg38_remapped with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam",
    conda:
        "../envs/truvari.yml"

rule retain_homozygous_hg38_remapped:
    # This rule takes in the BAM of unphased reads mapped to hg38 and selects reads that are mapped to both haplotypes.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unphased.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/hg38/homozygous.bam",
    conda:
        "../envs/truvari.yml"
    threads: 5
    shell:
        """
        samtools view -h -b -@{threads} -N {input.shared_reads} {input.bam} -o {output}
        """

use rule retain_homozygous_hg38_remapped as retain_homozygous_T2T_self_mapped with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/homozygous.bam"

use rule sniffles_standard as hg38_unphased_standard_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unphased.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/hg38/unphased.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unphased.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unphased.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unphased.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unphased.log"

rule T2T_unphased_standard_call:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unphased.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unphased.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unphased.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unphased.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unphased.log"
    shell:
        """
        sniffles --input {input.bam} \
        --vcf {output.vcf} \
        --snf {output.snf} \
        --reference {params.refgenome} \
        --threads {threads} \
        --mapq {params.mapq} \
        --output-rnames &> {log}
        """

use rule sniffles_standard as hg38_homozygous_standard_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/homozygous.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/hg38/homozygous.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/homozygous.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/homozygous.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/homozygous.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/homozygous.log"

use rule T2T_unphased_standard_call as T2T_homozygous_standard_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap}/homozygous.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/self/{hap}/homozygous.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap}/homozygous.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap}/homozygous.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap}/homozygous.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap}/homozygous.log"

use rule sniffles_standard as hg38_remapped_ambiguous_hap_standard_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of phased (but not unambiguously phased) reads mapped to hg38.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/all_{hap}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/all_{hap}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/all_{hap}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/all_{hap}.log"

use rule T2T_unphased_standard_call as T2T_ambiguous_hap_standard_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap1}/{hap2}_to_{hap1}.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/self/{hap1}/{hap2}_to_{hap1}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/{hap2}_to_{hap1}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/{hap2}_to_{hap1}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/{hap2}_to_{hap1}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/{hap2}_to_{hap1}.log"

use rule sniffles_standard as hg38_remapped_unambiguous_hap_standard_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of unambiguous {hap} phased reads mapped to hg38.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unambiguous_{hap}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unambiguous_{hap}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unambiguous_{hap}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/unambiguous_{hap}.log"

use rule T2T_unphased_standard_call as T2T_unambiguous_hap_standard_call with:
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap1}/unambiguous_{hap2}_to_{hap1}.bam",
        bai = "output/alignment/HG002/minimap2/standard/mapped/self/{hap1}/unambiguous_{hap2}_to_{hap1}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/unambiguous_{hap2}_to_{hap1}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/unambiguous_{hap2}_to_{hap1}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/unambiguous_{hap2}_to_{hap1}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        5
    resources:
        mem_mb=60000
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap1}/unambiguous_{hap2}_to_{hap1}.log"

# Section 4: Benchmarking mapped variants

rule hg38_full_CMRG_benchmark:
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.jl",
        benchmark = "benchmarks/HG002/unphased_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        benchmark_index = "benchmarks/HG002/unphased_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
        includebed = "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/hg38_full_CMRG_benchmark/{file}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        refgenome = config['reference']['fasta'],
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification?)
        truvari bench \
        -f {params.refgenome} \
        -b {input.benchmark} \
        -c {input.query} \
        -o {params.outdir}/bench \
        -r 1000 \
        --includebed {input.includebed} \
        --dup-to-ins \
        --passonly

        mv {params.outdir}/bench/* {params.outdir}/.
        rm -r {params.outdir}/bench
        """

use rule hg38_full_CMRG_benchmark as hg38_hap_CMRG_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.jl",
        benchmark = "benchmarks/HG002/{hap}_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        benchmark_index = "benchmarks/HG002/{hap}_HG002_GRCh38_CMRG_SV_v1.00.vcf.gz.tbi",
        includebed = "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/CMRG_{hap}_benchmark/{file}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as hg38_full_dipcall_benchmark with:
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.jl",
        benchmark = "benchmarks/HG002/unphased_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/unphased_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        includebed = "benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/full_dipcall_benchmark/{file}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as hg38_hap_dipcall_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{file}.jl",
        benchmark = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi",
        includebed = "benchmarks/HG002/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/dipcall_{hap}_benchmark/{file}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as T2T_MATERNAL_to_PATERNAL_CMRG_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.jl",
        benchmark = "benchmarks/HG002/liftover/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.vcf.gz",
        benchmark_index = "benchmarks/HG002/liftover/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.vcf.gz.tbi",
        includebed = "benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.bed"
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/MATERNAL/CMRG_MATERNAL_benchmark/unambiguous_MATERNAL_to_PATERNAL/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as T2T_PATERNAL_to_MATERNAL_CMRG_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.jl",
        benchmark = "benchmarks/HG002/liftover/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.vcf.gz",
        benchmark_index = "benchmarks/HG002/liftover/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.vcf.gz.tbi",
        includebed = "benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.bed"
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/PATERNAL/CMRG_PATERNAL_benchmark/unambiguous_PATERNAL_to_MATERNAL/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as T2T_MATERNAL_to_PATERNAL_dipcall_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_MATERNAL.jl",
        benchmark = "benchmarks/HG002/liftover/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.vcf.gz",
        benchmark_index = "benchmarks/HG002/liftover/MATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.vcf.gz.tbi",
        includebed = "benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_PATERNAL.bed"
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/MATERNAL/dipcall_MATERNAL_benchmark/unambiguous_MATERNAL_to_PATERNAL/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

use rule hg38_full_CMRG_benchmark as T2T_PATERNAL_to_MATERNAL_dipcall_benchmark with:
    # NOTE: The benchmark bed file is not split by hap * (does it need to be?)
    input:
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/unambiguous_PATERNAL.jl",
        benchmark = "benchmarks/HG002/liftover/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.vcf.gz",
        benchmark_index = "benchmarks/HG002/liftover/PATERNAL_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.vcf.gz.tbi",
        includebed = "benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_MATERNAL.bed"
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1.fasta.gz",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/PATERNAL/dipcall_PATERNAL_benchmark/unambiguous_PATERNAL_to_MATERNAL/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])

# Section 5: SV read extraction and sampling
rule extract_spanning_reads:
    input:
        script = "scripts/extract_spanning_per_sv.sh",
        vcf = "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/tp-comp.vcf.gz",
        bam = "output/alignment/HG002/minimap2/standard/mapped/{benchmark_ref}/unambiguous_{hap}.bam"
    output:
        outdir = directory("output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}"),
        var_rnames = "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}.txt"
    conda: "../envs/bcftools.yml"
    wildcard_constraints:
        hap = '[A-Za-z]+',
        chr = '[A-Za-z0-9]+'
    threads: 1
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    log:
        "logs/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}.log"
    shell:
        """
        bash {input.script} \
            {input.vcf} \
            {input.bam} \
            {wildcards.chr} \
            {output.var_rnames} \
            {output.outdir} \
            {threads} \
            {log}
        """

def get_var_fastqs(wildcards):
    pattern = f"output/alignment/HG002/minimap2/standard/variants/truvari/{wildcards.benchmark_ref}/{wildcards.benchmark}_{wildcards.hap}_benchmark/unambiguous_{wildcards.hap}/extracted_vars/{wildcards.chr}/*.fastq.gz"
    files = glob.glob(pattern)
    return files

rule aggregate_var_reads:
    input:
        get_var_fastqs
    output:
        "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}/all_var_reads.fastq.gz"
    shell:
        """
        cat {input} > {output}
        """

rule sample_n_var_reads:
    input:
        fastqs = get_var_fastqs
    output:
        gz = "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}/var_{n}_reads.fastq.gz"
    params:
        seed = 42,
        temp_fq = lambda wildcards, output: output.gz.strip('.gz')
    conda: "../envs/mapping.yml"
    log:
        "logs/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}/var_{n}_reads.log"
    threads: 1
    shell:
        """
        if [ -z "{input.fastqs}" ]; then
            echo "No variant FASTQ files found for {wildcards.hap} on {wildcards.chr}. Creating an empty output file." >> {log}
            touch {output.gz}
        else
            for fastq in {input.fastqs}; do
                # Extract var_id from the filename
                var_id=$(basename "$fastq" .fastq.gz)
                echo "Sampling {wildcards.n} from $var_id on {wildcards.chr}." >> {log}
                
                # Execute seqtk command
                seqtk sample -s{params.seed} $fastq {wildcards.n} >> {params.temp_fq}
                echo "Sampled {wildcards.n} from $var_id." >> {log}
            done
            gzip -c {params.temp_fq} > {output.gz}
            echo "Gzipped {params.temp_fq}." >> {log}
            rm {params.temp_fq}
        fi
        """

use rule aggregate_var_reads as merge_aggregated_reads with:
    input:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}/all_var_reads.fastq.gz", allow_missing = True, chr = chrs)
    output:
        "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/all_var_reads.fastq.gz"

use rule aggregate_var_reads as merge_n_sampled_reads with:
    input:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/extracted_vars/{chr}/var_{n}_reads.fastq.gz", allow_missing = True, chr = chrs)
    output:
        "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap}_benchmark/unambiguous_{hap}/simulation_spike/{n}_{hap}_spike.fastq.gz"

# Section 6: Simulation of mapping and spiking in SV reads

use rule minimap2 as baseline_self_assembly_mapping with:
    # This rule maps HiFi reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    # For example, this rule can map MATERNAL phased reads to either the MATERNAL assembly or the PATERNAL assembly.
    # We expect optimal mapping for matched haplotype mapping and worse mapping for cross-haplotype mapping.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap1}.fastq.gz"
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

use rule minimap2 as spiked_self_assembly_mapping with:
    # This rule maps HiFi reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    # For example, this rule can map MATERNAL phased reads to either the MATERNAL assembly or the PATERNAL assembly.
    # We expect optimal mapping for matched haplotype mapping and worse mapping for cross-haplotype mapping.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/variants/truvari/{benchmark_ref}/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/simulation_spike/{n}_{hap1}_spike.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{hap2}/unsorted/temp_{n}_{hap1}_spike_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{n}_{hap1}_spike_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 2

use rule minimap2 as aggregated_self_assembly_mapping with:
    # This rule maps HiFi reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    # For example, this rule can map MATERNAL phased reads to either the MATERNAL assembly or the PATERNAL assembly.
    # We expect optimal mapping for matched haplotype mapping and worse mapping for cross-haplotype mapping.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/extracted_vars/all_var_reads.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{hap2}/unsorted/temp_all_{hap1}_spike_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:all_{hap1}_spike_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 2

rule merge_spiked_self_assembly_mapping:
    input:
        baseline = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap2}_to_{hap2}.bam",
        spike = "output/alignment/HG002/minimap2/standard/mapped/self/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{hap2}/temp_{n}_{hap1}_spike_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/simulation/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{n}_{hap1}_spike_to_{hap2}.bam",
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 5
    shell:
        """
        samtools merge -@ {threads} {output} {input.baseline} {input.spike}
        """

use rule sniffles_standard as baseline_self_assembly_standard_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{hap1}_to_{hap2}/{hap1}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{hap1}_to_{hap2}/{hap1}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{hap1}_to_{hap2}/{hap1}_to_{hap2}.vcf.gz.tbi'
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
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.log"

use rule sniffles_mosaic as spiked_self_assembly_mosaic_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/simulation/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{n}_{hap1}_spike_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/simulation/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{n}_{hap1}_spike_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_mosaic/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_mosaic/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_mosaic/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz.tbi'
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
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_mosaic/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}.log"

use rule sniffles_standard as spiked_self_assembly_standard_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/simulation/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{n}_{hap1}_spike_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/simulation/{benchmark_ref}/{benchmark}_{hap1}_benchmark/{n}_{hap1}_spike_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz.tbi',
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
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        repeats = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.simpleRepeat.bed",
        mapq = config['sniffles']['mapq'],
    log:
        "logs/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_standard/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}.log"

rule check_multi_hap_SVs:
    input:
        script = "scripts/python/check_multi_hap_SVs.py",
        spiked = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz",
        hap1_rnames = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap1}_rnames.txt",
        hap2_rnames = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap2}_rnames.txt",
    output:
        mixed = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/mixed_hap.vcf.gz",
        hap1_only = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/hap1_only.vcf.gz",
        hap2_only = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/hap2_only.vcf.gz",
        report = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/multi_hap.report", # json formatted
        rnames = "output/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/multi_hap.rnames" # json formatted
    log:
        "logs/alignment/HG002/minimap2/standard/variants/simulation/{benchmark_ref}/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/multi_hap.log"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/sniffles.yml"
    threads: 1
    shell:
        """
        python {input.script} \
        {input.spiked} \
        {input.hap1_rnames} \
        {input.hap2_rnames} \
        {output.mixed} \
        {output.hap1_only} \
        {output.hap2_only} \
        {output.report} \
        {output.rnames} \
        {log}
        """

# Section 7: Evaluation of mapping and spiking simulation results

# Okay we are at the same problem where I can't compare cross map with hg38 sim so I need to work on the clever way to do this.
# The "baseline false positives" are found by matched hap mapping and calling.
# The "maximum expectation" is the difference between the matched hap calling and the cross map calling (all cross mapped).
# The "standard expectation is the difference between the cross map calling (all cross mapped) and SV-only cross map calling (all cross SV reads only).
# The "simulation recall" is the difference between the standard recall and where a progressively smaller number of reads are spiked in.


# rule hg38_simulation_hap_CMRG_benchmark:
#     # TODO: This currently benchmarks against the true positive set that was mapped against hg38 and benched against the hg38 CMRG benchmark, THEN lifted over.
#     # This is not really the right approach.
#     input:
#         query = "output/alignment/HG002/minimap2/standard/variants/simulation/hg38/sniffles_{setting}/CMRG_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz",
#         query_index = "output/alignment/HG002/minimap2/standard/variants/simulation/hg38/sniffles_{setting}/CMRG_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.vcf.gz.tbi",
#         jl = "output/alignment/HG002/minimap2/standard/variants/simulation/hg38/sniffles_{setting}/CMRG_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/callset.jl",
#         benchmark = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/CMRG_{hap1}_benchmark/unambiguous_{hap1}/tp-base.to_HG002_{hap2}.vcf.gz",
#         benchmark_index = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/CMRG_{hap1}_benchmark/unambiguous_{hap1}/tp-base.to_HG002_{hap2}.vcf.gz.tbi",
#         # includebed = "benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_{hap2}.bed"
#         # includebed = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/CMRG_{hap1}_benchmark/unambiguous_{hap1}/tp-base.to_HG002_{hap2}.refwidened.bed"
#     output:
#         expand("output/alignment/HG002/minimap2/standard/variants/simulation/hg38/sniffles_{setting}/CMRG_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/{outfiles}", allow_missing = True,
#                outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
#     params:
#         refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
#         outdir = lambda wildcards, output: os.path.dirname(output[0])
#     conda: "../envs/truvari.yml"
#     threads: 1
#     shell:
#         """
#         # --pctseq 0 required to analyze <DEL> (unresolved deletion, needs clarification?)
#         truvari bench \
#         -f {params.refgenome} \
#         -b {input.benchmark} \
#         -c {input.query} \
#         -o {params.outdir}/bench \
#         -r 1000 \
#         --dup-to-ins \
#         --passonly

#         mv {params.outdir}/bench/* {params.outdir}/.
#         rm -r {params.outdir}/bench
#         """

# Temporary register
rule benchmark_liftovervcf:
    input:
        vcf = "benchmarks/HG002/{benchmark}.vcf.gz",
        dict = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.dict"
    output:
        vcf = "benchmarks/HG002/liftover/{benchmark}.to_HG002_{hap}.vcf.gz",
        tbi = "benchmarks/HG002/liftover/{benchmark}.to_HG002_{hap}.vcf.gz.tbi",
        rejected = "benchmarks/HG002/liftover/{benchmark}.to_HG002_{hap}.rejected.vcf.gz",
    wildcard_constraints:
        hap = '[A-Za-z]+'
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap}.chain"
    threads: 1
    conda: "../envs/gatk.yml"
    log: "logs/benchmarks/HG002/liftover/{benchmark}.to_HG002_{hap}.log"
    shell:
        """
        picard -Xmx3G LiftoverVcf \
        -I {input.vcf} \
        -O {output.vcf} \
        -C {params.chain} \
        --REJECT {output.rejected} \
        --CREATE_INDEX true \
        -R {params.ref} \
        --VERBOSITY DEBUG 2> {log}
        """

use rule benchmark_liftovervcf as tp_liftovervcf with:
    input:
        vcf = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/tp-{set}.vcf.gz"
    output:
        vcf = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/tp-{set}.to_HG002_{hap2}.vcf.gz",
        tbi = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/tp-{set}.to_HG002_{hap2}.vcf.gz.tbi",
        rejected = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/tp-{set}.to_HG002_{hap2}.rejected.vcf.gz"
    params:
        ref = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta",
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap2}.chain"
    log: "logs/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap1}_benchmark/unambiguous_{hap1}/tp-{set}.to_HG002_{hap2}.log"

rule benchmark_liftoverbed:
    input:
        "benchmarks/HG002/{file}.bed",
    output:
        mapped = "benchmarks/HG002/liftover/{file}.to_HG002_{hap}.bed",
        unmapped = "benchmarks/HG002/liftover/{file}.to_HG002_{hap}.unmapped.bed"
    conda:
        "../envs/liftover.yml"
    params:
        chain = "/global/scratch/users/stacy-l/references/HG002/GRCh38_to_hg002v1.0.{hap}.chain"
    shell:
        """
        liftOver {input} {params.chain} {output.mapped} {output.unmapped}
        """
rule get_refwidened_CMRG_bed:
    input:
        "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
    output:
        "benchmarks/HG002/HG002_GRCh38_CMRG_SV_v1.00.refwidened.bed"
    shell:
        """
        echo "Creating refwidened BED file: ${output}"
        bcftools query -f '%INFO/REFWIDENED\n' {input} |
        awk -F':|-' '{
            chr = $1;
            start = $2;
            end = $3;
            print chr"\t"start"\t"end
        }' > {output}
        """

rule get_refwidened_vcf_bed:
    input:
        "{file}.vcf.gz"
    output:
        "{file}.refwidened.bed"
    params:
        padding = 1000
    conda:
        "../envs/truvari.yml"
    shell:
        """
        echo "Creating refwidened BED file: {output}"
        bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\n' {input} |
        awk -v padding={params.padding} -F'\t' '{{
            chr = $1;
            pos = $2;
            svlen = $3;
            start = pos - padding;
            end = pos + svlen + padding;
            print chr"\t"start"\t"end
        }}' > {output}
        """