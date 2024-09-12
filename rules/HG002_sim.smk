### generic rules for sorting, indexing, coverage, etc. ###

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

### read phasing
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

rule create_hap_read_name_lists:
    # This rule creates a list of read names mapped to each haplotype.
    input:
        maternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/MATERNAL.bam",
        paternal_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/PATERNAL.bam"
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

rule get_unambiguous_hap_bam:
    # TODO: Fix/unify this with the above– this should ideally be performed upstream directly after mapping to self and phasing + removing ambiguous
    # and before converting to the FASTQs that get mapped to hg38.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap}.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        filtered_bam = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.bam",
    conda:
        "../envs/truvari.yml"
    threads: 5
    shell:
        """
        samtools view -h -b -@{threads} -N ^{input.shared_reads} {input.bam} -o {output.filtered_bam}
        """

use rule get_unambiguous_hap_bam as get_unambiguous_hap_bam_hg38 with:
    # TODO: Fix/unify this with the above– this should ideally be performed upstream directly after mapping to self and phasing + removing ambiguous
    # and before converting to the FASTQs that get mapped to hg38.
    # Thus, "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.fastq.gz" should be fed to hg38 mapping.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/{hap}.bam",
        shared_reads = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/shared.rnames.txt"
    output:
        filtered_bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam",
    conda:
        "../envs/truvari.yml"

rule get_unambiguous_hap_fastq:
    # This rule converts the unambiguous hap-aligned reads in the BAM and converts them to a FASTQ format.
    # This step is used to phase the reads.
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

use rule sniffles_standard as hap_hg38_germline_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap} phased reads mapped to hg38.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam.bai"
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
        query_index = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.vcf.gz.tbi',
        jl = 'output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.jl',
        benchmark = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz",
        benchmark_index = "benchmarks/HG002/{hap}_GRCh38_HG2-T2TQ100-V1.0.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/all/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 1
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


### prepare reads for variant simulation ###
rule extract_spanning_reads:
    input:
        script = "scripts/extract_spanning_per_sv.sh",
        vcf = "output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/all/tp-comp.vcf.gz",
        bam = "output/alignment/HG002/minimap2/standard/mapped/hg38/unambiguous_{hap}.bam"
    output:
        var_rnames = "output/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}.txt",
        outdir = directory("output/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}")
    conda: "../envs/bcftools.yml"
    threads: 1
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
    log:
        "logs/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}.log"
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

### actual T2T hap-specific simulation ###

def get_var_fastqs(wildcards):
    pattern = f"output/alignment/HG002/minimap2/standard/variants/extracted_vars/{wildcards.hap}/{wildcards.chr}/*.fastq.gz"
    files = glob.glob(pattern)
    return files

rule sample_var_reads:
    input:
        fastqs = get_var_fastqs
    output:
        gz = "output/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}/var_{n}_reads.fastq.gz"
    params:
        seed = 42,
        temp_fq = lambda wildcards, output: output.gz.strip('.gz')
    conda: "../envs/mapping.yml"
    log:
        "logs/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}/var_{n}_reads.log"
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


rule merge_var_reads:
    input:
        expand("output/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap1}/{chr}/var_{n}_reads.fastq.gz", allow_missing = True, chr = chrs)
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{n}_{hap1}_spike.fastq.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule spike_var_reads:
    input:
        spike = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{n}_{hap1}_spike.fastq.gz",
        base = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap2}.fastq.gz"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{n}_{hap1}_spiked_{hap2}.fastq.gz"
    threads: 1
    shell:
        """
        cat {input.spike} {input.base} > {output}
        """

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

use rule minimap2 as spike_self_assembly_mapping with:
    # This rule maps HiFi reads to a haplotype assembly, allowing for self- and cross-mapping between haplotypes.
    # For example, this rule can map MATERNAL phased reads to either the MATERNAL assembly or the PATERNAL assembly.
    # We expect optimal mapping for matched haplotype mapping and worse mapping for cross-haplotype mapping.
    input:
        hifi = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{n}_{hap1}_spike.fastq.gz"
    output:
        temp("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/unsorted/{n}_{hap1}_to_{hap2}.bam")
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    params:
        refgenome = "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap2}.fasta.gz",
        readgroup = "@RG\\tID:HG002\\tDS:{n}_{hap1}_spike_to_{hap2}\\tPL:PACBIO",
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 1

rule merge_spiked_self_assembly_mapping:
    input:
        baseline = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.bam",
        spike = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{n}_{hap1}_to_{hap2}.bam"
    output:
        "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.bam",
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input.baseline} {input.spike}
        """

use rule sniffles_standard as baseline_self_assembly_germline_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
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

use rule sniffles_mosaic as spiked_self_assembly_mosaic_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.vcf.gz.tbi'
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
        minsupport = config['sniffles']['minsupport'],
        mapq = config['sniffles']['mapq'],
        mosaic_af_min = config['sniffles']['mosaic-af-min'],
        mosaic_af_max = config['sniffles']['mosaic-af-max'],
        mosaic_qc_strand = config['sniffles']['mosaic-qc-strand']
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/{n}_{hap1}_to_{hap2}.log"

use rule sniffles_standard as spiked_self_assembly_germline_call with:
    # This rule uses the standard Sniffles germline calling mode.
    # It call SVs from a BAM of all {hap1} phased reads mapped to {hap2}.
    input:
        bam = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.bam",
        index = "output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.bam.bai"
    output:
        vcf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.vcf.gz',
        snf='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.snf',
        tbi='output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.vcf.gz.tbi'
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
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/merged_{n}_{hap1}_to_{hap2}.log"

### Output vcf annotation

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
# Substitute all "motif" with "repeat" for truvari compatibility.
# zcat {file} | sed -e 's/\"motif\"/\"repeat\"/g' | bgzip - > {file}
    input:
        "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/{file}.vcf.gz"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf"),
        anno = "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf.gz",
        tbi = "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{subdirs}/annotated/{file}.trf.vcf.gz.tbi"
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

        bcftools sort -O z9 {output.temp} -o {output.anno} --write-index=tbi
        """

use rule truvari_anno_trf as sniffles_anno_trf with:
    input:
        "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/{file}.vcf.gz"
    output:
        temp = temp("output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/annotated/{file}.trf.vcf"),
        anno = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/annotated/{file}.trf.vcf.gz",
        tbi = "output/alignment/HG002/minimap2/standard/variants/{sniffles_setting}/self/{hap2}/annotated/{file}.trf.vcf.gz.tbi"
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

### Annotated vcf benchmarking against hg38 TP calls

use rule truvari_hg38_germline_all as truvari_baseline_self_assembly_benchmark with:
    # Broken atm due to incompatibility of hg38 variant calls against cross hap mapping
    # Chain file solution won't yield all variants
    input: 
        query = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.vcf.gz",
        query_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.vcf.gz.tbi",
        jl = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.jl",
        benchmark = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/annotated/{hap1}_to_{hap2}.trf.vcf.gz",
        benchmark_index = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/annotated/{hap1}_to_{hap2}.trf.vcf.gz.tbi"
    output:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{setting}/{n}_{hap1}_to_{hap2}/{outfiles}", allow_missing = True,
               outfiles = ["tp-base.vcf.gz", "tp-comp.vcf.gz", "fp.vcf.gz", "fn.vcf.gz", "summary.json", "params.json", "candidate.refine.bed", "log.txt"])
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 1
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])

### Annotated vcf intersection

rule truvari_consistency:
    input:
        baseline = "output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/annotated/{hap1}_to_{hap2}.trf.vcf.gz",
        spiked = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.vcf.gz"
    output:
        tsv = "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{setting}/{n}_{hap1}_to_{hap2}/spiked_consistency.tsv",
        json = "output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/{setting}/{n}_{hap1}_to_{hap2}/spiked_consistency.json"
    wildcard_constraints:
        hap1 = '[A-Za-z]+',
        hap2 = '[A-Za-z]+'
    conda: "../envs/truvari.yml"
    threads: 1
    shell:
        """
        truvari consistency --json {input.baseline} {input.spiked} --output {output.tsv} > {output.json}
        """

### Check for variants supported by multiple haplotype reads

rule check_multi_hap_SVs:
    input:
        script = "scripts/python/check_multi_hap_SVs.py",
        spiked = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.vcf.gz",
        hap1_rnames = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap1}_rnames.txt",
        hap2_rnames = "output/alignment/HG002/minimap2/standard/mapped/self/diploid/{hap2}_rnames.txt",
    output:
        mixed = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.mixed_hap.vcf.gz",
        hap1_only = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.hap1_only.vcf.gz",
        hap2_only = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.hap2_only.vcf.gz",
        report = "output/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.check_multi_hap_SVs.json"
    log:
        "logs/alignment/HG002/minimap2/standard/variants/sniffles_{setting}/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.check_multi_hap_SVs.log"
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
        {log}
        """