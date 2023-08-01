rule samtools_faidx:
    input: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.fa"
    output: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.fa.fai"
    conda: "../envs/omni-c.yaml"
    shell: "samtools faidx {input}"
    

rule genome_file:
    input: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.fa.fai"
    output: "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.genome"
    shell: "cut -f1,2 {input} > {output}"


rule bwa_index:
    input:
        "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.fa",
    output:
        idx=multiext("output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{species}/{ccs_opts}/{hifiasm_opts}/{species}.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.0.0/bio/bwa/index"

def get_omnic_reads(wildcards):
    path_trimmed = "output/trimmed-hic/{species}/{sample_name}-trimmed_{read}.fastq.gz"
    input_list = []
    samples = pd.read_table("rna_pepsamples.tsv", index_col=False)
    samples = samples[samples.species==wildcards.species][samples.type == "Hi-C"]
    samples_grouped = samples.groupby(samples.sample_name)
    for sample in set(samples["sample_name"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        if len(sample_subset) == 0:
            raise Exception("No files available for sample {}".format(sample))
        input_list.extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R1"].itertuples()])
        input_list.extend([path_trimmed.format(species = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R2"].itertuples()])
    return input_list


rule bwa_aln:
    input:
        fastq= get_omnic_reads,
        # Index can be a list of (all) files created by bwa, or one of them
        idx=multiext("output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.sam",
    params:
        extra="-5SP -T0",
    log:
        "logs/bwa_aln/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.log",
    threads: 32
    wrapper:
        "v1.0.0/bio/bwa/aln"


rule pairtools_parse:
    input: 
        sam = "output/Omni-C_BWAAligned/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.sam",
        genome = "output/hifiasm-fasta/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.fa"
    output: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.pairsam"
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        pairtools parse \
          --min-mapq 40 \
          --walks-policy 5unique \
          --max-inter-align-gap 30 \
          --nproc-in {threads} \
          --nproc-out {threads} \
          --chroms-path {input.genome} \
          {input.sam} > {output} 2> {log}
        """


rule pairtools_sort:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.aligned.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.pairsam")
    log: 
        "logs/pairtools_sort/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        mkdir /global/scratch2/mvazquez/tmp_pairtools
        pairtools sort \
          --tmpdir=/global/scratch2/mvazquez/tmp_pairtools \
          --nproc {threads} \
          {input} > {output} 2> {log}
        rm -r /global/scratch2/mvazquez/tmp_pairtools
        """


rule pairtools_dedup:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.pairsam"
    output: 
        temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.dedup.pairsam")
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.dedup.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        pairtools dedup \
         --nproc-in {threads} \
         --nproc-out {threads} \
         --mark-dups \
         --output-stats {output}.stats \
         --output {output} {input} 2> {log}
        """


rule pairtools_split:
    input: 
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.dedup.pairsam"
    output: 
        bam = temp("output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.bam"),
        pairs = "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.pairs",
    log: 
        "logs/pairtools_parse/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.split.log"
    conda: "../envs/omni-c.yaml"
    threads: 32
    shell: 
        """
        pairtools split \
        --nproc-in {threads} \
        --nproc-out {threads} \
        --output-pairs {output.pairs} \
        --output-sam {output.bam} \
        {input} 2> {log}
        """


rule samtools_sort:
    input:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.bam"
    output:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.bam"
    params:
        extra = "-m 100G"
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@.
    wrapper:
        "v1.0.0/bio/samtools/sort"


rule samtools_index:
    input:
        "output/Omni-C_pairsam/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.bam"
    output:
        "output/Omni-C_mapped/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.bam.bai"
    log:
        "output/Omni-C_mapped/{species}/{ccs_opts}/{hifiasm_opts}/{species}.p_ctg.sorted.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        32     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"
