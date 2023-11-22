rule survivor_create_error_profile:
    # Creates an error profile from a fractional subsample of the HPRC HG002 (Revio CCS) dataset.
    input:
        'output/mapping/hg38/winnowmap/standard/{control}.sorted.merged.bam'
    output:
        'output/mapping/hg38/simulations/in_silico/SURVIVOR/{control}_error_profile.txt'
    conda:
        "../envs/sniffles.yml"
    threads: 40
    params:
        minreadlen = config['simulations']['survivor']['error_profile']['minreadlen'],
        fraction = config['simulations']['survivor']['error_profile']['fraction']
    shell:
        """
        samtools view -@ {threads} -s {params.fraction} {input} | SURVIVOR scanreads {params.minreadlen} {output}
        """

rule sim_it_preproc_mapping:
    # Uses winnowmap to map a subset of control (HPRC HG002 Revio CCS) reads to a defined reference.
    input:
        reads = 'output/mapping/hg38/simulations/in_silico/Sim-it/{control}_{chr}_subset.fasta'
    output:
        # placeholder until I know what gets created here
        'output/mapping/hg38/simulations/in_silico/Sim-it/{control}_{chr}_status.success'
    conda:
        "../envs/Sim-it.yml"
    threads: 1
    shell:
        """
        perl code/Sim-it/train_error_profile.pl -ref {input.reference} -reads {input.reads}
        """

rule sim_it_create_error_profile:
    # Creates an error profile from a subset of reads from the control (HPRC HG002 Revio CCS) dataset.
    input:
        reads = 'output/mapping/hg38/simulations/in_silico/Sim-it/{control}_{chr}_subset.fasta'
    output:
        # placeholder until I know what gets created here
        'output/mapping/hg38/simulations/in_silico/Sim-it/{control}_{chr}_status.success'
    conda:
        "../envs/Sim-it.yml"
    params:
        reference = config['simulations']
    threads: 1
    shell:
        """
        perl code/Sim-it/train_error_profile.pl -ref {params.reference} -reads {input.reads}
        """

rule survivor_simgenome:
    # Using a custom error profile, reference genome, & config (params) file, generates a modified .fasta 
    # containing SVs defined by the params file. Simulated variants are described with the bed & vcf files,
    # and their identity is relative to the provided reference (unmodified) genome fasta ("option 0").
    # Note that this requires a .fa/.fasta reference file: compressed .fa.gz and etc formats will not work.
    input:
        error_profile = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/HG002_error_profile.txt',
        params = 'config/packages/SURVIVOR/{sim_id}.params'
    output:
        bed = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.bed',
        vcf = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.vcf',
        fasta = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.fasta',
        insertions = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.insertions.fa'
    conda:
        "../envs/sniffles.yml"
    threads: 1
    params:
        refgenome = config['simulations']['survivor']['refgenome'],
        snp_fraction = config['simulations']['survivor']['snp_fraction'],
        outdir = 'output/mapping/hg38/simulations/in_silico/SURVIVOR'
    shell:
        """
        SURVIVOR simSV {params.refgenome} {input.params} {params.snp_fraction} 0 {params.outdir}/{wildcards.sim_id}/{wildcards.sim_id}
        """

rule survivor_simreads:
    # Using a custom error profile, reference genome, & config (params) file, generates a modified .fasta 
    # containing SVs defined by the params file. Simulated variants are described with the bed & vcf files,
    # and their identity is relative to the provided reference (unmodified) genome fasta ("option 0").
    # Note that this requires a .fa/.fasta reference file: compressed .fa.gz and etc formats will not work.
    # This rule is triggered by the expand() call in merge_survivor_svsimreads for parallelization purposes.
    input:
        error_profile = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/HG002_error_profile.txt',
        fasta = 'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.fasta',
    output:
        temp('output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.temp_{n_sim}.fasta')
    conda:
        "../envs/sniffles.yml"
    threads: 1
    params:
        outdir = 'output/mapping/hg38/simulations/in_silico/SURVIVOR'
    shell:
        """
        echo "Generating simulation reads, process {wildcards.n_sim}..."
        SURVIVOR simreads {input.fasta} {input.error_profile} 1 {output}
        """

rule merge_survivor_svsimreads:
    # This rule triggers parallel generation of simulated reads by SURVIVOR w/ the survivor_simreads rule.
    # Note that target coverage is set in the configfile.
    input:
        expand('output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.temp_{n_sim}.fasta', allow_missing = True, 
                n_sim = np.arange(0, config['simulations']['survivor']['coverage']))
    output:
        'output/mapping/hg38/simulations/in_silico/SURVIVOR/{sim_id}/{sim_id}.simreads.fasta'
    conda:
        "../envs/sniffles.yml"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule sim_it_refreads:
    # Given a config file, runs Sim-it to create fastq files of simulated reference (no SV) reads.
    input:
        config = "config/packages/Sim-it/ref_reads.txt"
    output:
        bam = "output/mapping/hg38/simulations/in_silico/Sim-it/ref_reads/ref_reads_HAP12.bam",
        reads = expand("output/mapping/hg38/simulations/in_silico/Sim-it/ref_reads/Long_reads_ref_reads_{hap}.fasta", allow_missing = True, hap = ['HAP1', 'HAP12', 'HAP2']),
        log = "output/mapping/hg38/simulations/in_silico/Sim-it/ref_reads/log_ref_reads.txt"
    conda:
        "../envs/Sim-it.yml"
    threads: 1
    params:
        outdir = "output/mapping/hg38/simulations/in_silico/Sim-it"
    shell:
        """
        perl code/Sim-it/Sim-it1.3.4.pl -c {input.config} -o {params.outdir}/ref_reads
        """

rule sim_it_svsimreads:
    # Given a config file, runs Sim-it to create fastq files of simulated SV-containing reads.
    # Uses the default Sequel error profile for PacBio reads, unclear how to train a custom error profile right now :(
    input:
        config = "config/packages/Sim-it/{sim_id}.txt"
    output:
        reads = expand("output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/Long_reads_{sim_id}_{hap}.fasta", allow_missing = True, hap = ['HAP1', 'HAP12', 'HAP2']),
        log = "output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/log_{sim_id}.txt", 
        graphs = expand("output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/graph_{svtype}_{sim_id}.txt", allow_missing = True, svtype = ['DEL', 'INS', 'INV']),
        ref = "output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/{sim_id}.fasta",
        hap_refs = expand("output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/{sim_id}_haplotype{hap}.fasta", allow_missing = True, hap = [1, 2]),
        vcf = "output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/{sim_id}.vcf",
        # bam = "output/mapping/hg38/simulations/in_silico/Sim-it/{sim_id}/{sim_id}_HAP12.bam"
    conda:
        "../envs/Sim-it.yml"
    threads: 1
    params:
        outdir = "output/mapping/hg38/simulations/in_silico/Sim-it"
    shell:
        """
        perl code/Sim-it/Sim-it1.3.4.pl -c {input.config} -o {params.outdir}/{wildcards.sim_id}
        """