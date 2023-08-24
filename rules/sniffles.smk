rule sniffles_single_sample:
    input:
        bam = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.vcf.gz',
        snf='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.snf',
        tbi='output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats']
    log:
        "logs/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.bench.log"
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} \
        --reference {params.refgenome} --tandem-repeats {params.repeats} \
        --threads {threads}
        """

rule mosaic_single_sample:
    # Calls mosaic (somatic) SVs using the --mosaic option.
    input:
        bam = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam",
        index = "output/mapping/{refalias}/minimap2/standard/{specimen}.sorted.merged.bam.bai"
    output:
        vcf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz',
        snf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.snf',
        tbi='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz.tbi'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats']
    log:
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.bench.log"
    shell:
        """
        sniffles --input {input.bam} --vcf {output.vcf} --snf {output.snf} \
        --reference {params.refgenome} --tandem-repeats {params.repeats} \
        --threads {threads} --mosaic
        """

rule sniffles_multi_sample:
    # Note that this will call sniffles multi-input on *all* samples specified in the samples table.
    # If necessary, will build in support later for specifiying a certain subset of samples using a configfile.
    input: 
        expand('output/mapping/{refalias}/sniffles/standard/single_sample/{specimen}.snf', specimen = specimens, allow_missing = True)
    output: 
        'output/mapping/{refalias}/sniffles/standard/multi_sample/multi_sample.vcf'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats']
    log:
        "logs/mapping/{refalias}/sniffles/standard/multi_sample/multi_sample.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/standard/multi_sample/multi_sample.bench.log"
    shell:
        """
        sniffles --input {input} --vcf {output} \
        --reference {params.refgenome} --tandem-repeats {params.repeats} \
        --threads {threads}
        """

rule mosaic_multi_sample:
    # Same as the sniffles_multi_sample rule, but with mosaic (somatic) multi-call.
    input: 
        expand('output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.snf', specimen = specimens, allow_missing = True)
    output: 
        'output/mapping/{refalias}/sniffles/mosaic/multi_sample/multi_sample.vcf'
    conda:
        '../envs/sniffles.yml'
    threads:
        20
    resources:
        mem_mb=60000
    params:
        refgenome = config['reference']['fasta'],
        repeats = config['reference']['annotations']['repeats']
    log:
        "logs/mapping/{refalias}/sniffles/mosaic/multi_sample/multi_sample.log"
    benchmark:
        "logs/mapping/{refalias}/sniffles/mosaic/multi_sample/multi_sample.bench.log"
    shell:
        """
        sniffles --input {input} --vcf {output} \
        --reference {params.refgenome} --tandem-repeats {params.repeats} \
        --threads {threads} --mosaic
        """

rule filter_multi_vcf:
    input: 
        vcf='output/mapping/{refalias}/sniffles/{setting}/multi_sample/multi_sample.vcf'
    output: 
        filtered= 'output/mapping/{refalias}/sniffles/{setting}/multi_sample/multi_sample.filtered.vcf',
        fa = 'output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa'
    params: 
        repeatmasker_dir = "output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker"
    threads:
        10
    run:
        def write_fasta(vcf, outfile):
            # Takes the alt output of Sniffles2 and writes a fasta file for input through RepeatMasker.
            table = vcf[['ID', 'ALT']]
            with open(outfile, 'w') as f:
                for header, seq in table.to_records(index=False):
                    f.write(f'>{header}\n{seq}\n')

        vcf = pd.read_table(input.vcf, skiprows = 72) # skip vcf header
        filtered = vcf.copy()
        filtered = filtered[filtered['FILTER'] == 'PASS']
        filtered = filtered[~filtered['INFO'].str.contains('IMPRECISE')]
        filtered = filtered[~filtered['INFO'].str.contains('BND')]
        filtered = filtered[filtered['QUAL'].astype(int) >= 30]

        print("Writing filtered vcf to ", output.filtered)
        filtered.to_csv(output.filtered, index = False, sep ="\t")
        # # tack the header back on...
        # shell("head -72 {input.vcf} | cat - {output.filtered} > {output.filtered}")

        shell("mkdir -p {params.repeatmasker_dir}")
        print("Writing ALT sequences to ", output.fa)
        write_fasta(filtered, output.fa)

rule filtered_repeatmasker:
    input:
        'output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa'
    output:
        "output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa.tbl",
        "output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa.out.gff",
        "output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa.masked",
        "output/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.fa.cat"
    params:
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/mapping/{refalias}/sniffles/{setting}/multi_sample/repeatmasker/multi_sample.filtered.log"
    threads: 10
    conda:
        '../envs/sniffles.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} {input} &> {log}
        """