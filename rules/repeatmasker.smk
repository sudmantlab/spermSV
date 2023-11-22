rule alt_to_fasta:
    # Takes a (sniffles) vcf and writes INS alleles to a fasta file for RepeatMasker. 
    input: 
        vcf='output/mapping/{refalias}/sniffles/mosaic/single_sample/{specimen}.vcf.gz'
    output: 
        fa = 'analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa'
    threads:
        10
    run:
        def write_fasta(vcf, outfile):
            # Takes the alt output of Sniffles2 and writes a fasta file for input through RepeatMasker.
            table = vcf[['ID', 'ALT']]
            with open(outfile, 'w') as f:
                for header, seq in table.to_records(index=False):
                    f.write(f'>{header}\n{seq}\n')
        
        vcf = pd.read_table(input.vcf, skiprows = 241) # skip vcf header
        print("Writing ALT sequences to ", output.fa)
        write_fasta(vcf, output.fa)
        
rule annotate_alts:
    # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
    # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
    input:
        "analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa"
    output:
        "analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa.tbl",
        "analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa.out.gff",
        "analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa.masked",
        "analysis/{refalias}/{analysis}/repeatmasker/{specimen}.alt.fa.cat"
    params:
        outdir = "analysis/{refalias}/{analysis}/repeatmasker",
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/analysis/{refalias}/{analysis}/repeatmasker/{specimen}.log"
    threads: 10
    conda:
        '../envs/RepeatMasker.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} -dir {params.outdir} {input} &> {log}
        """