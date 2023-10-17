rule sniffles_to_repeatmasker:
    # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
    # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
    input:
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.filtered.fa"
    output:
        "output/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample/{specimen}.filtered.fa.tbl",
        "output/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample/{specimen}.filtered.fa.out.gff",
        "output/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample/{specimen}.filtered.fa.masked",
        "output/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample/{specimen}.filtered.fa.cat"
    params:
        outdir = "output/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample",
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/mapping/{refalias}/repeatmasker/sniffles/{setting}/single_sample/{specimen}.filtered.log"
    threads: 10
    conda:
        '../envs/sniffles.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} -dir {params.outdir} {input} &> {log}
        """