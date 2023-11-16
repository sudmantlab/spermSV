rule annotate_alts:
    # Takes a fasta file containing ALT allele sequences from Sniffles2 output, then analyzes the sequences using RepeatMasker.
    # Relevant outputs: a .tbl (difficult to parse) and a .gff (easier to parse) denoting identified repeat motifs.
    input:
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.alt.fa"
    output:
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.alt.fa.tbl",
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.alt.fa.out.gff",
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.alt.fa.masked",
        "output/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.alt.fa.cat"
    params:
        outdir = "output/mapping/{refalias}/sniffles/{setting}/single_sample",
        species = config['repeatmasker']['species'],
        engine = config['repeatmasker']['engine']
    log:
        "logs/mapping/{refalias}/sniffles/{setting}/single_sample/{specimen}.repeatmasker.log"
    threads: 10
    conda:
        '../envs/RepeatMasker.yml'
    shell:
        """
        RepeatMasker -pa {threads} -engine {params.engine} -nocut -gff \
        -species {params.species} -dir {params.outdir} {input} &> {log}
        """