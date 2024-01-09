rule minimap2:
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: temp("output/alignment/{refalias}/minimap2/standard/mapped/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10
    shell: 
        """
        minimap2 --version && minimap2 {params.refgenome} {input.hifi} -t {threads} -ax map-hifi -Y -y -L --eqx --cs --MD -R '{params.readgroup}' | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}
        """