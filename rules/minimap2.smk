rule minimap2_standard: 
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: temp("output/mapping/{refalias}/minimap2/standard/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/minimap2.yml"
    threads: 40
    shell: "minimap2 --version && minimap2 -R '{params.readgroup}' -t {threads} -ax map-hifi {params.refgenome} {input.hifi} | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}"

rule minimap2_softclip:
    # Runs minimap2 with -Y setting to soft-clip supplementary alignments.
    # Required by straglr.
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: temp("output/mapping/{refalias}/minimap2/softclip/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/minimap2.yml"
    threads: 40
    shell: "minimap2 --version && minimap2 -Y -R '{params.readgroup}' -t {threads} -ax map-hifi {params.refgenome} {input.hifi} | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}"