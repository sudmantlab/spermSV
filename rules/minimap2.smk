rule minimap2:
    # TODO: For some reason, the wildcard_constraints on refalias breaks the graph construction
    # for new specimens even when set to "^(hg38|CHM13)" and mapping *to* hg38.
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output: 
        temp("output/alignment/{refalias}/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
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

rule minimap2_to_diploid_self:
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz",
        hap1 = "output/assembly/hifiasm/{specimen}/{specimen}.hap1.scaffold.fasta",
        hap2 = "output/assembly/hifiasm/{specimen}/{specimen}.hap2.scaffold.fasta"
    output:
        temp("output/alignment/self_assembly/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' <(cat {input.hap1} {input.hap2}) {input.hifi} | samtools view -q {params.minQ} -b > {output}
        """

rule minimap2_to_haploid_self:
    input:
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz",
        hap = "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.fasta"
    output:
        temp("output/alignment/self_assembly/minimap2/standard/mapped/temp/{specimen}.{hap}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 10
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.hap} {input.hifi} | samtools view -q {params.minQ} -b > {output}
        """