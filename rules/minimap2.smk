rule minimap2:
    # TODO: For some reason, the wildcard_constraints on refalias breaks the graph construction
    # for new specimens even when set to "^(hg38|CHM13)" and mapping *to* hg38.
    input:
        # Skip HiFiAdapterFilt because Revio preprocesses the reads already
        hifi = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    output: 
        temp("output/alignment/{refalias}/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/mapping.yml"
    threads: 14
    shell: 
        """
        minimap2 --version && minimap2 {params.refgenome} {input.hifi} -t {threads} -ax map-hifi -Y -y -L --eqx --cs --MD -R '{params.readgroup}' | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}
        """

# TODO: Fix the pipe to minQ filter, save unfiltered because this might be messing with multi-mapping
rule minimap2_to_hg38_scaffolded:
    input:
        hifi = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz",
        fa = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta"
    output:
        temp("output/alignment/hg38_scaffolded/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        # minQ = config['samtools']['minQ'] # skip minQ for now
    conda: "../envs/mapping.yml"
    threads: 14
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.fa} {input.hifi} | samtools view -b > {output}
        """

rule minimap2_to_T2T_scaffolded:
    input:
        hifi = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz",
        fa = "output/assembly/hifiasm/{specimen}/T2T_scaffolded/{specimen}.diploid.fasta"
    output:
        temp("output/alignment/T2T_scaffolded/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        # minQ = config['samtools']['minQ'] # skip minQ for now
    conda: "../envs/mapping.yml"
    threads: 14
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.fa} {input.hifi} | samtools view -b > {output}
        """

# rule minimap2_to_haploid_self:
#     input:
#         hifi = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz",
#         hap = "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.scaffold.fasta"
#     output:
#         temp("output/alignment/hg38_scaffolded/minimap2/standard/mapped/temp/{specimen}.{hap}/{lane}/{smrtcell}.filt.bam")
#     params:
#         readgroup = config['minimap2']['readgroup'],
#         minQ = config['samtools']['minQ']
#     conda: "../envs/mapping.yml"
#     threads: 10
#     shell:
#         """
#         minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.hap} {input.hifi} | samtools view -q {params.minQ} -b > {output}
#         """