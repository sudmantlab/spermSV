rule meryl_kmer:
    input:
        refgenome = config['reference']['fasta']
    output:
        meryl_kmer = f"output/mapping/{{refalias}}/winnowmap/standard/meryl/repetitive_k{config['meryl']['k']}_{{refalias}}.txt"
    conda: "../envs/winnowmap.yml"
    params:
        DB = "output/mapping/{refalias}/winnowmap/standard/meryl/merylDB",
        k = config['meryl']['k'],
        distinct = config['meryl']['distinct'],
    shell:
        """
        meryl count k={params.k} output {params.DB} {input.refgenome}
        meryl print greater-than distinct={params.distinct} {params.DB} > {output.meryl_kmer}
        """

rule winnowmap:
    input: 
        meryl_kmer = f"output/mapping/{{refalias}}/winnowmap/standard/meryl/repetitive_k{config['meryl']['k']}_{{refalias}}.txt",
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output:
        temp("output/mapping/{refalias}/winnowmap/standard/{specimen}/{lane}/{smrtcell}.filt.bam")
    conda: "../envs/winnowmap.yml"
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['winnowmap']['readgroup'],
        minQ = config['samtools']['minQ']
    shell:
        """ 
        winnowmap -W {input.meryl_kmer} -x map-pb -a -Y -L --eqx  --cs {params.refgenome} {input.hifi} -R {params.readgroup}  | samtools view -q {params.minQ} -hbT {params.refgenome} -o {output}
        """