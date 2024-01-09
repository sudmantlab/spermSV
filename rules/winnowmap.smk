rule meryl_kmer:
    input:
        refgenome = config['reference']['fasta']
    output:
        meryl_kmer = f"output/alignment/{{refalias}}/winnowmap/standard/mapped/meryl/repetitive_k{config['meryl']['k']}_{{refalias}}.txt"
    conda: "../envs/winnowmap.yml"
    threads: 20
    params:
        DB = "output/alignment/{{refalias}}/winnowmap/standard/mapped/meryl/merylDB",
        k = config['meryl']['k'],
        distinct = config['meryl']['distinct'],
    shell:
        """
        meryl count k={params.k} output {params.DB} {input.refgenome}
        meryl print greater-than distinct={params.distinct} {params.DB} > {output.meryl_kmer}
        """

rule winnowmap:
    input: 
        meryl_kmer = f"output/alignment/{{refalias}}/winnowmap/standard/mapped/meryl/repetitive_k{config['meryl']['k']}_{{refalias}}.txt",
        hifi = "output/preprocessing/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    output:
        temp("output/alignment/{refalias}/winnowmap/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    conda: "../envs/winnowmap.yml"
    threads: 10
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['winnowmap']['readgroup'],
        minQ = config['samtools']['minQ']
    shell:
        """ 
        winnowmap -t {threads} -W {input.meryl_kmer} -x map-pb -a -Y -y -L --eqx --cs --MD {params.refgenome} {input.hifi} -R {params.readgroup}  | samtools view -q {params.minQ} -hbT {params.refgenome} -o {output}
        """