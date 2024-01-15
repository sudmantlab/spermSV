svision_path = "/global/scratch/users/stacy-l/spermSV/code/SVision"

if os.system("echo $PATH | grep SVision") == 256:
    os.environ["PATH"] += os.pathsep + svision_path

rule svision:
    input:
        bam = "output/hg38_no_alts/minimap2/standard/{specimen}/{specimen}.sorted.merged.bam",
        # make sure modelfiles have been downloaded from google drive (yeah)
        # and placed in a new model directory in the repo
        modelfiles = expand("code/SVision/model/svision-cnn-model.ckpt.{suffix}", suffix = ["data-00000-of-00001", "index", "meta"])
    output:
        "output/hg38_no_alts/svision/{specimen}.svision.s5.graph.vcf",
        "output/hg38_no_alts/svision/{specimen}.graph_exactly_match.txt",
        "output/hg38_no_alts/svision/{specimen}.graph_symmetry_match.txt",
        # not including the graph CSV.gfa output right now because I have no idea where that's coming out?
        # praying for no overwrite issues (SCREAM)
    threads: 20
    conda: "../envs/svision.yml" 
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),,
        refgenome = "/global/scratch/users/stacy-l/references/hg38_ucsc/hg38.fa",
        script_path = "code/SVision/SVision",
        modeldir = "code/SVision/model"
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.script_path} -o {params.outdir} \
        -b {input.bam} \
        -m {params.modeldir}/svision-cnn-model.ckpt \
        -g {params.refgenome} \
        -n {wildcards.specimen} -s 5 --graph --qname
        """