import pandas as pd
configfile: "config/snakemake/hprc.config.yml"
workdir: config['workdir']
localrules: hprc_merge_gams

rule hprc_vg_index:
    input:
        # download from HPRC
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz"
    output:
        "/global/scratch/users/stacy-l/references/HPRC/{prefix}.longread.zipcodes", 
        "/global/scratch/users/stacy-l/references/HPRC/{prefix}.longread.withzip.min"
    log:
        "/global/scratch/users/stacy-l/spermSV/logs/references/HPRC/{prefix}.index.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 40
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home $(dirname {input.gbz}) \
        {params.container} \
        vg autoindex --workflow lr-giraffe \
        --prefix {wildcards.prefix} \
        --gbz {input.gbz}) 2> {log}
        """

rule hprc_vg_giraffe_align:
    input:
        fastq = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz",
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz",
        indices = ["/global/scratch/users/stacy-l/references/HPRC/{prefix}.longread.zipcodes", "/global/scratch/users/stacy-l/references/HPRC/{prefix}.longread.withzip.min"]
    output:
        gam = temp("output/alignment/{prefix}/giraffe/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.gam"),
        report = "logs/alignment/{prefix}/giraffe/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.report.tsv"
    log:
        "logs/alignment/{prefix}/giraffe/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg giraffe -b hifi \
        -Z {input.gbz} \
        -f {input.fastq} \
        --report-name {output.report} \
        -p > {output.gam}) 2> {log}
        """

def get_temp_gams_per_sample(wildcards, sample_table = config['sample_table']) -> List[str]:
    """
    A helper function for creating a list of expected per-smrtcell gam outputs per sample to be merged into a single gam.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    basepath = "output/alignment/{prefix}/giraffe/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.mapped.gam"
    input_samples = [basepath.format(prefix=wildcards.prefix, specimen=s[0], lane=s[2], smrtcell=s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    return input_samples

rule merge_gams:
    input:
        get_temp_gams_per_sample
    output:
        "output/alignment/{prefix}/giraffe/mapped/{specimen}.mapped.gam"
    threads: 1
    shell:
        "cat {input} > {output}"

rule hprc_vg_stats:
    input:
        gam = "output/alignment/{prefix}/giraffe/mapped/{specimen}.mapped.gam"
    output:
        stats = "output/alignment/{prefix}/giraffe/mapped/{specimen}.mapped.stats"
    log:
        "logs/alignment/{prefix}/giraffe/mapped/{specimen}.mapped.stats.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 8
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg stats -a {input.gam} > {output.stats}) 2> {log}
        """

rule hprc_vg_surject:
    input:
        gam = "output/alignment/{prefix}/giraffe/mapped/{specimen}.mapped.gam",
        gbz = "/global/scratch/users/stacy-l/references/HPRC/{prefix}.gbz"
    output:
        bam = "output/alignment/{prefix}/giraffe/mapped/{specimen}.surjected.bam"
    log:
        "logs/alignment/{prefix}/giraffe/mapped/{specimen}.surjected.log"
    params:
        home = config["workdir"],
        container = config["singularity"]["vg"]["container"],
        bind_paths = " ".join([f"-B {path}:{path}" for path in config["singularity"]["vg"]["bind_paths"]]),
        singularity_args = config["singularity"]["vg"]["args"]
    threads: 16
    shell:
        """
        (singularity exec {params.singularity_args} \
        {params.bind_paths} \
        --home {params.home} \
        {params.container} \
        vg surject -b \
        -x {input.gbz} \
        {input.gam} \
        --prune-low-cplx > {output.bam}) 2> {log}
        """