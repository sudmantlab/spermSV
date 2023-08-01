import os

rule ccs_consensus_minPasses_minRQ_chunk:
    # version: subprocess.check_output("ccs --version | head -1", shell=True)
    input: "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.subreads.bam"
    output: 
        ccs="output/HiFi-CCS/{specimen}/{lane}/{smrtcell}.ccs.{chunk}_40.bam",
        report="output/HiFi-CCS/{specimen}/{lane}/{smrtcell}.ccs_report.{chunk}_40.txt"
    log: "logs/HiFi-CCS/{specimen}/{lane}/{smrtcell}.{chunk}.logs"
    params:
        loglevel = "DEBUG",
        minPasses = 3,
        minRQ = 0.99
    wildcard_constraints:
        chunk="[0-9]?[0-9]",
        species="[A-Z]_[a-z]+"
    threads: 4
    conda: "../envs/CCS.yml"
    shell: "ccs -j {threads} --min-passes {params.minPasses} --min-rq {params.minRQ} --report-file {output.report} --log-level {params.loglevel} --log-file {log} --chunk {wildcards.chunk}/40 {input} {output.ccs}"

rule ccs_chunk_merge:
    input: expand("output/HiFi-CCS/{specimen}/{lane}/{smrtcell}.ccs.{chunk}_40.bam", chunk=range(1,41), allow_missing = True)
    output: "output/HiFi-CCS/{specimen}/{lane}/{smrtcell}.ccs.bam"
    wildcard_constraints:
        isDefault="(defaults-)?",
        chunk="\d{1,2}",
        species="[A-Z]_[a-z]+"
    threads: 32
    conda: "../envs/CCS.yml"
    shell: "samtools merge -@{threads} {output} {input}"