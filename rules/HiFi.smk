import os
import pandas as pd

def get_hifiasm_inputs(wildcards):
    hifi_path = "output/HiFiAdapterFilt/{specimen}/{lane}/{smrtcell}.ccs.filt.fastq.gz"
    samples = pd.read_table("samples.tsv", index_col=False, dtype=str)
    samples = samples[samples["Specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(specimen=s[0], lane=s[1], smrtcell = s[2]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

rule hifiasm_noopts:
    version: subprocess.run(["hifiasm --version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').strip('\n')
    input: get_hifiasm_inputs
    output: 
        p_ctg_hap1 = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.r_utg.noseq.gfa",
        ec = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.ec.bin",
        ovlp_reverse = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm/{specimen}/no_opts/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm/{specimen}/no_opts/{specimen}.asm"
    threads: 52
    log: "output/hifiasm/{specimen}/no_opts/{specimen}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    # For reference, -l2 purges haplotigs, but we want to retain those.
    # shell: "hifiasm -l2 -o {params.prefix} -t {threads} {input} > {log} 2>&1"
    shell: "hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1" 

rule gfaToFa:
    version: subprocess.run(["gfatools version"], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split('\n')[0]
    input: "output/hifiasm/{specimen}/no_opts/{specimen}.asm.bp.{genometype}.gfa"
    output: "output/hifiasm-fasta/{specimen}/no_opts/{specimen}.{genometype}.fa"
    log: "logs/gfaToFa/{specimen}/no_opts/{specimen}.{genometype}.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "gfatools gfa2fa {input} > {output} 2> {log}"

