rule hifiasm:
    input: get_fastqs_per_sample
    output:
        p_ctg_hap1 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.r_utg.noseq.gfa",
        ec = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ec.bin",
        ovlp_reverse = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/assembly/hifiasm/{specimen}/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/assembly/hifiasm/{specimen}/{specimen}.asm"
    threads: 40
    log: 
        "logs/assembly/hifiasm/{specimen}/{specimen}.asm.log"
    conda: 
        "../envs/HiFiAssembly.yml"
    shell: 
        """
        hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1
        """

rule gfaToFa:
    input: "output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.{hap}.p_ctg.gfa"
    output: "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    wildcard_constraints:
        hap = ['hap1', 'hap2']
    log: "logs/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: 
        "gfatools gfa2fa {input} > {output} 2> {log}"

