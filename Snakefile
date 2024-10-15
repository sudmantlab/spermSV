from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
import glob
import os
import numpy as np
import pandas as pd

configfile: "config/snakemake/hg38.config.yml"
workdir: config['workdir']

# include helper functions
include: "rules/common.smk"
include: "rules/make_symlinks.smk"

# preprocessing steps
include: "rules/preprocessing.smk"

# Assembly
include: "rules/assembly.smk"

# Mapping
include: "rules/minimap2.smk"
include: "rules/winnowmap.smk"
include: "rules/duplomap.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"
include: "rules/pmdv.smk"
include: "rules/straglr.smk"
include: "rules/trgt.smk"

# Simulations & Benchmarking
include: "rules/spike_in.smk"
include: "rules/spike_in_svtier1.smk"
include: "rules/spike_in_cmrg.smk"
include: "rules/in_silico_simulations.smk"
include: "rules/HG002_sim.smk"
# include: "rules/HG002_hg38_sim.smk"
include: "rules/HG002_graphsim.smk"

ruleorder: T2T_self_mapping > minimap2
ruleorder: disambiguate_T2T_self_mapped > generic_sort
ruleorder: disambiguate_hg38_remapped > generic_sort

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{benchmark}_{hap}_benchmark/{file}/{output}", 
        benchmark = ["dipcall", "CMRG"], hap = ["PATERNAL", "MATERNAL"], file = ["unphased", "homozygous","all_PATERNAL", "all_MATERNAL", "unambiguous_PATERNAL", "unambiguous_MATERNAL"], output = ["log.txt", "tp-comp.vcf.gz"]),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{benchmark}_{hap}_benchmark/{file}/{output}", 
        benchmark = ["CMRG"], hap = ["PATERNAL", "MATERNAL"], file = ["unphased", "homozygous","all_PATERNAL", "all_MATERNAL", "unambiguous_PATERNAL", "unambiguous_MATERNAL"], output = ["log.txt", "tp-comp.vcf.gz"]),
        expand("output/alignment/HG002/minimap2/standard/variants/simulation/sniffles_mosaic/CMRG_{hap1}_benchmark/{hap1}_spike_to_{hap2}/1_{hap1}_spike_to_{hap2}/log.txt",
        benchmark = ["CMRG"], hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),
        # # benchmarks here
        # expand("benchmarks/HG002/liftover/HG002_GRCh38_CMRG_SV_v1.00.to_HG002_{hap}.{ext}", hap = ["PATERNAL", "MATERNAL"], ext = ["vcf.gz", "bed"]),
        # expand("benchmarks/HG002/liftover/{hap1}_HG002_GRCh38_CMRG_SV_v1.00.to_HG002_{hap2}.{ext}", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"], ext = ["vcf.gz"]),
        # # this one appears to fail due to a "malformed allele? see logs. this is why there's no dipcall simulations rn
        # # expand("benchmarks/HG002/liftover/GRCh38_HG2-T2TQ100-V1.0.to_HG002_{hap}.{ext}", hap = ["PATERNAL", "MATERNAL"], ext = ["vcf.gz"]),
        # # expand("benchmarks/HG002/liftover/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.to_HG002_{hap}.bed", hap = ["PATERNAL", "MATERNAL"]),
        # # expand("benchmarks/HG002/liftover/{hap1}_GRCh38_HG2-T2TQ100-V1.0.to_HG002_{hap2}.{ext}", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"], ext = ["vcf.gz"])

# obtain rulegraph with:
# snakemake --rulegraph \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/MATERNAL.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/MATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/PATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/PATERNAL/annotated/merged_1_MATERNAL_spike_to_PATERNAL.trf.check_multi_hap_SVs.report \
# | dot -Tpng > simulation_rulegraph.png