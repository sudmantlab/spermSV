from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
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
include: "rules/HG002_hg38_sim.smk"
include: "rules/HG002_graphsim.smk"
ruleorder: diploid_self_mapping > minimap2

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        # expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{hap1}_to_{hap2}.coverage.txt", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),
        # # expand("output/alignment/HG002/minimap2/standard/mapped/reports/gatk/1.0_{hap1}_SVs_to_{hap2}", zip, hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/{hap1}_to_{hap2}.vcf.gz", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{hap1}_to_{hap2}/summary.json", hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/1.0_{hap1}_SVs_to_{hap2}/summary.json", hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{spike}_{hap1}_SVs_to_{hap2}/summary.json", spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{spike}_{hap1}_SVs_to_{hap2}/summary.json", spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{hap1}_to_{hap2}/{files}.jl", files = ['tp-base', 'fp', 'fn'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/1.0_{hap1}_SVs_to_{hap2}/{files}.jl", files = ['tp-base', 'fp', 'fn'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{spike}_{hap1}_SVs_to_{hap2}/{files}.jl", files = ['tp-base', 'fp', 'fn'], spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{spike}_{hap1}_SVs_to_{hap2}/{files}.jl", files = ['tp-base', 'fp', 'fn'], spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/1.0_{hap1}_SVs_to_{hap2}/annotated/{files}.trf.jl", files = ['tp-base', 'fp', 'fn'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{spike}_{hap1}_SVs_to_{hap2}/annotated/{files}.trf.jl", files = ['tp-base', 'fp', 'fn'], spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{spike}_{hap1}_SVs_to_{hap2}/annotated/{files}.trf.jl", files = ['tp-base', 'fp', 'fn'], spike = ['0.3', '0.2', '0.18', '0.15', '0.12', '0.1', '0.08', '0.05', '0.03', '0.01'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['PATERNAL', 'MATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/germline/{hap1}_to_{hap2}/annotated/{files}.trf.jl", files = ['tp-base', 'fp', 'fn'], hap1 = ['PATERNAL', 'MATERNAL'], hap2 = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/annotated/{hap1}_to_{hap2}.trf.jl", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),
        # new stuff for new variant spike implementation
        expand("output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/{hap}.jl", hap = ['MATERNAL', 'PATERNAL']),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/hg38/{hap}/germline/{hap}/{set}/summary.json", hap = ['MATERNAL', 'PATERNAL'], set = ['all', 'sub10kB'])