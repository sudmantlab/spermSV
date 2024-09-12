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
include: "rules/HG002_hg38_sim.smk"
include: "rules/HG002_graphsim.smk"

ruleorder: diploid_self_mapping > minimap2
ruleorder: get_unambiguous_hap_bam > generic_sort
ruleorder: get_unambiguous_hap_bam_hg38 > generic_sort

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        expand("output/alignment/HG002/minimap2/standard/mapped/self/diploid/unambiguous_{hap}.fastq.gz", hap=["PATERNAL", "MATERNAL"]),
        expand("output/alignment/HG002/minimap2/standard/variants/extracted_vars/{hap}/{chr}", hap=["PATERNAL", "MATERNAL"], chr = chrs),
        expand("output/alignment/HG002/minimap2/standard/mapped/self/{hap2}/{n}_{hap1}_to_{hap2}.bam", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"], n = [str(i) for i in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{n}_{hap1}_to_{hap2}/log.txt", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/mosaic/{n}_{hap1}_to_{hap2}/spiked_consistency.json", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/sniffles_mosaic/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.check_multi_hap_SVs.json", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)]),
        # germline
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/standard/{n}_{hap1}_to_{hap2}/log.txt", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/truvari/self/{hap2}/standard/{n}_{hap1}_to_{hap2}/spiked_consistency.json", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/{hap2}/annotated/merged_{n}_{hap1}_to_{hap2}.trf.check_multi_hap_SVs.json", hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"],n = [str(i) for i in range(1, 11)])