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
include: "rules/somrit.smk"
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
        expand('output/alignment/HG002/minimap2/standard/variants/simulation/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/agg_{hap1}_spike_to_{hap2}/callset.vcf.gz', setting = ["mosaic", "standard"], benchmark = ["CMRG", "dipcall"], hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),
#         expand("output/alignment/hg38/minimap2/somrit/mapped/{specimen}_realigned_classified_filtered.tsv", specimen=specimens),
        expand("output/alignment/HG002/minimap2/standard/variants/simulation/sniffles_mosaic/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/{n}_{hap1}_spike_to_{hap2}/multi_hap.report", benchmark = ["CMRG", "dipcall"], hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"], n = [str(x) for x in range(1, 11)]),
        expand("output/alignment/HG002/minimap2/standard/variants/simulation/sniffles_{setting}/{benchmark}_{hap1}_benchmark/{hap1}_spike_to_{hap2}/agg_{hap1}_spike_to_{hap2}/spiked_consistency.json", setting = ["mosaic", "standard"], benchmark = ["CMRG", "dipcall"], hap1 = ["PATERNAL", "MATERNAL"], hap2 = ["PATERNAL", "MATERNAL"]),

# obtain rulegraph with:
# snakemake --rulegraph \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/MATERNAL.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/MATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/PATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/PATERNAL/annotated/merged_1_MATERNAL_spike_to_PATERNAL.trf.check_multi_hap_SVs.report \
# | dot -Tpng > simulation_rulegraph.png