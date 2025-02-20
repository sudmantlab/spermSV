from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "config/snakemake/hg38.config.yml"
workdir: config['workdir']
localrules: make_bam_symlink, make_fastq_symlink

# include helper functions
include: "rules/common.smk"
include: "rules/make_symlinks.smk"

# preprocessing steps
include: "rules/preprocessing.smk"

# Assembly
include: "rules/assembly.smk"
include: "rules/get_assembly_stats.smk"
include: "rules/flagger.smk"

# Mapping
include: "rules/minimap2.smk"
include: "rules/winnowmap.smk"
include: "rules/duplomap.smk"
include: "rules/somrit.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Graph methods
include: "rules/minigraph-cactus.smk"

# Variant calling
include: "rules/sniffles.smk"
# include: "rules/pmdv.smk"
# include: "rules/straglr.smk"
# include: "rules/trgt.smk"
include: "rules/multimap_analysis.smk"

# Simulations & Benchmarking
include: "rules/spike_in.smk"
include: "rules/spike_in_svtier1.smk"
include: "rules/spike_in_cmrg.smk"
include: "rules/in_silico_simulations.smk"
include: "rules/HG002_sim.smk"
# include: "rules/HG002_hg38_sim.smk"
include: "rules/HG002_graphsim.smk"

ruleorder: hifiasm_900 > hifiasm
ruleorder: minimap2_to_hg38_scaffolded > minimap2
ruleorder: minimap2_to_T2T_scaffolded > minimap2
ruleorder: disambiguate_T2T_self_mapped > generic_sort
ruleorder: disambiguate_hg38_remapped > generic_sort
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: preprocess_scaffolded_variants > preprocess_variants

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        # assembly
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = [x for x in specimens if x != '900']),
        expand("output/assembly/hifiasm/{specimen}/quast/T2T_scaffolded/report.html", specimen = [x for x in specimens if x != '900']),
        expand("output/assembly/hifiasm/{specimen}/quast/hg38_scaffolded/report.html", specimen = [x for x in specimens if x != '900']),
        expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = specimens, hap = ['hap1', 'hap2']),
        expand("output/alignment/hprc_personalized/mapped/hprc-v1.1-mc-chm13.d9/{specimen}.surjected.bam", specimen = specimens),
        'output/assembly/assembly_stats/all.tsv',
        # coverage stats
        expand("output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}.coverage.html", ref = ['hg38', 'T2T'], specimen = [x for x in specimens if x != '900']),
        # variant calling
        # expand("analysis/{refalias}/denovo/files/minimap2/standard/variants/all.active.single_type.low_div.repeatmasker_insertions.tsv", refalias = ['hg38', 'hg38_scaffolded', 'T2T_scaffolded']),
        # expand("analysis/{refalias}/denovo/files/minimap2/standard/variants/all.qc_all.active.single_type.low_div.repeatmasker_insertions.tsv", refalias = ['hg38', 'hg38_scaffolded', 'T2T_scaffolded']),
        # "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.qc_all.active.single_type.low_div.repeatmasker_insertions.tsv",
        # "analysis/hprc_personalized/denovo/files/giraffe/longread/variants/hprc_personalized-all.active.single_type.low_div.repeatmasker_insertions.tsv",
        # "analysis/T2T_scaffolded/denovo/files/minimap2/standard/variants/T2T_scaffolded-all.qc_all.active.single_type.low_div.repeatmasker_insertions.tsv",
        # "analysis/T2T_scaffolded/denovo/files/minimap2/standard/variants/T2T_scaffolded-all.active.single_type.low_div.repeatmasker_insertions.tsv"

# obtain rulegraph with:
# snakemake --rulegraph \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/MATERNAL.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/MATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/PATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/PATERNAL/annotated/merged_1_MATERNAL_spike_to_PATERNAL.trf.check_multi_hap_SVs.report \
# | dot -Tpng > simulation_rulegraph.png

