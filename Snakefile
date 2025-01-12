from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
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
ruleorder: sniffles_mosaic_hg38_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_hg38_scaffolded > sniffles_standard
ruleorder: sniffles_mosaic_T2T_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_T2T_scaffolded > sniffles_standard
ruleorder: preprocess_hg38_scaffolded_variants > preprocess_variants
ruleorder: preprocess_T2T_scaffolded_variants > preprocess_variants

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        # expand("output/alignment/hg38/minimap2/standard/mapped/{specimen}.sorted.merged.bam", specimen = specimens),
        # expand("output/alignment/hg38/minimap2/standard/variants/sniffles_mosaic/{specimen}.vcf.gz", specimen = specimens),
        # expand("output/alignment/hg38/minimap2/duplomap/variants/sniffles_mosaic/{specimen}.vcf.gz", specimen = specimens),
        # # expand("output/alignment/hg38/minimap2/duplomap/coverage_stats/{specimen}/{chr}.coverage.pdf", specimen = specimens, chr = chrs)
        expand("output/alignment/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.filt.vcf.gz", specimen = [x for x in specimens if x not in ['HG002', '900']]),
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz", specimen = [x for x in specimens if x not in ['900']]),
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.filt.vcf.gz", specimen = [x for x in specimens if x != '900']),
        expand("output/alignment/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz", specimen = [x for x in specimens if x not in ['HG002', '900']]),
        expand("output/assembly/hifiasm/{specimen}/quast/T2T_scaffolded/report.html", specimen = [x for x in specimens if x not in ['HG002', '900']]),
        expand("output/assembly/hifiasm/{specimen}/quast/hg38_scaffolded/report.html", specimen = [x for x in specimens if x != '900']),
        expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = [x for x in specimens if x != '900'], hap = ['hap1', 'hap2']),
        'output/assembly/assembly_stats/all.tsv',
        expand("output/alignment/hg38_scaffolded/minimap2/standard/mapped/{specimen}.{suffix}.bam", specimen = [x for x in specimens if x != '900'], suffix = ['sorted.merged', 'sorted.merged.mapQ_modified']),
        expand("output/alignment/T2T_scaffolded/minimap2/standard/mapped/{specimen}.{suffix}.bam", specimen = [x for x in specimens if x not in ['900', 'HG002']], suffix = ['sorted.merged', 'sorted.merged.mapQ_modified']),
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.vcf.gz", specimen = [x for x in specimens if x not in ['900']]),
        expand("output/alignment/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.mapQ_modified.qc_all.vcf.gz", specimen = [x for x in specimens if x not in ['HG002', '900']]),

# obtain rulegraph with:
# snakemake --rulegraph \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/hg38/MATERNAL.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/MATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/truvari/hg38/PATERNAL/standard/all/tp-comp.vcf.gz \
# output/alignment/HG002/minimap2/standard/variants/sniffles_standard/self/PATERNAL/annotated/merged_1_MATERNAL_spike_to_PATERNAL.trf.check_multi_hap_SVs.report \
# | dot -Tpng > simulation_rulegraph.png

