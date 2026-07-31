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

# Assembly and QC
include: "rules/assembly.smk"
# include: "rules/get_assembly_stats.smk"

# Alignment (and realignment)
include: "rules/minimap2.smk"
include: "rules/duplomap.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"

# Preliminary graph variant calling
include: "rules/minigraph-cactus.smk"

ruleorder: minimap2_to_hg38_scaffolded > minimap2
ruleorder: minimap2_to_T2T_scaffolded > minimap2
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: preprocess_scaffolded_variants > preprocess_variants

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        # # assembly
        # expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = [x for x in specimens if x != '900']),
        # expand("output/assembly/hifiasm/{specimen}/quast/T2T_scaffolded/report.html", specimen = [x for x in specimens if x != '900']),
        # expand("output/assembly/hifiasm/{specimen}/quast/hg38_scaffolded/report.html", specimen = [x for x in specimens if x != '900']),
        # expand("output/assembly/hifiasm/{specimen}/quast/raw/report.html", specimen = specimens),
        # expand("output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}.coverage.html", ref = ['hg38', 'T2T'], specimen = [x for x in specimens if x != '900']),
        # expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = specimens, hap = ['hap1', 'hap2']),
        # 'output/assembly/assembly_stats/all.tsv',
        # expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/repeatmasker/{specimen}.{hap}.repeatmasker.bed.gz", specimen = [x for x in specimens if x != '900'], hap = ['hap1', 'hap2']),
        # liftover
        # expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/GRCh38_repeatmasker.bed.gz", specimen = [x for x in specimens if x != '900'], hap = ['hap1', 'hap2']),
        # expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/GRCh38_genomicSuperDups.bed.gz", specimen = [x for x in specimens if x != '900'], hap = ['hap1', 'hap2']),
        # expand("output/assembly/hifiasm/{specimen}/hg38_scaffolded/{hap}/GRCh38_gencodeV44.bed.gz", specimen = [x for x in specimens if x != '900'], hap = ['hap1', 'hap2']),
        # expand("output/alignment/hg38_scaffolded/minimap2/standard/alu_positions/all_intersected_{file}.bed",  file = ['GRCh38_repeatmasker', 'GRCh38_genomicSuperDups', 'GRCh38_gencodeV44']),
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.bed",
        expand("output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.reverse_GRCh38.{file}.bed",  file = ['gencodeV48', 'genomicSuperDups', 'centromeres', 'microsat']),
        "output/alignment/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/alu_positions/all_alu_positions.de_novo_repeatmasker.bed"
