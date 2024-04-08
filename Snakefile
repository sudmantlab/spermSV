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

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        expand("output/assembly/hifiasm/{specimen}/{specimen}.asm.bp.{hap}.p_ctg.gfa", hap = ['hap1', 'hap2'], specimen = specimens),
        expand("output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa", specimen = specimens, hap = ['hap1', 'hap2'])
        # tandem repeat genotyping
        # expand('output/alignment/hg38/winnowmap/standard/variants/trgt/repeat_catalog/{specimen}.sorted.vcf.gz',  specimen = specimens),
        # expand('output/alignment/hg38/winnowmap/standard/variants/trgt/pathogenic/{specimen}.sorted.vcf.gz',  specimen = specimens),
        # alu analysis pipeline rebuild
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/all.filt.gff", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap']),
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/all.overlap_repetitive.gff", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap']),
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/all.annotate_{source}.gff", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap'], source = ['repeatMasker', 'simpleRepeat', 'genomicSuperDups', 'centromeres', 'microsat']),
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/all.annotate_{database}.gff", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap'], database = ['gencode', 'refseq']),
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/repeatmasker/all.filt.alt.fa", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap']),
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/all.annotate_repeatmasker.tsv", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap']),
        # edit distance calcs
        # expand("analysis/hg38/denovo/files/{mapper}/{setting}/variants/{specimen}.annotate_edit_distance.tsv", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap'], specimen = specimens),
        # pepper-margin-deepvariant calls; chrs is a list in common.smk
        # expand("output/alignment/hg38/{mapper}/{setting}/variants/pmdv/{specimen}/{region}/{specimen}.{region}.vcf.gz", mapper = ['minimap2', 'winnowmap'], setting = ['standard', 'duplomap'], specimen = specimens, region = chrs)