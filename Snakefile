# Create data symlinks (if necessary)
include: "rules/make_symlinks.smk"

# CCS Generation
# include: "rules/CCS_split.smk"

# HiFi QC
include: "rules/HiFiAdapterFilt.smk"

# Genome Assembly: hifiasm
include: "rules/HiFi.smk"

# Mapping: minimap2
include: "rules/minimap2.smk"

# Read-mapped SV calling
include: "rules/sniffles.smk"

# DeepVariant – small variant calling
include: "rules/DeepVariant.smk"

# Tandem repeat expansion detection & genotyping
include: "rules/straglr.smk"
include: "rules/trgt.smk"

# Genome QC
# include: "rules/genomeQC.smk"

# Trinity
#include: "rules/QC_RNASeq.smk"
#include: "rules/trimmomatic.smk"
#include: "rules/trinity.smk"

# Annotate
#include: "rules/mitoHiFi.smk"

# Omni-C mapping to genome
#include: "rules/omniC.smk"