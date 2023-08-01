Rules
=====


This folder contains the Snakemake rules used during assembly and QC. 

Generally, each Snakefile contains rules pertaining to a certain aspect of the assembly. All rules were derived and/or modified from `sudmantlab/panpan` unless otherwise stated.

- `Hifi.smk` contains rules associated with generating and assembling HiFi reads;
- `genomeQC.smk` contains rules associated with generating QC metrics at each step.
- `minimap2.smk` contains rules associated with running minimap2, modified from `sudmantlab/MyotisGenomeAssembly`.
- `make_symlinks.smk` contains rules associated with creating symlinks in `data` folder from the `globus_data` transfer, maintaining the canonical file structure associated with rules from `sudmantlab/panpan` and `sudmantlab/MyotisGenomeAssembly`.