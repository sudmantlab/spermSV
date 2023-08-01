# This rule creates symlinks of the Globus direct data transfer files in globus_data, to the "data" folder 
# specified in downstream Snakemake rules. This preserves the structure for analysis AND allows for Globus
# file transfer sync without requiring manual transfers to/from either directory.

# Options are provided below to symlink either bam or fastq: default should be bam.
# ***Don't do both***, the HiFiAdapterFilt script doesn't like how that plays with the existing file/directory
# structure, and it'll just assume you want *everything* ending in .bam, .fastq, .fastq.gz, etc. processed.
# The filtering step will thus take 2-4X as long.

import os
import pandas as pd

def make_bam_symlinks(wildcards):
    symlink_path_format = "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.ccs.bam"
    samples = pd.read_table("samples.tsv", index_col=False)
    samples = samples.to_records(index=False)
    symlink_paths = [symlink_path_format.format(specimen=s[0], lane=s[1], smrtcell = s[2]) for s in samples]
    return symlink_paths

# def make_fastq_symlinks(wildcards):
#     symlink_path_format = "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
#     samples = pd.read_table("samples.tsv", index_col=False)
#     samples = samples.to_records(index=False)
#     symlink_paths = [symlink_path_format.format(specimen=s[0], lane=s[1], smrtcell = s[2]) for s in samples]
#     return symlink_paths

rule make_bam_symlink:
    input:
        "/global/scratch/users/stacy-l/spermSV/globus_data/{lane}/{smrtcell}.bam"
    output:
        "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.ccs.bam"
    shell:
        """
        ln -s {input} {output} 
        """

# rule make_fastq_symlink:
#     input:
#         "/global/scratch/users/stacy-l/spermSV/globus_data/{lane}/{smrtcell}.fastq.gz"
#     output:
#         "data/PacBio-HiFi/homo_sapiens/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
#     shell:
#         """
#         ln -s {input} {output} 
#         """