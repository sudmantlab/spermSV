# This "helper" rule creates symlinks of the Globus direct data transfer files in globus_data, to the "data" folder 
# specified in downstream Snakemake rules. This preserves the structure for analysis AND allows for Globus
# file transfer sync without requiring manual transfers to/from either directory.

# Helper functions exist in rules/common.smk to symlink either bam or fastq: default should be bam.
# ***Don't do both***, the HiFiAdapterFilt script doesn't like how that plays with the existing file/directory
# structure, and it'll just assume you want *everything* ending in .bam, .fastq, .fastq.gz, etc. processed.
# The filtering step will thus take 2-4X as long.

rule make_bam_symlink:
    input:
        "/global/scratch/users/stacy-l/spermSV/globus_data/{lane}/{smrtcell}.bam"
    output:
        "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.bam"
    shell:
        """
        ln -s {input} {output} 
        """

rule make_fastq_symlink:
    input:
        "/global/scratch/users/stacy-l/spermSV/globus_data/{lane}/{smrtcell}.fastq.gz"
    output:
        "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    shell:
        """
        ln -s {input} {output} 
        """