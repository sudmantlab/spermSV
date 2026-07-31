### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()
specimens_by_group = sample_table.groupby('group')['specimen'].unique().apply(list).to_dict()

chrs = ['chr' + str(n) for n in np.arange(1, 23).tolist()+['X', 'Y']]

### rules/make_symlinks.smk ###

def make_bam_symlinks(wildcards):
    symlink_path_format = "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.bam"
    samples = pd.read_table("samples.tsv", index_col=False)
    samples = samples.to_records(index=False)
    symlink_paths = [symlink_path_format.format(specimen=s[0], lane=s[2], smrtcell = s[3]) for s in samples]
    return symlink_paths

def make_fastq_symlinks(wildcards):
    symlink_path_format = "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.fastq.gz"
    samples = pd.read_table("samples.tsv", index_col=False)
    samples = samples.to_records(index=False)
    symlink_paths = [symlink_path_format.format(specimen=s[0], lane=s[2], smrtcell = s[3]) for s in samples]
    return symlink_paths

### rules/samtools_utils.smk ###

def get_temp_bams_per_sample(wildcards, sample_table = config['sample_table']) -> List[str]:
    """
    A helper function for creating a list of expected per-smrtcell bam outputs per sample to be merged into a single bam.
    Adjusted to optionally include {hap} in the basepath if wildcards.hap exists.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    
    # Check if '{hap}' is in the basepath and format accordingly
    if hasattr(wildcards, 'hap'):
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}.{hap}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], hap=wildcards.hap, lane=s[2], smrtcell=s[3]) for s in samples]
    else:
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], lane=s[2], smrtcell=s[3]) for s in samples]
    
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

def get_fastqs_per_sample(wildcards, sample_table = config['sample_table']) -> List:
    """
    A helper function for creating a list of fastq.gz paths for each sample ahead of hifiasm assembly.
    TODO: The base path should maybe be set in the configfile.
    """
    basepath = "output/preprocessing/uBAMtoFastq/{specimen}/{lane}/{smrtcell}.ccs.fastq.gz"
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [basepath.format(specimen=s[0], group = s[1], lane=s[2], smrtcell = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples