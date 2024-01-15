### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()
specimens_by_group = sample_table.groupby('group')['specimen'].unique().apply(list).to_dict()

# workaround right now for not being able to predefine chromosomes/contigs for parallelization of DeepVariant
chrs = ['chr' + str(n) for n in np.arange(1, 2).tolist()+['X', 'Y', 'M']]

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

def get_bams_per_sample(wildcards, sample_table = config['sample_table']):
    # A helper function for creating a list of expected bam outputs per sample post mapping + sorting, to be fed into a sample collation/merging rule.
    # Does not rely on checkpoint output: instead, relies on the metadata table to create expected paths. 
    bam_path = "output/mapping/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [bam_path.format(refalias=wildcards.refalias, mapper = wildcards.mapper, specimen=s[0], group = s[1], lane=s[2], smrtcell = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

### rules/sniffles.smk ###

def filter_vcf(vcf, outfile, write_fasta = True):
    # Reads in a vcf file as a pandas dataframe, then applies several filtering steps to remove likely inaccurate (FILTER fail or QUAL cutoff)
    # or imprecise (INFO field flags) SV calls.
    # If write_fasta is enabled, writes fasta output of alt allele sequences for downstream analyses (ex. repeatmasker).

    def write_fasta(vcf, outfile):
        # Takes the alt allele sequence output from VCF files and writes a fasta file.
        table = vcf[['ID', 'ALT']]
        print("Writing ALT sequences to ", outfile)
        with open(outfile, 'w') as f:
            for header, seq in table.to_records(index=False):
                f.write(f'>{header}\n{seq}\n')

    unfiltered = pd.read_table(vcf, skiprows = 72) # skip vcf header
    filtered = unfiltered[unfiltered['FILTER'] == 'PASS']
    filtered = filtered[~filtered['INFO'].str.contains('IMPRECISE')]
    filtered = filtered[~filtered['INFO'].str.contains('BND')]
    filtered = filtered[filtered['QUAL'].astype(int) >= 30]

    print("Writing filtered vcf to ", outfile)
    filtered.to_csv(outfile, index = False, sep ="\t")

    if write_fasta:
        # writes fasta output to the same path basename as outfile, but with .fa ext
        write_fasta(filtered, outfile.replace('unheadered.filtered.vcf', 'filtered.fa'))


# # tack the header back on...
        # shell("head -72 {input.vcf} | cat - {output.filtered} > {output.filtered}")