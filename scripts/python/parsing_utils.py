import pandas as pd

class SnifflesVCF(pd.DataFrame):
    def __init__(self):
        # inherit all methods and props from DataFrame
        super().__init__()

    def create(self, vcf_path, header_len = 241):
        # Takes in a sniffles (mosaic) vcf (.vcf.gz) and parses it into a pandas df.
        # The header length will vary based on the reference genome and other operations applied to the vcf.
        # 241 is the default value for the hg38 (no alts) reference genome + default sniffles mosaic output.

        vcf = pd.read_table(vcf_path, skiprows = header_len)
        vcf['RNAMES'] = vcf['INFO'].str.extract("(?<=RNAMES=)(.+)(?=;C)")

        # Adds additional columns for filtering downstream
        vcf['SVTYPE'] = vcf['INFO'].str.extract("(?<=SVTYPE=)(\\w+)(?=;)")
        vcf['SVLEN'] = vcf['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)")

        self.__dict__.update(vcf.__dict__)
        return self

class AnnotatedGFF(pd.DataFrame):
    def __init__(self):
        # inherit all methods and props from DataFrame
        super().__init__()

        # specific attributes
        self.origin_vcf = None
        self.origin_gff = None
        self.motif_filter = None
        self.read_filter = None 
    
    def create(self, vcf_path, gff_path, header_len = 241):
        # Creates a variant table based on the VCF standard, annotated with information from RepeatMasker (GFF).
        # Returns a DataFrame-like table.

        # load VCF, extract info values to independent columns for downstream parsing
        vcf = pd.read_table(vcf_path, skiprows = header_len) # 241 default value is for hg38 no alts x sniffles, will vary on the reference
        vcf['RNAMES'] = vcf['INFO'].str.extract("(?<=RNAMES=)(.+)(?=;C)")
        vcf['SVLEN'] = vcf['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)")
        
        # read in base repeatmasker GFF, extract attribute values to independent columns for 1downstream parsing
        merged = pd.read_table(gff_path, skiprows = 3, header = None)

        # from gff format spec
        # modified colnames to make values clear relative to repeatmasker
        # attribute contains motif name + start/end of the match in the consensus seq
        merged.columns = ['ID', 'source', 'feature', 'motif_start', 'motif_end', 'divergence', 'strand', 'frame', 'attribute']
        merged['motif'] = merged['attribute'].str.extract('(?<=\\:)(.*)(?=\")')
        merged['match_len'] = merged['motif_end'] - merged['motif_start'] # gives total length matching motif consensus, distinct from VCF's variant length
        merged = merged[['ID', 'motif', 'motif_start', 'motif_end', 'divergence', 'strand']]
        merged = merged.merge(vcf, on = 'ID')
        merged['motif_start'] = merged['POS'] + merged['motif_start']
        merged['motif_end'] = merged['POS'] + merged['motif_end']
        merged[['ref_reads', 'alt_reads']] = merged['SAMPLE'].str.split(':', expand = True)[[2,3]].astype(int)

        self.__dict__.update(merged.__dict__)
        self.origin_vcf = vcf_path
        self.origin_gff = gff_path
        return self
    
    def as_aGFF(self, df):
        # Sets a copy of a DataFrame with the given class.
        self.__dict__update(df.__dict__)
        return self
    
### utility functions that can operate on aGFF or DataFrame

def filter_motif(df, motif, inplace = False):
    # Filters the annotated table to only variants of a certain motif, using a simple substring regex match.
    
    if inplace:
        df.motif_filter = motif
        df = df[df['motif'].str.contains(motif)].reset_index(drop = True)
    else:
        return df[df['motif'].str.contains(motif)].reset_index(drop = True)     

def filter_reads(df, ref_max, alt_max, inplace = False):
    # Filters the annotated table by restricting to <= ref, <= alt read support.
    if inplace:
        df.read_filter = {"ref_max": ref_max, "alt_max": alt_max}
        df = df.query(f"ref_reads <= {ref_max} & alt_reads <= {alt_max}")
    else:
        return df.query(f"ref_reads <= {ref_max} & alt_reads <= {alt_max}")

def export_rnames(df, filepath):
    # Writes a .txt file containing QNAMES for each support read in the table.
    # This file can be used in a call to samtools view with the -N flag to subset a BAM file containing only these reads.

    flat_rnames = pd.Series([read for sublist in [x.split(',') for x in df['RNAMES'].to_list()] for read in sublist]).drop_duplicates()

    print(f'Saving .txt file to {filepath} ...')
    flat_rnames.to_csv(filepath, index = False, header = None)

def write_fasta(vcf, outfile):
    # Takes the ALT sequences and writes to a fasta file.
    table = vcf[['ID', 'ALT']]
    with open(outfile, 'w') as f:
        for header, seq in table.to_records(index=False):
            f.write(f'>{header}\n{seq}\n')

def as_gff3(df, save = False, filepath = None):
    # Creates a GFF3 formatted version containing the sniffles variant ID in the attribute field.
    # This can be used as a key to merge with other annotated tables.

    gff3 = df.copy()
    gff3.rename(columns = {'#CHROM': 'seqname', 'POS': 'start'}, inplace = True)
    gff3['source'] = 'sniffles_mosaic'
    gff3['feature'] = 'variation'
    gff3['start'] = gff3['start'].astype('int') # ensure int type for gff3 compatibility
    gff3['end'] = gff3['INFO'].str.extract("(?<=END=)(\\w+)(?=;)")
    gff3['score'] = '.'
    gff3['strand'] =  gff3['INFO'].str.extract("(?<=STRAND=)(.)(?=;)").fillna('.')
    gff3['frame'] = '.'
    gff3['attribute'] = gff3['ID']
    
    # fix BND (translocation) values with NaN for end by assigning end = start 
    bnd_index = gff3.index[gff3['end'].isna()]
    gff3.loc[bnd_index, 'end'] = gff3.loc[bnd_index, 'start']

    gff3 = gff3[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

    if save:
        if filepath == None:
            raise SyntaxError("Please provide a file path to save the GFF3 file.")
        else:
            print(f'Saving GFF3 to {filepath} ...')
            gff3.to_csv(filepath, index = False, header = None, sep = '\t')
    else:
        return gff3