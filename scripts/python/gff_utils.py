import pandas as pd
import os

class AnnotatedGFF:
    def __init__(self):
        # define accessible attributes
        self.origin_vcf = None
        self.origin_gff = None
        self.motif = None
        self.table = None
        self.read_filter = None 
    
    def create(self, vcf_path, gff_path, motif, header_len = 241):
        # Instantiates an annotated RepeatMasker GFF with additional variant info from the origin VCF file.
        # Additionally, filters the GFF to a desired subset of motifs, using a simple substring regex match.
        # Returns a dataframe of the annotated GFF output. 

        # load VCF, extract info values to independent columns for downstream parsing
        vcf = pd.read_table(vcf_path, skiprows = header_len) # 241 default value is for hg38 no alts x sniffles, will vary on the reference
        vcf['SVLEN'] = vcf['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)").astype('float')
        vcf['RNAMES'] = vcf['INFO'].str.extract("(?<=RNAMES=)(.+)(?=;C)")
        vcf[['ref_reads', 'alt_reads']] = vcf['SAMPLE'].str.split(':', expand = True)[[2,3]].astype(int)
        
        # read in base GFF, extract attribute values to independent columns for downstream parsing
        aGFF = pd.read_table(gff_path, skiprows = 3, header = None)
        aGFF.columns = ['ID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] # from gff format spec
        aGFF['motif'] = aGFF['attribute'].str.extract('(?<=\:)(.*)(?=\")')
        aGFF['match_len'] = aGFF['end'] - aGFF['start'] # gives total length matching motif consensus, distinct from VCF's variant length
        aGFF = aGFF.merge(vcf[['SVLEN', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'ref_reads', 'alt_reads', 'RNAMES']], on = 'ID')
        aGFF['start'] = aGFF['POS'] + aGFF['start']
        aGFF['end'] = aGFF['POS'] + aGFF['end']
        aGFF = aGFF[aGFF['motif'].str.contains(motif)].reset_index(drop = True) # filtering step for motif

        # define accessible attributes
        self.origin_vcf = vcf_path
        self.origin_gff = gff_path
        self.motif = motif
        self.table = aGFF
        self.read_filter = None
        return self
    
    def filter_reads(self, ref_max, alt_max, inplace = False):
        # Filters the annotated table by restricting to <= ref, <= alt read support.
        
        self.read_filter = {'ref_max': ref_max, 'alt_max': alt_max}
        if inplace:
            self.table = self.table.query(f"ref_reads <= {ref_max} & alt_reads <= {alt_max}")

        else:
            return self.table.query(f"ref_reads <= {ref_max} & alt_reads <= {alt_max}")

    def export_rnames(self, filename):
        # Writes a .txt file containing QNAMES for each support read in the table.
        # This file can be used in a call to samtools view with the -N flag to subset a BAM file containing only these reads.

        flat_rnames = pd.Series([read for sublist in [x.split(',') for x in self.table['RNAMES'].to_list()] for read in sublist]).drop_duplicates()

        print(f'Saving .txt file to {filename} ...')
        flat_rnames.to_csv(filename, index = False, header = None)

    def to_gff3(self, filename, return_df = True):
        # Writes a GFF3 formatted version containing the RepeatMasker motif and sniffles variant ID in the attribute field.
        # By default, returns a copy of the GFF3 for inspection.
        export = self.table.copy()
        export['seqname'] = export['#CHROM']
        export['attribute'] = ['motif "{a}"; sniffles_ID "{b}"; support "{c}"'.format(a = motif, b = ID, c = str((ref, alt))) for (index, motif, ID, ref, alt) in export[['motif', 'ID', 'ref_reads', 'alt_reads']].itertuples()]
        export = export[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

        print(f'Saving export file to {filename} ...')
        export.to_csv(filename, index = False, header = None, sep = '\t')

        if return_df:
            return export