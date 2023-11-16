import pandas as pd

def annotate_gff(vcf_path, gff_path, motif, header_len = 241):
    # Creates a GFF that annotates the RepeatMasker GFF output with additional variant info from the origin VCF file.
    # Additionally, filters the GFF to a desired subset of motifs, using a simple substring regex match.
    # Returns a dataframe of the annotated GFF output. 

    # modifications to VCF for parsing later
    vcf = pd.read_table(vcf_path, skiprows = header_len) # 241 default value is for hg38 no alts x sniffles, will vary on the reference
    vcf['SVLEN'] = vcf['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)").astype('float')
    vcf['RNAMES'] = vcf['INFO'].str.extract("(?<=RNAMES=)(.+)(?=;C)")
    vcf[['ref_reads', 'alt_reads']] = vcf['SAMPLE'].str.split(':', expand = True)[[2,3]].astype(int)
    
    # read in base GFF
    ann_gff = pd.read_table(gff_path, skiprows = 3, header = None)
    ann_gff.columns = ['ID', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] # from gff format
    ann_gff['motif'] = ann_gff['attribute'].str.extract('(?<=\:)(.*)(?=\")')
    ann_gff['match_len'] = ann_gff['end'] - ann_gff['start'] # gives total length matching motif consensus
    ann_gff = ann_gff.merge(vcf[['SVLEN', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'ref_reads', 'alt_reads', 'RNAMES']], on = 'ID')
    ann_gff['start'] = ann_gff['POS'] + ann_gff['start']
    ann_gff['end'] = ann_gff['POS'] + ann_gff['end']
    ann_gff = ann_gff[ann_gff['motif'].str.contains(motif)].reset_index(drop = True) # filtering step for motif

    return ann_gff.sort_values('match_len', ascending = False)

def rnames_to_txt(filename, ann_gff):
    # Takes in an annotated GFF from annotate_gff and writes a .txt file containing QNAMES for each support read in the table.
    # This file can be used in a call to samtools view with the -N flag to subset a BAM file containing only these reads.

    flat_rnames = pd.Series([read for sublist in [x.split(',') for x in ann_gff['RNAMES'].to_list()] for read in sublist]).drop_duplicates()

    print(f'Saving .txt file to {filename} ...')
    flat_rnames.to_csv(filename, index = False, header = None)

def save_gff3(filename, gff_df, return_df = True):
    # Takes an input annotated GFF and returns a GFF3 formatted version containing the motif and variant ID in the attribute field.
    gff3 = gff_df.copy()
    gff3['seqname'] = gff3['#CHROM']
    gff3['attribute'] = ['motif "{a}"; sniffles_ID "{b}"; support "{c}"'.format(a = motif, b = ID, c = str((ref, alt))) for (index, motif, ID, ref, alt) in gff3[['motif', 'ID', 'ref_reads', 'alt_reads']].itertuples()]
    gff3 = gff3[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

    print(f'Saving GFF3 file to {filename} ...')
    gff3.to_csv(filename, index = False, header = None, sep = '\t')

    if return_df:
        return gff3