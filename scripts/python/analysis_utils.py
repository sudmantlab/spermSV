import pandas as pd
import numpy as np
from parsing_utils import *

class calc:
    def edit_distance(contig, position, alt_rnames, bamfile):
        # Takes a variant of interest (id), and a DataFrame of variant information (variants).
        # At minimum, variants must contain #CHROM, POS, ID, and RNAMES.
        # Returns a dictionary containing the average & stdev of edit distance for variant 
        # support reads and non-support reads.

        tags = parse.samtags(bamfile, contig, position)
        
        # alt_reads = variants.query(f"ID == '{id}'")
        # alt_rnames, position = alt_reads['RNAMES'].values, alt_reads['POS'].unique()
        alt_tags, ref_tags = tags.loc[alt_rnames], tags.loc[~tags.index.isin(alt_rnames)]
        
        # number reads, mean edit distance, stdev edit distance

        return {'alt_n': alt_tags.shape[0],
                'alt_median': alt_tags['NM'].median(),
                'alt_mean': alt_tags['NM'].mean(),
                'alt_std': alt_tags['NM'].std(),
                'ref_n': ref_tags.shape[0],
                'ref_median': ref_tags['NM'].median(),
                'ref_mean': ref_tags['NM'].mean(),
                'ref_std': ref_tags['NM'].std(),
                'delta_median': alt_tags['NM'].median() - ref_tags['NM'].median(),
                'delta_mean': alt_tags['NM'].mean() - ref_tags['NM'].mean()}