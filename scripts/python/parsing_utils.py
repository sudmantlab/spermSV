import pandas as pd
import numpy as np
import pysam

class parse:
    def sniffles_vcf(vcf_path, header_len = 241):
        # Takes in a sniffles (mosaic) vcf (.vcf.gz) and parses it a DataFrame.
        # The header length will vary based on the reference genome and other operations applied to the vcf.
        # 241 is the default value for the hg38 (no alts) reference genome + default sniffles mosaic output.
        vcf = pd.read_table(vcf_path, skiprows = header_len)

        # Extracts specific fields to columns for downstream use
        vcf['RNAMES'] = vcf['INFO'].str.extract("(?<=RNAMES=)(.+)(?=;C)")
        vcf['SVTYPE'] = vcf['INFO'].str.extract("(?<=SVTYPE=)(\\w+)(?=;)")
        vcf['SVLEN'] = vcf['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)").astype('int')
        vcf[['ref_reads', 'alt_reads']] = vcf['SAMPLE'].str.split(':', expand = True)[[2,3]].astype(int)
        vcf['origin'] = vcf['ID'].str.split('_', expand = True)[0]
        
        return vcf
    
    def repeatmasker(rm_path):
        # Takes in a RepeatMasker fixed-width outfile and parses it into a DataFrame.
        # RepeatMasker outfiles are not quite true fixed-width.
        # The original outfile column names are relabeled for clarity.

        out = pd.read_fwf(rm_path, header = 1)
        out = out[['score', 'div.', 'del.', 'ins.', 'sequence', 'begin end', '(left)', 'repeat', 'class/family', 'begin  end', '(left).1']]
        out[['query_start', 'query_end']] = out['begin end'].str.split(expand = True).astype(int)

        # the consensus sequence start/end/left coordinates are hard to parse cleanly
        # TODO: consider a cleaner way to parse this, in the meantime leave as str
        # if the start coordinates are encapsulated by (), it indicates a start position on the reverse complement of the consensus
        # if the left coordinates are encapsulated by (), it indicates the remaining length of the consensus past the match
        out[['consensus_start', 'consensus_end']] = out['begin  end'].str.split(expand = True)
        out.rename(columns = {'score': 'SW_score', 'div.': 'divergence', 'del.': 'deleted', 'ins.': 'inserted', 
                            'sequence': 'ID', '(left)': 'query_left', 'repeat': 'motif', 'class/family': 'family', 
                            '(left).1': 'consensus_left'}, inplace = True)
        
        # gives total length of insertion sequence that matched consensus sequence
        out['consensus_len'] = out['query_end'] - out['query_start']

        return out[['ID', 'motif', 'family', 'SW_score', 'divergence', 'deleted', 'inserted', 'query_start', 'query_end', 'query_left',
                    'consensus_start', 'consensus_end', 'consensus_len', 'consensus_left']]

    def samtags(bamfile, contig, position, interval = 500):
        # Given an AlignmentFile from pysam, fetches reads within a given region and returns a 
        # DataFrame of read tags for all reads within the given +/- interval (default 500 bp).
        iter = bamfile.fetch(contig, position-interval, position+interval)
        tags = dict()
        for x in iter:
            tags[x.query_name] = dict(x.tags)
        return pd.DataFrame.from_dict(tags, orient = 'index')


class annotate:
    def repeatmasker(vcf_path, rm_path, header_len = 241, save = False, filepath = None):
        # Reads a VCF and RepeatMasker outfile, merging the contents of the outfile with the VCF.

        vcf = parse.sniffles_vcf(vcf_path, header_len)
        rm = parse.repeatmasker(rm_path)

        merged = vcf.merge(rm, on = 'ID')
        # position query start/end coords relative to variant POS
        merged['query_start'] = merged['POS'] + merged['query_start']
        merged['query_end'] = merged['POS'] + merged['query_end']
        merged['match_len'] = np.abs(merged['query_end'] - merged['query_start'])
        merged['match_len_svlen_perc'] = (merged['match_len']/merged['SVLEN'])*100

        if save:
            write.table(merged, filepath)
        else:
            return merged
    
class transform:
    def gff3(df, save = False, filepath = None):
        # Creates a GFF3 formatted version containing the sniffles variant ID in the attribute field.
        # This can be used as a key to merge with other annotated tables.
        # Export using write.table() with a .gff suffix.

        df.rename(columns = {'#CHROM': 'seqname', 'POS': 'start'}, inplace = True)
        df['source'] = 'sniffles_mosaic'
        df['feature'] = 'variation'
        df['start'] = df['start'].astype('int') # ensure int type for df compatibility
        df['end'] = np.maximum(df['start'], df['start'] + df['SVLEN'])
        df['score'] = '.'
        df['strand'] =  df['INFO'].str.extract("(?<=STRAND=)(.)(?=;)").fillna('.')
        df['frame'] = '.'
        df['attribute'] = df['ID']
        
        # fix BND (translocation) values with NaN or lesser value for end by assigning end = start 
        bnd_index = df.index[df['end'].isna()]
        df.loc[bnd_index, 'end'] = df.loc[bnd_index, 'start']

        # subset to GFF3 standard columns
        df = df[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]

        if save:
            df.to_csv(filepath, sep = '\t', index = False, header = False) # header breaks bedtools parsing
        else:
            return df
    
class filter:
    def motif(df, motif, inplace = False):
        # Wrapper to filter the annotated table to only variants of a certain motif.
        return df[df['motif'].str.contains(motif)].reset_index(drop = True)

    def origin(df, origin, inplace = False):
        # Wrapper to filter the annotated table to only variants from a specific specimen.
        return df[df['origin'] == origin].reset_index(drop = True)
        
    def reads(df, ref_max = 0, alt_max = 0, inplace = False):
        # Wrapper to filter the annotated table to only variants with less than hard-set ref/alt support.
        # Defaults to 0 for both.
        return df.query(f"ref_reads <= {ref_max} & alt_reads <= {alt_max}")

class write:
    def rnames(df, filepath):
        # Writes a .txt file containing RNAMES for the given table.
        # This file can be used in a call to samtools view with the -N flag to subset a BAM file containing only these reads.
        # TODO: Assert .txt extension before publish time

        candidates = df[['ID', 'origin', 'RNAMES']]
        candidates.loc[:, 'RNAMES'] = candidates['RNAMES'].str.split(',')
        candidates = candidates.explode('RNAMES')

        # flattens a list of lists
        rnames = [read for sublist in [x.split(',') for x in df['RNAMES'].to_list()] for read in sublist]

        print(f'Saving .txt file to {filepath} ...')
        rnames.to_csv(filepath, index = False, header = None)

    def alt_to_fasta(df, filepath):
        # Takes the ALT sequences from a VCF file and writes to a fasta file.
        # TODO: Assert .fa extension before publish time
        table = df[['ID', 'ALT']]

        print(f'Writing ALT sequences to {filepath} ...')
        with open(filepath, 'w') as f:
            for header, seq in table.to_records(index=False):
                f.write(f'>{header}\n{seq}\n')
    
    def table(df, filepath, index = False):
        # Wrapper for .to_csv with presets in mind for genomic table outputs.
        if filepath == None:
                raise SyntaxError("Please provide a file path.")
        else:
            print(f'Saving df to {filepath} ...')
            if filepath.split('.')[-1] == 'tsv':
                df.to_csv(filepath, index = index, sep = '\t')
            elif filepath.split('.')[-1] == 'csv':
                df.to_csv(filepath, index = index, sep = ',')
            else:
                # save tab delimited as default if other extension
                # TODO: make less messy later
                df.to_csv(filepath, index = index, sep = '\t')