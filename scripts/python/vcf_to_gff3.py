#!/usr/bin/env python3
import pandas as pd
import numpy as np
import gzip
import io
import time

def read_vcf_gz(vcf_file):
    """Read compressed VCF file, strip headers, and return as DataFrame"""
    # Create string buffer to store non-header lines
    string_buffer = io.StringIO()
    
    # Read compressed file and filter header lines
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                string_buffer.write(line)
    
    # Reset buffer position
    string_buffer.seek(0)
    
    # Read into DataFrame
    df = pd.read_csv(string_buffer, sep='\t', header=None,
                    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                          'INFO', 'FORMAT', 'SAMPLE'])
    
    # Close buffer
    string_buffer.close()
    
    return df

def process_vcf(vcf_file):
    """Process VCF file into DataFrame"""
    # Read VCF file
    print(f"[{pd.Timestamp.now()}] Reading and parsing VCF file...")
    if vcf_file.endswith('.vcf.gz'):
        vcf_df = read_vcf_gz(vcf_file)
        
    elif vcf_file.endswith('.vcf'):
        vcf_df = pd.read_csv(vcf_file, sep='\t', header=None,
                    names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 
                          'INFO', 'FORMAT', 'SAMPLE'])
    else:
        raise ValueError("VCF file must end with either '.vcf.gz' or '.vcf'")

    return vcf_df

def vcf_to_gff3(df):
    """Convert VCF DataFrame to GFF3 format"""
    df = df.copy()  # Create a copy to avoid modifying the original DataFrame
    df.rename(columns={'CHROM': 'seqname', 'POS': 'start'}, inplace=True)
    df['SVLEN'] = df['INFO'].str.extract("(?<=SVLEN=)(.+)(?=;END)").astype('int')
    df['source'] = 'sniffles_mosaic'
    df['feature'] = 'variation'
    df['start'] = df['start'].astype('int')  # ensure int type for df compatibility
    df['end'] = np.maximum(df['start'], df['start'] + df['SVLEN'])
    df['score'] = '.'
    df['strand'] = df['INFO'].str.extract("(?<=STRAND=)(.)(?=;)").fillna('.')
    df['frame'] = '.'
    df['attribute'] = df['ID']
    
    # fix BND (translocation) values with NaN or lesser value for end by assigning end = start 
    bnd_index = df.index[df['end'].isna()]
    df.loc[bnd_index, 'end'] = df.loc[bnd_index, 'start']

    # subset to GFF3 standard columns
    df = df[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]
    return df


def main(vcf_file, output_file):
    """Main function to process vcf"""
    # Print input and output file paths with timestamps
    print(f"[{pd.Timestamp.now()}] Processing VCF file: {vcf_file}")

    # Process input files
    vcf = process_vcf(vcf_file)
    gff3 = vcf_to_gff3(vcf)
    
    # Save to file
    print(f"[{pd.Timestamp.now()}] Saving output...")
    gff3.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"[{pd.Timestamp.now()}] Output saved to: {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: vcf_to_gff3.py <file.vcf.gz> <file.gff>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])