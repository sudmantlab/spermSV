#!/usr/bin/env python3
import pandas as pd
import numpy as np
import gzip
import io
import time

def parse_info_field(info_str):
    """Extract specific values from INFO field"""
    info_dict = dict(item.split('=') if '=' in item else (item, True) 
                    for item in info_str.split(';'))
    return info_dict

def parse_sample_field(sample_str):
    """Extract read counts from sample field"""
    try:
        parts = sample_str.split(':')
        return int(parts[2]), int(parts[3])  # DR (ref_reads), DV (alt_reads)
    except (IndexError, ValueError):
        return 0, 0

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
    
    # Extract INFO field values
    print(f"[{pd.Timestamp.now()}] Extracting INFO fields...")
    info_fields = vcf_df['INFO'].apply(parse_info_field)
    vcf_df['RNAMES'] = info_fields.apply(lambda x: x.get('RNAMES', ''))
    vcf_df['SVTYPE'] = info_fields.apply(lambda x: x.get('SVTYPE', ''))
    vcf_df['SVLEN'] = info_fields.apply(lambda x: int(x.get('SVLEN', 0)))
    
    # Extract sample field values
    print(f"[{pd.Timestamp.now()}] Processing sample fields...")
    sample_fields = vcf_df['SAMPLE'].apply(parse_sample_field)
    vcf_df['ref_support'] = sample_fields.apply(lambda x: x[0])
    vcf_df['alt_support'] = sample_fields.apply(lambda x: x[1])
    
    # Extract origin
    vcf_df['origin'] = vcf_df['ID'].str.split('_', expand=True)[0]
    
    return vcf_df

def process_repeatmasker(rm_file):
    """Process RepeatMasker file into DataFrame"""
    print(f"[{pd.Timestamp.now()}] Processing RepeatMasker file...")
    rm_df = pd.read_csv(rm_file, sep='\t')
    
    # Calculate additional columns
    rm_df['consensus_len'] = np.abs(rm_df['consensus_end'] - rm_df['consensus_start'])
    rm_df['match_len'] = np.abs(rm_df['query_end'] - rm_df['query_start'])
    
    return rm_df

def merge_data(vcf_df, rm_df):
    """Merge VCF and RepeatMasker data"""
    print(f"[{pd.Timestamp.now()}] Merging VCF and RepeatMasker data...")
    # Merge on ID fields
    merged_df = pd.merge(vcf_df, rm_df, 
                        left_on='ID', 
                        right_on='ID',
                        how='left')
    
    # Calculate percentage of match_len relative to SVLEN
    merged_df['match_len_svlen_perc'] = (merged_df['match_len'] / 
                                        merged_df['SVLEN'].abs() * 100).round(2)
    
    return merged_df

def main(vcf_file, rm_file):
    """Main function to process and merge files"""
    # Print input and output file paths with timestamps
    print(f"[{pd.Timestamp.now()}] Processing VCF file: {vcf_file}")
    print(f"[{pd.Timestamp.now()}] Processing RepeatMasker TSV file: {rm_file}")

    # Process input files
    vcf_df = process_vcf(vcf_file)
    rm_df = process_repeatmasker(rm_file)
    
    # Generate output file paths
    vcf_path = vcf_file
    base_path = vcf_path.rsplit('.', 3)[0]  # Remove '.filt.vcf.gz'
    
    all_annotations = f"{base_path}.all.repeatmasker_insertions.tsv"
    active_annotations = f"{base_path}.active.repeatmasker_insertions.tsv"
    single_type_annotations = f"{base_path}.active.single_type.repeatmasker_insertions.tsv"
    multi_type_annotations = f"{base_path}.active.multi_type.repeatmasker_insertions.tsv"
    multi_low_div_annotations = f"{base_path}.active.multi_type.low_div.repeatmasker_insertions.tsv"
    single_low_div_annotations = f"{base_path}.active.single_type.low_div.repeatmasker_insertions.tsv"
    
    # Merge & filter data
    all = merge_data(vcf_df, rm_df)
    active = all.query("family in ['SINE/Alu', 'LINE/L1', 'Retroposon/SVA']")
    single = active[active['ID'].isin(active.value_counts('ID').loc[lambda x: x == 1].index)]
    multiple = active[~active['ID'].isin(active.value_counts('ID').loc[lambda x: x == 1].index)]
    multi_low_div = multiple.query(f"divergence <= 3 & match_len_svlen_perc >= 25").reset_index(drop = True)
    single_low_div = single.query(f"divergence <= 3 & match_len_svlen_perc >= 25").reset_index(drop = True)
    
    # Save to file
    print(f"[{pd.Timestamp.now()}] Saving files...")
    # Match data frames with output paths and save as TSVs without index
    all.to_csv(all_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved all annotations to: {all_annotations}")

    active.to_csv(active_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved active annotations to: {active_annotations}")

    single.to_csv(single_type_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved single type annotations to: {single_type_annotations}")

    multiple.to_csv(multi_type_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved multi type annotations to: {multi_type_annotations}")

    multi_low_div.to_csv(multi_low_div_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved multi type, low divergence annotations to: {multi_low_div_annotations}")

    single_low_div.to_csv(single_low_div_annotations, sep='\t', index=False)
    print(f"[{pd.Timestamp.now()}] Saved single type, low divergence annotations to: {single_low_div_annotations}")
    

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: annotate_vcf_repeatmasker.py <vcf_file.gz> <repeatmasker_file>")
        print("Output files will be saved to the same directory as the input vcf.")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2])