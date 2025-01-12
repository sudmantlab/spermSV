import cyvcf2
from collections import Counter
import json
import sys
import os
import argparse
import numpy as np

def summarize_vcf(filename):
    vcf = cyvcf2.VCF(filename)
    # Initialize counters
    filter_counts = Counter()
    svtype_counts = Counter()
    svlengths_by_type = {}
    rnames_lengths_by_type = {}
    # Single pass through the file
    for variant in vcf:
        # FILTER values
        if variant.FILTER == None:
            filter_counts['PASS'] += 1
        else:
            filter_counts[variant.FILTER] += 1
        # SVTYPE values
        svtype = variant.INFO.get('SVTYPE')
        svtype_counts[svtype] += 1
        # SVLEN values
        svlen = variant.INFO.get('SVLEN')
        if svlen is None:
            svlen = 0
        # RNAMES by SVTYPE
        rnames = variant.INFO.get('RNAMES')
        if isinstance(rnames, str):
            rnames = [rnames] if rnames else []
        elif rnames is None:
            rnames = []
        if svtype not in rnames_lengths_by_type:
            rnames_lengths_by_type[svtype] = []
        rnames_lengths_by_type[svtype].append(len(rnames))
        
        # SVLEN by SVTYPE
        if svtype not in svlengths_by_type:
            svlengths_by_type[svtype] = []
        svlengths_by_type[svtype].append(svlen)
    
    # Calculate stats for each SV type
    svlen_stats = {}
    for svtype, lengths in svlengths_by_type.items():
        lengths = np.array(lengths)
        svlen_stats[svtype] = {
            "min": float(round(np.min(lengths), 2)),
            "max": float(round(np.max(lengths), 2)),
            "mean": float(round(np.mean(lengths), 2)),
            "median": float(round(np.median(lengths), 2)),
            "stddev": float(round(np.std(lengths), 2)),
        }
    rnames_lengths_stats = {}
    for svtype, lengths in rnames_lengths_by_type.items():
        lengths = np.array(lengths)
        rnames_lengths_stats[svtype] = {
            "min": float(round(np.min(lengths), 2)),
            "max": float(round(np.max(lengths), 2)),
            "mean": float(round(np.mean(lengths), 2)),
            "median": float(round(np.median(lengths), 2)),
            "stddev": float(round(np.std(lengths), 2)),
        }
    
    return {
        "filename": filename,
        "total_variants": sum(filter_counts.values()),
        "filter_counts": dict(filter_counts),
        "svtype_counts": dict(svtype_counts),
        "svlen_stats_by_type": svlen_stats,
        "rnames_lengths_stats_by_type": rnames_lengths_stats
    }

def main():
    parser = argparse.ArgumentParser(description='Analyze VCF file and output summary statistics')
    parser.add_argument('vcf_path', help='Path to the input VCF file')
    args = parser.parse_args()

    # Generate output path
    output_path = os.path.splitext(args.vcf_path)[0]
    if output_path.endswith('.vcf'):  # Handle double extension for .vcf.gz
        output_path = os.path.splitext(output_path)[0]
    output_path += '.summary.json'

    # Analyze VCF and write results
    try:
        results = summarize_vcf(args.vcf_path)
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Summary written to: {output_path}")
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()