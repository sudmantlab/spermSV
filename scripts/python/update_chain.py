#!/usr/bin/env python

import sys
from collections import defaultdict

def parse_agp(agp_file):
    """
    Parse the AGP file and create a mapping from contig names to reference chromosome names.
    """
    contig_to_chrom = {}
    
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip incomplete lines
            
            obj = fields[0]  # Reference chromosome (e.g., chr1_RagTag)
            comp_type = fields[4]  # Component type (W for successful scaffolding)
            comp_id = fields[5]  # Component ID (e.g., h1tg000176l)
            
            if comp_type == 'W':
                # Map the contig name to the reference chromosome
                contig_to_chrom[comp_id] = obj
    
    return contig_to_chrom

def update_chain_file(chain_file, contig_to_chrom, output_file):
    """
    Update the .chain file by replacing contig names with reference chromosome names.
    """
    with open(chain_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('chain'):
                fields = line.split('\t')
                
                # Replace contig names in column 8
                if fields[7] in contig_to_chrom:
                    fields[7] = contig_to_chrom[fields[7]]  # Replace contig name in column 8
                
                # Write the updated line to the output file
                out.write('\t'.join(fields) + '\n')
            else:
                # Write non-chain lines (e.g., alignment scores) as-is
                out.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python update_chain.py <agp> <chain> <updated_chain>")
        sys.exit(1)
    
    agp_file = sys.argv[1]
    chain_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Parse the AGP file to create the contig-to-chromosome mapping
    contig_to_chrom = parse_agp(agp_file)
    
    # Update the .chain file with the mapping
    update_chain_file(chain_file, contig_to_chrom, output_file)
    
    print(f"Updated chain file saved to {output_file}")

if __name__ == "__main__":
    main()