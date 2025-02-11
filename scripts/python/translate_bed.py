#!/usr/bin/env python

import sys
from collections import defaultdict
from intervaltree import IntervalTree

def parse_agp(agp_file):
    """
    Parse the AGP file and build an interval tree for each chromosome.
    """
    trans = defaultdict(IntervalTree)
    
    with open(agp_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip header lines
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue  # Skip incomplete lines
            if fields[4] != 'W':
                continue  # Skip gap lines
            
            obj = fields[0]  # Reference chromosome (e.g., chr1_RagTag)
            obj_beg = int(fields[1])  # Start position in reference
            obj_end = int(fields[2])  # End position in reference
            part_num = int(fields[3])  # Part number (not used here)
            comp_type = fields[4]  # Component type (W for successful scaffolding)
            comp_id = fields[5]  # Component ID (e.g., h1tg000176l)
            comp_beg = int(fields[6])  # Start position in component
            comp_end = int(fields[7])  # End position in component
            orientation = fields[8]  # Orientation (+, -, etc.)

            trans[obj][obj_beg:obj_end] = (comp_id, comp_beg, comp_end, orientation)
    
    return trans

def translate_bed_coordinates(bed_file, trans, hap):
    """
    Translate BED file coordinates using the interval tree.
    """
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('track') or line.startswith('browser'):
                print(line.strip())  # Print track/browser lines as-is
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                print(line.strip())  # Print malformed lines as-is
                continue
            
            chrom = fields[0]  # Chromosome in reference genome
            start = int(fields[1])  # Start position in reference
            end = int(fields[2])  # End position in reference
            
            # Append "_RagTag" to match the AGP file format
            chrom_ragtag = f"{chrom}_RagTag"
            
            if chrom_ragtag not in trans:
                # Chromosome not in AGP, skip
                continue
            
            # Find overlapping intervals in the interval tree
            overlaps = trans[chrom_ragtag][start:end]
            if not overlaps:
                print(line.strip())  # No overlap found, print as-is
                continue
            
            # Assuming there's only one overlapping interval (as per AGP spec)
            overlap = list(overlaps)[0]
            comp_id, comp_beg, comp_end, orientation = overlap.data
            
            # Calculate new coordinates in the de novo assembly
            new_start = comp_beg + (start - overlap.begin)
            new_end = comp_beg + (end - overlap.begin)
            
            # Update the BED line with the new coordinates
            fields[0] = f"{chrom}_RagTag_{hap}" # match to the hap-tagged scaffolded fasta
            fields[1] = str(new_start)
            fields[2] = str(new_end)
            
            print('\t'.join(fields))

def main():
    if len(sys.argv) != 4:
        print("Usage: python translate_bed.py <agp_file.tsv> <hap> <bed_file.bed>")
        sys.exit(1)
    
    agp_file = sys.argv[1]
    hap = sys.argv[2]
    bed_file = sys.argv[3]
    
    # Parse the AGP file and build the interval tree
    trans = parse_agp(agp_file)
    
    # Translate the BED file coordinates
    translate_bed_coordinates(bed_file, trans, hap)

if __name__ == "__main__":
    main()