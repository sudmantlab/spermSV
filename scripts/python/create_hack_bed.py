import random
import pandas as pd
import time
from Bio import SeqIO
import argparse
import os
import re

def get_random_position(chrom, seq_length, chrom_info):
    chrom_length = chrom_info[chrom_info['chrom'] == chrom]['length'].values[0]
    return random.randint(1, chrom_length - seq_length)

def create_hack_bed(chrom_info_file, alu_count, l1_count, alu_fasta, L1_fasta, output_file, exclude_random, verbose):
    # Load chromosome information
    chrom_info = pd.read_csv(chrom_info_file, sep='\t', names=['chrom', 'length'])
    chrom_info = chrom_info[chrom_info['chrom'].str.startswith('chr')]

    if exclude_random:
        if verbose:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Filtering out random/unknown contigs from {chrom_info_file}")
        # Filter chromosomes to include only chr1-22, X, and Y
        chrom_pattern = re.compile(r'^chr(\d{1,2}|X|Y)$')
        chrom_info = chrom_info[chrom_info['chrom'].str.match(chrom_pattern)]

    # Load AluY and L1 sequences
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Loading AluY sequences from {alu_fasta}")
    aluy_seqs = list(SeqIO.parse(alu_fasta, 'fasta'))
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Loading L1 sequences from {L1_fasta}")
    l1_seqs = list(SeqIO.parse(L1_fasta, 'fasta'))
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Finished loading sequences")

    # Generate insertions
    insertions = []

    # Random AluY insertions
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Generating {alu_count} random AluY insertions")
    for i in range(alu_count):
        seq = random.choice(aluy_seqs).seq
        seq = ''.join(random.choice('ACTG') if base not in 'ACTG' else base for base in str(seq))
        chrom = random.choice(chrom_info['chrom'].tolist())
        pos = get_random_position(chrom, len(seq), chrom_info)
        insertions.append((chrom, pos, pos+1, 'insertion', str(seq), 0))
        if verbose and (i+1) % 1000 == 0:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Generated {i+1} AluY insertions")

    # Random L1 insertions
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Generating {l1_count} random L1 insertions")
    for i in range(l1_count):
        seq = random.choice(l1_seqs).seq
        seq = ''.join(random.choice('ACTG') if base not in 'ACTG' else base for base in str(seq))
        chrom = random.choice(chrom_info['chrom'].tolist())
        pos = get_random_position(chrom, len(seq), chrom_info)
        insertions.append((chrom, pos, pos+1, 'insertion', str(seq), 0))
        if verbose and (i+1) % 1000 == 0:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Generated {i+1} L1 insertions")

    # Sort insertions by chromosome and position
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Sorting insertions")
    insertions.sort(key=lambda x: (x[0], x[1]))

    # Create path if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write to BED file
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Writing insertions to BED file")
    with open(output_file, 'w') as f:
        for ins in insertions:
            f.write('\t'.join(map(str, ins)) + '\n')
    if verbose:
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Finished writing {len(insertions)} insertions to BED file")

    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Created HACk BED file with {len(insertions)} insertions")
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - Saved to {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a HACk BED file with random Alu and L1 insertions. Example usage: python create_hack_bed.py --chrom_info chrom_info.tsv --alu_count 1000 --l1_count 500 --alu_fasta alu.fasta --L1_fasta L1.fasta --output output.bed --verbose > log.txt")
    parser.add_argument("--chrom_info", help="Path to the chromosome info TSV file", required=True)
    parser.add_argument("--alu_count", type=int, help="Number of Alu insertions", required=True)
    parser.add_argument("--L1_count", type=int, help="Number of L1 insertions", required=True)
    parser.add_argument("--alu_fasta", help="Path to the Alu fasta file", required=True)
    parser.add_argument("--L1_fasta", help="Path to the L1 fasta file", required=True)
    parser.add_argument("--output", help="Output BED file name", required=True)
    parser.add_argument("--exclude_random", help="Exclude random/unknown contigs (default True)", action="store_false")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output. To save to a log file, redirect stdout: python create_hack_bed.py ... --verbose > log.txt")
    
    args = parser.parse_args()

    create_hack_bed(args.chrom_info, args.alu_count, args.L1_count, args.alu_fasta, args.L1_fasta, args.output, args.exclude_random, args.verbose)