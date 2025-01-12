#!/usr/bin/env python3
import pysam
import json
from collections import defaultdict
from dataclasses import dataclass, asdict, field
import argparse
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
import subprocess

@dataclass
class ReadStats:
    """Data class to store alignment statistics for a single read"""
    total_alignments: int = 0
    primary_alignments: int = 0
    secondary_alignments: int = 0
    is_unmapped: bool = False
    primary_contig: str = ""  # Changed default from "ref" to empty string
    secondary_contigs: list = field(default_factory=list)  # Store secondary alignment contigs

    def to_dict(self):
        """Convert ReadStats object to dictionary for JSON serialization"""
        return {
            'total_alignments': self.total_alignments,
            'primary_alignments': self.primary_alignments,
            'secondary_alignments': self.secondary_alignments,
            'is_unmapped': self.is_unmapped,
            'primary_contig': self.primary_contig,
            'secondary_contigs': self.secondary_contigs
        }

def get_contig_chunks(bam_path, chunk_size=1000000):
    """
    Generate chunks of genomic regions for each contig.
    
    Args:
        bam_path: Path to BAM file
        chunk_size: Size of each genomic chunk in base pairs
        
    Returns:
        List of tuples (contig, start, end) for each chunk
    """
    chunks = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for contig in bam.references:
            contig_length = bam.get_reference_length(contig)
            for start in range(0, contig_length, chunk_size):
                end = min(start + chunk_size, contig_length)
                chunks.append((contig, start, end))
    return chunks

def process_chunk(args) -> tuple[dict[str, ReadStats], dict]:
    """
    Process a genomic chunk of BAM file and return read stats, modified reads, and chunk summary.

    Reads are modified to mapQ = 254 where a primary and a single secondary exists on the matched opposing haplotype scaffold.
    """
    bam_file, contig, start, end = args
    read_stats = defaultdict(ReadStats)
    current_read = None
    current_alignments = []
    
    # Initialize chunk summary
    chunk_summary = {
        "single_alignment_reads": 0,
        "secondary_alignment_reads": {
            "count": 0,
            "histogram": defaultdict(int),
            "single_opposing_secondary": 0
        },
        "unmapped_reads": 0
    }
    
    # logging.info(f"Processing chunk {contig}:{start:,}-{end:,}.")
    
    with pysam.AlignmentFile(bam_file, "rb", threads=1) as bam:
        for read in bam.fetch(contig, start, end, until_eof=True):
            if read.query_name is None:
                continue
            # Update statistics for each read
            if read.is_unmapped:
                read_stats[read.query_name].is_unmapped = read.is_unmapped
                chunk_summary["unmapped_reads"] += 1
            elif read.is_secondary:
                read_stats[read.query_name].secondary_alignments += 1
                read_stats[read.query_name].total_alignments += 1
                read_stats[read.query_name].secondary_contigs.append(read.reference_name)
            else: # thus primary
                read_stats[read.query_name].primary_alignments += 1
                read_stats[read.query_name].total_alignments += 1
                read_stats[read.query_name].primary_contig = read.reference_name
        
    # Collate summary stats for this chunk
    for read_name, stats in read_stats.items():
        if stats.secondary_alignments == 0:
            chunk_summary["single_alignment_reads"] += 1
        else:
            chunk_summary["secondary_alignment_reads"]["count"] += 1
            chunk_summary["secondary_alignment_reads"]["histogram"][stats.secondary_alignments] += 1
                
    return dict(read_stats), chunk_summary

def is_opposing_haplotype(contig1, contig2):
    """
    Check if contigs are same chromosome but different haplotypes.
    Args:
        contig1: String containing contig name (e.g. 'chr1_RagTag_hap2')
        contig2: List containing single contig name (e.g. ['chr1_RagTag_hap1'])
    """
    # Obtain contig + hap basenames
    contig1 = contig1.split('_')
    base1, hap1 = contig1[0], contig1[-1]
    contig2 = contig2.split('_')
    base2, hap2 = contig2[0], contig2[-1]
    
    # Check if base names match and haplotypes differ
    return (base1 == base2 and
            hap1 != hap2)

def merge_stats(chunk_readstats, chunk_summaries):
    """
    Merge read statistics and summaries from multiple chunks.
    
    Args:
        chunk_readstats: List of dictionaries containing read statistics from each chunk
        chunk_summaries: List of dictionaries containing summary statistics from each chunk
        
    Returns:
        Tuple of (merged read statistics, merged summary statistics)
    """
    merged_readstats = dict()
    reads_to_modify = set()
    
    # Iterate through read stats from each chunk
    logging.info(f"Merging results from {len(chunk_readstats)} chunks.")
    for chunk_dict in chunk_readstats:
        for read_name, stats in chunk_dict.items():
            if read_name in merged_readstats:
                # Only update stats if we've seen this read before
                merged_readstats[read_name].total_alignments += stats.total_alignments
                merged_readstats[read_name].primary_alignments += stats.primary_alignments
                merged_readstats[read_name].secondary_alignments += stats.secondary_alignments
                merged_readstats[read_name].is_unmapped |= stats.is_unmapped
                merged_readstats[read_name].secondary_contigs.extend(stats.secondary_contigs)
                # Update primary contig if this chunk has one
                if stats.primary_contig:
                    merged_readstats[read_name].primary_contig = stats.primary_contig
            else:
                # For new reads, just copy the stats object directly
                merged_readstats[read_name] = ReadStats(
                    total_alignments=stats.total_alignments,
                    primary_alignments=stats.primary_alignments,
                    secondary_alignments=stats.secondary_alignments,
                    is_unmapped=stats.is_unmapped,
                    primary_contig=stats.primary_contig,
                    secondary_contigs=stats.secondary_contigs.copy()
                )
    # Check for reads with single secondaries on matched scaffold, opposing haplotype
    for read_name, stats in merged_readstats.items():
        if stats.primary_alignments == 1 and stats.secondary_alignments == 1:
            if is_opposing_haplotype(stats.primary_contig, stats.secondary_contigs[0]):
                reads_to_modify.add(read_name)
            
    # Create placeholder for merged summary statistics
    merged_summary = {
        "single_alignment_reads": {
            "count": sum(s["single_alignment_reads"] for s in chunk_summaries)
        },
        "secondary_alignment_reads": {
            "count": sum(s["secondary_alignment_reads"]["count"] for s in chunk_summaries),
            "histogram": defaultdict(int),
            "single_opposing_secondary": len(reads_to_modify)
        },
        "unmapped_reads": {
            "count": sum(s["unmapped_reads"] for s in chunk_summaries)
        }
    }
    
    # Merge histograms
    for summary in chunk_summaries:
        for count, freq in summary["secondary_alignment_reads"]["histogram"].items():
            merged_summary["secondary_alignment_reads"]["histogram"][count] += freq
            
    # Convert histogram to string keys
    merged_summary["secondary_alignment_reads"]["histogram"] = {
        str(k): v for k, v in sorted(merged_summary["secondary_alignment_reads"]["histogram"].items())
    }
    
    merged_summary["total_reads"] = (
        merged_summary["single_alignment_reads"]["count"] +
        merged_summary["secondary_alignment_reads"]["count"] +
        merged_summary["unmapped_reads"]["count"]
    )
    
    return merged_readstats, merged_summary, reads_to_modify

def process_bam(bam_path, threads=None, chunk_size=1_000_000):
    """
    Process BAM file and generate alignment statistics per read and aggregated across the BAM.
    After processing, returns a set of reads that have single matched scaffold, opposing haplotype secondaries for
    modification.
    
    Args:
        bam_path: Path to input BAM file
        threads: Number of threads to use (defaults to CPU count)
        chunk_size: Size of genomic chunks in base pairs
    
    Returns:
        merged_readstats: Statistics on mapping for each read (dict)
        merged_summary: Statistics on read mapping aggregated across the whole BAM (dict)
        reads_to_modify: Read names for single matched scaffold, opposing haplotype secondaries (set)
    """
    if threads is None:
        threads = multiprocessing.cpu_count()
    
    # Get genomic chunks
    chunks = get_contig_chunks(bam_path, chunk_size)
    chunk_args = [(bam_path, contig, start, end) for contig, start, end in chunks]
    
    logging.info(f"Processing {len(chunks)} genomic chunks of {chunk_size:,} bp with {threads} threads")
    
    # Process chunks in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        chunk_results = list(executor.map(process_chunk, chunk_args))
    
    # Combine results
    merged_readstats, merged_summary, reads_to_modify = merge_stats(
        [stats for stats, _ in chunk_results],
        [summary for _, summary in chunk_results]
    )
    
    return merged_readstats, merged_summary, reads_to_modify

def modify_bam(bam_path, output_bam, reads_to_modify):
    """
    Given the list of reads to modify, set mapQ = 254 for each of the reads in the bamfile.
    """
    logging.info("Writing modified BAM file...")
    with pysam.AlignmentFile(bam_path, "rb") as input_bam:
        with pysam.AlignmentFile(output_bam, "wb", template=input_bam) as output_bam:
            for read in input_bam.fetch(until_eof=True):
                if read.query_name in reads_to_modify:
                    read.mapping_quality = 254
                output_bam.write(read)

def main():
    parser = argparse.ArgumentParser(description='Process BAM file and generate alignment statistics')
    parser.add_argument('bam_file', help='Input BAM file path')
    parser.add_argument('--output', '-o', help='Output BAM file path', required=True)
    parser.add_argument('--threads', '-t', type=int, help='Number of threads to use')
    parser.add_argument('--chunk-size', type=int, default=1_000_000,
                       help='Size of genomic chunks to process in base pairs')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Generate output paths
    output_base = Path(args.output)
    readstats_json = output_base.with_suffix('.readstats.json')
    summary_json = output_base.with_suffix('.summary.json')
    modified_txt = output_base.with_suffix('.reads.txt')
    
    logging.info(f"Processing {args.bam_file} with {args.threads or multiprocessing.cpu_count()} threads")
    
    # Process BAM file statistics and flag reads for modification
    read_stats, summary_stats, reads_to_modify = process_bam(args.bam_file, args.threads, args.chunk_size)

    # Convert ReadStats objects to dictionaries for JSON serialization
    read_stats_dict = {read_name: stats.to_dict() for read_name, stats in read_stats.items()}

    # Write statistics and modified reads
    with open(readstats_json, 'w') as f:
        json.dump(read_stats_dict, f, indent=2)
        logging.info(f"Detailed results written to {readstats_json}")
    with open(summary_json, 'w') as f:
        json.dump(summary_stats, f, indent=2)
        logging.info(f"Summary results written to {summary_json}")
    with open(modified_txt, 'w') as f:
        for read_name in reads_to_modify:
            f.write(f"{read_name}\n")
        logging.info(f"Modified read names written to {modified_txt}")

    # Modify flagged reads in output BAM
    modify_bam(args.bam_file, args.output, reads_to_modify)
    logging.info(f"Modified BAM written to {args.output}")

    # Index output BAM
    logging.info("Indexing output BAM file...")
    pysam.index(args.output)

if __name__ == '__main__':
    main()