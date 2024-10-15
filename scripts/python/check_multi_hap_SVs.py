#!/usr/bin/env python3

import sys
import gzip
import json
import logging
import statistics

def setup_logging(log_file):
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def read_hap_rnames(filename):
    with open(filename, 'r') as f:
        return set(line.strip() for line in f)

def process_vcf(input_vcf, mixed_output, spike_hap_only_output, base_hap_only_output, hap1_rnames, hap2_rnames, report_json, rnames_json):
    """
    Input:
        input_vcf: (Sniffles) VCF file path with multiple haplotype variants 
        mixed_output: Output VCF path for mixed haplotype variants
        spike_hap_only_output: Output VCF path for haplotype 1 (spiked) only variants
        base_hap_only_output: Output VCF path for haplotype 2 (base) only variants
        hap2_rnames: txt file path with haplotype 2 read names
        report_json: json file path for output report
    """
    total_variants = 0
    mixed_variants = 0
    spike_hap_only_variants = 0
    base_hap_only_variants = 0
    rnames_counts = []
    rnames_log = dict()

    with gzip.open(input_vcf, 'rt') as vcf, \
         gzip.open(mixed_output, 'wt') as mixed, \
         gzip.open(spike_hap_only_output, 'wt') as spike_hap_only, \
         gzip.open(base_hap_only_output, 'wt') as base_hap_only:
        
        # Write headers to all output files
        for line in vcf:
            if line.startswith('#'):
                mixed.write(line)
                spike_hap_only.write(line)
                base_hap_only.write(line)
            else:
                break

        # Process variant records
        while True:
            line = vcf.readline()
            if not line:
                break

            total_variants += 1
            fields = line.split('\t')
            info = fields[7]
            rnames = None
            for item in info.split(';'):
                if item.startswith('RNAMES='):
                    rnames = item[7:].split(',')
                    if type(rnames) == str: # force list even if single rname (str)
                        rnames = [rnames]
                    rnames_counts.append(len(rnames))
                    break

            if rnames:
                hap1_reads = set(rnames) & hap1_rnames
                hap2_reads = set(rnames) & hap2_rnames
                rnames_log[fields[2]] = {
                    'count': len(rnames), 
                    'hap1_rnames': list(hap1_reads), 
                    'hap2_rnames': list(hap2_reads)} # variant ID to rnames even if None
                if hap1_reads and hap2_reads:
                    mixed.write(line)
                    mixed_variants += 1
                elif hap1_reads:
                    spike_hap_only.write(line)
                    spike_hap_only_variants += 1
                elif hap2_reads:
                    base_hap_only.write(line)
                    base_hap_only_variants += 1

            else:
                logging.warning(f"No RNAMES found for variant with ID {fields[2]}:\n{line}")

            if total_variants % 100000 == 0:
                logging.info(f"Processed {total_variants} variants")

    # Generate report
    report = {
        "total_variants": total_variants,
        "mixed_variants": mixed_variants,
        "spike_hap_only_variants": spike_hap_only_variants,
        "base_hap_only_variants": base_hap_only_variants,
        "mixed_percentage": (mixed_variants / total_variants) * 100 if total_variants > 0 else 0,
        "spike_hap_only_percentage": (spike_hap_only_variants / total_variants) * 100 if total_variants > 0 else 0,
        "base_hap_only_percentage": (base_hap_only_variants / total_variants) * 100 if total_variants > 0 else 0,
        "rnames_mean": statistics.mean(rnames_counts) if rnames_counts else 0,
        "rnames_median": statistics.median(rnames_counts) if rnames_counts else 0,
        "rnames_mode": statistics.mode(rnames_counts) if rnames_counts else 0,
        "rnames_stdev": statistics.stdev(rnames_counts) if len(rnames_counts) > 1 else 0,
    }

    with open(report_json, 'w') as f:
        json.dump(report, f, indent=2)

    with open(rnames_json, 'w') as f:
        json.dump(rnames_log, f, indent=2)

    logging.info(f"Total variants: {total_variants}")
    logging.info(f"Mixed haplotype variants: {mixed_variants}")
    logging.info(f"Spiked haplotype only variants: {spike_hap_only_variants}")
    logging.info(f"Base haplotype only variants: {base_hap_only_variants}")
    logging.info(f"Mean RNAMES count: {report['rnames_mean']:.2f}")
    logging.info(f"Median RNAMES count: {report['rnames_median']:.2f}")
    logging.info(f"Mode RNAMES count: {report['rnames_mode']:.2f}")
    logging.info(f"Standard deviation of RNAMES count: {report['rnames_stdev']:.2f}")

if __name__ == "__main__":
    input_vcf = sys.argv[1]
    hap1_rnames_file = sys.argv[2]
    hap2_rnames_file = sys.argv[3]
    mixed_output = sys.argv[4]
    spike_hap_only_output = sys.argv[5]
    base_hap_only_output = sys.argv[6]
    report_json = sys.argv[7]
    rnames_json = sys.argv[8]
    log_file = sys.argv[9]

    setup_logging(log_file)

    logging.info("Starting VCF processing")
    hap1_rnames, hap2_rnames = read_hap_rnames(hap1_rnames_file), read_hap_rnames(hap2_rnames_file)
    process_vcf(input_vcf, mixed_output, spike_hap_only_output, base_hap_only_output, hap1_rnames, hap2_rnames, report_json, rnames_json)
    logging.info("VCF processing completed")