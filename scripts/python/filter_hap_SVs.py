#!/usr/bin/env python3

import sys
import logging
import subprocess
import json
from datetime import datetime

def setup_logging(log_file):
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def read_hap_reads(filename):
    hap_reads = set()
    with open(filename, 'r') as f:
        for line in f:
            hap_reads.add(line.strip())
    logging.info(f"Loaded {len(hap_reads)} hap read names")
    return hap_reads

def process_vcf(vcf_filename, output_filename, hap_reads):
    total_variants = 0
    filtered_variants = 0
    passed_variants = 0
    
    with subprocess.Popen(['bgzip', '-c', '-d', vcf_filename], stdout=subprocess.PIPE, universal_newlines=True) as vcf, \
         subprocess.Popen(['bgzip', '-c'], stdin=subprocess.PIPE, stdout=open(output_filename, 'wb'), universal_newlines=True) as out:
        for line in vcf.stdout:
            if line.startswith('#'):
                out.stdin.write(line)
                continue
            
            total_variants += 1
            
            fields = line.split('\t')
            info = fields[7]
            rnames = None
            for item in info.split(';'):
                if item.startswith('RNAMES='):
                    rnames = item[7:].split(',')
                    break
            
            if rnames is None or not any(read in hap_reads for read in rnames):
                out.stdin.write(line)
                passed_variants += 1
            else:
                filtered_variants += 1
            
            if total_variants % 100000 == 0:
                logging.info(f"Processed {total_variants} variants")

    logging.info(f"Total variants: {total_variants}")
    logging.info(f"Filtered variants: {filtered_variants}")
    logging.info(f"Passed variants: {passed_variants}")
    
    return total_variants, filtered_variants, passed_variants

def write_json_output(json_filename, total, filtered, passed):
    data = {
        "total_variants": total,
        "filtered_variants": filtered,
        "passed_variants": passed,
        "filtered_percentage": filtered / total * 100 if total > 0 else 0,
        "passed_percentage": passed / total * 100 if total > 0 else 0
    }
    
    with open(json_filename, 'w') as f:
        json.dump(data, f, indent=2)
    
    logging.info(f"Written JSON output to {json_filename}")

if __name__ == "__main__":
    hap_reads_file = sys.argv[1]
    input_vcf = sys.argv[2]
    output_vcf = sys.argv[3]
    json_file = sys.argv[4]
    log_file = sys.argv[5]

    setup_logging(log_file)
    
    logging.info("Starting VCF filtering process")
    hap_reads = read_hap_reads(hap_reads_file)
    total, filtered, passed = process_vcf(input_vcf, output_vcf, hap_reads)
    
    logging.info("VCF filtering process completed")
    logging.info(f"Filtered {filtered} out of {total} variants ({filtered/total:.2%})")
    logging.info(f"Retained {passed} variants ({passed/total:.2%})")

    # Create tabix index
    subprocess.run(['tabix', '-p', 'vcf', output_vcf])
    logging.info(f"Created tabix index for {output_vcf}")

    # Write JSON output
    write_json_output(json_file, total, filtered, passed)