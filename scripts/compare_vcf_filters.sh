#!/bin/bash

# Initialize variables
qc_all=false
tsv_output=false
input_file=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --qc-all)
            qc_all=true
            shift
            ;;
        --tsv)
            tsv_output=true
            shift
            ;;
        *)
            input_file="$1"
            shift
            ;;
    esac
done

# Check if TSV file is provided
if [ -z "$input_file" ]; then
    echo "Usage: $0 [--qc-all] [--tsv] config/snakemake/samples.tsv"
    exit 1
fi

# Define references and base path
base_path="output/alignment"

# Set file suffix based on qc_all flag
suffix=".vcf.gz"
if [ "$qc_all" = true ]; then
    suffix=".qc_all.vcf.gz"
fi

# Print header
if [ "$tsv_output" = true ]; then
    printf "Specimen\tFilter\tCount(hg38_ref)\tCount(hg38_scaf)\tCount(T2T_scaf)\n"
else
    printf "%-16s %-16s %-16s %-16s %-16s\n" "Specimen" "Filter" "Count(hg38_ref)" "Count(hg38_scaf)" "Count(T2T_scaf)"
fi

# Read unique specimen names from TSV file (skip header if present)
tail -n +2 "$input_file" | cut -f1 | sort -u | while read specimen; do
    # Create temporary files for all references
    hg38_tmp=$(mktemp)
    hg38_scaf_tmp=$(mktemp)
    t2t_tmp=$(mktemp)
    
    # Process hg38 file
    zcat "${base_path}/hg38/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | cut -f7 | awk '{print ($1=="" || $1==".") ? "PASS" : $1}' | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$hg38_tmp"
    
    # Process hg38 scaffolded file
    zcat "${base_path}/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | cut -f7 | awk '{print ($1=="" || $1==".") ? "PASS" : $1}' | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$hg38_scaf_tmp"
    
    # Process T2T file
    zcat "${base_path}/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | cut -f7 | awk '{print ($1=="" || $1==".") ? "PASS" : $1}' | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$t2t_tmp"
    
    # Join the results and print
    join -a1 -a2 -e"0" -o "1.1 1.2 2.2" <(sort -k1,1 "$hg38_tmp") <(sort -k1,1 "$hg38_scaf_tmp") | \
    join -a1 -a2 -e"0" -o "1.1 1.2 1.3 2.2" - <(sort -k1,1 "$t2t_tmp") | \
    while read -r filter hg38_count hg38_scaf_count t2t_count; do
        if [ "$tsv_output" = true ]; then
            printf "%s\t%s\t%s\t%s\t%s\n" "${specimen}" "${filter}" "${hg38_count}" "${hg38_scaf_count}" "${t2t_count}"
        else
            printf "%-16s %-16s %-16s %-16s %-16s\n" "${specimen}" "${filter}" "${hg38_count}" "${hg38_scaf_count}" "${t2t_count}"
        fi
    done
    
    # Clean up temporary files
    rm "$hg38_tmp" "$hg38_scaf_tmp" "$t2t_tmp"
done