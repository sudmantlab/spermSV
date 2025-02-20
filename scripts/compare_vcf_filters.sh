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


# Define all possible filters
read -r -d '' FILTERS << EOM
PASS
GT
SUPPORT_MIN
STDEV_POS
STDEV_LEN
COV_MIN
COV_MIN_GT
COV_CHANGE_DEL
COV_CHANGE_DUP
COV_CHANGE_INS
COV_CHANGE_FRAC_US
COV_CHANGE_FRAC_SC
COV_CHANGE_FRAC_CE
COV_CHANGE_FRAC_ED
MOSAIC_VAF
NOT_MOSAIC_VAF
ALN_NM
STRAND_BND
STRAND
STRAND_MOSAIC
SVLEN_MIN
SVLEN_MIN_MOSAIC
EOM

# Read unique specimen names from TSV file (skip header if present)
tail -n +2 "$input_file" | cut -f1 | sort -u | while read specimen; do
    # Create temporary files for all references
    hg38_tmp=$(mktemp)
    hg38_scaf_tmp=$(mktemp)
    t2t_tmp=$(mktemp)
    
    # For all files, exclude BND lines
    # Process hg38 file
    gzip -dc "${base_path}/hg38/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | grep -v "BND" | cut -f7 | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$hg38_tmp"
    
    # Process hg38 scaffolded file
    gzip -dc "${base_path}/hg38_scaffolded/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | grep -v "BND" | cut -f7 | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$hg38_scaf_tmp"
    
    # Process T2T file
    gzip -dc "${base_path}/T2T_scaffolded/minimap2/standard/variants/sniffles_mosaic/${specimen}${suffix}" 2>/dev/null | \
    grep -v "^#" | grep -v "BND" | cut -f7 | sort | uniq -c | sort -nr | \
    awk '{print $2" "$1}' > "$t2t_tmp"

    # Print results for each filter, defaulting to 0 if not found
    while read -r filter; do
        hg38_count=$(awk -v filter="$filter" '$1 == filter {print $2}' "$hg38_tmp" || echo "0")
        hg38_scaf_count=$(awk -v filter="$filter" '$1 == filter {print $2}' "$hg38_scaf_tmp" || echo "0")
        t2t_count=$(awk -v filter="$filter" '$1 == filter {print $2}' "$t2t_tmp" || echo "0")
        
        # If no match was found, awk will output nothing, so we need to convert empty strings to 0
        hg38_count=${hg38_count:-0}
        hg38_scaf_count=${hg38_scaf_count:-0}
        t2t_count=${t2t_count:-0}
        
        if [ "$tsv_output" = true ]; then
            printf "%s\t%s\t%s\t%s\t%s\n" "${specimen}" "${filter}" "${hg38_count}" "${hg38_scaf_count}" "${t2t_count}"
        else
            printf "%-16s %-16s %-16s %-16s %-16s\n" "${specimen}" "${filter}" "${hg38_count}" "${hg38_scaf_count}" "${t2t_count}"
        fi
    done <<< "$FILTERS"

    # Clean up temporary files
    rm "$hg38_tmp" "$hg38_scaf_tmp" "$t2t_tmp"
done