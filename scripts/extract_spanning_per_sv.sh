#!/bin/bash

set -e

# Check if we have the correct number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <sniffles_vcf> <bam_file> <chrom> <output_var_rnames> <output_dir> <threads> <report>"
    exit 1
fi

sniffles_vcf=$1
bam_file=$2
chrom=$3
output_var_rnames=$4
output_dir=$5
threads=$6
report=$7

# Check if input files exist
for file in "$sniffles_vcf" "$bam_file"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist."
        exit 1
    fi
done

echo "Processing Sniffles VCF file: $sniffles_vcf"
echo "Processing chromosome: $chrom"
echo "Processing BAM file: $bam_file"
echo "Output variant read names file: $output_var_rnames"
echo "Output directory: $output_dir"
echo "Using $threads threads"
echo "Report: $report"

# Create output directories if they don't exist
mkdir -p "$output_dir"

# Extract variant IDs, positions, and read names from Sniffles VCF for the specified chromosome
echo "Extracting variant information for chromosome $chrom..."
bcftools filter -r "$chrom" "$sniffles_vcf" | \
bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/SVLEN\t%INFO/SVTYPE\t%INFO/RNAMES\n' > "$output_var_rnames"

# Initialize log file
echo -e "Variant_ID\tReads_at_locus\tRNAMES_count\tSpanning_reads" > "$report"

# Process each variant
while IFS=$'\t' read -r var_id chrom pos svlen svtype rnames; do
    # Create a temporary file with read names for this variant
    temp_rnames=$(mktemp)
    echo "$rnames" | tr ',' '\n' > "$temp_rnames"
    rnames_count=$(wc -l < "$temp_rnames")
    
    # Calculate SV end position
    if [[ "$svtype" == "DEL" || "$svtype" == "INV" ]]; then
        end=$((pos + 1))
    else
        end=$((pos + svlen))
    fi
    
    # Extract reads for this variant, check if they span the SV, and create BAM file
    samtools view -@ "$threads" -h "$bam_file" "$chrom:$pos-$end" | \
    awk -v OFS='\t' -v var_start="$pos" -v var_end="$end" -v var_len="$svlen" -v temp_file="$temp_rnames" '
    BEGIN {
        while (getline < temp_file) {
            variant_reads[$1] = 1
        }
        close(temp_file)
        reads_at_locus = 0
        spanning_reads = 0
    }
    $1 ~ /^@/ {print; next}
    {
        reads_at_locus++
        if ($1 in variant_reads) {
            cigar = $6
            span = 0
            while(match(cigar, /[0-9]+[MDNX=]/)) {
                span += substr(cigar, RSTART, RLENGTH-1)
                cigar = substr(cigar, RSTART+RLENGTH)
            }
            start = $4
            end = start + span
            if(start <= var_start && end >= var_end && span >= var_len) {
                print
                spanning_reads++
            }
        }
    }
    END {
        print reads_at_locus "\t" length(variant_reads) "\t" spanning_reads > "/dev/stderr"
    }
    ' 2> >(read locus_reads_count rnames_count spanning_reads_count; \
           echo -e "$var_id\t$locus_reads_count\t$rnames_count\t$spanning_reads_count" >> "$report") | \
    samtools view -@ "$threads" -b > "$output_dir/${var_id}.bam"
    
    # Convert BAM to FASTQ
    samtools fastq -@ "$threads" "$output_dir/${var_id}.bam" | gzip > "$output_dir/${var_id}.fastq.gz"
    
    # Remove temporary file
    rm "$temp_rnames"
done < "$output_var_rnames"

echo "All variants processed."