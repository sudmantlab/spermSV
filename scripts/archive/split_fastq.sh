#!/bin/bash

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 input.fastq.gz num_chunks num_cores total_mem_gb output_prefix format"
    echo "Example: $0 input.fastq.gz 4 40 256 ASD6463_chunk fastq"
    echo "Format options: fastq or fasta"
    exit 1
fi

input_file=$1
num_chunks=$2
num_cores=$3
total_mem=$4
output_prefix=$5
format=$6

# Validate format argument
if [[ "$format" != "fastq" && "$format" != "fasta" ]]; then
    echo "Error: format must be either 'fastq' or 'fasta'"
    exit 1
fi

# Calculate cores per chunk while leaving some headroom
cores_per_chunk=$((num_cores / num_chunks - 1))
# Calculate memory per chunk (in GB)
mem_per_chunk=$((total_mem / num_chunks))
# Calculate block size as 2% of memory per chunk, converted to MB
block_size=$((mem_per_chunk * 20))

# Use pigz for parallel compression with multiple cores per chunk
export PIGZ="-p $cores_per_chunk"

echo "Using settings:"
echo "Cores per chunk: $cores_per_chunk"
echo "Memory per chunk: ${mem_per_chunk}GB"
echo "Block size: ${block_size}M"
echo "Output format: $format"

# Create RAM disk if we have enough memory
if [ -d "/dev/shm" ]; then
    TMPDIR="/dev/shm/fastq_split_$$"
    mkdir -p "$TMPDIR"
    trap 'rm -rf "$TMPDIR"' EXIT
else
    TMPDIR=$(mktemp -d)
    trap 'rm -rf "$TMPDIR"' EXIT
fi

# Create chunk number file
seq 1 $num_chunks > "$TMPDIR/chunk_nums.txt"

# Use parallel to split and compress simultaneously
# Block size is now calculated based on available memory
if [ "$format" = "fastq" ]; then
    zcat -f "$input_file" | \
        parallel --pipe -N4 --round-robin -j"$num_chunks" --block ${block_size}M --files \
        "cat > >(pigz $PIGZ > ${output_prefix}_{#}.fastq.gz)" :::: "$TMPDIR/chunk_nums.txt"
else
    # For fasta, convert fastq to fasta then split
    zcat -f "$input_file" | \
        awk 'NR%4==1{printf ">%s\n",substr($0,2)}NR%4==2{print}' | \
        parallel --pipe --round-robin -j"$num_chunks" --block ${block_size}M --files \
        "cat > >(pigz $PIGZ > ${output_prefix}_{#}.fasta.gz)" :::: "$TMPDIR/chunk_nums.txt"
fi

# Print file paths
echo "Files saved at:"
for f in ${output_prefix}_*.${format}.gz; do 
    echo "$f"
done

# Clean up temporary files
rm -rf "$TMPDIR"
trap - EXIT
