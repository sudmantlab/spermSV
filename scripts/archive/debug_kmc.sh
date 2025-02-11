#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 <input.fastq.gz> <kmc_path> <output_dir>"
    echo "Example: $0 ASD6463_chunk_1.fastq.gz /path/to/kmc output_tests"
    exit 1
}

# Check arguments
if [ "$#" -ne 3 ]; then
    usage
fi

INPUT_FILE=$1
KMC_PATH=$2
OUTPUT_DIR=$3

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/tmp"

# Log file for results
RESULTS_LOG="$OUTPUT_DIR/results_summary.txt"
# Log the command used to launch the script
echo "Command: $0 $*" >> "$RESULTS_LOG"
echo "KMC Testing Results" > "$RESULTS_LOG"
echo "===================" >> "$RESULTS_LOG"
echo "Input file: $INPUT_FILE" >> "$RESULTS_LOG"
echo "Date: $(date)" >> "$RESULTS_LOG"
echo "" >> "$RESULTS_LOG"

# Function to run KMC with specific parameters and log results
run_kmc_test() {
    local mem=$1
    local threads=$2
    local subset_size=$3
    local test_name="mem${mem}_t${threads}_s${subset_size}"
    local log_file="$OUTPUT_DIR/logs/${test_name}.log"
    local tmp_dir="$OUTPUT_DIR/tmp/${test_name}"
    local input_file="$INPUT_FILE"
    
    mkdir -p "$tmp_dir"
    
    echo "Running test: $test_name"
    echo "===================" >> "$RESULTS_LOG"
    echo "Test: $test_name" >> "$RESULTS_LOG"
    echo "Memory: ${mem}GB" >> "$RESULTS_LOG"
    echo "Threads: $threads" >> "$RESULTS_LOG"
    echo "Subset size: $subset_size reads" >> "$RESULTS_LOG"
    
    # Create subset if needed
    if [ "$subset_size" != "full" ]; then
        input_file="$tmp_dir/subset.fastq.gz"
        echo "Creating subset with $subset_size reads..."
        zcat "$INPUT_FILE" | head -n $((subset_size * 4)) | gzip > "$input_file"
    fi
    
    # Run KMC with time and memory monitoring
    {
        echo "Command output:" >> "$RESULTS_LOG"
        /usr/bin/time -v $KMC_PATH \
            -k29 -m$mem -sm -t$threads -hp -v \
            "$input_file" \
            "$tmp_dir/output" \
            "$tmp_dir" \
            2>&1 | tee -a "$log_file" "$RESULTS_LOG"
            
        echo "" >> "$RESULTS_LOG"
        echo "Resource usage:" >> "$RESULTS_LOG"
        grep "Maximum resident set size" "$log_file" >> "$RESULTS_LOG"
        grep "User time" "$log_file" >> "$RESULTS_LOG"
        grep "System time" "$log_file" >> "$RESULTS_LOG"
    } 2>&1
    
    # Check if KMC produced output files
    if [ -f "$tmp_dir/output.kff" ]; then
        echo "Output file created successfully" >> "$RESULTS_LOG"
        # Get file size
        size=$(ls -lh "$tmp_dir/output.kff" | awk '{print $5}')
        echo "Output file size: $size" >> "$RESULTS_LOG"
    else
        echo "No output file created" >> "$RESULTS_LOG"
    fi
    
    echo "" >> "$RESULTS_LOG"
    echo "-------------------" >> "$RESULTS_LOG"
    echo "" >> "$RESULTS_LOG"
}

# Test different memory configurations
echo "Testing memory impact..."
for mem in 100 200 300; do
    run_kmc_test $mem 40 "full"
done

# Test different input sizes with optimal memory and threads
echo "Testing input size impact..."
for size in 100000 500000 1000000 5000000; do
    run_kmc_test 300 40 $size
done

echo "Testing complete. Results are in $RESULTS_LOG"

# Create summary visualization of results
echo "Creating results summary..."
{
    echo "Summary of KMC Tests"
    echo "==================="
    echo ""
    echo "Memory Impact (with 40 threads):"
    grep "Memory:" "$RESULTS_LOG" | cut -d' ' -f2
    grep "Output file size:" "$RESULTS_LOG" | cut -d' ' -f4
    
    echo ""
    echo "Input Size Impact:"
    grep "Subset size:" "$RESULTS_LOG" | cut -d' ' -f3
    grep "Maximum resident set size" "$RESULTS_LOG" | awk '{print $6}'
} > "$OUTPUT_DIR/summary.txt"