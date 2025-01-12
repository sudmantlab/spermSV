#!/bin/bash

# Function to join header lines and create proper column names
process_header() {
    echo -e "SW_score\tdivergence\tdeleted\tinserted\tID\t"\
"query_start\tquery_end\tquery_left\tstrand\tmotif\t"\
"family\tconsensus_start\tconsensus_end\tconsensus_left"
}

# Function to process data lines
process_data() {
    # Skip first 3 lines (headers and blank line), then process data
    tail -n +4 "$1" | awk '
    NF > 0 {  # Skip empty lines
        # Remove leading tab if present
        sub(/^\t/, "")
        
        # Get the basic fields
        sw_score=$1
        divergence=$2
        deleted=$3
        inserted=$4
        query_sequence=$5
        query_start=$6
        query_end=$7
        query_left=$8
        
        # Split the strand and motif from the 9th field
        strand=$9
        motif=$10
        for(i=11; i<=NF; i++) {
            if($i ~ /^(Simple_repeat|SINE|LINE|LTR|Retroposon|Satellite)/) {
                break
            }
            motif = motif " " $i
        }
        
        # Get family and remaining fields
        family=$i
        if($(i+1) != "") {
            if($(i+1) ~ /^[0-9]+$/ || $(i+1) ~ /^\([0-9]+\)$/) {
                # Next field is a number or (number), so no family part
                consensus_start=$(i+1)
                consensus_end=$(i+2)
                consensus_left=$(i+3)
                id=($(i+4) != "*" ? $(i+4) : "")
            } else {
                # Include family part
                family = family "/" $(i+1)
                consensus_start=$(i+2)
                consensus_end=$(i+3)
                consensus_left=$(i+4)
            }
        }
        
        # Clean up parentheses
        gsub(/[()]/, "", query_left)
        gsub(/[()]/, "", consensus_start)
        gsub(/[()]/, "", consensus_end)
        gsub(/[()]/, "", consensus_left)
        
        # Print fields in TSV format
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            sw_score, divergence, deleted, inserted, query_sequence, query_start, query_end, query_left,
            strand, motif, family,
            consensus_start, consensus_end, consensus_left
    }'
}

# Main execution
if [ -z "$1" ]; then
    echo "Usage: $0 <repeatmasker_output_file>"
    echo "The tsv will be saved in the same directory with a .tsv suffix."
    exit 1
fi

# Process header and data
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing input from $1"
{
    process_header
    process_data "$1"
} > "${1%.out}.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processed $(wc -l < "${1%.out}.tsv") rows."
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Output saved to ${1%.out}.tsv"