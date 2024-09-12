#!/bin/bash

set -e

vcf_file=$1
bam_file=$2
haplotype_reads=$3
output_fastq=$4
output_bed=${output_fastq/.fastq.gz/.bed}
threads=$5
haplotype=$6

echo "Processing VCF file: $vcf_file"
echo "Processing BAM file: $bam_file"
echo "Using haplotype reads file: $haplotype_reads"
echo "Output FASTQ file: $output_fastq"
echo "Using $threads threads"

# Extract SV regions
echo "Creating BED file: $output_bed"
bcftools query -f '%CHROM\t%POS\t%INFO/SVLEN\t%INFO/SVTYPE\n' $vcf_file | \
awk '{
    start=$2; 
    svlen=$3; 
    svtype=$4;
    if (svtype == "DEL" || svtype == "INV") 
        print $1"\t"start"\t"start+svlen"\t"svlen; 
    else 
        print $1"\t"start"\t"start+1"\t"svlen;
}' > $output_bed

echo "Number of SV regions found:"
wc -l $output_bed

echo "Number of haplotype-specific reads:"
wc -l $haplotype_reads

echo "Extracting spanning reads from $bam_file..."

# Process BAM file using regions from output_bed
samtools view -@ $threads --region-file $output_bed -h $bam_file | \
awk -v OFS='\t' -v haplotype_file="$haplotype_reads" '
    BEGIN {
        while (getline < haplotype_file) {
            haplotype_reads[$1] = 1
        }
        close(haplotype_file)
        reads_processed = 0
        reads_output = 0
    }
    NR==FNR {
        sv[$1][++s[$1]] = $2 "\t" $3 "\t" $4
        next
    }
    $1 ~ /^@/ {print; next}
    {
        reads_processed++
        if (reads_processed % 100000 == 0) {
            print "Processed " reads_processed " reads. Output " reads_output " reads." > "/dev/stderr"
        }
    }
    $1 in haplotype_reads {
        cigar = $6
        span = 0
        while(match(cigar, /[0-9]+[MDNX=]/)) {
            span += substr(cigar, RSTART, RLENGTH-1)
            cigar = substr(cigar, RSTART+RLENGTH)
        }
        start = $4
        end = start + span
        for(i=1; i<=s[$3]; i++) {
            split(sv[$3][i], a)
            if(start <= a[1] && end >= a[2] && span >= a[3]) {
                print  # This line prints the matching read
                reads_output++
                break
            }
        }
    }
    END {
        print "Total reads processed: " reads_processed > "/dev/stderr"
        print "Total reads output: " reads_output > "/dev/stderr"
    }
' $output_bed - | \
samtools fastq -@ $threads -c 6 -T '*' -0 $output_fastq

echo "Spanning reads extracted. Output FASTQ file info:"
ls -l $output_fastq