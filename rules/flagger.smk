rule prepare_flagger_inputs:
    input:
        fasta = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta",
        fai = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta.fai"
    output:
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.bed})
        
        # Create whole genome bed from fai
        cat {input.fai} | awk '{{print $1"\t0\t"$2}}' > {output.bed}
        
        # Create annotations JSON
        echo '{{"whole_genome": "{output.bed}"}}' > {output.json}
        """

use rule prepare_flagger_inputs as HG002_unscaffolded_prepare with:
    input:
        fasta = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa",
        fai = "output/assembly/hifiasm/HG002/HG002.unscaffolded.diploid.fa.fai"
    output:
        bed = "output/assembly/flagger/HG002/unscaffolded/HG002_whole_genome.bed",
        json = "output/assembly/flagger/HG002/unscaffolded/HG002_annotations_path.json"

# Convert BAM to coverage file
rule bam_to_coverage:
    input:
        bam = "output/alignment/hg38_scaffolded/minimap2/standard/mapped/{specimen}.sorted.merged.bam",
        bed = "output/assembly/flagger/{specimen}/whole_genome.bed",
        json = "output/assembly/flagger/{specimen}/annotations_path.json"
    output:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            bam2cov \
                --bam {input.bam} \
                --output {output.cov} \
                --annotationJson {input.json} \
                --threads {threads} \
                --baselineAnnotation whole_genome
        """

# Run HMM-Flagger
rule run_hmm_flagger:
    input:
        cov = "output/assembly/flagger/{specimen}/coverage_file.cov.gz",
        alpha_tsv = "config/packages/flagger/alpha_optimum_trunc_exp_gaussian_w_16000_n_50.HiFi_DC_1.2_DEC_2024.v1.1.0.tsv"
    output:
        "output/assembly/flagger/{specimen}/prediction_summary_final.tsv"
    threads: 16
    resources:
        mem_mb = 32000
    params:
        base_dir = config['workdir']
    shell:
        """        
        singularity exec \
            --bind {params.base_dir}:{params.base_dir} \
            docker://mobinasri/flagger:v1.1.0 \
            hmm_flagger \
                --input {input.cov} \
                --outputDir $(dirname {output}) \
                --alphaTsv {input.alpha_tsv} \
                --labelNames Err,Dup,Hap,Col \
                --threads {threads}
        """

rule summarize_predictions:
    input:
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = [x for x in specimens if x != '900'])
    output:
        "output/assembly/flagger/all_prediction_summary_final.tsv"
    threads: 1
    shell:
        """
        echo -e "specimen\tErr\tDup\tHap\tCol\tUnk" > {output}
        find output/assembly/flagger -name "prediction_summary_final.tsv" | while read file; do
            specimen=$(echo "$file" | awk -F'/' '{{print $(NF-1)}}')
            awk -v specimen="$specimen" '$1=="PREDICTION" && $2=="base_level" && $3=="percentage" && $4=="annotation" && $5=="whole_genome" {{print specimen "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}}' "$file" >> {output}
        done
        """