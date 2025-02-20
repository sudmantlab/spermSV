# HAPLOTYPES = ["hap1", "hap2"]
# GENOME_SIZES = [3100000000]

rule run_assembly_stats:
    input:
        "output/assembly/hifiasm/{specimen}/{specimen}.{hap}.fa"
    output:
        'output/assembly/assembly_stats/{specimen}.{hap}.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename {wildcards.specimen}.{wildcards.hap} \
        --genomesize {params.genomesize}
        """

rule run_reference_stats:
    input:
        "/global/scratch/users/stacy-l/references/hg38_HGSVC/hg38.no_alt.fa"
    output:
        'output/assembly/assembly_stats/hg38.ref.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename hg38.ref \
        --genomesize {params.genomesize}
        """

rule run_reference_hap_stats:
    input:
        "/global/scratch/users/stacy-l/references/HG002/hg002v1.0.1_{hap}.fasta"
    output:
        'output/assembly/assembly_stats/T2T-HG002.{hap}.tsv'
    params:
        genomesize = 3100000000
    shell:
        """
        /global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats \
        {input} {output} \
        --genomename T2T-HG002.{wildcards.hap} \
        --genomesize {params.genomesize}
        """

rule combine_all_stats:
    input:
        expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = [x for x in specimens if x != '901'], hap = ['hap1', 'hap2'])
        # expand("output/assembly/assembly_stats/{specimen}.{hap}.tsv", specimen = glob_wildcards("output/assembly/assembly_stats/{specimen}.hap1.tsv").specimen, hap = ['hap1', 'hap2'])
    output:
        'output/assembly/assembly_stats/all.tsv'
    run:
        import pandas as pd
        import glob
        import os

        # Get all matching tsv files at runtime
        input_pattern = "output/assembly/assembly_stats/*.*.tsv"
        all_files = glob.glob(input_pattern)
        
        if not all_files:
            raise ValueError(f"No files found matching the pattern: {input_pattern}")

        # Read and concatenate all files
        df_list = []
        for file in all_files:
            try:
                df = pd.read_csv(file, sep='\t')
                df_list.append(df)
            except pd.errors.EmptyDataError:
                print(f"Warning: {file} is empty and will be skipped.")
        
        if not df_list:
            raise ValueError("No valid data found in any of the input files.")
        
        combined_df = pd.concat(df_list, ignore_index=True)
        
        # Write the combined dataframe to the output file
        combined_df.to_csv(output[0], sep='\t', index=False)