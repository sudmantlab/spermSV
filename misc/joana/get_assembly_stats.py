import numpy as np
import pandas as pd
import os 
gsizes = [3018593705] # Clint contigs
#gsizes = [3073728508] # bonobo contigs
#gsizes = [3054815472] # t2t

sample_table = pd.read_csv('/global/scratch/users/joana_rocha/PANPAN/pepsamples.tsv', '\t')
print(sample_table)
sample_table = sample_table.set_index("SampleID", drop=False)
print(sample_table)
gnames = sample_table['Specimen'].unique()
print(gnames)
#gnames = ["GCF_002880755.1_Clint_PTRv2_genomic", "GCA_013052645.3_Mhudiblu_PPA_v2", "GCF_009914755.1_T2T-CHM13v2.0_genomic", "mGorGor1", "mPanPan1", "mPanTro3"]
haps = ['.hap1.', '.hap2.', '.' ]
#haps = ['.']

rule all:
    input:
        expand('assembly_stats/{gname}{hap}p_ctg_{gsize}.csv', gname=gnames,hap=haps, gsize=gsizes),
        'assembly_stats/Pan_updated_Final.tsv'
       
rule run_assembly_stats_hifiasm_fastas:
    input: 'Hifiasm-fasta_shortcut//{gname}{hap}p_ctg.fa'
    output: 'assembly_stats/{gname}{hap}p_ctg_{gsize}.csv'
    shell: '/global/scratch/users/joana_rocha/software/assemblystats/target/release/assemblystats {input} {output} --genomename {wildcards.gname} --genomesize {wildcards.gsize}'

rule concat_all_csvs:
    input: expand('assembly_stats/{gname}{hap}p_ctg_{gsize}.csv', gname=gnames, hap=haps, gsize=gsizes)
    output: 'assembly_stats/Pan_updated_Final.tsv'
    run: 
       df=pd.concat([pd.read_csv(path, sep= '\t')  for path in input], ignore_index=True)
       df['type'] = 'primary'
       df.loc[df.genomeName.str.contains('.hap1'), 'type'] = 'hap1'
       df.loc[df.genomeName.str.contains('.hap2'), 'type'] = 'hap2'
       df.genomeName = df.genomeName.str.split('.').str[0]
       print(df)
       df.to_csv(output[0], sep="\t", index=False)

    
